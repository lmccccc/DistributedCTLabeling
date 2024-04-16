#include "../Blogel/utils/Combiner.h"
#include "../Blogel/blogel/BVertex.h"
#include "../Blogel/blogel/Block.h"
#include "../Blogel/blogel/BWorker.h"
#include "../Blogel/blogel/BGlobal.h"
#include "../Blogel/blogel/Heap.h"

#include <iostream>
#include <vector>
#include <float.h>

using namespace std;

int src = 1;

struct SPEdge
{
    double len;
    int nb;
    int block;
    int worker;
};

ibinstream& operator<<(ibinstream& m, const SPEdge& v)
{
    m << v.len;
    m << v.nb;
    m << v.block;
    m << v.worker;
    return m;
}

obinstream& operator>>(obinstream& m, SPEdge& v)
{
    m >> v.len;
    m >> v.nb;
    m >> v.block;
    m >> v.worker;
    return m;
}

//====================================

struct SPValue
{
    double dist;
    int from;
    vector<SPEdge> edges;
    int split; //v.edges[0, ..., inSplit] are local to block
};

ibinstream& operator<<(ibinstream& m, const SPValue& v)
{
    m << v.dist;
    m << v.from;
    m << v.edges;
    return m;
}

obinstream& operator>>(obinstream& m, SPValue& v)
{
    m >> v.dist;
    m >> v.from;
    m >> v.edges;
    return m;
}

//====================================

struct SPMsg
{
    double dist;
    int from;
};

ibinstream& operator<<(ibinstream& m, const SPMsg& v)
{
    m << v.dist;
    m << v.from;

    return m;
}

obinstream& operator>>(obinstream& m, SPMsg& v)
{
    m >> v.dist;
    m >> v.from;
    return m;
}

//====================================

//set field "active" if v.dist changes
//set back after block-computing
class SPVertex : public BVertex<VertexID, SPValue, SPMsg>
{
public:
    virtual void compute(MessageContainer& messages)
    {
        if (step_num() == 1)//first time
        {
            if (id == src)
            {
                value().dist = 0;
                value().from = -1;
            } //remain active, to be processed by block-computing
            else
            {
                value().dist = DBL_MAX;
                value().from = -1;
                vote_to_halt();
            }
        }
        else
        {
            SPMsg min;
            min.dist = DBL_MAX;
            for (int i = 0; i < messages.size(); i++)
            {
                SPMsg msg = messages[i];
                if (min.dist > msg.dist)//new dist
                {
                    min = msg;
                }
            }
            if (min.dist < value().dist)
            {
                value().dist = min.dist;
                value().from = min.from;
            } //remain active, to be processed by block-computing
            else
                vote_to_halt();//not updated
        }
    }
};

//====================================

class SPBlock : public Block<char, SPVertex, char>
{
public:
    typedef qelem<double, int> elemT;

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes) //multi-source Dijkstra (assume a super src node)
    { //heap is better than queue, since each vertex is enheaped only once
        //collect active seeds
        //cout << bid <<" " << _my_rank << endl;
        heap<double, int> hp;

        int n = size;
        vector<bool> tag(n, false);
        vector<elemT*> eles(n);
        for (int i = 0; i < n; i++)
            eles[i] = NULL; //init eles[]

        for (int i = begin; i < begin + size; i++)
        {
            SPVertex& vertex = *(vertexes[i]);
            if (vertex.is_active())
            {
                double key = vertexes[i]->value().dist; //distance
                int val = i - begin; //logic ID
                eles[val] = new elemT(key, val);
                hp.add(*eles[val]); //1. add active vertex to minheap
                vertex.vote_to_halt(); //2. deactivate the vertex
            }
        }
        //run dijkstra's alg
        while (hp.size() > 0)//for each active vertex
        {
            elemT& u = hp.peek();//ele[1]?
            if (u.key == DBL_MAX)
                break;
            hp.remove();
            tag[u.val] = true;
            int index = begin + u.val;
            SPVertex& uVertex = *(vertexes[index]);
            vector<SPEdge>& edges = uVertex.value().edges;
            int split = uVertex.value().split;
            double udist = uVertex.value().dist;
            //in-block processing
            for (int i = 0; i <= split; i++)//for each local edge
            {
                SPEdge& v = edges[i];
                int logID = v.worker - begin;
                if (tag[logID] == false)
                {
                    double alt = udist + v.len;
                    SPVertex& vVertex = *(vertexes[v.worker]);
                    double& vdist = vVertex.value().dist;
                    if (alt < vdist)//smaller distance
                    {
                        if (eles[logID] == NULL)
                        {
                            eles[logID] = new elemT(alt, logID);
                            hp.add(*eles[logID]);
                        }
                        else
                        {
                            eles[logID]->key = alt;
                            hp.fix(*eles[logID]);
                        }
                        vdist = alt;//update val
                        vVertex.value().from = uVertex.id;//update route
                    }
                }
            }
            //out-block msg passing
            for (int i = split + 1; i < edges.size(); i++)//send msg according to cut edges
            {
                SPEdge& v = edges[i];
                SPMsg msg;
                msg.dist = udist + v.len;
                msg.from = uVertex.id;
                uVertex.send_message(v.nb, v.worker, msg);//send dis to its neighbors
            }
        }
        //free elemT objects
        for (int i = 0; i < n; i++)
        {
            if (eles[i] != NULL)
                delete eles[i];
        }
        vote_to_halt();
    }
};

//====================================

class SPBlockWorker : public BWorker<SPBlock>
{
    char buf[1000];

public:
    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)// sort edges
    {//b0 : v1 v2 v3 v4// b1: v5 v6 // b2: v7
        hash_map<int, int> map;
        for (int i = 0; i < vertexes.size(); i++)//vertexes that the block have
            map[vertexes[i]->id] = i;//id to v's position
        //////
        if (_my_rank == MASTER_RANK)
            cout << "Splitting in/out-block edges ..." << endl;
        for (BlockIter it = blocks.begin(); it != blocks.end(); it++)
        {
            SPBlock* block = *it;
            for (int i = block->begin; i < block->begin + block->size; i++)//begin: v offset, size: v cnt
            {
                SPVertex* vertex = vertexes[i];
                vector<SPEdge>& edges = vertex->value().edges;
                vector<SPEdge> tmp;
                vector<SPEdge> tmp1;
                for (int j = 0; j < edges.size(); j++)
                {
                    if (edges[j].block == block->bid)//edge belongs to the same block
                    {
                        edges[j].worker = map[edges[j].nb]; //workerID->array index
                        tmp.push_back(edges[j]);
                    }
                    else
                        tmp1.push_back(edges[j]);//belongs to different block
                }
                edges.swap(tmp);
                vertex->value().split = edges.size() - 1;
                edges.insert(edges.end(), tmp1.begin(), tmp1.end());//internal edges at the front of the list, splitted by "split"
                /*
                if(vertex->id == 1)
                {
                    for(int j = 0 ; j <= vertex->value().split; j ++)
                    {
                        cout << edges[j].nb << " " << edges[j].block   << endl;
                    }
                    for(int j = vertex->value().split + 1; j < edges.size(); j ++)
                    {
                        cout << edges[j].nb << " " << edges[j].block << " " << edges[j].worker << endl;
                    }
                }
                */
#ifdef _DEBUG
                cout << _my_rank << " vertex " << vertex->id << " edge size=" << edges.size() << " split=" << vertex->value().split << endl;
#endif
            }
        }
        if (_my_rank == MASTER_RANK)
        {
            cout << "In/out-block edges split" << endl;
        }
    }

    //input line format: vid blockID workerID \t nb1 nb2 ...
    //nbi format: vid edgeLength blockID workerID
    virtual SPVertex* toVertex(char* line)
    {
        //id bid \t wid nb1 block1 worker1 nb2 block2 worker2 ...
        char* pch;
        pch = strtok(line, " ");
        SPVertex* v = new SPVertex;
        v->id = atoi(pch);
        pch = strtok(NULL, " ");
        v->bid = atoi(pch);
        pch = strtok(NULL, " ");
        v->wid = atoi(pch);
        vector<SPEdge>& edges = v->value().edges;
		while (pch = strtok(NULL, " "))
        {
            SPEdge trip;
            trip.nb = atoi(pch);
//            pch = strtok(NULL, " ");
//            trip.len = atof(pch);
            trip.len = 1;
            pch = strtok(NULL, " ");
            trip.block = atoi(pch);
            pch = strtok(NULL, " ");
            trip.worker = atoi(pch);
            edges.push_back(trip);
        }
#ifdef _DEBUG
        cout << _my_rank << " get vertex id=" << v->id << " block=" << v->bid << " worker=" << v->wid <<
            " edge size=" << edges.size() << endl;
#endif
        ////////
        if (v->id == src)// source
        {
            v->value().dist = 0;
            v->value().from = -1;
        }
        else// halt at first?
        {
            v->value().dist = DBL_MAX;
            v->value().from = -1;
            v->vote_to_halt();
        }
        return v;
    }

    virtual void toline(SPBlock* b, SPVertex* v, BufferedWriter& writer)
    {
        if (v->value().dist != DBL_MAX)
        {
            sprintf(buf, "%d\t%f %d", v->id, v->value().dist, v->value().from);
            writer.write(buf);
        }
        writer.write("\n");
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
};

class SPCombiner : public Combiner<SPMsg>
{
public:
    virtual void combine(SPMsg& old, const SPMsg& new_msg)
    {

#ifdef _DEBUG
        cout <<  " self defined combine" << endl;
#endif
        if (old.dist > new_msg.dist)
            old = new_msg;
    }
};

void blogel_app_sssp(string in_path, string out_path)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    SPBlockWorker worker;
    worker.set_compute_mode(SPBlockWorker::VB_COMP);
    SPCombiner combiner;
    worker.setCombiner(&combiner);//maybe it's v combiner
    worker.run(param);
}

//override: Bvertex->Block->Bworker->Combiner
//vertex: implement value, message, override compute
//Block:compute
//Worker: BlockInit, toVertex, toLine
//Combiner: combine, assined as a parameter