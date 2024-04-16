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

struct CTbinEdge
{
    VertexID id;
    int block;// neighbor block id = worker id
    int order;// id order
};

ibinstream& operator<<(ibinstream& m, const CTbinEdge& v)
{
    m << v.id;
    m << v.block;
    m << v.order;
    return m;
}

obinstream& operator>>(obinstream& m, CTbinEdge& v)
{
    m >> v.id;
    m >> v.block;
    m >> v.order;
    return m;
}

struct CTbinMsg
{
    VertexID id;
    int block;
    int order;
};

ibinstream& operator<<(ibinstream& m, const CTbinMsg& v)
{
    m << v.id;
    m << v.block;
    m << v.order;

    return m;
}

obinstream& operator>>(obinstream& m, CTbinMsg& v)
{
    m >> v.id;
    m >> v.block;
    m >> v.order;
    return m;
}

//====================================

struct CTbinValue
{
    vector<CTbinEdge> edges;
    int order;
    int split;
};

ibinstream& operator<<(ibinstream& m, const CTbinValue& v)
{
    m << v.edges;
    m << v.order;
    m << v.split;
    return m;
}

obinstream& operator>>(obinstream& m, CTbinValue& v)
{
    m >> v.edges;
    m >> v.order;
    m >> v.split;
    return m;
}

map<VertexID, int> id2od;//id to order

class CTbinVertex : public BVertex<VertexID, CTbinValue, CTbinMsg> {//<id, vertex value, message info>
public:

    virtual void compute(MessageContainer& messages)
    {

        if (step_num() == 1)//first time
        {
        }
        else if(step_num() == 2)
        {
            for (int i = 0; i < messages.size(); i++)
            {
                int start = value().split;
                int end = value().edges.size();
                for (int j = start; j < end; ++j) {
                    if(value().edges[j].id == messages[i].id){
                        value().edges[j].block = messages[i].block;
                        value().edges[j].order = messages[i].order;
                    }
                }
            }

            //reorder for all nbrs
            sort(value().edges.begin(), value().edges.begin()+value().split, [](CTbinEdge& x, CTbinEdge& y){
                return x.order < y.order;
            });
            sort(value().edges.begin()+value().split, value().edges.end(), [](CTbinEdge& x, CTbinEdge& y){
                if(x.block == y.block) return x.order < y.order;
                return x.block < y.block;
            });

            vote_to_halt();
        }
        else{
            cout << "error v compute step " << step_num() << endl;
            exit(-1);
        }
    }

};


class CTbinBlock : public Block<CTbinValue, CTbinVertex, CTbinMsg> {//<vertex value, vertex(compute), msgtype>
public:

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes)
    {

        if (step_num() == 1)//
        {
            //update order
            for (int i = 0; i < vertexes.size(); ++i) {//for each vertex
                //in-block processing
                auto &uVertex = vertexes[i];
                auto& edges = uVertex->value().edges;
                int local_size = uVertex->value().split;

                for (int j = 0; j < local_size; j++) // for each local edge
                {
                    vertexes[i]->value().edges[j].order = id2od[vertexes[i]->value().edges[j].id];
                }

                //out-block msg passing
                int start = uVertex->value().split;
                int end = uVertex->value().edges.size();
                for (int j = start; j < end; j++)//send msg according to cut edges
                {
                    CTbinEdge &v = edges[j];
                    CTbinMsg msg;
                    msg.id = uVertex->id;
                    msg.block = get_worker_id();
                    msg.order = uVertex->value().order;
                    uVertex->send_message(v.id, v.block, msg);//send dis to its neighbors
                }
            }
            vote_to_halt();
        }
        else if (step_num() == 2)
        {
            vote_to_halt();
            return;
        }
        else{
            cout << "error step " <<step_num() << endl;
            exit(-1);
        }
    }
};

class CTbinBlockWorker : public BWorker<CTbinBlock> {
    char buf[1000];
    vector<int> vertex_dist;

public:

    void load_dist(string dist_path){
        FILE* in = getRHandle(dist_path.c_str());
        LineReader reader(in);
        bool first = true;
        vertex_dist.clear();
        vertex_dist.push_back(-1);

        while (true)
        {
            reader.readLine();
            if (!reader.eof()){
                auto line = reader.getLine();
                int block = atoi(line);
                vertex_dist.push_back(block);
            }
            else
                break;
        }
        fclose(in);
    }


    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)
    {
        //get m, n
        int n = vertexes.size();
        int m = 0;
        for(auto& v : vertexes) m += v->value().edges.size();

        //sort vertexes by deg
        if(get_worker_id() == MASTER_RANK) cout << _my_rank << " reordering ..." << endl;

        auto cmp = [](const CTbinVertex* x, const CTbinVertex* y)
        {
            auto cut_cnt1 = x->value().edges.size() - x->value().split;
            auto cut_cnt2 = y->value().edges.size() - y->value().split;
            if(cut_cnt1 == cut_cnt2)
                return x->value().edges.size() > y->value().edges.size();// no cut edge, compare edge cnt
            return cut_cnt1 > cut_cnt2; // at least one of them has cut edge, compare cut edge cnt
        };

        sort(vertexes.begin(), vertexes.end(), cmp);// sort by deg
        for( int i = 0; i < vertexes.size(); ++i ) {
            vertexes[i]->value().order = i;
            id2od[vertexes[i]->id] = i;
        }
#ifdef _DEBUG
        cout << _my_rank << " id2od=" ;
        for( int i = 0; i <= 12; ++i ) {
            cout << id2od[i] << " ";
        }
        cout << endl;
#endif
    }


    void load_graph(const char* inpath)
    {
//        tolinecnt = 0;
        string instr(inpath);
        int pos = instr.rfind('/');
        instr = instr.substr(0, pos);

        FILE* in = getRHandle(instr.c_str());
        LineReader reader(in);
        reader.readLine();// m n
        char* line = reader.getLine();
        char* pch;
        pch = strtok(line, " ");
        int n = atoi(pch);
        pch = strtok(NULL, " ");
        int m = atoi(pch);

        for (int idx = 1; idx <= n; ++idx) {
            reader.readLine();
            if (!reader.eof()){
                if(vertex_dist[idx] == get_worker_id()) {
                    auto v = toVertex(reader.getLine(), idx);
                    load_vertex(v);
                }
            }
            else
                break;
        }
        fclose(in);
    }

    virtual CTbinVertex* toVertex(char* line)
    {
        return NULL;
    }

    //C version
    virtual CTbinVertex* toVertex(char* line, int id)
    {
        //nb1 nb2 ....
        CTbinVertex* v = new CTbinVertex;
        v->id = id;
        v->bid = get_worker_id();
        v->wid = get_worker_id();

        auto& edges = v->value().edges;
        char* pch;
        pch = strtok(line, " ");//id
        while (pch)
        {
            CTbinEdge t;
            t.id = atoi(pch);//id id
            t.block = vertex_dist[t.id];
            t.order = -1;
            edges.push_back(t);

            pch = strtok(NULL, " ");
        }
//        sort edges, cut edge shall be at the end of edge list
        vector<CTbinEdge> tmp;
        vector<CTbinEdge> tmp1;
        for (int j = 0; j < v->value().edges.size(); j++)
        {
            if (edges[j].block == get_worker_id())//edge belongs to the same block
            {
                tmp.push_back(edges[j]);
            }
            else
                tmp1.push_back(edges[j]);//belongs to different block
        }
        edges.swap(tmp);
        v->value().split = edges.size();
        edges.insert(edges.end(), tmp1.begin(), tmp1.end());//internal edges at the front of the list, splitted by "split"

        return v;
    }

//    int tolinecnt;

    virtual void toline(CTbinBlock* b, CTbinVertex* v, BufferedWriter& writer)
    {
        //id bid order size csize bid1 id1 bid2 id2...
        sprintf(buf, "%d %d %d %d %d ", v->id, v->bid, v->value().order, v->value().edges.size(), v->value().split);
        writer.write(buf);
        auto& e = v->value().edges;
        for (int i = 0; i < e.size(); ++i) {
            sprintf(buf,"%d %d ", e[i].block, e[i].order);
            writer.write(buf);
        }
//        if(get_worker_id() == MASTER_RANK && tolinecnt < 10) {
//            cout << "save line " << tolinecnt << " id=" << v->id << " order=" << v->value().order << endl;
//        }
//        tolinecnt++;

        writer.write("\n");

    }
};


class CTbinCombiner : public Combiner<CTbinMsg>
{
public:
    virtual void combine(CTbinMsg& old, const CTbinMsg& new_msg)
    {
    }
};

void blogel_CTlabeling_create_bin(string in_path, string out_path, string dist_path)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    CTbinBlockWorker worker;
    worker.set_compute_mode(CTbinBlockWorker::VB_COMP);
    CTbinCombiner combiner;
    worker.setCombiner(&combiner);//maybe it's v combiner

    worker.load_dist(dist_path);
    worker.run(param);
#ifdef _DEBUG
    cout << _my_rank << " create bin finfished" << endl;
#endif
}