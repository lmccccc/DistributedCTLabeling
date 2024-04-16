#include "../Blogel/utils/Combiner.h"
#include "../Blogel/blogel/BVertex.h"
#include "../Blogel/blogel/Block.h"
#include "../Blogel/blogel/BWorker.h"
#include "../Blogel/blogel/BGlobal.h"
#include "../Blogel/blogel/Heap.h"
#include <iostream>
#include <vector>
#include <float.h>
#include <omp.h>
#include <sys/time.h>

using namespace std;
#define MAXD 120

typedef unsigned short tint;
typedef char dint;

string to_string(dint v){
    return to_string((int)v);
}

struct trip{
    int id;
    int cost;
};


struct CTValue//value that each edge contains, in CT this means core label and tree label
{
    int order;
    int rank;

    //tree
    int f, h, rid, rsize, w;
    vector<int> nbrs;// ->order
    vector<int> anc; //of size h
    vector<dint> dis;


    //core

    vector<int> cutnbrs;
    vector<int> cutblock;

    //query
    vector<trip> to;//<worker, id, order, cost>
};

ibinstream& operator<<(ibinstream& m, const CTValue& v)
{
    m << v.order;
    m << v.rank;
    m << v.f;
    m << v.h;
    m << v.rid;
    m << v.rsize;
    m << v.w;
    m << v.nbrs;
    m << v.anc;
    m << v.dis;
    m << v.cutnbrs;
    m << v.cutblock;
    return m;
}

obinstream& operator>>(obinstream& m, CTValue& v)
{
    m >> v.order;
    m >> v.rank;
    m >> v.f;
    m >> v.h;
    m >> v.rid;
    m >> v.rsize;
    m >> v.w;
    m >> v.nbrs;
    m >> v.anc;
    m >> v.dis;
    m >> v.cutnbrs;
    m >> v.cutblock;
    return m;
}

hash_map<VertexID, int> id2od;//id -> order
vector<int> id2wk;//id -> worker id


class CTVertex : public BVertex<VertexID, CTValue, char> {//<id, vertex value, message info>, no need of message
public:

    virtual void compute(MessageContainer& messages)
    {
    }
};

class CTBlock : public Block<CTValue, CTVertex, char> {//<vertex value, vertex(compute), msgtype>
public:

    int query_by_tree(VertexContainer& vertexes, int u, int v) {
        auto &tu = vertexes[u]->value(), &tv = vertexes[v]->value();
        int len = min(tu.h,tv.h);
        int d = INT_MAX;
        //iterate all tree dis that u and v can access, and find the smallest dis sum
        //dis is the tree label
        for(int i = tu.rsize, j = 0; i < len && tu.anc[j] == tv.anc[j]; ++i,++j) {//the same anc
            if (get_worker_id() == MASTER_RANK) {
                cout << " tu.dis[" << i << "]=" << (int)tu.dis[i] << " tv.dis[" << i << "]=" << (int)tv.dis[i] << endl;
            }
            d = min(d, tu.dis[i] + tv.dis[i]);
        }
        return d;
    }

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes)
    {

        for (int i = 0; i < vertexes.size(); ++i) {
            auto& vi = vertexes[i]->value();
            for (int j = 0; j < vi.to.size(); ++j) {
                int to_id = vi.to[j].id;
                if(id2wk[to_id] == get_worker_id()){
                    //tree query
                    int to_order = id2od[to_id];
                    auto& vj = vertexes[to_order]->value();
                    if(vi.rid == vj.rid) {
                        if(get_worker_id() == MASTER_RANK) {
                            cout << " query by tree";
                            cout << " vi=" << vertexes[i]->id << " vj=" << vertexes[to_order]->id << endl;
                        }
                        int res = query_by_tree(vertexes, i, to_order);
                        vi.to[j].cost = res;
                    }

                    // others unimplenmented
                }
            }
        }
        vote_to_halt();
    }
};

class CTBlockWorker : public BWorker<CTBlock> {
    char buf[1000];
    string dist_path;
    int max_w, n_threads;

    vector<int> res;
    int global_max_n;

    vector<char> vertex_dist;
public:
    typedef char dint;

    void set_max_w(int w){
        max_w = w;
    }

    void set_n_threads(int nt){
        n_threads = nt;
    }

    void set_dist_path(string p){
        dist_path = p;
    }


    void load_dist(string dist_path){
        FILE* in = getRHandle(dist_path.c_str());
        LineReader reader(in);
        bool first = true;
        id2wk.clear();
        id2wk.push_back(-1);
        while (true)
        {
            reader.readLine();
            if (!reader.eof()){
                auto line = reader.getLine();
                int block = atoi(line);
                id2wk.push_back(block);
            }
            else
                break;
        }
        fclose(in);
    }

    void construct_query(VertexContainer& vertexes){
        int query_cnt = 0;
        int max_query = 10;
        for (int i = 100; i <= global_max_n; ++i) {

            if(i == 103 && get_worker_id() == MASTER_RANK){
                cout << "i=" << i << " id2wk=" << id2wk[i];
                if(id2od.find(i) == id2od.end()){

                    cout << " error no such i ";
                    cout << "  belongs to " << id2wk[i] << endl;
                }
                else{
                    cout << " order=" << id2od[i] << endl;
                }
            }

            if(id2wk[i] != get_worker_id()) continue;
            for (int j = 1; j <= global_max_n; ++j) {
                //for tree query only
                if(i == j || id2wk[j] != get_worker_id()) continue;

                auto& vi = vertexes[id2od[i]];
                auto& vj = vertexes[id2od[j]];
                if(vi->value().rank != -1 && vj->value().rank != -1 && vi->value().rid == vj->value().rid) {
                    if(get_worker_id() == MASTER_RANK){
                        cout << "query v id=" << i << " order=" << id2od[i] << " j id=" << j << " order=" << id2od[j] << endl;
                    }
                    vi->value().to.push_back({j, -1});
                    query_cnt++;
                    if(query_cnt == max_query) return;
                }
            }
        }
    }


    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)//
    {
        //load node dist file
        int c = 0;
        for (int i = 0; i < vertexes.size(); ++i) {
            if(c < 10 && get_worker_id() == MASTER_RANK) {
                cout << get_worker_id() << " get id " << vertexes[i]->id << " idx=" << i << " order=" << vertexes[i]->value().order << endl;
            }
            ++c;
            id2od[vertexes[i]->id] = i;
        }
        load_dist(dist_path);
        //get tree query
        global_max_n = id2wk.size();
        construct_query(vertexes);
    }


    virtual CTVertex* toVertex(char* line)
    {//id order rank rid rsize h w (order1 block1)*w dis*h anc*(h-w)

        CTVertex* v = new CTVertex;
        char* pch;

        pch = strtok(line, " ");//id
        v->id = atoi(pch);

        pch = strtok(NULL, " ");//bid
        v->wid = v->bid = atoi(pch);

        pch = strtok(NULL, " ");//order
        v->value().order = atoi(pch);

        pch = strtok(NULL, " ");//rank
        v->value().rank = atoi(pch);

        pch = strtok(NULL, " ");//rid
        v->value().rid = atoi(pch);

        pch = strtok(NULL, " ");//rsize
        v->value().rsize = atoi(pch);

        pch = strtok(NULL, " ");//h
        v->value().h = atoi(pch);

        pch = strtok(NULL, " ");//w
        v->value().w = atoi(pch);

        if(v->value().rank != -1) {
            for (int i = 0; i < v->value().w; ++i) {
                pch = strtok(NULL, " ");
                v->value().nbrs.push_back(atoi(pch));

            }
//        pch = strtok(NULL, "\t");

            for (int i = 0; i < v->value().h; ++i) {
                pch = strtok(NULL, " ");
                v->value().dis.push_back(atoi(pch));

            }
//        pch = strtok(NULL, "\t");


            for (int i = 0; i < v->value().h - v->value().w; ++i) {
                pch = strtok(NULL, " ");
                v->value().anc.push_back(atoi(pch));
            }
        }
        return v;
    }

    virtual void toline(CTBlock* b, CTVertex* v, BufferedWriter& writer)
    {
        for (int i = 0; i < v->value().to.size(); ++i) {
            //fromid \t toid \t cost
            auto& info = v->value().to[i];
            sprintf(buf, "%d %d %d \n", v->id, info.id, info.cost);
            writer.write(buf);
//            string s;
//            s = s + to_string(v->id) + "\t";
//            s = s + to_string(info.id) + "\t";
//            s = s + to_string(info.cost) + "\t";
//            writer.write(s.c_str());
//            writer.write("\n");
        }
    }
};



void blogel_ctlabelingquery(string in_path, string out_path, string dist_path, int max_w, int n_threads)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    CTBlockWorker worker;
    worker.set_max_w(max_w);
    worker.set_n_threads(n_threads);
    worker.set_dist_path(dist_path);
    worker.run(param);
}