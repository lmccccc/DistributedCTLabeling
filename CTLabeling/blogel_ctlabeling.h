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


struct CTValue//value that each edge contains, in CT this means core label and tree label
{
    int order;
//    int split; //v.edges[0, ..., inSplit] are local to block
    int rank;


    int f, h, rid, rsize, w;
    vector<int> nbrs;// ->order
    vector<int> cost;
    vector<int> anc; //of size h
    vector<dint> dis;
    vector<int> ch;

    vector<int> cutnbrs;
    vector<int> cutblock;
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
    m << v.cost;
    m << v.anc;
    m << v.dis;
    m << v.ch;
    m << v.cutnbrs;
    m << v.cutblock;
    return m;
}

obinstream& operator>>(obinstream& m, CTValue& v)
{
    m >> v.order;
    m >> v.rank;
    m >> v.rank;
    m >> v.f;
    m >> v.h;
    m >> v.rid;
    m >> v.rsize;
    m >> v.w;
    m >> v.nbrs;
    m >> v.cost;
    m >> v.anc;
    m >> v.dis;
    m >> v.ch;
    m >> v.cutnbrs;
    m >> v.cutblock;
    return m;
}



class CTVertex : public BVertex<VertexID, CTValue, char> {//<id, vertex value, message info>, no need of message
public:

    virtual void compute(MessageContainer& messages)
    {
        vote_to_halt();
    }
};

class CTBlock : public Block<CTValue, CTVertex, char> {//<vertex value, vertex(compute), msgtype>
public:

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes)
    {
        vote_to_halt();
    }
};

class CTBlockWorker : public BWorker<CTBlock> {
    char buf[1000];
    int max_w, n_threads;
    int n, m, n_core, m_core;
    vector<int> score, _score, ord;
    vector<vector<int>> nbr, cost;
    vector<vector<pair<int, int>>> E;//<order, dis>

public:
    typedef char dint;

    class IntVal{
    public:
        int x;
    public:
        IntVal(){
            x=-1;
        }

        IntVal(int x){
            this->x=x;
        }

        //Lambda设置map排序规则,同样适合于set

//        bool operator<(const IntVal &v) const {
//            if(score[x] == score[v.x])
//                return x<v.x;
//            return score[x] < score[v.x];
//        }
    };

    void set_max_w(int w){
        max_w = w;
    }

    void set_n_threads(int nt){
        n_threads = nt;
    }


    //reduce: decompose the graph into core and trees.
    //input: graph
    //output: Core Tree
    void reduce(VertexContainer& vertexes, int max_w, int n_threads) {// reduce: decomposition core step
        omp_set_num_threads(n_threads);
        double t = omp_get_wtime();
        this->max_w = max_w;
        n = vertexes.size();

        cout << get_worker_id() << " n= " << n << endl;
        score.resize(n);
        _score.resize(n);
        nbr.resize(n);
        cost.resize(n);


        vector<bool> changed(n,false);
        for(int i = 0; i < n; ++i)
            score[i] = _score[i] = vertexes[i]->value().nbrs.size();



        auto cmp = [&](const IntVal& left, const IntVal& right) {
            if(score[left.x] == score[right.x])
                return left.x<right.x;
            return score[left.x] < score[right.x];
        };

        set<IntVal, decltype(cmp)> q(cmp);// ordered quque?

//        nbr.resize(n); cost.resize(n);

        int r = 0;

        for(int i = 0; i < n; ++i) score[i] = _score[i] = vertexes[i]->value().nbrs.size();

//        printf( "%d, Initializing q...", get_worker_id() );
        cout << get_worker_id() << " Initializing q... " << endl;
        vector<bool> active(n,false);
        for(int u = 0; u < n; ++u) {
            int deg = vertexes[u]->value().nbrs.size();
            if (deg < max_w && vertexes[u]->value().cutnbrs.size() == 0) {
                q.insert(IntVal(u));

                active[u] = true;
            } // nodes whose degree less than max_w are active
//            if(deg-1 != split) cutV.push_back(u);
        }
        cout << get_worker_id() << " t=" << omp_get_wtime()-t << " secs" << endl;

        cout << get_worker_id() << " Initializing E..." << endl;
//             printf( ", t=%0.3lf secs\nInitializing E...", omp_get_wtime()-t );
        E.resize(n);//nodes and its neighbors

        for(int u = 0; u < n; ++u) {
            auto& unbrs = vertexes[u]->value().nbrs;
            int deg = unbrs.size();
            for (int i = 0; i < deg; ++i) {
                int nbod = unbrs[i];
                E[u].push_back(make_pair(unbrs[i],1));// all neighbors, only internal edges
            }
        }
        cout << get_worker_id() << " t=" << omp_get_wtime()-t << " secs" << endl;
        cout << get_worker_id() << " Reducing Graph..." << endl;
//             printf( ", t=%0.3lf secs\nReducing Graph...\n", omp_get_wtime()-t );
        int cnt = 0;
        vector<pair<int,int>> tmp;
        //每次选deg最小的node,加入nbr,E删除与该node连接的边,但是添加到自己的nbr中. rank记录标记成tree node的顺序
#ifdef _DEBUG
        for (int l = 0; l < get_num_workers(); ++l) {
            if(l == get_worker_id()) {
                cout << get_worker_id() << " q=";
                set<IntVal>::iterator it;
                for (it = q.begin(); it != q.end(); it++) {
                    cout << it->x << " ";
                }
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        int circle = 0;
#endif
        while(!q.empty()) {
            int x = q.begin()->x;// x is the core of this iteration, from q which is assigned as active
            if(vertexes[x]->value().cutnbrs.size() > 0){
                cout << "error, invalid node with cut edge, x=" << x<< endl;
                exit(-1);
            }
#ifdef _DEBUG
            for (int l = 0; l < get_num_workers(); ++l) {
                if(l == get_worker_id()){
                    cout << get_worker_id() << " circle=" << circle <<
                        " removing edge " << x << " order=" << vertexes[x]->value().order <<
                        " nbr size=" << nbr[x].size() << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            circle++;
#endif
            while(changed[x]) { // get the first that not changed
                q.erase(x);
                score[x] = _score[x];
                q.insert(x);// move top to the tail if changed, and update score
                changed[x] = false;
                x = q.begin()->x;
            }

            //may never used
            if(score[x] >= max_w) break; // key step, break if any of which scores higher than max, means CT decomposition is finished

            ord.push_back(x);
            q.erase(x);// pop x

            // x is no longer core node, order that assigned as tree
            vertexes[x]->value().rank = r;
            r++;

            for(auto &it:E[x]){nbr[x].push_back(it.first); cost[x].push_back(it.second);}// collect active nodes

            for(auto &y:nbr[x]) {// for each neighbor of all nodes
                if(E[y].size() >= max_w * 2) {active[y] = false; q.erase(y);} // remove some nodes
                if(!active[y]) continue;
                for(size_t i=0;i<E[y].size();++i)
                    if(E[y][i].first == x) {E[y].erase(E[y].begin()+i); break;} // X is y's neighbor, then remove E[y][x]
                _score[y] = (int) E[y].size();// update size
                changed[y] = true;// mark changed
            }

            //update E, remove tree neighbor, add correspond edges, update its score
            for(size_t i = 0; i < nbr[x].size(); ++i) {
                int u = nbr[x][i];

                if(!active[u]) {
                    E[u].push_back(make_pair(x,-cost[x][i]));// negative means not active, is core?
                    continue;
                }
                tmp.clear();
                size_t j=0, k=0;
                while(j<nbr[x].size()&&k<E[u].size())
                    if(j==i) ++j;
                    else if(nbr[x][j]<E[u][k].first) {tmp.push_back(make_pair(nbr[x][j], cost[x][i]+cost[x][j])); ++j;}
                    else if(nbr[x][j]>E[u][k].first) {tmp.push_back(E[u][k]); ++k;}
                    else {
                        if(E[u][k].second < cost[x][i]+cost[x][j]) tmp.push_back(E[u][k]);
                        else tmp.push_back(make_pair(nbr[x][j], cost[x][i]+cost[x][j]));
                        ++j; ++k;
                    }
                for(;j<nbr[x].size();++j) if(j!=i) tmp.push_back(make_pair(nbr[x][j], cost[x][i]+cost[x][j]));
                for(;k<E[u].size();++k) tmp.push_back(E[u][k]);
                E[u] = tmp;// update E
                if(_score[u] != (int) E[u].size()) {
                    changed[u] = true;
                    _score[u] = (int) E[u].size();
                }
            }

            if((++cnt) * score[x] > 1000000) {
                printf( "%d nodes reduced, score[x]=%d, remaining size=%0.3lf%% t=%0.3lf secs\n",
                        r, (n-r)*100.0/n, score[x], omp_get_wtime()-t);
                cnt = 0;
            }
#ifdef _DEBUG
            for (int l = 0; l < get_num_workers(); ++l) {
                if(l == get_worker_id()){
                    cout << get_worker_id() << " new score=";
                    for (int i = 0; i < n; ++i) {
                        cout << score[i] << " ";
                    }
                    cout << endl;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
            circle++;
#endif
        }
#ifdef _DEBUG
        for(int i = 0; i < get_num_workers(); ++i) {
            if(_my_rank == i) {
                cout << _my_rank << " after reduce E=" << endl;
                for (int i = 0; i < E.size(); ++i) {
                    cout << "E[" << i << "] = ";
                    for (int j = 0; j < E[i].size(); ++j) {
                        cout << " (" << E[i][j].first << "," << E[i][j].second << ") ";
                    }
                    cout << endl;
                }
                cout << endl;

                cout << "ord=";
                for (int j = 0; j < ord.size(); ++j) {
                    cout << ord[j] << " ";
                }
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif

        printf( "Reordering edges...\n" );

#pragma omp parallel
        {
            vector<int> ve;
            vector<int> buf(n,-1);

#pragma omp for schedule(dynamic)
            for(int u = 0; u < n; ++u)    //for each node with more than 1 degree
                if(!active[u] && E[u].size()>0) { // if not active and
                    auto &e = E[u];
                    ve.clear();
                    for(size_t i = 0; i < e.size(); ++i) {
                        if(e[i].second >= 0) {
                            int v = e[i].first, w = e[i].second;
                            int rankv = vertexes[v]->value().rank;
                            if(rankv>=0) continue;
                            else if(buf[v] == -1) {buf[v] = w; ve.push_back(v);}
                            else if(w < buf[v]) buf[v] = w;
                        } else {
                            auto &s = E[e[i].first];
                            for(size_t j=0; j<s.size(); ++j) {
                                int v = s[j].first, w = s[j].second - e[i].second;
                                int rankv = vertexes[v]->value().rank;
                                if(v == u || rankv>=0) continue;
                                else if(buf[v] == -1) {buf[v] = w; ve.push_back(v);}
                                else if(w < buf[v]) buf[v] = w;
                            }
                        }
                    }
                    e.resize(ve.size());
                    for(size_t i = 0; i < ve.size(); ++i) {
                        e[i]=make_pair(ve[i], buf[ve[i]]);
                        buf[ve[i]] = -1;
                    }
                    sort(e.begin(), e.end());
                }
        }

        n_core = 0;
        m_core = 0;
        for(int u = 0; u < n; ++u)
            if( vertexes[u]->value().rank == -1 ) {//-1 means core node
                ++n_core;
                m_core += (int) E[u].size();
            }

        printf( "%d Reducing finished, t=%0.3lf secs\nn_core=%d,m_core=%lld,node_rate=%0.3lf,edge_rate=%0.3lf\n",
                _my_rank, omp_get_wtime()-t, n_core, m_core, n_core*1.0/n, m_core*1.0/m );

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void create_tree(VertexContainer& vertexes) {
        printf( "Creating Tree...\n" );
        double t = omp_get_wtime();
        for(int u = 0; u < n; ++u) {
            vertexes[u]->value().order = u;//?seems to be useless
            vertexes[u]->value().nbrs.clear();
            vertexes[u]->value().cost.clear();
        }
        vector<pair<int,int>> v_pair;
        int maxh = 0, cnt_root = 0, maxdep = 0, max_sub_tree = 1;
        vector<int> tcnt(n,0);
        double tw = 0;
#ifdef _DEBUG
        cout << "ord = ";
        for (int i = 0; i < ord.size(); ++i) {
            cout << ord[i] << " ";
        }
        cout << endl;
#endif
        for(int i = (int) ord.size()-1; i >= 0; --i) {//
            int x = ord[i];//get tree node
#ifdef _DEBUG
        if(_my_rank == MASTER_RANK)
            cout << "x=" << x << endl;
#endif
            CTValue &tn = vertexes[x]->value();
            v_pair.clear();
            for(int j = 0; j < (int) nbr[x].size(); ++j) {// nbr used
                int y = nbr[x][j];
                int ranky = vertexes[y]->value().rank;
                if(ranky == -1) v_pair.push_back(make_pair(n,j));// core neighbor =n
                else v_pair.push_back(make_pair(ranky,j));// tree neighbor =rank
            }
            sort(v_pair.begin(),v_pair.end());// core in the tail
            reverse(v_pair.begin(), v_pair.end());//decrease order
            int w = (int) nbr[x].size();
            tn.nbrs.resize(w);
            tn.cost.resize(w);
            for(int j=0; j<w; ++j) {
                tn.nbrs[j] = nbr[x][v_pair[j].second];
                tn.cost[j] = cost[x][v_pair[j].second];
            }
            tn.w = w;//width of the tree node
            tn.order = x;//id
            tn.f = -1;//father?
            for(auto &u:nbr[x]) {
                int ranku = vertexes[u]->value().rank;
                if (ranku != -1 && (tn.f == -1 || ranku < vertexes[tn.f]->value().rank))// tree node and (not changed or smaller rank)
                    tn.f = u;//neighbors' smallest rank
            }
            if(tn.f == -1) {//no father
                tn.h = tn.w + 1;
                ++cnt_root;
                ++tcnt[x];
                tn.rid = x;
                tn.rsize = tn.w;
                tn.anc.push_back(x);//itself id
            } else {
                CTValue& treef = vertexes[tn.f]->value();
                tn.h = treef.h+1;
//                tn.h = treef.h - treef.w + tn.w + 1;
                treef.ch.push_back(x);//children
                tn.rid = treef.rid;
                ++tcnt[tn.rid];
                max_sub_tree = max(max_sub_tree, tcnt[tn.rid]);
                tn.rsize = treef.rsize;
                tn.anc = treef.anc;
                tn.anc.push_back(x);//father's anc and together with the id itself
            }
            tw += tn.rsize;
            maxh = max(maxh, tn.h);
            maxdep = max(maxdep, (int)tn.anc.size());
        }

#ifdef _DEBUG
        for(int k = 0; k < _num_workers; ++k) {
            if(k == _my_rank) {
                cout << _my_rank << " tree = " << endl;
                for (int i = 0; i < vertexes.size(); ++i) {
                    cout << i << " ";
                    cout << "id=" << vertexes[i]->id << " ";
                    cout << "ori order=" << vertexes[i]->value().order << " ";
                    cout << "new order=" << i << " ";
                    cout << "w=" << vertexes[i]->value().w << " ";// width
                    cout << "rank=" << vertexes[i]->value().rank << " ";//rank
                    cout << "f=" << vertexes[i]->value().f << " ";//father node
                    cout << "h=" << vertexes[i]->value().h << " ";//node height
                    cout << "rid=" << vertexes[i]->value().rid << " ";//id that is core and connected to this tree
                    cout << "rsize=" << vertexes[i]->value().rsize
                         << " ";//width that is a core and connected to this tree
                    cout << endl;

                    cout << "   neighbor size " << vertexes[i]->value().nbrs.size() << " " << vertexes[i]->value().cost.size() << " ";
                    for (int j = 0; j < vertexes[i]->value().nbrs.size(); ++j) {
                        cout << " (" << vertexes[i]->value().nbrs[j] << "," << vertexes[i]->value().cost[j] << ") ";
                    }
                    cout << endl << endl;
                }
                cout << endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif



        printf( "Core tree constructed, maxh=%d, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%0.3lf, t=%0.3lf secs\n",
                maxh, maxdep, cnt_root, max_sub_tree, tw/(n-n_core), omp_get_wtime()-t);
    }

    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)//
    {
        n = vertexes.size();
        reduce(vertexes, max_w, n_threads);
        create_tree(vertexes);
        compute_tree_label(vertexes);

        //each node: rank, rid, rsize, h, w, nbr*w, dis*h, anc*(h-w)
        if(_my_rank == MASTER_RANK) {
            sleep(1);
            cout << "tee label=" << endl;
            for(int i = 0; i < vertexes.size(); ++i){
                auto& tn = vertexes[i];
                cout << "i=" << i << " "
                    << tn->id << " "
                    << tn->value().order << " "
                    << tn->value().rank << " "
                    << tn->value().rid << " "
                    << tn->value().rsize << " "
                    << tn->value().h << " "
                    << tn->value().w << " "
                    << tn->value().f <<  "\t";
                for (int j = 0; j < tn->value().w; ++j) cout << tn->value().nbrs[j] << " ";
                cout << "\t";
                for (int j = 0; j < tn->value().h - tn->value().w; ++j) cout << tn->value().anc[j] << " ";
                cout << "\t";
                for (int j = 0; j < tn->value().h; ++j) cout << (int)tn->value().dis[j] << " ";
                cout << "\t";
                cout << endl;
            }
            cout << endl;
        }
    }



    void compute_tree_label(VertexContainer& vertexes, int x, int rsize, vector<CTVertex*> &s) {//
        s.push_back(vertexes[x]); // node itself
        auto &tn = vertexes[x]->value();
        tn.dis.resize(tn.h);// to each layer's distance?
        int pos = 0;
        vector<int> p(tn.w);
        for(int j = 0; j < tn.w; ++j) { //
            while(s[pos]->value().order != tn.nbrs[j]) ++pos;//find neighbor position in s
            p[j] = pos;
        }

        for(int i = 0; i < tn.h-1; ++i) {
            tn.dis[i] = -1;//tn到i的距离?
            for(int j = 0; j < tn.w; ++j) {
                //dis to this layer is the smallest dis to all neighbors of this node
                int w = tn.cost[j], k = p[j], nowdis = -1;//w: each neighbor cost, k: neighbor position in s, nowdis:
                if(k<=i) {// neighbor position <= height
                    if(i>=rsize)//height >= connected core node size
                        nowdis = s[i]->value().dis[k];
                    else if(k==i)
                        nowdis=0;
                }
                else if(k>=rsize) nowdis = s[k]->value().dis[i];
                if(nowdis>=0 && (tn.dis[i]==-1 || nowdis+w<tn.dis[i])) tn.dis[i]=min(nowdis+w,MAXD);//MAXD=120
            }
        }
        tn.dis[tn.h-1] = 0;//to last layer dis=0
        for(int &u:vertexes[x]->value().ch) // for each children
            compute_tree_label(vertexes, u, rsize, s); // compute label. what is label? dis!
        s.pop_back();
    }
    
    void compute_tree_label(VertexContainer& vertexes) {// tree label
        printf( "Computing Tree Label...\n" );
        double t = omp_get_wtime();
        vector<CTVertex*> s;
        for(int v=0; v<n; ++v){//each v
            if(vertexes[v]->value().rank >= 0 && vertexes[v]->value().f == -1) {// not core
                s.clear();
                for(int i = 0; i < vertexes[v]->value().w; ++i) {
                    auto& nb = vertexes[v]->value().nbrs[i];
                    s.push_back(vertexes[nb]);//all node elements
                }
                compute_tree_label(vertexes, v, vertexes[v]->value().rsize, s);//construct dis for tree[v] and its children
            }
        }
        double t_size = 0;
        int maxdis = 0;
        for(int v=0; v<n; ++v) {
            if(vertexes[v]->value().rank >= 0) {
                t_size += vertexes[v]->value().dis.size() * 1.0 * (sizeof(int)+sizeof(dint));//for mem cnt
                vector<pair<int, int>>().swap(E[v]);//clear E
                for(auto &d:vertexes[v]->value().dis) maxdis = max(maxdis, (int)d);//max distance
            } else vector<pair<int,int>>(E[v]).swap(E[v]);//clear E
        }

#ifdef _DEBUG
        if(_my_rank == MASTER_RANK) {
            cout << "tee dis=" << endl;
            for (int i = 0; i < n; ++i) {
                cout << i << " dis=";
                for (int j = 0; j < vertexes[i]->value().dis.size(); ++j) {
                    cout << (int) vertexes[i]->value().dis[j] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
#endif

        printf( "Tree Label Computed, t=%0.3lf secs, maxdis=%d, tree label size=%0.3lf MB\n", omp_get_wtime()-t, maxdis, t_size/(1024.0*1024.0));
    }

    virtual CTVertex* toVertex(char* line)
    {
        //id bid split \t nbid nbblock nborder \t nbid2 nbblock2 nborder2 \t ...
        char* pch;
        pch = strtok(line, " ");//id
        CTVertex* v = new CTVertex;
        v->value().rank = -1;
        v->id = atoi(pch);
        pch = strtok(NULL, " ");//nid/bid
        v->bid = atoi(pch);
        v->wid = v->bid;

        pch = strtok(NULL, " ");//order
        v->value().order = atoi(pch);
        pch = strtok(NULL, "\t");//split
        int split = atoi(pch);
//        v->value().split = atoi(pch);

        //order is empty
        for (int i = 0; i <= split; ++i) {
            pch = strtok(NULL, " ");//id

            pch = strtok(NULL, " ");//bid/nid

            pch = strtok(NULL, "\t");//order
            v->value().nbrs.push_back(atoi(pch));
        }

        while (pch = strtok(NULL, " ")) // neighbor id
        {
            pch = strtok(NULL, " ");//bid/nid
            int block = atoi(pch);
            pch = strtok(NULL, "\t");//order
            int order = atoi(pch);
            v->value().cutnbrs.push_back(order);
            v->value().cutblock.push_back(block);
        }
//
//#ifdef _DEBUG
//        cout << _my_rank;
//        cout << " get vertex id=" << v->id;
//        cout << " order=" << v->value().order << " block=" << v->bid;
//        cout << " local edge size=" << v->value().nbrs.size() << " cut edge size=" << v->value().cutnbrs.size() << endl;
//#endif
        return v;
    }

    virtual void toline(CTBlock* b, CTVertex* v, BufferedWriter& writer)
    {//each node: rank, rid, rsize, h, w, nbr*w, dis*h, anc*(h-w)
        
        string s;//id order rank rid rsize h w (order1 block1)*w dis*h anc*(h-w)
        s = s + to_string(v->id) + " ";
        s = s + to_string(v->bid) + " ";
        s = s + to_string(v->value().order) + " ";
        s = s + to_string(v->value().rank) + " ";
        s = s + to_string(v->value().rid) + " ";
        s = s + to_string(v->value().rsize) + " ";
        s = s + to_string(v->value().h) + " ";
        s = s + to_string(v->value().w) + " ";
        for (int i = 0; i < v->value().w; ++i) {
            s = s + to_string(v->value().nbrs[i]) + " ";
        }
        for (int i = 0; i < v->value().h; ++i) {
            s = s + to_string(v->value().dis[i]) + " ";
        }
        for (int i = 0; i < (v->value().h - v->value().w); ++i) {
            s = s + to_string(v->value().anc[i]) + " ";
        }
        writer.write(s.c_str());
        writer.write("\n");
    }
};

//class CTCombiner : public Combiner<char>
//{
//public:
//    virtual void combine(char& old, const char& new_msg)
//    {
//    }
//};

void blogel_ctlabeling(string in_path, string out_path, int max_w, int n_threads)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    CTBlockWorker worker;
//    CTCombiner combiner;
//    worker.setCombiner(&combiner);
    worker.set_max_w(max_w);
    worker.set_n_threads(n_threads);
    worker.run(param);
}