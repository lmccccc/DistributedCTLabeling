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
#include <deque>
#include <mutex>
#include <shared_mutex>
#include <random>
#include <functional>

using namespace std;
#define MAXD 120
#define PRUNEHOP 2
#define N_ROOTS 0
#define MAX_BP_THREADS 8

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

ibinstream& operator<<(ibinstream& mstr, const CTValue& v)
{
    mstr << v.order;
    mstr << v.rank;
    mstr << v.f;
    mstr << v.h;
    mstr << v.rid;
    mstr << v.rsize;
    mstr << v.w;
    mstr << v.nbrs;
    mstr << v.cost;
    mstr << v.anc;
    mstr << v.dis;
    mstr << v.ch;
    mstr << v.cutnbrs;
    mstr << v.cutblock;
    return mstr;
}

obinstream& operator>>(obinstream& mstr, CTValue& v)
{
    mstr >> v.order;
    mstr >> v.rank;
    mstr >> v.rank;
    mstr >> v.f;
    mstr >> v.h;
    mstr >> v.rid;
    mstr >> v.rsize;
    mstr >> v.w;
    mstr >> v.nbrs;
    mstr >> v.cost;
    mstr >> v.anc;
    mstr >> v.dis;
    mstr >> v.ch;
    mstr >> v.cutnbrs;
    mstr >> v.cutblock;
    return mstr;
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
    int global_n, n, n_core;
    long long global_m, m, m_core;
    vector<int> score, _score, ord;
    vector<vector<int>> nbr, cost;
    vector<vector<pair<int, int>>> E;//<order, dis>

public:
    typedef char dint;
    int removed_edge;

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
    // pruning: remove edges that longer than 2-4 hops
    void reduce(VertexContainer& vertexes, int max_w, int n_threads) {// reduce: decomposition core step
        omp_set_num_threads(n_threads);
        double t = omp_get_wtime();
        this->max_w = max_w;

        score.resize(n);
        _score.resize(n);
        nbr.resize(n);
        cost.resize(n);


        vector<bool> changed(n,false);
        for(int i = 0; i < n; ++i)
            score[i] = _score[i] = vertexes[i]->value().nbrs.size();

        auto cmp = [&](const IntVal& left, const IntVal& right) {
            if(score[left.x] == score[right.x])
                return left.x < right.x;
            return score[left.x] < score[right.x];
        };

        set<IntVal, decltype(cmp)> q(cmp);// ordered quque?

//        nbr.resize(n); cost.resize(n);

        int r = 0;

        for(int i = 0; i < n; ++i) score[i] = _score[i] = vertexes[i]->value().nbrs.size();

//        printf( "%d, Initializing q...", get_worker_id() );
        if(get_worker_id() == MASTER_RANK) {
            cout << get_worker_id() << " Initializing q... " << endl;
        }

        vector<bool> active(n,false);
        for(int u = 0; u < n; ++u) {
            int deg = vertexes[u]->value().nbrs.size();
            if (deg < max_w && vertexes[u]->value().cutnbrs.size() == 0) {
                q.insert(IntVal(u));

                active[u] = true;
            } // nodes whose degree less than max_w are active
//            if(deg-1 != split) cutV.push_back(u);
        }
        if(get_worker_id() == MASTER_RANK) {
            cout << get_worker_id() << " t=" << omp_get_wtime() - t << " secs" << endl;
            cout << get_worker_id() << " Initializing E..." << endl;
        }
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
        if(get_worker_id() == MASTER_RANK) {
            cout << get_worker_id() << " t=" << omp_get_wtime() - t << " secs" << endl;
            cout << get_worker_id() << " Reducing Graph..." << endl;
        }
//             printf( ", t=%0.3lf secs\nReducing Graph...\n", omp_get_wtime()-t );
        int cnt = 0;
        vector<pair<int,int>> tmp;
        //每次选deg最小的node,加入nbr,E删除与该node连接的边,但是添加到自己的nbr中. rank记录标记成tree node的顺序

        removed_edge = 0;
        while(!q.empty()) {

            int x = q.begin()->x;// x is the core of this iteration, from q which is assigned as active
            if(vertexes[x]->value().cutnbrs.size() > 0){
                cout << "error, invalid node with cut edge, x=" << x<< endl;
                exit(-1);
            }
            while(changed[x]) { // get the first that not changed
                q.erase(x);
                score[x] = _score[x];
                q.insert(x);// move top to the tail if changed, and update score
                changed[x] = false;
                x = q.begin()->x;
            }

            //may never used
            if(score[x] >= max_w) break; // break if any of which scores higher than max, means CT decomposition is finished

            ord.push_back(x);
            q.erase(x);// pop x

            // x is no longer core node, order that assigned as tree
            vertexes[x]->value().rank = r;
            r++;

            //nbr empty at first
            //E[x] is neighbor and n-hop neighbor of x
            //we add x's neighbor in E to
            for(auto &it:E[x]){nbr[x].push_back(it.first); cost[x].push_back(it.second);}// collect active nodes

            for(auto &y:nbr[x]) { // for each neighbor of all nodes
                if(E[y].size() >= max_w * 2) {active[y] = false; q.erase(y);} // remove some nodes
                if(!active[y]) continue;
                for(size_t i=0;i<E[y].size();++i) {
                    if (E[y][i].first == x) {
                        E[y].erase(E[y].begin() + i);
                        break;
                    } // X is y's neighbor, then remove E[y][x]
                }
                _score[y] = (int) E[y].size();// update size
                changed[y] = true;// mark changed
            }

//            cout << "hop_label size=" << hop_label.size() << " ";
//            for (auto& itr: hop_label) {
//                cout << "(" << itr.first << "," << itr.second.size() << ") ";
//            }
//            cout << endl;

            //update E, remove tree neighbor, add correspond edges, update its score
            for(size_t i = 0; i < nbr[x].size(); ++i) {
                int u = nbr[x][i];

                if(!active[u]) {
                    E[u].push_back(make_pair(x,-cost[x][i]));// negative means not active, is core?
                    continue;
                }
                tmp.clear();
                size_t j=0, k=0;
                //------------------------------------------
                //optimization
                //check if new edge is the clothest amoung PRUNEDHOPS
                while(j<nbr[x].size()&&k<E[u].size()) {
                    if (j == i) ++j;//if u == nbr[x][j], ignore
                    else if (nbr[x][j] < E[u][k].first) {//x have neighbor [x][j] but u dont have
                        tmp.push_back(make_pair(nbr[x][j], cost[x][i] + cost[x][j]));
                        ++j;
                    }
                    else if (nbr[x][j] > E[u][k].first) {//x dont have but u have
                        tmp.push_back(E[u][k]);
                        ++k;
                    }
                    else {// both have, save the smallest
                        if (E[u][k].second < cost[x][i] + cost[x][j]) {
                            tmp.push_back(E[u][k]);
                        }
                        else {
                            tmp.push_back(make_pair(nbr[x][j], cost[x][i] + cost[x][j]));
                        }
                        ++j;
                        ++k;
                    }
                }
                for(;j<nbr[x].size();++j){
                    if(j!=i) {
                        tmp.push_back(make_pair(nbr[x][j], cost[x][i]+cost[x][j]));
                    }
                }
                for(;k<E[u].size();++k) {
                    tmp.push_back(E[u][k]);
                }
                E[u] = tmp;// update E
                if(_score[u] != (int) E[u].size()) {
                    changed[u] = true;
                    _score[u] = (int) E[u].size();
                }
            }

            if((++cnt) * score[x] > 1000000) {
                if(get_worker_id() == MASTER_RANK)
                printf( "%d nodes reduced, score[x]=%d, remaining size=%0.3lf%% t=%0.3lf secs\n",
                        r, (n-r)*100.0/n, score[x], omp_get_wtime()-t);
                cnt = 0;
            }
        }

        cout << _my_rank << " reduce finished, time=" << omp_get_wtime() - t << "s" << endl;

        if(get_worker_id() == MASTER_RANK) printf( "Reordering edges...\n" );

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
        for(int u = 0; u < n; ++u) {
            if (vertexes[u]->value().rank == -1) {//-1 means core node
                ++n_core;
                m_core += E[u].size() + vertexes[u]->value().cutnbrs.size();
            }
        }
        auto global_n_core = all_sum(n_core);
        auto global_m_core = all_sum_LL(m_core);

        for(int i = 0; i < _num_workers; ++i){
            if( i == _my_rank) 
                printf("%d Reducing finished, t=%0.3lf secs\nn_core=%d,m_core=%lld,node_rate=%0.3lf,edge_rate=%0.3lf\n",
                       _my_rank, omp_get_wtime() - t, n_core, m_core, n_core * 1.0 / n, m_core * 1.0 / m);
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(0.1);
        }
        if(_my_rank == MASTER_RANK){
            cout << "global_n_core = " << global_n_core << ", global_m_core=" << global_m_core 
                << ", node_rate=" << global_n_core * 1.0 / global_n 
                << ", edge_rate=" << global_m_core * 1.0 / global_m << endl;
        }
        

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void create_tree(VertexContainer& vertexes) {
        //nbr, vertexes.rank
        if(get_worker_id() == MASTER_RANK) printf( "Creating Tree...\n" );
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

        for(int i = (int) ord.size()-1; i >= 0; --i) {//
            int x = ord[i];//get tree node

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



        if(get_worker_id() == MASTER_RANK)
            printf( "Core tree constructed, maxh=%d, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%0.3lf, t=%0.3lf secs\n",
                maxh, maxdep, cnt_root, max_sub_tree, tw/(n-n_core), omp_get_wtime()-t);
    }

    struct queue_info{
        int idx;
        int f_idx;
        int cost;
        uint8_t depth;

        queue_info(int _idx, int _f, int _cost, uint8_t _dep):idx(_idx), f_idx(_f), cost(_cost), depth(_dep){}
    };

    void construct_hops(int u, vector<pair<int, int>>& hops, vector<pair<int, int>>& bucket, deque<queue_info>& q){

        q.push_back({u, -1, 0, 0});
        int nb_greatest = 0;

        while(!q.empty()) {
            int cur = q.front().idx;
            int f_idx = q.front().f_idx;
            int depth = q.front().depth;
            int _cost = q.front().cost;
            q.pop_front();

            if(depth > 0 && _cost >= nb_greatest) continue;

            int start = hops.size();
            hops.resize(start + nbr[cur].size());
            //push current layer
            for (int i = 0; i < nbr[cur].size(); ++i) {
                if (nbr[cur][i] != f_idx) {
                    hops[start+i].first = nbr[cur][i];
                    hops[start+i].second =  _cost + cost[cur][i];
                    if(depth == 0 && nb_greatest < hops[start+i].second)  nb_greatest = hops[start+i].second;
                    if(depth < PRUNEHOP) q.push_back({nbr[cur][i], cur, hops[start+i].second , depth+1});
                }
            }
            int end = hops.size();
            bucket.push_back(make_pair(start, end));//current layer's position in hop list
        }
    }


    void remove_redundant_edges_nbr(VertexContainer& vertexes){//remove edges that are longer than other route
        if(get_worker_id() == MASTER_RANK) cout << "removing redundant edges..." << endl;
        long long ori_core_edge = 0;
        long long removed_cnt = 0;

        //start
        double construct_t = 0;
        double cacu_t = 0;
        double update_t = 0;
        double t = omp_get_wtime();
        vector<pair<int, int>> hop_label;//idx, cost
        vector<pair<int, int>> bucket;//start, end
        deque<queue_info> q;

        //for debug usage
        int max_hop_nb_cnt = 0;
        int max_bucket_cnt = 0;
        int max_circles = 0;
        int max_cnt = 0;
        int max_mem = 0;
        long long one_edge = 0;
        long long edgec = 0;
        long long cutc = 0;
        for (auto i = 0; i < nbr.size(); ++i) {
            edgec +=  nbr[i].size();
            for (auto j = 0; j < nbr[i].size(); ++j) {
                if(cost[i][j] == 1) one_edge++;
            }
        }
        if(get_worker_id() == MASTER_RANK)
            cout << get_worker_id() << " total edge count=" << edgec << " one dis edge count=" << one_edge << " percentage=" << one_edge * 100.0 / edgec << "%" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        one_edge = 0;
        //--------------



//#pragma omp parallel
//        {
//            int pid = omp_get_thread_num(), np = omp_get_num_threads();
//
//
//#pragma omp for schedule(dynamic) reduction(max: max_hop_nb_cnt, max_bucket_cnt, max_circles, max_cnt) reduction(+: one_edge)
            for (int i = 0; i < nbr.size(); ++i) {//each vertex

                //construct N-hop list for all neighbors of x
                double construct_t0 = omp_get_wtime();
                hop_label.clear();
                bucket.clear();
                q.clear();
                construct_hops(i, hop_label, bucket, q);//source vertex, stored vector
                if (max_mem < hop_label.size()) max_mem = hop_label.size();
                construct_t += omp_get_wtime() - construct_t0;

                double cacu_t0 = omp_get_wtime();
                vector<bool> valid(nbr[i].size(), true);

                int cacu_cnt = 0;
                for (int nbstart = bucket[0].first; nbstart < bucket[0].second; ++nbstart) {
                    if (hop_label[nbstart].second == 1) {
                        one_edge++;
                        continue;
                    }
                    for (int j = 1; j < bucket.size(); ++j) {
                        int &start = bucket[j].first;
                        int &end = bucket[j].second;
                        while (start < end && hop_label[start].first < hop_label[nbstart].first) {
                            cacu_cnt++;
                            start++;
                        }
                        if (start < end && hop_label[start].first == hop_label[nbstart].first) {
                            //equals to the neighbor && closer
                            if (hop_label[start].second <= hop_label[nbstart].second) {
                                valid[nbstart] = false;
                                break;
                            }
                            start++;
                            cacu_cnt++;
                        }
                    }
                }
                //---debug---
                if (max_hop_nb_cnt < hop_label.size()) {
                    max_hop_nb_cnt = hop_label.size();
                    max_bucket_cnt = bucket.size();
                }
                if (max_cnt < cacu_cnt) max_cnt = cacu_cnt;
                //-----------

                cacu_t += omp_get_wtime() - cacu_t0;

                double update_t0 = omp_get_wtime();
                ori_core_edge += nbr[i].size();
                vector<int> new_nbr;
                vector<int> new_cost;
                for (int k = 0; k < valid.size(); ++k) {
                    if (valid[k]) {
                        new_nbr.push_back(nbr[i][k]);
                        new_cost.push_back(cost[i][k]);
                    } else removed_cnt++;
                }

                nbr[i] = new_nbr;
                cost[i] = new_cost;
                update_t += omp_get_wtime() - update_t0;


                if (i != 0 && i % 100000 == 0 && get_worker_id() == MASTER_RANK) {
                    cout << get_worker_id() << " " << i << " vertexes updated, "
                         << omp_get_wtime() - t << "s consumed, "
                         << removed_cnt << " edges removed, about "
                         << removed_cnt * 100.0 / ori_core_edge << "% "
                         << " construct " << construct_t << "s, "
                         << " cacu " << cacu_t << "s, "
                         << " update " << update_t << "s, "
                         << " max memory cost " << max_mem * sizeof(pair<int, int>) / (1024 * 1024) << "MB. " << endl;
                    cout << " max nb cnt=" << max_hop_nb_cnt << " max buckets=" << max_bucket_cnt
                         << " max caculation cnt=" << max_cnt << endl;
                    max_hop_nb_cnt = 0;
                    max_bucket_cnt = 0;
                }
            }
//        }
        removed_cnt = all_sum_LL(removed_cnt);
        ori_core_edge = all_sum_LL(ori_core_edge);
        one_edge = all_sum_LL(one_edge);
        if(get_worker_id() == MASTER_RANK) {
            cout << get_worker_id() << " Removing redundant edge finished, t=" << omp_get_wtime() - t
                << " secs. Removed edge cnt " << removed_cnt
                << ", rest " << (ori_core_edge - removed_cnt) * 100.0 / ori_core_edge << "%" << endl;
            cout << "rest edge pnt " <<  (global_m - removed_cnt) * 100.0 / global_m << endl;
            cout << "original edge cnt = " << global_m << " reduced edge cnt=" << ori_core_edge << endl;
            cout << "one distance edge cnt = " << one_edge << " percentage=" << one_edge * 100.0 / ori_core_edge << "%" << endl;
        }
    }

    void construct_hops_E(int u, vector<pair<int, int>>& hops, vector<pair<int, int>>& bucket, deque<queue_info>& q){
        int ori_depth = 1;
        q.push_back(queue_info(u, -1, 0, ori_depth));
        int nb_greatest = 0;

        while(!q.empty()) {
            int cur = q.front().idx;
            int f_idx = q.front().f_idx;
            int depth = q.front().depth;
            int _cost = q.front().cost;
            q.pop_front();

            if(depth > ori_depth && _cost >= nb_greatest) continue;

            int start = hops.size();
            hops.resize(start + E[cur].size());
            //push current layer
            for (int i = 0; i < E[cur].size(); ++i) {
                if (E[cur][i].first != f_idx) {
                    hops[start+i].first = E[cur][i].first;
                    hops[start+i].second =  _cost + E[cur][i].second;
                    if(depth == ori_depth && nb_greatest < hops[start+i].second)  nb_greatest = hops[start+i].second;
                    if(depth < PRUNEHOP) q.push_back({E[cur][i].first, cur, hops[start+i].second , depth+1});
                }
            }
            int end = hops.size();
            bucket.push_back(make_pair(start, end));//current layer's position in hop list
        }
    }

    void construct_flags_E(int u, 
                           vector<pair<int, int>>& hops, 
                           vector<bool>& flags, 
                           deque<queue_info>& q, 
                           uniform_real_distribution<>& dis, 
                           mt19937& gen, 
                           double& prob)
    {
        const int ori_depth = 1;
        int nb_greatest = 0;

        hops.resize(E[u].size());
        for (int i = 0; i < E[u].size(); ++i) {
            hops[i].first = E[u][i].first;
            hops[i].second = E[u][i].second;
            if(nb_greatest < hops[i].second)  nb_greatest = hops[i].second;
            q.push_back({hops[i].first, u, hops[i].second , ori_depth + 1});
        }

        int cur, f_idx, _cost;
        uint8_t depth;
        int itr, first, second;
        double randomValue;
        while(!q.empty()) {
            randomValue = dis(gen);
            if(randomValue > prob) {
                q.pop_front();
                continue;
            }

            cur = q.front().idx;
            f_idx = q.front().f_idx;
            depth = q.front().depth;
            _cost = q.front().cost;
            q.pop_front();

            //check current layer
            itr = 0;
            for (int i = 0; i < E[cur].size(); ++i) {
                if (E[cur][i].first != f_idx) {
                    first = E[cur][i].first;
                    second = _cost + E[cur][i].second;

                    while(itr < hops.size() && (hops[itr].first < first || !flags[itr])) itr++;
                    if(itr == hops.size()) break;
                    if(hops[itr].first == first && hops[itr].second >= second){//closer route
                        flags[itr] = false;//mark updated
                        if(hops[itr].second == nb_greatest){//update the longest neighbor
                            nb_greatest = -1;
                            for (int j = 0; j < hops.size(); ++j)
                                if(flags[itr] && nb_greatest < hops[j].second)
                                    nb_greatest = hops[j].second;
                        }
                    }
                    if(depth < PRUNEHOP && second < nb_greatest){ //push next hop nbrs
                        q.push_back({E[cur][i].first, cur, second , depth+1});
                    }
                }
            }
        }
    }

    long long all_min_LL(long long my_copy)
    {
        long long tmp = 0;
        MPI_Allreduce(&my_copy, &tmp, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
        return tmp;
    }

    void remove_redundant_edges_E(VertexContainer& vertexes){//remove edges that are longer than other route
        if(get_worker_id() == MASTER_RANK) cout << "removing redundant E edges..." << endl;
        long long ori_core_edge = 0;
        long long removed_cnt = 0;

        //start
        double construct_t = 0;
        double cacu_t = 0;
        double update_t = 0;
        double t = omp_get_wtime();

        //for debug usage
        int max_hop_nb_cnt = 0;
        int max_bucket_cnt = 0;
        int max_circles = 0;
        int max_cnt = 0;
        int max_mem = 0;
        long long one_edge = 0;
        long long new_core_edge = 0;
        long long cut_edge = 0;
        //--------------

        vector<vector<bool>> valid;
        valid.resize(E.size());
        for (size_t i = 0; i < E.size(); i++)
        {
            valid[i].resize(E[i].size());
            fill(valid[i].begin(), valid[i].end(), true);
        }
        

        //debug---
        long long op_cnt = 0;
        for (int i = 0; i < E.size(); ++i) {
            if(usd_bp[i] || vertexes[i]->value().rank >= 0) continue;
            for(int j = 0; j < E[i].size(); ++j){
                op_cnt += E[E[i][j].first].size();
            }
        }

        auto min_op = all_min_LL(op_cnt);
        double sel_rate = min_op * 1.0 / op_cnt; //select rate for current worker, to ensure work balance
        
        for(int i = 0; i < _num_workers; ++i){
            if( i == _my_rank) 
                cout << _my_rank << "E size=" << E.size() 
                    << ", 2 hop dectation need " << op_cnt 
                    << " operations, select rate=" << sel_rate
                    << " , about " << min_op / 100000000.0 << "seconds." <<endl;
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(0.1);
        }




#pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();
            vector<pair<int, int>> hop_label;//idx, cost
            // vector<pair<int, int>> bucket;//start, end
            deque<queue_info> q;
            double construct_t0, cacu_t0, update_t0;
            double print_time = omp_get_wtime();

            //minum operation seed
            random_device rd;  // 用于获取随机数种子
            mt19937 gen(rd()); // 使用种子初始化Mersenne Twister生成器
            uniform_real_distribution<> dis(0.0, 1.0); // 定义分布范围[0.0, 1.0)

            #pragma omp for schedule(dynamic)
            for (int i = 0; i < E.size(); ++i) {//each vertex
                if (vertexes[i]->value().rank >= 0 || usd_bp[i]) {
                    continue;
                }
                //construct N-hop list for all neighbors of x
                construct_t0 = omp_get_wtime();
                hop_label.clear();
                // bucket.clear();
                q.clear();
                
//                construct_hops_E(i, hop_label, bucket, q);//source vertex, stored vector
                if (pid == 0 && max_mem < hop_label.size()) max_mem = hop_label.size();
                if (pid == 0) construct_t += omp_get_wtime() - construct_t0;

                cacu_t0 = omp_get_wtime();
                construct_flags_E(i, hop_label, valid[i], q, dis, gen, sel_rate);

                //---debug---
                if (max_hop_nb_cnt < hop_label.size()) {
                    max_hop_nb_cnt = hop_label.size();
                    // max_bucket_cnt = bucket.size();
                }
                // if (max_cnt < cacu_cnt) max_cnt = cacu_cnt;
                //-----------

                if(pid == 0) cacu_t += omp_get_wtime() - cacu_t0;

                if (omp_get_wtime() - print_time > 60 && get_worker_id() == MASTER_RANK) {
                    print_time = omp_get_wtime();
                    cout << "tid " << pid << ": " << i << " E updated, "
                         << omp_get_wtime() - t << "s consumed, "
                         << " construct " << construct_t << "s, "
                         << " cacu " << cacu_t << "s, "
                         << " update " << update_t << "s, "
                         << " max memory cost " << max_mem * sizeof(pair<int, int>) / (1024 * 1024) << "MB. "
                         << " max nb cnt=" << max_hop_nb_cnt << " max buckets=" << max_bucket_cnt
                         << " max caculation cnt=" << max_cnt << endl;
                    max_hop_nb_cnt = 0;
                    max_bucket_cnt = 0;
                }
            }
        }
        cout << _my_rank << " multi-thread E reduction finished, cost " << omp_get_wtime() - t << "sec." << endl ;
        
        vector<pair<int, int>> new_E;
        for(int i = 0; i < E.size(); ++i){
            if (vertexes[i]->value().rank >= 0 || usd_bp[i]){
                E[i].clear();//update
                continue;
            }
            ori_core_edge += E[i].size() + vertexes[i]->value().cutnbrs.size();

            new_E.clear();
            for(int j = 0; j < valid[i].size(); ++j){
                if (valid[i][j]) new_E.push_back(E[i][j]);
                if (E[i][j].second == 1) one_edge++;
            }

            new_core_edge += new_E.size() + vertexes[i]->value().cutnbrs.size();
            cut_edge += vertexes[i]->value().cutnbrs.size();
            removed_cnt += E[i].size() - new_E.size();
            E[i] = new_E;//update
        }

        removed_cnt = all_sum_LL(removed_cnt);
        ori_core_edge = all_sum_LL(ori_core_edge);
        one_edge = all_sum_LL(one_edge);
        new_core_edge = all_sum_LL(new_core_edge);
        cut_edge = all_sum_LL(cut_edge);
        if(get_worker_id() == MASTER_RANK) {
            cout << "----------------------------------------------------"<< endl;
            cout << "Removing redundant edge finished, t=" << omp_get_wtime() - t << "secs." << endl;
//            cout << ", rest " << (ori_core_edge - removed_cnt) * 100.0 / ori_core_edge << "%" << endl;
            cout << "original core edge cnt " << ori_core_edge << " percentage " << ori_core_edge * 100.0 / global_m << "%" << endl;
            cout << "new core edge cnt " << new_core_edge << " percentage " << new_core_edge * 100.0 / global_m << "%" << endl;
            cout << "Removed edge cnt " << removed_cnt << " removed edge percentage " <<  removed_cnt * 100.0 / global_m << "%" << endl;
            cout << "original edge cnt = " << global_m << " reduced edge cnt " << ori_core_edge << endl;
            cout << "one distance edge cnt = " << one_edge << " percentage " << one_edge * 100.0 / global_m << "%" << endl;
            cout << "cut edge cnt = " << cut_edge << " percentage " << cut_edge * 100.0 / global_m << "%" << endl;
            cout << "----------------------------------------------------"<< endl;
        }
    }

    string tmp_path;
    void set_tmp_path(string str){
        tmp_path = str;
    }

    struct BPLabel {
        uint8_t bpspt_d[N_ROOTS];//for better memory performance?
        uint64_t bpspt_s[N_ROOTS][2];
    };
    BPLabel *label_bp;
    bool *usd_bp;

    void create_bp(){
        label_bp = new BPLabel[n];
        usd_bp = new bool[n];
        memset( usd_bp, 0, sizeof(bool) * n );
    }

    void compute_bp_label(){
        if(_my_rank == MASTER_RANK) ( "Constructing BP Label...\n" );
        double t = omp_get_wtime();
        vector<int> v_vs[N_ROOTS];

        int r = 0;
        for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt) {
            while (r < n && usd_bp[r]) ++r;
            if (r == n) {
                for (int v = 0; v < n; ++v) label_bp[v].bpspt_d[i_bpspt] = MAXD;
                continue;
            }
            usd_bp[r] = true;
            v_vs[i_bpspt].push_back(r);
            int ns = 0;
            for (int i = 0; i < vertexes[r]->value().nbrs.size(); ++i) {
                int v = vertexes[r]->value().nbrs[i];
                if (!usd_bp[v]) {
                    usd_bp[v] = true;
                    v_vs[i_bpspt].push_back(v);
                    if (++ns == 64) break;
                }
            }
        }

        omp_set_num_threads(min(min(n_threads, N_ROOTS),MAX_BP_THREADS));
    #pragma omp parallel
        {
            int pid = omp_get_thread_num(), np = omp_get_num_threads();
            if(_my_rank == MASTER_RANK && pid == 0 ) printf( "n_threads_bp = %d\n", np );
            vector<uint8_t> tmp_d(n);
            vector<pair<uint64_t, uint64_t> > tmp_s(n);
            vector<int> que(n);
            vector<pair<int, int> > child_es(m/2);

    #pragma omp for schedule(dynamic)
            for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt) {
                if(_my_rank == MASTER_RANK) printf( "[%d]", i_bpspt );

                if( v_vs[i_bpspt].size() == 0 ) continue;
                fill(tmp_d.begin(), tmp_d.end(), MAXD);
                fill(tmp_s.begin(), tmp_s.end(), make_pair(0, 0));

                r = v_vs[i_bpspt][0];
                int que_t0 = 0, que_t1 = 0, que_h = 0;
                que[que_h++] = r;
                tmp_d[r] = 0;
                que_t1 = que_h;

                for( size_t i = 1; i < v_vs[i_bpspt].size(); ++i) {
                    int v = v_vs[i_bpspt][i];
                    que[que_h++] = v;
                    tmp_d[v] = 1;
                    tmp_s[v].first = 1ULL << (i-1);
                }

                for (int d = 0; que_t0 < que_h; ++d) {
                    int num_child_es = 0;

                    for (int que_i = que_t0; que_i < que_t1; ++que_i) {
                        int v = que[que_i];

                        for (int i = 0; i < vertexes[v]->value().nbrs.size(); ++i) {
                            int tv = vertexes[v]->value().nbrs[i];
                            int td = d + 1;

                            if (d == tmp_d[tv]) {
                                if (v < tv) {
                                    tmp_s[v].second |= tmp_s[tv].first;
                                    tmp_s[tv].second |= tmp_s[v].first;
                                }
                            } else if( d < tmp_d[tv]) {
                                if (tmp_d[tv] == MAXD) {
                                    que[que_h++] = tv;
                                    tmp_d[tv] = td;
                                }
                                child_es[num_child_es].first  = v;
                                child_es[num_child_es].second = tv;
                                ++num_child_es;
                            }
                        }
                    }

                    for (int i = 0; i < num_child_es; ++i) {
                        int v = child_es[i].first, c = child_es[i].second;
                        tmp_s[c].first  |= tmp_s[v].first;
                        tmp_s[c].second |= tmp_s[v].second;
                    }

                    que_t0 = que_t1;
                    que_t1 = que_h;
                }

                for (int v = 0; v < n; ++v) {
                    label_bp[v].bpspt_d[i_bpspt] = tmp_d[v];
                    label_bp[v].bpspt_s[i_bpspt][0] = tmp_s[v].first;
                    label_bp[v].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
                }
            }
        }
        printf( "\n%d BP Label Constructed, bp_size=%0.3lf MB, t = %0.3lf secs\n", _my_rank, sizeof(BPLabel)*n/(1024.0*1024.0), omp_get_wtime() - t );
    }

    string bp_path;
    void set_bp_path(string path){
        bp_path = path;
    }
    void save_bp_label(){
        char tmp[5];
        sprintf(tmp, "%d", _my_rank);
        string p = bp_path + "/part_" + tmp;

        BufferedWriter* writer = new BufferedWriter(p.c_str());

        writer->check();

        sprintf(buf, "%d\n", E.size());//size
        writer->write(buf);
        for (int i = 0; i < n; i++)
        {
            sprintf(buf, " ", E.size());//size
            writer->write(buf);
        }
        for (int i = 0; i < E.size(); i++)
        {
            writer->check();
            sprintf(buf, "%d ", E[i].size());//E[i] size
            writer->write(buf);
            for (int j = 0; j < E[i].size(); ++j) {
                sprintf(buf, "%d %d ", E[i][j].first, E[i][j].second);//E[i][j]
                writer->write(buf);
            }
            writer->write("\n");
        }
        delete writer;
    }

    void delete_bp(){
        if(label_bp) delete[] label_bp; 
        if(usd_bp) delete[] usd_bp;
    }

    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)//
    {
        tolinecnt = 0;
        n = vertexes.size();
        if(_my_rank == MASTER_RANK) cout << "n=" << n << endl;
        global_n = all_sum(n);

        m = 0;
        global_m = 0;
        for (int i = 0; i < n; ++i) {
            m += vertexes[i]->value().nbrs.size();
            global_m += vertexes[i]->value().nbrs.size() + vertexes[i]->value().cutnbrs.size();
        }
        global_m = all_sum_LL(global_m);

        for(int i = 0; i < _num_workers; ++i){
            if( i == _my_rank) cout << _my_rank << " local node=" << n << ", local edge cnt=" << m << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(0.1);
        }
        
        if (_my_rank == MASTER_RANK) {
            check_tmp_graph(tmp_path);//check E dir
            check_tmp_graph(bp_path);//check bp dir
        }
        //key mem usage: vertexes.nbrs, nbr, cost, E

        reduce(vertexes, max_w, n_threads);
        


        create_bp();
        compute_bp_label();
        save_bp_label();
//        remove_redundant_edges_nbr(vertexes);
        remove_redundant_edges_E(vertexes);
        delete_bp();
        
        save_E(tmp_path);//save E
        create_tree(vertexes);
        compute_tree_label(vertexes);
    }



    void compute_tree_label(VertexContainer& vertexes, int x, int rsize, vector<CTVertex*> &s) {//
        s.push_back(vertexes[x]); // node itself
        auto &tn = vertexes[x]->value();
        tn.dis.resize(tn.h);// to each layer's distance?
        int pos = 0;
        vector<int> p(tn.w);
        for(int j = 0; j < tn.w; ++j) { //each neighbor of tree[x] is in s?
            //s中必然有一个是x的nbr? x是什么? 遍历到的tree node,
            //s是什么? 第一个是x, 后来是x的children
            while(s[pos]->value().order != tn.nbrs[j]) {//
                ++pos;//find neighbor position in s
//                cout << "nbr=" << tn.nbrs[j] << " s order=" << s[pos]->value().order << endl;
            }
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
        for(int &u:vertexes[x]->value().ch) {// for each children
            compute_tree_label(vertexes, u, rsize, s); // compute label. what is label? dis!
        }
        s.pop_back();
    }
    
    void compute_tree_label(VertexContainer& vertexes) {// tree label
        if(get_worker_id() == MASTER_RANK) printf( "Computing Tree Label...\n" );
        double t = omp_get_wtime();
        vector<CTVertex*> s;
        for(int v=0; v<n; ++v){//each v
            if(vertexes[v]->value().rank >= 0 && vertexes[v]->value().f == -1) {// not core and root node
                s.clear();
                for(int i = 0; i < vertexes[v]->value().w; ++i) {//push back each neighbor of v
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
            } else vector<pair<int,int>>(E[v]).swap(E[v]);//
        }


        printf( "%d Tree Label Computed, t=%0.3lf secs, maxdis=%d, tree label size=%0.3lf MB\n",
                get_worker_id(), omp_get_wtime()-t, maxdis, t_size/(1024.0*1024.0));
    }

    virtual void load_graph(const char* inpath)
    {
        FILE* in = getRHandle(inpath);
        LineReader reader(in);
        while (true)
        {
            reader.readLine();
            if (!reader.eof()){
                load_vertex(toVertex(reader.getLine()));
            }
            else
                break;
        }
        fclose(in);
    }

    virtual CTVertex* toVertex(char* line)
    {
        //id bid order size csize nbblock nborder2 ... cutbid1 cutid1 cutbid2 cutid2...
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
        pch = strtok(NULL, " ");//size
        int edge_size = atoi(pch);
        pch = strtok(NULL, " ");//split
        int split = atoi(pch);

        //order is empty
        v->value().nbrs.resize(split);
        for (int i = 0; i < split; ++i) {
            pch = strtok(NULL, " ");//bid

            pch = strtok(NULL, " ");//id
            v->value().nbrs[i] = atoi(pch);
        }

        int cut_size = edge_size-split;
        v->value().cutnbrs.resize(cut_size);
        v->value().cutblock.resize(cut_size);
        for (int i = 0; i < cut_size; ++i) {
            pch = strtok(NULL, " ");//bid/nid
            int block = atoi(pch);
            pch = strtok(NULL, " ");//id
            int order = atoi(pch);
            v->value().cutnbrs[i] = order;
            v->value().cutblock[i] = block;
        }



        return v;
    }

    int tolinecnt;

    virtual void toline(CTBlock* b, CTVertex* v, BufferedWriter& writer)
    {//each node: rank, rid, rsize, h, w, nbr*w, dis*h, anc*(h-w)

        sprintf(buf, "%d %d %d %d %d %d %d %d ",
                v->id, v->bid, v->value().order, v->value().rank, v->value().rid, v->value().rsize, v->value().h, v->value().w);
        writer.write(buf);
        if(v->value().rank != -1) {
            for (int i = 0; i < v->value().w; ++i) {
                sprintf(buf, "%d ", v->value().nbrs[i]);
                writer.write(buf);
            }

            for (int i = 0; i < v->value().h; ++i) {
                sprintf(buf, "%d ", v->value().dis[i]);
                writer.write(buf);
            }

            while (v->value().anc.size() < v->value().h - v->value().w) v->value().anc.push_back(-1);
            for (int i = 0; i < v->value().h - v->value().w; ++i) {
                sprintf(buf, "%d ", v->value().anc[i]);
                writer.write(buf);
            }

        }
        writer.write("\n");
    }

    void save_E(string tmp_outpath){
        char tmp[5];
        sprintf(tmp, "%d", _my_rank);
        tmp_outpath = tmp_outpath + "/part_" + tmp;
        BufferedWriter* writer = new BufferedWriter(tmp_outpath.c_str());

        writer->check();

        sprintf(buf, "%d\n", E.size());//size
        writer->write(buf);
        for (int i = 0; i < E.size(); i++)
        {
            writer->check();
            sprintf(buf, "%d ", E[i].size());//E[i] size
            writer->write(buf);
            for (int j = 0; j < E[i].size(); ++j) {
                sprintf(buf, "%d %d ", E[i][j].first, E[i][j].second);//E[i][j]
                writer->write(buf);
            }
            writer->write("\n");
        }
        delete writer;


    }

    bool check_tmp_graph(string p){
        if (access(p.c_str(), F_OK) == 0) {
            if (rm_dir(p.c_str()) == -1) {
                fprintf(stderr, "Error deleting %s!\n", p.c_str());
                exit(-1);
            }
            int created = mkdir(p.c_str(), MODE);
            if (created == -1) {
                fprintf(stderr, "Failed to create folder %s! after removing \n", p.c_str());
                exit(-1);
            }
        } else {
            int created = mkdir(p.c_str(), MODE);
            if (created == -1) {
                fprintf(stderr, "Failed to create folder %s!\n", p.c_str());
                exit(-1);
            }
        }
        return true;
    }
};



void blogel_ctlabeling(string in_path, string out_path, string tmp_path, string bp_path, int max_w, int n_threads)
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
    worker.set_tmp_path(tmp_path);//tmp graph path(E)
    worker.set_bp_path(bp_path);//bp label
    worker.run(param);
}