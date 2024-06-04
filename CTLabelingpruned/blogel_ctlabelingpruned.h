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

struct BPLabel {
    uint8_t bpspt_d[N_ROOTS];//for better memory performance?
    uint64_t bpspt_s[N_ROOTS][2];
};

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

    //bp
    bool usd_bp;
    BPLabel label_bp;
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

struct rootMsg
{
    int order;
    int block;
    int deg;
    vector<int> nbr;
    vector<int> cutnbr;
    vector<int> cutblock;
};

ibinstream& operator<<(ibinstream& m, const rootMsg& v)
{
    m << v.order;
    m << v.block;
    m << v.deg;
    m << v.nbr;
    m << v.cutnbr;
    m << v.cutblock;
    return m;
}

obinstream& operator>>(obinstream& m, rootMsg& v)
{
    m >> v.order;
    m >> v.block;
    m >> v.deg;
    m >> v.nbr;
    m >> v.cutnbr;
    m >> v.cutblock;
    return m;
}



class CTVertex : public BVertex<VertexID, CTValue, rootMsg> {//<id, vertex value, message info>, no need of message
public:

    virtual void compute(MessageContainer& messages)
    {
        vote_to_halt();
    }
};

class CTBlock : public Block<CTValue, CTVertex, rootMsg> {//<vertex value, vertex(compute), msgtype>
public:

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes)
    {
        vote_to_halt();
    }

    // struct queue_info{
    //     int idx;
    //     int f_idx;
    //     int cost;
    //     uint8_t depth;

    //     queue_info(int _idx, int _f, int _cost, uint8_t _dep):idx(_idx), f_idx(_f), cost(_cost), depth(_dep){}
    // };


    // void construct_flags_E(int u, 
    //                        vector<pair<int, int>>& hops, 
    //                        vector<bool>& flags, 
    //                        deque<queue_info>& q, 
    //                        uniform_real_distribution<>& dis, 
    //                        mt19937& gen, 
    //                        double& prob)
    // {
    //     const int ori_depth = 1;
    //     int nb_greatest = 0;

    //     hops.resize(E[u].size());
    //     for (int i = 0; i < E[u].size(); ++i) {
    //         hops[i].first = E[u][i].first;
    //         hops[i].second = E[u][i].second;
    //         if(nb_greatest < hops[i].second)  nb_greatest = hops[i].second;
    //         q.push_back({hops[i].first, u, hops[i].second , ori_depth + 1});
    //     }

    //     int cur, f_idx, _cost;
    //     uint8_t depth;
    //     int itr, first, second;
    //     double randomValue;
    //     while(!q.empty()) {
    //         randomValue = dis(gen);
    //         if(randomValue > prob) {
    //             q.pop_front();
    //             continue;
    //         }

    //         cur = q.front().idx;
    //         f_idx = q.front().f_idx;
    //         depth = q.front().depth;
    //         _cost = q.front().cost;
    //         q.pop_front();

    //         //check current layer
    //         itr = 0;
    //         for (int i = 0; i < E[cur].size(); ++i) {
    //             if (E[cur][i].first != f_idx) {
    //                 first = E[cur][i].first;
    //                 second = _cost + E[cur][i].second;

    //                 while(itr < hops.size() && (hops[itr].first < first || !flags[itr])) itr++;
    //                 if(itr == hops.size()) break;
    //                 if(hops[itr].first == first && hops[itr].second >= second){//closer route
    //                     flags[itr] = false;//mark updated
    //                     if(hops[itr].second == nb_greatest){//update the longest neighbor
    //                         nb_greatest = -1;
    //                         for (int j = 0; j < hops.size(); ++j)
    //                             if(flags[itr] && nb_greatest < hops[j].second)
    //                                 nb_greatest = hops[j].second;
    //                     }
    //                 }
    //                 if(depth < PRUNEHOP && second < nb_greatest){ //push next hop nbrs
    //                     q.push_back({E[cur][i].first, cur, second , depth+1});
    //                 }
    //             }
    //         }
    //     }
    // }

    // long long all_min_LL(long long my_copy)
    // {
    //     long long tmp = 0;
    //     MPI_Allreduce(&my_copy, &tmp, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
    //     return tmp;
    // }

//     void remove_redundant_edges_E(VertexContainer& vertexes){//remove edges that are longer than other route
//         if(get_worker_id() == MASTER_RANK) cout << "removing redundant E edges..." << endl;
//         long long ori_core_edge = 0;
//         long long removed_cnt = 0;

//         //start
//         double construct_t = 0;
//         double cacu_t = 0;
//         double update_t = 0;
//         double t = omp_get_wtime();

//         //for debug usage
//         int max_hop_nb_cnt = 0;
//         int max_bucket_cnt = 0;
//         int max_circles = 0;
//         int max_cnt = 0;
//         int max_mem = 0;
//         long long one_edge = 0;
//         long long new_core_edge = 0;
//         long long cut_edge = 0;
//         //--------------

//         vector<vector<bool>> valid;
//         valid.resize(E.size());
//         for (size_t i = 0; i < E.size(); i++)
//         {
//             valid[i].resize(E[i].size());
//             fill(valid[i].begin(), valid[i].end(), true);
//         }
        

//         //debug---
//         long long op_cnt = 0;
//         for (int i = 0; i < E.size(); ++i) {
//             if(usd_bp[i] || vertexes[i]->value().rank >= 0) continue;
//             for(int j = 0; j < E[i].size(); ++j){
//                 op_cnt += E[E[i][j].first].size();
//             }
//         }

//         auto min_op = all_min_LL(op_cnt);
//         double sel_rate = min_op * 1.0 / op_cnt; //select rate for current worker, to ensure work balance
        
//         for(int i = 0; i < _num_workers; ++i){
//             if( i == _my_rank) 
//                 cout << _my_rank << "E size=" << E.size() 
//                     << ", 2 hop dectation need " << op_cnt 
//                     << " operations, select rate=" << sel_rate
//                     << " , about " << min_op / 100000000.0 << "seconds." <<endl;
//             MPI_Barrier(MPI_COMM_WORLD);
//             sleep(0.1);
//         }




// #pragma omp parallel
//         {
//             int pid = omp_get_thread_num(), np = omp_get_num_threads();
//             vector<pair<int, int>> hop_label;//idx, cost
//             // vector<pair<int, int>> bucket;//start, end
//             deque<queue_info> q;
//             double construct_t0, cacu_t0, update_t0;
//             double print_time = omp_get_wtime();

//             //minum operation seed
//             random_device rd;  // 用于获取随机数种子
//             mt19937 gen(rd()); // 使用种子初始化Mersenne Twister生成器
//             uniform_real_distribution<> dis(0.0, 1.0); // 定义分布范围[0.0, 1.0)

//             #pragma omp for schedule(dynamic)
//             for (int i = 0; i < E.size(); ++i) {//each vertex
//                 if (vertexes[i]->value().rank >= 0 || usd_bp[i]) {
//                     continue;
//                 }
//                 //construct N-hop list for all neighbors of x
//                 construct_t0 = omp_get_wtime();
//                 hop_label.clear();
//                 // bucket.clear();
//                 q.clear();
                
//                 if (pid == 0 && max_mem < hop_label.size()) max_mem = hop_label.size();
//                 if (pid == 0) construct_t += omp_get_wtime() - construct_t0;

//                 cacu_t0 = omp_get_wtime();
//                 construct_flags_E(i, hop_label, valid[i], q, dis, gen, sel_rate);

//                 //---debug---
//                 if (max_hop_nb_cnt < hop_label.size()) {
//                     max_hop_nb_cnt = hop_label.size();
//                     // max_bucket_cnt = bucket.size();
//                 }
//                 // if (max_cnt < cacu_cnt) max_cnt = cacu_cnt;
//                 //-----------

//                 if(pid == 0) cacu_t += omp_get_wtime() - cacu_t0;

//                 if (omp_get_wtime() - print_time > 60 && get_worker_id() == MASTER_RANK) {
//                     print_time = omp_get_wtime();
//                     cout << "tid " << pid << ": " << i << " E updated, "
//                          << omp_get_wtime() - t << "s consumed, "
//                          << " construct " << construct_t << "s, "
//                          << " cacu " << cacu_t << "s, "
//                          << " update " << update_t << "s, "
//                          << " max memory cost " << max_mem * sizeof(pair<int, int>) / (1024 * 1024) << "MB. "
//                          << " max nb cnt=" << max_hop_nb_cnt << " max buckets=" << max_bucket_cnt
//                          << " max caculation cnt=" << max_cnt << endl;
//                     max_hop_nb_cnt = 0;
//                     max_bucket_cnt = 0;
//                 }
//             }
//         }
//         cout << _my_rank << " multi-thread E reduction finished, cost " << omp_get_wtime() - t << "sec." << endl ;
        
//         vector<pair<int, int>> new_E;
//         for(int i = 0; i < E.size(); ++i){
//             if (vertexes[i]->value().rank >= 0 || usd_bp[i]){
//                 E[i].clear();//update
//                 continue;
//             }
//             ori_core_edge += E[i].size() + vertexes[i]->value().cutnbrs.size();

//             new_E.clear();
//             for(int j = 0; j < valid[i].size(); ++j){
//                 if (valid[i][j]) new_E.push_back(E[i][j]);
//                 if (E[i][j].second == 1) one_edge++;
//             }

//             new_core_edge += new_E.size() + vertexes[i]->value().cutnbrs.size();
//             cut_edge += vertexes[i]->value().cutnbrs.size();
//             removed_cnt += E[i].size() - new_E.size();
//             E[i] = new_E;//update
//         }

//         removed_cnt = all_sum_LL(removed_cnt);
//         ori_core_edge = all_sum_LL(ori_core_edge);
//         one_edge = all_sum_LL(one_edge);
//         new_core_edge = all_sum_LL(new_core_edge);
//         cut_edge = all_sum_LL(cut_edge);
//         if(get_worker_id() == MASTER_RANK) {
//             cout << "----------------------------------------------------"<< endl;
//             cout << "Removing redundant edge finished, t=" << omp_get_wtime() - t << "secs." << endl;
// //            cout << ", rest " << (ori_core_edge - removed_cnt) * 100.0 / ori_core_edge << "%" << endl;
//             cout << "original core edge cnt " << ori_core_edge << " percentage " << ori_core_edge * 100.0 / global_m << "%" << endl;
//             cout << "new core edge cnt " << new_core_edge << " percentage " << new_core_edge * 100.0 / global_m << "%" << endl;
//             cout << "Removed edge cnt " << removed_cnt << " removed edge percentage " <<  removed_cnt * 100.0 / global_m << "%" << endl;
//             cout << "original edge cnt = " << global_m << " reduced edge cnt " << ori_core_edge << endl;
//             cout << "one distance edge cnt = " << one_edge << " percentage " << one_edge * 100.0 / global_m << "%" << endl;
//             cout << "cut edge cnt = " << cut_edge << " percentage " << cut_edge * 100.0 / global_m << "%" << endl;
//             cout << "----------------------------------------------------"<< endl;
//         }
//     }

    // struct BPLabel {
    //     uint8_t bpspt_d[N_ROOTS];//for better memory performance?
    //     uint64_t bpspt_s[N_ROOTS][2];
    // };
    // BPLabel *label_bp;
    // bool *usd_bp;

    // void create_bp(){
    //     label_bp = new BPLabel[n];
    //     usd_bp = new bool[n];
    //     memset( usd_bp, 0, sizeof(bool) * n );
    // }

    // struct cutv{
    //     int v;
    //     int block;
    // };

    // void compute_bp_label(VertexContainer& vertexes, vector<rootMsg> top_n){
    //     if(_my_rank == MASTER_RANK) ( "Constructing BP Label...\n" );
    //     double t = omp_get_wtime();
    //     vector<int> v_vs[N_ROOTS];
    //     vector<int> v_vs_block[N_ROOTS];
    //     vector<bool> usd_top(top_n.size(), false);
    //     int top_ind = 0;
    //     for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt){
    //         while(usd_top[top_ind]) ++top_ind;
    //         if (top_ind == top_n.size()) {
    //             for (int v = 0; v < n; ++v) label_bp[v].bpspt_d[i_bpspt] = MAXD;
    //             continue;
    //         }
    //         usd_top[top_ind] = true;
    //         if(top_n[i_bpspt].block == _my_rank) usd_bp[top_n[i_bpspt].order] = true;
    //         v_vs[i_bpspt].push_back(top_n[i_bpspt].order);
    //         v_vs_block[i_bpspt].push_back(top_n[i_bpspt].block);
            
    //         int ns = 0;
    //         if(top_n[i_bpspt].block == _my_rank){
    //             for (int i = 0; i < top_n[top_ind].nbr.size(); ++i) {
    //                 int v = top_n[top_ind].nbr[i];
    //                 if (!usd_bp[v]) {
    //                     usd_bp[v] = true;
    //                     v_vs[i_bpspt].push_back(v);
    //                     v_vs_block[i_bpspt].push_back(_my_rank);
    //                     if (++ns == 64) break;
    //                 }
    //             }
    //         }
    //         else{
    //             for (int i = 0; i < top_n[top_ind].cutnbr.size(); ++i) {
    //                 int v = top_n[top_ind].cutnbr[i];
    //                 int b = top_n[top_ind].cutblock[i];
    //                 if (b == _my_rank && !usd_bp[v]) {
    //                     usd_bp[v] = true;
    //                     v_vs[i_bpspt].push_back(v);
    //                     v_vs_block[i_bpspt].push_back(b);
    //                     if (++ns == 64) break;
    //                 }
    //             }
    //         }
    //     }

    //     omp_set_num_threads(min(min(n_threads, N_ROOTS),MAX_BP_THREADS));
    // #pragma omp parallel
    //     {
    //         int pid = omp_get_thread_num(), np = omp_get_num_threads();
    //         if(_my_rank == MASTER_RANK && pid == 0 ) printf( "n_threads_bp = %d\n", np );
    //         vector<uint8_t> tmp_d(n);
    //         vector<pair<uint64_t, uint64_t> > tmp_s(n);
    //         vector<int> que(n);
    //         vector<pair<int, int> > child_es(m/2);
    //         int r, rb;

    // #pragma omp for schedule(dynamic)
    //         for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt) {
    //             if(_my_rank == MASTER_RANK) printf( "[%d]", i_bpspt );

    //             if( v_vs[i_bpspt].size() == 0 ) continue;
    //             fill(tmp_d.begin(), tmp_d.end(), MAXD);
    //             fill(tmp_s.begin(), tmp_s.end(), make_pair(0, 0));

    //             r = v_vs[i_bpspt][0];
    //             rb = v_vs_block[i_bpspt][0];
    //             int que_t0 = 0, que_t1 = 0, que_h = 0;
    //             que[que_h++] = r;
    //             tmp_d[r] = 0;
    //             que_t1 = que_h;

    //             for( size_t i = 1; i < v_vs[i_bpspt].size(); ++i) {
    //                 int v = v_vs[i_bpspt][i];
    //                 que[que_h++] = v;
    //                 tmp_d[v] = 1;
    //                 tmp_s[v].first = 1ULL << (i-1);
    //             }

    //             for (int d = 0; que_t0 < que_h; ++d) {
    //                 int num_child_es = 0;

    //                 for (int que_i = que_t0; que_i < que_t1; ++que_i) {
    //                     int v = que[que_i];

    //                     for (int i = 0; i < vertexes[v]->value().nbrs.size(); ++i) {
    //                         int tv = vertexes[v]->value().nbrs[i];
    //                         int td = d + 1;

    //                         if (d == tmp_d[tv]) {
    //                             if (v < tv) {
    //                                 tmp_s[v].second |= tmp_s[tv].first;
    //                                 tmp_s[tv].second |= tmp_s[v].first;
    //                             }
    //                         } else if( d < tmp_d[tv]) {
    //                             if (tmp_d[tv] == MAXD) {
    //                                 que[que_h++] = tv;
    //                                 tmp_d[tv] = td;
    //                             }
    //                             child_es[num_child_es].first  = v;
    //                             child_es[num_child_es].second = tv;
    //                             ++num_child_es;
    //                         }
    //                     }
    //                 }

    //                 for (int i = 0; i < num_child_es; ++i) {
    //                     int v = child_es[i].first, c = child_es[i].second;
    //                     tmp_s[c].first  |= tmp_s[v].first;
    //                     tmp_s[c].second |= tmp_s[v].second;
    //                 }

    //                 que_t0 = que_t1;
    //                 que_t1 = que_h;
    //             }

    //             for (int v = 0; v < n; ++v) {
    //                 label_bp[v].bpspt_d[i_bpspt] = tmp_d[v];
    //                 label_bp[v].bpspt_s[i_bpspt][0] = tmp_s[v].first;
    //                 label_bp[v].bpspt_s[i_bpspt][1] = tmp_s[v].second & ~tmp_s[v].first;
    //             }
    //         }
    //     }
    //     printf( "\n%d BP Label Constructed, bp_size=%0.3lf MB, t = %0.3lf secs\n", _my_rank, sizeof(BPLabel)*n/(1024.0*1024.0), omp_get_wtime() - t );
    // }
    
    // string tmp_path;
    // string bp_path;

    // void save_bp_label(){
    //     char tmp[5];
    //     sprintf(tmp, "%d", _my_rank);
    //     string p = bp_path + "/part_" + tmp;

    //     BufferedWriter* writer = new BufferedWriter(p.c_str());

    //     writer->check();

    //     sprintf(buf, "%d\n", E.size());//size
    //     writer->write(buf);
    //     for (int i = 0; i < n; i++)
    //     {
    //         sprintf(buf, " ", E.size());//size
    //         writer->write(buf);
    //     }
    //     for (int i = 0; i < E.size(); i++)
    //     {
    //         writer->check();
    //         sprintf(buf, "%d ", E[i].size());//E[i] size
    //         writer->write(buf);
    //         for (int j = 0; j < E[i].size(); ++j) {
    //             sprintf(buf, "%d %d ", E[i][j].first, E[i][j].second);//E[i][j]
    //             writer->write(buf);
    //         }
    //         writer->write("\n");
    //     }
    //     delete writer;
    // }

    // void delete_bp(){
    //     if(label_bp) delete[] label_bp; 
    //     if(usd_bp) delete[] usd_bp;
    // }






};

class CTBlockWorker : public BWorker<CTBlock> {
    char buf[1000];
    typedef char dint;
    int max_w, n_threads;
    string tmp_path;
    string bp_path;


    int global_n, n, n_core;
    long long global_m, m, m_core;
    vector<int> score, _score, ord;
    vector<vector<int>> nbr, cost;
    vector<vector<pair<int, int>>> E;//<order, dis>

public:

    
    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)//
    {
        // for(int i = 0; i < blocks.size(); ++i){
        //     blocks[i]->max_w = max_w;
        //     blocks[i]->n_threads = n_threads;
        //     blocks[i]->tmp_path = tmp_path;
        // }
        if(_my_rank == MASTER_RANK) cout << " loading bp " << endl;
        load_bp(bp_path);
        if(_my_rank == MASTER_RANK) cout << " load bp suc " << endl;

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
        }

        reduce(vertexes, max_w, n_threads);
        
        save_E(tmp_path);//save E
        create_tree(vertexes);
        compute_tree_label(vertexes);
    }

    
    int query_by_bp( int u, int v) {
        BPLabel &idx_u = vertexes[u]->value().label_bp, &idx_v = vertexes[v]->value().label_bp;
        int d = MAXD;//120
        for (int i = 0; i < N_ROOTS; ++i) {//i < 4, what's n_roots?
            int td = idx_u.bpspt_d[i] + (int) idx_v.bpspt_d[i];
            if (td - 2 <= d)
                td += (idx_u.bpspt_s[i][0] & idx_v.bpspt_s[i][0]) ? -2 :
                        ((idx_u.bpspt_s[i][0] & idx_v.bpspt_s[i][1]) | (idx_u.bpspt_s[i][1] & idx_v.bpspt_s[i][0])) ? -1 : 0;
            if (td < d) d = td;
        }
        return d == MAXD ? INT_MAX : d;
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
        int bp_dis, tmp_dis;
        double bp_t = 0, bp_tstart = 0;
        long long removed_edge = 0;

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
        int cnt = 0;
        vector<pair<int,int>> tmp;
        //每次选deg最小的node,加入nbr,E删除与该node连接的边,但是添加到自己的nbr中. rank记录标记成tree node的顺序

        //for record

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
                    E[u].push_back(make_pair(x,-cost[x][i]));// not active, is core?
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
                printf( "%d nodes reduced, score[x]=%d, remaining size=%0.3lf%% t=%0.3lf secs, bp query cost %0.3lf secs\n",
                        r, (n-r)*100.0/n, score[x], omp_get_wtime()-t, bp_t);
                cnt = 0;
            }
        }

        cout << _my_rank << " reduce finished, time=" << omp_get_wtime() - t << "s" << endl;

        if(get_worker_id() == MASTER_RANK) printf( "Reordering edges...\n" );

#pragma omp parallel reduction(+:removed_edge)
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
                            if(rankv>=0) continue;//remove tree
                            else if(query_by_bp(u, v) <= w){removed_edge++; continue;}
                            else if(buf[v] == -1) {buf[v] = w; ve.push_back(v);}
                            else if(w < buf[v]) buf[v] = w;
                        } else {
                            auto &s = E[e[i].first];
                            for(size_t j=0; j<s.size(); ++j) {//check if its neighbor have route to itself
                                int v = s[j].first, w = s[j].second - e[i].second;
                                int rankv = vertexes[v]->value().rank;
                                if(v == u || rankv>=0) continue;
                                else if(query_by_bp(u, v) <= w){removed_edge++; continue;}
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
        auto global_m_removed = all_sum_LL(removed_edge);

        for(int i = 0; i < _num_workers; ++i){
            if( i == _my_rank) 
                printf("%d Reduce finished, t=%0.3lf secs\nn_core=%d,m_core=%lld,node_rate=%0.3lf,edge_rate=%0.3lf\n",
                       _my_rank, omp_get_wtime() - t, n_core, m_core, n_core * 1.0 / n, m_core * 1.0 / m);
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(0.1);
        }
        if(_my_rank == MASTER_RANK){
            cout 
                << "global n=" << global_n << ", global m=" << global_m
                << ", global_n_core = " << global_n_core << ", global_m_core=" << global_m_core 
                << ", removed=" << global_m_removed 
                << ", node_rate=" << global_n_core * 1.0 / global_n 
                << ", edge_rate=" << global_m_core * 1.0 / global_m 
                << ", removed_rate=" << global_m_removed * 1.0 / global_m << endl;
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



        if(get_worker_id() == MASTER_RANK)
            printf( "Core tree constructed, maxh=%d, maxdep=%d, cnt_root=%d, max_stree=%d, avg_rsize=%0.3lf, t=%0.3lf secs\n",
                maxh, maxdep, cnt_root, max_sub_tree, tw/(n-n_core), omp_get_wtime()-t);
    }

    void set_max_w(int w){
        max_w = w;
    }

    void set_n_threads(int nt){
        n_threads = nt;
    }

    void set_tmp_path(string str){
        tmp_path = str;
    }
    void set_bp_path(string path){
        bp_path = path;
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

    void load_bp(string bppath)
    {
        bppath += "/part_" + to_string(_my_rank); 
        FILE* in = getRHandle(bppath.c_str());
        LineReader reader(in);
        int ind = 0;
        while (true)
        {
            reader.readLine();
            if (!reader.eof()){
                load_vertex_bp(ind, reader.getLine());
                ind++;
            }
            else
                break;
        }
        fclose(in);
    }

    void load_vertex_bp(int ind, char* line){
        //id bid order size csize nbblock nborder2 ... cutbid1 cutid1 cutbid2 cutid2...
        char* pch;
        pch = strtok(line, " ");//usdbp
        vertexes[ind]->value().usd_bp = atoi(pch);

        for (size_t i = 0; i < N_ROOTS; i++)
        {
            pch = strtok(NULL, " ");//d
            vertexes[ind]->value().label_bp.bpspt_d[i] = atoi(pch);
            pch = strtok(NULL, " ");//s[0]
            vertexes[ind]->value().label_bp.bpspt_s[i][0] = atoll(pch);
            pch = strtok(NULL, " ");//s[1]
            vertexes[ind]->value().label_bp.bpspt_s[i][1] = atoll(pch);
        }
        return;
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




};



void blogel_ctlabeling(string in_path, string bp_path, string out_path, string tmp_path, int max_w, int n_threads)
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
    worker.set_bp_path(bp_path);//bp label
    worker.set_tmp_path(tmp_path);//tmp graph path(E)
    worker.run(param);
}