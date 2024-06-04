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
#include <queue>

using namespace std;
#define MAXD 120
#define PRUNEHOP 2
#define N_ROOTS 10
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
    BPLabel label_bp;
    bool usd_bp;

    vector<int> nbrs;
    vector<int> cutnbrs;
    vector<int> cutblock;
};

ibinstream& operator<<(ibinstream& mstr, const CTValue& v)
{
    mstr << v.order;
    mstr << v.rank;
    mstr << v.usd_bp;
    mstr << v.nbrs;
    mstr << v.cutnbrs;
    mstr << v.cutblock;
    return mstr;
}

obinstream& operator>>(obinstream& mstr, CTValue& v)
{
    mstr >> v.order;
    mstr >> v.rank;
    mstr >> v.usd_bp;
    mstr >> v.nbrs;
    mstr >> v.cutnbrs;
    mstr >> v.cutblock;
    return mstr;
}


struct bpMsg
{
    uint8_t i;
    uint8_t td;
    int v;
    int tv;
    pair<uint64_t, uint64_t> tmp_s_v;
};

ibinstream& operator<<(ibinstream& m, const bpMsg& v)
{
    m << v.i;
    m << v.td;
    m << v.v;
    m << v.tv;
    m << v.tmp_s_v.first;
    m << v.tmp_s_v.second;
    return m;
}

obinstream& operator>>(obinstream& m, bpMsg& v)
{
    m >> v.i;
    m >> v.td;
    m >> v.v;
    m >> v.tv;
    m >> v.tmp_s_v.first;
    m >> v.tmp_s_v.second;
    return m;
}


struct maxmsg{
    CTValue val;
    int wid;
};

ibinstream& operator<<(ibinstream& mstr, const maxmsg& v)
{
    mstr << v.val;
    mstr << v.wid;
    return mstr;
}

obinstream& operator>>(obinstream& mstr, maxmsg& v)
{
    mstr >> v.val;
    mstr >> v.wid;
    return mstr;
}



class CTVertex : public BVertex<VertexID, CTValue, bpMsg> {//<id, vertex value, message info>, no need of message
public:

    virtual void compute(MessageContainer& messages)
    {
        vote_to_halt();
    }
};

class CTBlock : public Block<CTValue, CTVertex, bpMsg> {//<vertex value, vertex(compute), msgtype>
public:

    vector<int> v_vs[N_ROOTS];
    vector<int> v_vs_block[N_ROOTS];

    int n_threads;
    int global_n, n, n_core;
    long long global_m, m, m_core;
    vector<uint8_t> tmp_d[N_ROOTS];
    vector<pair<uint64_t, uint64_t> > tmp_s[N_ROOTS];

    virtual void compute(MessageContainer& messages, VertexContainer& vertexes)
    {
        if(_my_rank == MASTER_RANK) cout << _my_rank << " block circle " << step_num() << " msg cnt=" << messages.size() << endl;
        vector<bpMsg> remote_msg[N_ROOTS];
        vector<bpMsg> msg_order[N_ROOTS];
        int msgcnt = 0;
        if(step_num() > 1){
            if(messages.size() == 0){ // no new msg, halt to stop
                vote_to_halt();
                return;
            }

            for (size_t i = 0; i < messages.size(); i++) // push msg to different n_roots
            {
                remote_msg[messages[i].i].push_back(messages[i]);
            }
            
            for (uint8_t i = 0; i < N_ROOTS; i++) // sort msg
            {
                sort(remote_msg[i].begin(), remote_msg[i].end(), [](const bpMsg& a, const bpMsg& b){
                    if(a.tv < b.tv) return true;//sort by id, and then td
                    else if(a.tv == b.tv){
                        return a.td < b.td;
                    }
                    return false;
                });

                //combine
                int ptr = 0;
                for(int j = 1; j < remote_msg[i].size(); ++j){
                    if(remote_msg[i][ptr].tv == remote_msg[i][j].tv){//same v
                        if(remote_msg[i][ptr].td == remote_msg[i][j].td){//same d
                            remote_msg[i][ptr].tmp_s_v.first |= remote_msg[i][j].tmp_s_v.first;
                            remote_msg[i][ptr].tmp_s_v.second |= remote_msg[i][j].tmp_s_v.second;
                        }
                        else if(remote_msg[i][ptr].td < remote_msg[i][j].td){//closer d
                            remote_msg[i][ptr] =  remote_msg[i][j];
                        }
                    }
                    else{ // different v
                        ptr++;
                    }
                }
                remote_msg[i].resize(ptr);

                sort(remote_msg[i].begin(), remote_msg[i].end(), [](const bpMsg& a, const bpMsg& b){//reorder by td
                    return a.td < b.td;
                });
            }

            for(int i = 0; i < _num_workers; ++i){
                if( i == _my_rank){
                    cout << _my_rank << " msg cnt = " ;
                    for (size_t j = 0; j < N_ROOTS; j++)
                    {
                        cout << remote_msg[j].size() << " ";
                    }
                    cout << endl;
                }
                // MPI_Barrier(MPI_COMM_WORLD);
                sleep(0.1);
            }
        }
        

    
        vector<pair<int, int> > child_es(m/2);
        vector<int> que(n);
        double t = omp_get_wtime();
        double local_t = 0, remote_t = 0, msg_t = 0, child_t = 0, local_start, remote_start, msg_start, child_start;
        int closer = 0, newtmp = 0;
        for (uint8_t i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt) { // for each root
            if(_my_rank == MASTER_RANK) printf( "[%d]", i_bpspt );
            int r, rb;
            int que_t0 = 0, que_t1 = 0, que_h = 0;
            int msgItr = 0;
            bool has_root = false;
            vector<bool> visited(n);
            fill(visited.begin(), visited.end(), false);
            if(step_num() == 1){//first time, no remote msg, start from v_vs, push them to queue
                if( v_vs[i_bpspt].size() == 0 ) continue;

                r = v_vs[i_bpspt][0];
                rb = v_vs_block[i_bpspt][0];

                if(rb == _my_rank){
                    que[que_h++] = r;
                    tmp_d[i_bpspt][r] = 0;
                    has_root = true;
                    visited[r] = true;
                }
                
                for( size_t i = 1; i < v_vs[i_bpspt].size(); ++i) {
                    int v = v_vs[i_bpspt][i];
                    int vb = v_vs_block[i_bpspt][i];
                    if(vb == _my_rank){
                        que[que_h++] = v;
                        tmp_d[i_bpspt][v] = 1;
                        tmp_s[i_bpspt][v].first = 1ULL << (i-1);
                        visited[v] = true;
                    }
                }
                if(has_root) que_t1 = 1;
                else que_t1 = que_h;
            }

            //first time: update locally, easy and fast
            //other circles: each remote info contains different d, which means they update more complicated
            // if(_my_rank == MASTER_RANK) cout << _my_rank << " initially quet0=" << que_t0 << " que_h=" << que_h << " time cost " << omp_get_wtime()-t << endl;
            int d;
            t = omp_get_wtime();
            if(step_num() == 1) d = has_root ? 0 : 1;
            else d = 0;
            for (; que_t0 < que_h || msgItr < remote_msg[i_bpspt].size(); ++d) {//increasing d
                int num_child_es = 0;
                // if(_my_rank == MASTER_RANK) cout << _my_rank << " d=" << d << " que_t0=" << que_t0 << " que_t1=" << que_t1 << " msgcnt=" << msgcnt << " time cost " << omp_get_wtime()-t << endl;
                for (int que_i = que_t0; que_i < que_t1; ++que_i) {// iterate queue
                    int v = que[que_i];

                    uint8_t td = d + 1;
                    local_start = omp_get_wtime();
                    for (int i = 0; i < vertexes[v]->value().nbrs.size(); ++i) {//each nbr of vertexes
                        int tv = vertexes[v]->value().nbrs[i];
                        //1. tmpd < d: do nothing
                        //2. tmpd == d: switch tmps.second
                        //3. tmpd == td: add tmps to tv
                        //4. tmpd > td: recover tmpd and tmps
                        if (d == tmp_d[i_bpspt][tv]){// same level instead of child
                            if (v < tv) {//update each other's route
                                tmp_s[i_bpspt][v].second |= tmp_s[i_bpspt][tv].first;
                                tmp_s[i_bpspt][tv].second |= tmp_s[i_bpspt][v].first;
                            }
                        } else if(td < tmp_d[i_bpspt][tv]){
                            if (!visited[tv]) {
                                que[que_h++] = tv;
                                visited[tv] = true;
                            }
                            
                            tmp_d[i_bpspt][tv] = td;
                            tmp_s[i_bpspt][tv].first  = tmp_s[i_bpspt][v].first;
                            tmp_s[i_bpspt][tv].second = tmp_s[i_bpspt][v].second;
                        } else if(td == tmp_d[i_bpspt][tv]){//td(d+1) == tmp_d or tmp_d == MAX_D
                            //original: only never used tv is added into queue
                            //new: any closer tv is added. note that we must ignore cross neighbor

							//new vertex
                            if (!visited[tv]) {
                                que[que_h++] = tv;
                                tmp_d[i_bpspt][tv] = td;
                                visited[tv] = true;
                            }

                            child_es[num_child_es].first  = v; // add to child
                            child_es[num_child_es].second = tv;
                            ++num_child_es;
                        }
                    }
                    local_t += omp_get_wtime() - local_start;

                    // for remote nbrs
                    remote_start = omp_get_wtime();
                    for (int i = 0; i < vertexes[v]->value().cutnbrs.size(); ++i) {//each remote nbr of vertexes
                        int tv = vertexes[v]->value().cutnbrs[i];
                        int tb = vertexes[v]->value().cutblock[i];
                        send_message(tb, tb, {i_bpspt, td, v, tv, tmp_s[i_bpspt][v]});
                        msgcnt++;
                    }
                    remote_t += omp_get_wtime() - remote_start;
                }
                
                child_start = omp_get_wtime();
                for (int i = 0; i < num_child_es; ++i) {//update tmps
                    int v = child_es[i].first, c = child_es[i].second;
                    tmp_s[i_bpspt][c].first  |= tmp_s[i_bpspt][v].first;
                    tmp_s[i_bpspt][c].second |= tmp_s[i_bpspt][v].second;
                }
                child_t += omp_get_wtime() - child_start;

                msg_start = omp_get_wtime();
                //check if there are remote msg that matches d
                num_child_es = 0;
                while(msgItr < remote_msg[i_bpspt].size() && remote_msg[i_bpspt][msgItr].td == d){
                    uint8_t td = remote_msg[i_bpspt][msgItr].td;
                    int tv = remote_msg[i_bpspt][msgItr].tv;
                
                    //1. td < tmpd: recover tv
                    //2. td == tmpd: check route, add if not exist
                    //3. td > tmpd: do nothing
                    if(td < tmp_d[i_bpspt][tv]){ // update dis and route
                        if(!visited[tv]){
                            que[que_h++] = tv;
                            visited[tv] = true;
                        }
                        
                        
                        tmp_d[i_bpspt][tv] = td;
                        tmp_s[i_bpspt][tv].first  = remote_msg[i_bpspt][msgItr].tmp_s_v.first;
                        tmp_s[i_bpspt][tv].second = remote_msg[i_bpspt][msgItr].tmp_s_v.second;
                        closer++;
                    }
                    else if(td == tmp_d[i_bpspt][tv]){ // join remote route
                        
                        if((~tmp_s[i_bpspt][tv].first & remote_msg[i_bpspt][msgItr].tmp_s_v.first) 
                            || (~tmp_s[i_bpspt][tv].second & remote_msg[i_bpspt][msgItr].tmp_s_v.second)){//local route can't cover remote route
                            
                            if(!visited[tv]){
                                que[que_h++] = tv;
                                visited[tv] = true;
                            }
                            tmp_s[i_bpspt][tv].first  |= remote_msg[i_bpspt][msgItr].tmp_s_v.first;
                            tmp_s[i_bpspt][tv].second |= remote_msg[i_bpspt][msgItr].tmp_s_v.second;
                            newtmp++;
                        }
                    }
                    msgItr++;
                }
                // if(_my_rank == MASTER_RANK) cout << " d=" << d <<  " get remote msg que=" << que_h << endl;
                msg_t = omp_get_wtime() - msg_start;
                        
                que_t0 = que_t1;
                que_t1 = que_h;
            }
        }
        
        if(_my_rank == MASTER_RANK) cout << _my_rank << " local proc cost " << local_t 
                << " remote proc cost " << remote_t 
                << " msg proc cost " << msg_t 
                << " child proc cost " << child_t << endl;
        
        // for(int i = 0; i < _num_workers; ++i){
        //     if( i == _my_rank){
        //         cout << _my_rank << " from msg, closer " << closer << " equal and new route " << newtmp << endl;
        //     }
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     sleep(0.1);
        // }
    }

    
    int query_by_bp(VertexContainer& vertexes, int u, int v) {
        BPLabel &idx_u = vertexes[u]->value().label_bp, &idx_v = vertexes[u]->value().label_bp;
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


    long long all_min_LL(long long my_copy)
    {
        long long tmp = 0;
        MPI_Allreduce(&my_copy, &tmp, 1, MPI_LONG_LONG_INT, MPI_MIN, MPI_COMM_WORLD);
        return tmp;
    }

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



    void set_n_threads(int nt){
        n_threads = nt;
    }


};

class CTBlockWorker : public BWorker<CTBlock> {
    char buf[1000];
    int n, global_n;
    long long m, global_m;
public:
    typedef char dint;
    int n_threads;

    
    virtual void blockInit(VertexContainer& vertexes, BlockContainer& blocks)// for bp init
    {
        blocks[0]->n_threads = n_threads;
        auto& v_vs = blocks[0]->v_vs;
        auto& v_vs_block = blocks[0]->v_vs_block;

        tolinecnt = 0;
        n = vertexes.size();
        blocks[0]->n = n;
        if(_my_rank == MASTER_RANK) cout << "n=" << n << endl;
        global_n = all_sum(n);

        m = 0;
        global_m = 0;
        long long int cut_m = 0;
        for (int i = 0; i < n; ++i) {
            m += vertexes[i]->value().nbrs.size();
            global_m += vertexes[i]->value().nbrs.size() + vertexes[i]->value().cutnbrs.size();
            cut_m += vertexes[i]->value().cutnbrs.size();
        }
        blocks[0]->m = m;
        global_m = all_sum_LL(global_m);

        for(int i = 0; i < _num_workers; ++i){
            if( i == _my_rank) cout << "rank id:" << _my_rank 
                                    << ", bock id:" << blocks[0]->bid 
                                    << ", local node=" << n 
                                    << ", local edge cnt=" << m  
                                    << ", cut edge cnt=" << cut_m << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            sleep(0.1);
        }
        
        //root vertexes choice shall be implemented before blogel computation
        
        int cut_pos = 0;
        for (size_t i = 0; i < vertexes.size(); i++) if(vertexes[i]->value().cutnbrs.size() == 0){ cut_pos=i; break;}
        int itr1 = 0, itr2 = cut_pos;
        for (int i_bpspt = 0; i_bpspt < N_ROOTS; ++i_bpspt){
            int deg1 = -1;
            if(itr1 < cut_pos) {
                deg1 = vertexes[itr1]->value().nbrs.size() + vertexes[itr1]->value().cutnbrs.size();
            }
            int deg2 = -1;
            if(itr2 < vertexes.size()) deg2 = vertexes[itr2]->value().nbrs.size();

            CTValue* max_v;
            if(deg1 == -1 && deg2 == -1){
                //do nothing, keep nbr and cutnbr list empty
            }
            else if(deg1 > deg2) max_v = &vertexes[itr1]->value();
            else max_v = &vertexes[itr2]->value();
            // get all vertex
            maxmsg global_max;
            if(_my_rank == MASTER_RANK){
                vector<maxmsg> gather_max;
                gather_max.resize(_num_workers);
                gather_max[MASTER_RANK] = {*max_v, _my_rank};
                masterGather(gather_max);
                
                int global_max_deg = -1;
                int global_max_src = -1;
                

                for (size_t i = 0; i < gather_max.size(); i++)//choose root
                {
                    int deg = gather_max[i].val.nbrs.size() + gather_max[i].val.cutnbrs.size();
                    if(global_max_deg < deg) {
                        global_max_deg = deg;
                        global_max_src = i;
                    }
                }
                int order = gather_max[global_max_src].val.order;

                if(global_max_deg == 0){//no valid root
                    for (int v = 0; v < n; ++v) vertexes[v]->value().label_bp.bpspt_d[i_bpspt] = MAXD;
                    continue;
                }
                cout << _my_rank << " initialise root, circle " << i_bpspt 
                    << " finally choose one from " << global_max_src 
                    << " with size " << global_max_deg << endl;
                global_max = gather_max[global_max_src];
                masterBcast(global_max);
            }
            else{
                maxmsg slamsg = {*max_v, _my_rank};
                slaveGather(slamsg);
                slaveBcast(global_max);
            }
            
            //v_vs init
            int ns = 0;
            vector<vector<int>> global_usd(_num_workers);
            v_vs[i_bpspt].push_back(global_max.val.order);
            v_vs_block[i_bpspt].push_back(global_max.wid);
            if(global_max.wid == _my_rank) {
                if(deg1 > deg2) itr1++;//update itr
                else itr2++;
                vertexes[global_max.val.order]->value().usd_bp = true;
            }
            global_usd[global_max.wid].push_back(global_max.val.order);


            for (int i = 0; i <  global_max.val.nbrs.size(); ++i) {
                int v = global_max.val.nbrs[i];
                bool usd;
                if(global_max.wid == _my_rank) usd = vertexes[v]->value().usd_bp;
                MPI_Bcast(&usd, 1, MPI_INT, global_max.wid, MPI_COMM_WORLD);//broadcast if usd
                if (!usd) {
                    if(global_max.wid == _my_rank) vertexes[v]->value().usd_bp = true;
                    v_vs[i_bpspt].push_back(v);
                    v_vs_block[i_bpspt].push_back(global_max.wid);
                    if (++ns == 64) break;
                }
            }
            for (int i = 0; i <  global_max.val.cutnbrs.size(); ++i) {
                int v = global_max.val.cutnbrs[i];
                int vb = global_max.val.cutblock[i];
                bool usd;
                if(vb == _my_rank) usd = vertexes[v]->value().usd_bp;
                MPI_Bcast(&usd, 1, MPI_INT, vb, MPI_COMM_WORLD);//broadcast if usd
                if (!usd) {//local neighbors of remote vertexes
                    if(vb == _my_rank) vertexes[v]->value().usd_bp = true;
                    v_vs[i_bpspt].push_back(v);
                    v_vs_block[i_bpspt].push_back(vb);
                    if (++ns == 64) break;
                }
            }
            
            blocks[0]->tmp_d[i_bpspt].resize(n);
            blocks[0]->tmp_s[i_bpspt].resize(n);
            fill(blocks[0]->tmp_d[i_bpspt].begin(), blocks[0]->tmp_d[i_bpspt].end(), MAXD);
            fill(blocks[0]->tmp_s[i_bpspt].begin(), blocks[0]->tmp_s[i_bpspt].end(), make_pair(0, 0));
        }
        
    }

    void set_n_threads(int nt){
        n_threads = nt;
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
    {//each node: usd d s s d s s d s s 
    
        sprintf(buf, "%d ", (int)v->value().usd_bp);
        writer.write(buf);
        for (size_t i_bpspt = 0; i_bpspt < N_ROOTS; i_bpspt++)
        {
            int order = v->value().order;
            v->value().label_bp.bpspt_d[i_bpspt] = b->tmp_d[i_bpspt][order];
            v->value().label_bp.bpspt_s[i_bpspt][0] = b->tmp_s[i_bpspt][order].first;
            v->value().label_bp.bpspt_s[i_bpspt][1] = b->tmp_s[i_bpspt][order].second & ~blocks[0]->tmp_s[i_bpspt][order].first;
            // if(_my_rank == 0 && v->value().order == 0 && i_bpspt == 9){
            //     cout << _my_rank << " d[9]=" << (int)v->value().label_bp.bpspt_d[i_bpspt] 
            //          << " s[9][0]=" << v->value().label_bp.bpspt_s[i_bpspt][0] 
            //          << " s[9][1]=" << v->value().label_bp.bpspt_s[i_bpspt][1] << endl;
            // }
            // sprintf(buf, "%d %llu %llu ",
            //         (int)v->value().label_bp.bpspt_d[i_bpspt], v->value().label_bp.bpspt_s[i_bpspt][0], v->value().label_bp.bpspt_s[i_bpspt][1]);
            sprintf(buf, "%d %d ",
                    v->id, (int)v->value().label_bp.bpspt_d[i_bpspt]);
            
            
            writer.write(buf);
        }
        writer.write("\n");
        
    }
};


void blogel_ctlabeling(string in_path, string out_path, int n_threads)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    CTBlockWorker worker;
    // CTCombiner combiner;
    // worker.setCombiner(&combiner);
    worker.set_n_threads(n_threads);
    worker.run(param);
}