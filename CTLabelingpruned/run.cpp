#include "blogel_ctlabelingpruned.h"
int main(int argc, char* argv[])
{
    init_workers();
    char *inpath = argv[1];
    char *outpath = argv[2];
    char *tmppath = argv[3];
    char *bppath = argv[4];
    if(argc < 7 || inpath == NULL || outpath == NULL){
        cout << "invalid input path" << endl;
        exit(-1);
    }
    int maxw = atoi(argv[5]);
    int threads = atoi(argv[6]);

    if(get_worker_id() == MASTER_RANK) {
        cout << "in dir=" << inpath << endl;
        cout << "out dir=" << outpath << endl;
        cout << "core graph dir=" << tmppath << endl;
        cout << "bp dir=" << bppath << endl;
        cout << "max width=" << maxw << endl;
        cout << "threads=" << threads << endl;
    }
    blogel_ctlabeling(string(inpath), string(outpath), string(tmppath), string(bppath), maxw, threads);
    worker_finalize();
    return 0;
}
