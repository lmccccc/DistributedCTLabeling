#include "blogel_ctlabelingbp.h"
int main(int argc, char* argv[])
{
    init_workers();
    char *inpath = argv[1];
    char *outpath = argv[2];
    if(argc < 4 || inpath == NULL || outpath == NULL){
        cout << "invalid input path" << endl;
        exit(-1);
    }
    int threads = atoi(argv[3]);

    if(get_worker_id() == MASTER_RANK) {
        cout << "in dir=" << inpath << endl;
        cout << "out dir=" << outpath << endl;
        cout << "threads=" << threads << endl;
    }
    blogel_ctlabeling(string(inpath), string(outpath), threads);
    worker_finalize();
    return 0;
}
