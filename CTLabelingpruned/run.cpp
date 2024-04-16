#include "blogel_ctlabelingpruned.h"
int main(int argc, char* argv[])
{
    init_workers();
    char *inpath = argv[1];
    char *outpath = argv[2];
    if(argc < 3 || inpath == NULL || outpath == NULL){
        cout << "invalid input path" << endl;
        exit(-1);
    }

    if(get_worker_id() == MASTER_RANK) {
        cout << "in dir=" << inpath << endl;
        cout << "out dir=" << outpath << endl;
    }
    blogel_ctlabeling(string(inpath), string(outpath), 100, 4);
    worker_finalize();
    return 0;
}
