#include "blogel_ctlabeling.h"
int main(int argc, char* argv[])
{
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
    init_workers();
    blogel_ctlabeling("../exp/ctbin", "../exp/ctlabel", 4, 4);
    worker_finalize();
    return 0;
}
