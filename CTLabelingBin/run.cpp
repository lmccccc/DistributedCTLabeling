#include "ctlabeling_create_bin.h"
int main(int argc, char* argv[])
{
    init_workers();
    char *inpath = argv[1];
    char *outpath = argv[2];
    char *distpath = argv[3];
    if(argc < 4 || inpath == NULL || outpath == NULL){
        cout << "invalid input path" << endl;
        exit(-1);
    }

    if(get_worker_id() == MASTER_RANK) {
        cout << "in dir=" << inpath << endl;
        cout << "out dir=" << outpath << endl;
        cout << "dist file=" << distpath << endl;
    }
    blogel_CTlabeling_create_bin(string(inpath), string(outpath), string(distpath));
    worker_finalize();
    return 0;
}
