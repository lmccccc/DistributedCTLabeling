#include "blogel_app_sssp.h"
int main(int argc, char* argv[])
{
    init_workers();
    blogel_app_sssp("../vor/mydis", "../exp/sssp");
    worker_finalize();
    return 0;
}
