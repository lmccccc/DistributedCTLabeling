#include "../Blogel/blogel/Voronoi.h"
#include <iostream>
#include <sstream>
#include "../Blogel/blogel/BGlobal.h"
using namespace std;

class vorPart : public BPartWorker {
    char buf[1000];

public:

    virtual BPartVertex* toVertex(char* line)//id x \t y cnt neighbor1 len neighbor2 len ...
    {
        char* pch;
//        BPartVertex* v = new BPartVertex;
//        v->value().content = line; //first set content!!! line will change later due to "strtok"
//        pch = strtok(line, " ");
//        v->id = atoi(pch);
//        pch = strtok(NULL, " "); //filter x
//        pch = strtok(NULL, "\t"); //filter y
//        pch = strtok(NULL, " ");
//        int num = atoi(pch);
//
//        while (num --) {
//            pch = strtok(NULL, " ");
//            int nb = atoi(pch);
//            v->value().neighbors.push_back(nb);
//            strtok(NULL, " "); //edge length
//        }

        pch = strtok(line, "\t");
        BPartVertex* v = new BPartVertex;
        v->id = atoi(pch);
        pch = strtok(NULL, " ");
        int num = atoi(pch);
        //v->value().color=-1;//default is -1
        while (num--) {
            pch = strtok(NULL, " ");
            v->value().neighbors.push_back(atoi(pch));
        }

        return v;
    }

    virtual void toline(BPartVertex* v, BufferedWriter& writer) //key: "vertexID blockID slaveID"
    { //val: list of "vid block slave "
        sprintf(buf, "%d %d %d\t", v->id, v->value().color, _my_rank);
        writer.write(buf);

        vector<triplet>& vec = v->value().nbsInfo;
        hash_map<int, triplet> map;
        for (int i = 0; i < vec.size(); i++) {
            map[vec[i].vid] = vec[i];
        }
        ////////
        stringstream ss(v->value().content);
        string token;
        ss >> token; //vid
        ss >> token; //x
        ss >> token; //y
        ss >> token; //number
        int num = v->value().neighbors.size();
        for (int i = 0; i < num; i++) {
            ss >> token;
            int vid = atoi(token.c_str());
            ss >> token;
            double elen = atof(token.c_str());
            triplet trip = map[vid];
//            sprintf(buf, "%d %f %d %d ", vid, elen, trip.bid, trip.wid);

            sprintf(buf, "%d %d %d ", vec[i].vid, vec[i].bid, vec[i].wid);
            writer.write(buf);
        }
        writer.write("\n");
    }
};

int blogel_sssp_vorPart(string in_path, string out_path)
{
    WorkerParams param;
    param.input_path = in_path;
    param.output_path = out_path;
    param.force_write = true;
    param.native_dispatcher = false;
    bool to_undirected = false;
    //////
    set_sampRate(0.01);
    set_maxHop(50);
    set_maxVCSize(1000000000);
    set_factor(2);
    set_stopRatio(0.9);
    set_maxRate(0.1);
    //////
    vorPart worker;
    worker.run(param, to_undirected);
}
