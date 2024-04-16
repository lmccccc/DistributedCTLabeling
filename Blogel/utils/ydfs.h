#ifndef YDHDFS_H
#define YDHDFS_H

//#include "hdfs.h"
#include <string.h> //memcpy, memchr
#include <stdlib.h> //realloc
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>

#include "global.h"
using namespace std;


#define FS_BUF_SIZE 65536
#define LINE_DEFAULT_SIZE 4096
#define FS_BLOCK_SIZE 8388608 //8M
#define MODE (S_IRWXU | S_IRWXG | S_IRWXO)

typedef uint tOffset;
typedef int32_t tSize;
unsigned char kObjectKindFile = DT_REG;

//====== get File System ======

//fsFS getFS()
//{
//	hdfsFS fs = hdfsConnect("default", 0);
//	if(!fs) {
//		fprintf(stderr, "Failed to connect to HDFS!\n");
//		exit(-1);
//	}
//	return fs;
//}

//hdfsFS getlocalFS()
//{
//    hdfsFS lfs = hdfsConnect(NULL, 0);
//    if (!lfs) {
//        fprintf(stderr, "Failed to connect to 'local' FS!\n");
//        exit(-1);
//    }
//    return lfs;
//}

//====== get File Handle ======

FILE* getRHandle(const char* path)
{
    FILE* f = fopen(path, "r");
    if (f == NULL) {
        fprintf(stderr, "Failed to open %s for reading!\n", path);
        exit(-1);
    }
    return f;
}

FILE* getWHandle(const char* path)
{
    FILE* f = fopen(path, "w");
    if (f == NULL) {
        fprintf(stderr, "Failed to open %s for writing!\n", path);
        exit(-1);
    }
    return f;
}

FILE* getRWHandle(const char* path)
{
    FILE* f = fopen(path, "w+");
    if (f == NULL) {
        fprintf(stderr, "Failed to open %s!\n", path);
        exit(-1);
    }
    return f;
}

//====== Read line ======

//logic:
//buf[] is for batch reading from HDFS file
//line[] is a line buffer, the string length is "length", the buffer size is "size"
//after each readLine(), need to check eof(), if it's true, no line is read due to EOF
struct LineReader {
    //static fields
    char buf[FS_BUF_SIZE];
    tSize bufPos;
    tSize bufSize;
//    hdfsFS fs;
    FILE* handle;
    bool fileEnd;

    //dynamic fields
    char* line;
    int length;
    int size;

    LineReader(FILE* handle)
        : bufPos(0)
        , length(0)
        , size(LINE_DEFAULT_SIZE)
    {
//        this->fs = fs;
        this->handle = handle;
        fileEnd = false;
        fill();
        line = (char*)malloc(LINE_DEFAULT_SIZE * sizeof(char));
    }

    ~LineReader()
    {
        free(line);
    }

    //internal use only!
    void doubleLineBuf()
    {
        size *= 2;
        line = (char*)realloc(line, size * sizeof(char));
    }

    //internal use only!
    void lineAppend(char* first, int num)
    {
        while (length + num + 1 > size)
            doubleLineBuf();
        memcpy(line + length, first, num);
        length += num;
    }

    //internal use only!
    void fill()
    {
//        bufSize = hdfsRead(fs, handle, buf, HDFS_BUF_SIZE);
        bufSize = fread(buf, sizeof(char), FS_BUF_SIZE, handle);
        if (bufSize == -1) {
            fprintf(stderr, "Read Failure!\n");
            exit(-1);
        }
        bufPos = 0;
        if (bufSize < FS_BUF_SIZE)
            fileEnd = true;
    }

    //user interface
    //the line starts at "line", with "length" chars
    void readLine()
    {
        length = 0;
        if (bufPos == bufSize)
            return;
        char* pch = (char*)memchr(buf + bufPos, '\n', bufSize - bufPos);
        if (pch == NULL) {
            lineAppend(buf + bufPos, bufSize - bufPos);
            bufPos = bufSize;
            if (!fileEnd)
                fill();
            else
                return;
            pch = (char*)memchr(buf, '\n', bufSize);
            while (pch == NULL) {
                lineAppend(buf, bufSize);
                if (!fileEnd)
                    fill();
                else
                    return;
                pch = (char*)memchr(buf, '\n', bufSize);
            }
        }
        int validLen = pch - buf - bufPos;
        lineAppend(buf + bufPos, validLen);
        bufPos += validLen + 1; //+1 to skip '\n'
        if (bufPos == bufSize) {
            if (!fileEnd)
                fill();
            else
                return;
        }
    }

    char* getLine()
    {
        line[length] = '\0';
        return line;
    }

    bool eof()
    {
        return length == 0 && fileEnd;
    }
};

int rm_dir(std::string dir_full_path)
{
    DIR* dirp = opendir(dir_full_path.c_str());
    if(!dirp)
    {
        return -1;
    }
    struct dirent *dir;
    struct stat st;
    while((dir = readdir(dirp)) != NULL)
    {
        if(strcmp(dir->d_name,".") == 0
           || strcmp(dir->d_name,"..") == 0)
        {
            continue;
        }
        std::string sub_path = dir_full_path + '/' + dir->d_name;
        if(lstat(sub_path.c_str(),&st) == -1)
        {
            cout << "rm_dir:lstat " << sub_path << " error" << endl;
            continue;
        }
        if(S_ISDIR(st.st_mode))
        {
            if(rm_dir(sub_path) == -1) // 如果是目录文件，递归删除
            {
                closedir(dirp);
                return -1;
            }
            rmdir(sub_path.c_str());
        }
        else if(S_ISREG(st.st_mode))
        {
            unlink(sub_path.c_str());     // 如果是普通文件，则unlink
        }
        else
        {
            cout << "rm_dir:st_mode " << sub_path << " error" << endl;
            continue;
        }
    }

    if(rmdir(dir_full_path.c_str()) == -1)//delete dir itself.
    {
        cout << "rmdir fail" << endl;
        closedir(dirp);
        return -1;
    }

    closedir(dirp);
    return 0;
}

//====== Dir Check ======
int dirCheck(const char* indir, const char* outdir, bool print, bool force) //returns -1 if fail, 0 if succeed
{
//    hdfsFS fs = getHdfsFS();
    if (access(indir, F_OK) != 0) {
        if (print)
            fprintf(stderr, "Input path \"%s\" does not exist!\n", indir);
//        hdfsDisconnect(fs);
        return -1;
    }

    if (access(outdir, F_OK) == 0) {
        if (force) {
            if (rm_dir(outdir) == -1) {
                if (print)
                    fprintf(stderr, "Error deleting %s!\n", outdir);
                exit(-1);
            }
            int created = mkdir(outdir, MODE);
            if (created == -1) {
                if (print)
                    fprintf(stderr, "Failed to create folder %s!\n", outdir);
                exit(-1);
            }
        } else {
            if (print)
                fprintf(stderr, "Output path \"%s\" already exists!\n", outdir);
//            hdfsDisconnect(fs);
            return -1;
        }
    } else {

        int created = mkdir(outdir, MODE);
        if (created == -1) {
            if (print)
                fprintf(stderr, "Failed to create folder %s!\n", outdir);
            exit(-1);
        }
    }
//    hdfsDisconnect(fs);
    return 0;
}

int dirCheck(vector<string> indirs, const char* outdir, bool print, bool force) //returns -1 if fail, 0 if succeed
{
//    hdfsFS fs = getHdfsFS();
    for (int i = 0; i < indirs.size(); i++) {
        const char* indir = indirs[i].c_str();
        if (access(indir, F_OK) != 0) {
            if (print)
                fprintf(stderr, "Input path \"%s\" does not exist!\n", indir);
//            hdfsDisconnect(fs);
            return -1;
        }
    }
    if (access(outdir, F_OK) == 0) {
        if (force) {
            if (rm_dir(outdir) == -1) {
                if (print)
                    fprintf(stderr, "Error deleting %s!\n", outdir);
                exit(-1);
            }
            int created = mkdir(outdir, MODE);
            if (created == -1) {
                if (print)
                    fprintf(stderr, "Failed to create folder %s!\n", outdir);
                exit(-1);
            }
        } else {
            if (print)
                fprintf(stderr, "Output path \"%s\" already exists!\n", outdir);
//            hdfsDisconnect(fs);
            return -1;
        }
    } else {
        int created = mkdir(outdir, MODE);
        if (created == -1) {
            if (print)
                fprintf(stderr, "Failed to create folder %s!\n", outdir);
            exit(-1);
        }
    }
//    hdfsDisconnect(fs);
    return 0;
}

int dirCheck(const char* outdir, bool force) //returns -1 if fail, 0 if succeed
{
//    hdfsFS fs = getHdfsFS();

    if (access(outdir, F_OK) == 0) {
        if (force) {
            if (rm_dir(outdir) == -1) {
                fprintf(stderr, "Error deleting %s!\n", outdir);
                exit(-1);
            }
            int created = mkdir(outdir, MODE);
            if (created == -1) {
                fprintf(stderr, "Failed to create folder %s!\n", outdir);
                exit(-1);
            }
        } else {
            fprintf(stderr, "Output path \"%s\" already exists!\n", outdir);
//            hdfsDisconnect(fs);
            return -1;
        }
    } else {
        int created = mkdir(outdir, MODE);
        if (created == -1) {
            fprintf(stderr, "Failed to create folder %s!\n", outdir);
            exit(-1);
        }
    }
//    hdfsDisconnect(fs);
    return 0;
}

//====== Write line ======

const char* newLine = "\n";

struct LineWriter {
//    hdfsFS fs;
    const char* path;
    int me; //-1 if there's no concept of machines (like: hadoop fs -put)
    int nxtPart;
    int curSize;

    FILE* curHdl;

    LineWriter(const char* path, int me)
        : nxtPart(0)
        , curSize(0)
    {
        this->path = path;
//        this->fs = fs;
        this->me = me;
        curHdl = NULL;
        //===============================
        //if(overwrite==true) readDirForce();
        //else readDirCheck();
        //===============================
        //1. cannot use above, otherwise multiple dir check/delete will be done during parallel writing
        //2. before calling the constructor, make sure "path" does not exist
        nextHdl();
    }

    ~LineWriter()
    {
//        if (hdfsFlush(fs, curHdl)) {
//            fprintf(stderr, "Failed to 'flush' %s\n", path);
//            exit(-1);
//        }
//        hdfsCloseFile(fs, curHdl);
        fclose(curHdl);
    }

    /*//================== not for parallel writing =====================
    //internal use only!
    void readDirCheck()
{
    	if(access(path, F_OK)==0)
    	{
    		fprintf(stderr, "%s already exists!\n", path);
    		exit(-1);
    	}
    	int created=mkdir(path, MODE);
    	if(created==-1)
    	{
    		fprintf(stderr, "Failed to create folder %s!\n", path);
    		exit(-1);
    	}
}

    //internal use only!
    void readDirForce()
{
    	if(access(path, F_OK)==0)
    	{
    		if(hdfsDelete(fs, path)==-1)
    		{
    			fprintf(stderr, "Error deleting %s!\n", path);
    			exit(-1);
    		}
    	}
    	int created=mkdir(path, MODE);
    	if(created==-1)
    	{
    		fprintf(stderr, "Failed to create folder %s!\n", path);
    		exit(-1);
    	}
}
    */ //================== not for parallel writing =====================

    //internal use only!
    void nextHdl()
    {
        //set fileName
        char fname[20];
        strcpy(fname, "part_");
        char buffer[10];
        if (me >= 0) {
            sprintf(buffer, "%d", me);
            strcat(fname, buffer);
            strcat(fname, "_");
        }
        sprintf(buffer, "%d", nxtPart);
        strcat(fname, buffer);
        //flush old file
        if (nxtPart > 0) {
//            if (hdfsFlush(fs, curHdl)) {
//                fprintf(stderr, "Failed to 'flush' %s\n", path);
//                exit(-1);
//            }
            fclose(curHdl);
        }
        //open new file
        nxtPart++;
        curSize = 0;
        char* filePath = new char[strlen(path) + strlen(fname) + 2];
        strcpy(filePath, path);
        strcat(filePath, "/");
        strcat(filePath, fname);
        curHdl = getWHandle(filePath);
        delete[] filePath;
    }

    void writeLine(char* line, int num)
    {
        if (curSize + num + 1 > FS_BLOCK_SIZE) //+1 because of '\n'
        {
            nextHdl();
        }
        tSize numWritten = fwrite(line, sizeof(char), num, curHdl);
        if (numWritten == -1) {
            fprintf(stderr, "Failed to write file!\n");
            exit(-1);
        }
        curSize += numWritten;
        numWritten = fwrite(newLine, sizeof(char), 1, curHdl);
        if (numWritten == -1) {
            fprintf(stderr, "Failed to create a new line!\n");
            exit(-1);
        }
        curSize += 1;
    }
};

//====== Put: local->HDFS ======

void put(char* localpath, char* fspath)
{
    if (dirCheck(fspath, false) == -1)
        return;
//    hdfsFS fs = getHdfsFS();
//    hdfsFS lfs = getlocalFS();

    FILE* in = getRHandle(localpath);
    LineReader* reader = new LineReader(in);
    LineWriter* writer = new LineWriter(fspath, -1);
    while (true) {
        reader->readLine();
        if (!reader->eof()) {
            writer->writeLine(reader->line, reader->length);
        } else
            break;
    }
//    hdfsCloseFile(lfs, in);
    fclose(in);
    delete reader;
    delete writer;

//    hdfsDisconnect(lfs);
//    hdfsDisconnect(fs);
}

void putf(char* localpath, char* fspath) //force put, overwrites target
{
    dirCheck(fspath, true);
//    hdfsFS fs = getHdfsFS();
//    hdfsFS lfs = getlocalFS();

    FILE* in = getRHandle(localpath);
    LineReader* reader = new LineReader(in);
    LineWriter* writer = new LineWriter(fspath, -1);
    while (true) {
        reader->readLine();
        if (!reader->eof()) {
            writer->writeLine(reader->line, reader->length);
        } else
            break;
    }
//    hdfsCloseFile(lfs, in);
    fclose(in);
    delete reader;
    delete writer;

//    hdfsDisconnect(lfs);
//    hdfsDisconnect(fs);
}

//====== BufferedWriter ======
struct BufferedWriter {
//    hdfsFS fs;
    const char* path;
    int me; //-1 if there's no concept of machines (like: hadoop fs -put)
    int nxtPart;
    vector<char> buf;
    FILE* curHdl;

    BufferedWriter(const char* path)
    {
        this->path = path;
        this->me = -1;
        this->curHdl = getWHandle(this->path);
    }
    BufferedWriter(const char* path, int me)
        : nxtPart(0)
    {
        this->path = path;
        this->me = me;
        curHdl = NULL;
        nextHdl();
    }

    ~BufferedWriter()
    {
        tSize numWritten = fwrite(&buf[0], sizeof(char), buf.size(), curHdl);
        if (numWritten == -1) {
            fprintf(stderr, "Failed to write file!\n");
            exit(-1);
        }
        buf.clear();

//        if (fsFlush(fs, curHdl)) {
//            fprintf(stderr, "Failed to 'flush' %s\n", path);
//            exit(-1);
//        }
//        hdfsCloseFile(fs, curHdl);
        fclose(curHdl);
    }

    //internal use only!
    void nextHdl()
    {
        //set fileName
        char fname[20];

        if (me >= 0) {
            sprintf(fname, "part_%d_%d", me, nxtPart);
        } else {
            sprintf(fname, "part_%d", nxtPart);
        }

        //flush old file
        if (nxtPart > 0) {
//            if (hdfsFlush(curHdl)) {
//                fprintf(stderr, "Failed to 'flush' %s\n", path);
//                exit(-1);
//            }
//            hdfsCloseFile(fs, curHdl);
            fclose(curHdl);
        }
        //open new file
        nxtPart++;
        char* filePath = new char[strlen(path) + strlen(fname) + 2];
        sprintf(filePath, "%s/%s", path, fname);
        curHdl = getWHandle(filePath);
        delete[] filePath;
    }

    void check()
    {
        if (buf.size() >= FS_BLOCK_SIZE) {
            tSize numWritten = fwrite(&buf[0], sizeof(char), buf.size(), curHdl);
            if (numWritten == -1) {
                fprintf(stderr, "Failed to write file!\n");
                exit(-1);
            }
            buf.clear();
            if (me != -1) // -1 means "output in the specified file only"
            {
                nextHdl();
            }
        }
    }

    void write(const char* content)
    {
        int len = strlen(content);
        buf.insert(buf.end(), content, content + len);
    }
};

//====== Dispatcher ======

struct sizedFName {
    char* fname;
    tOffset size;

    bool operator<(const sizedFName& o) const
    {
        return size > o.size; //large file goes first
    }
};

struct sizedFString {
    string fname;
    tOffset size;

    bool operator<(const sizedFString& o) const
    {
        return size > o.size; //large file goes first
    }
};

const char* rfind(const char* str, char delim)
{
    int len = strlen(str);
    int pos = 0;
    for (int i = len - 1; i >= 0; i--) {
        if (str[i] == delim) {
            pos = i;
            break;
        }
    }
    return str + pos;
}


typedef struct fileInfo{
    unsigned char mKind;
    off_t mSize;
    char mName[256];

    fileInfo(const unsigned short kind, char* name, const off_t size){
        this->mKind = kind;
        this->mSize = size;
        memcpy(this->mName, name, 256);
    }
};

typedef struct direntInfo{
    unsigned char type;
//    char path[256];
    char name[256];

    direntInfo(const unsigned char type, const char* name){
        this->type = type;
        memcpy(this->name, name, 256);
    }
};

int getFiles(const string& path, vector<direntInfo> &files, vector<string> suffixs = {})
{
    int FileCnt = 0;
    DIR *dirp;
    struct dirent *dp;

    dirp = opendir(path.c_str());
    if (dirp == NULL) {
        printf("opendir %s failed\n",path.c_str());
        return -1;
    }

    while ((dp = readdir(dirp)) != NULL) {
        string curpath(path);
        if (path.back() != '/') {
            curpath += string("/") += dp->d_name;
        } else {
            curpath += dp->d_name;
        }
        //如果是目录，递归查找
        if(dp->d_type == DT_DIR) {
//            if(0 != strcmp(dp->d_name,".") && 0 != strcmp(dp->d_name,"..")){
//                FileCnt += getFiles(curpath, files, suffixs);
//            }
        }
            //判断是否为文件以及文件后缀名
        else if (dp->d_type == DT_REG) {
            if (suffixs.size() <=0 ) {
                direntInfo tempf(dp->d_type, curpath.c_str());
                files.push_back(tempf);
                FileCnt++;
            } else {
                for (auto suffix : suffixs) {
                    if (string(strrchr(dp->d_name,'.')) == suffix) {
                        direntInfo tempf = {dp->d_type, curpath.c_str()};
                        files.push_back(tempf);
                        FileCnt++;
                        break;
                    }
                }
            }
        }
    }

    closedir(dirp);
    return FileCnt;
}

//从指定路径获取文件路径（有序）
//int getFilesOrdered(const string& path, vector<direntInfo> &files, vector<string> suffixs = {}){
//    int FileCnt = 0;
//    FileCnt = getFiles(path, files, suffixs);
//    sort(files.begin(), files.end());
//    return FileCnt;
//}


vector<fileInfo> fsListDirectory(const char* path, int* numFiles){
    vector<direntInfo> files;
    vector<fileInfo> infoList;
    int FileCnt = getFiles(path, files);
    struct stat st;//定义结构体变量，保存所获取的文件属性
    for(auto file : files){

        int res = stat(file.name, &st);
        if(res == -1)//获取文件属性失败，errno设置为合适的值
        {
            perror("stat fail");
            exit(1);
        }
        fileInfo t = {file.type, file.name, st.st_size};
        infoList.push_back(t);
    }
    *numFiles = files.size();
    return infoList;
}

vector<string>* dispatchRan(const char* inDir, int numSlaves) //remember to "delete[] assignment" after used
{ //locality is not considered for simplicity
    vector<string>* assignment = new vector<string>[numSlaves];
//    hdfsFS fs = getHdfsFS();
    int numFiles;
    vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);
    if (fileinfo.empty()) {
        fprintf(stderr, "Failed to list directory %s!\n", inDir);
        exit(-1);
    }
    tOffset* assigned = new tOffset[numSlaves];
    for (int i = 0; i < numSlaves; i++)
        assigned[i] = 0;
    //sort files by size
    vector<sizedFName> sizedfile;
    for (int i = 0; i < numFiles; i++) {
        if (fileinfo[i].mKind == kObjectKindFile) {
            sizedFName cur = { fileinfo[i].mName, fileinfo[i].mSize };
            sizedfile.push_back(cur);
        }
    }
    sort(sizedfile.begin(), sizedfile.end());
    //allocate files to slaves
    vector<sizedFName>::iterator it;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < numSlaves; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
    }
    delete[] assigned;
//    hdfsFreeFileInfo(fileinfo, numFiles);
    return assignment;
}

//considers locality
//1. compute avg size, define it as quota
//2. sort files by size
//3. for each file, if its slave has quota, assign it to the slave
//4. for the rest, run the greedy assignment
//(libhdfs do not have location info, but we can check slaveID from fileName)
//*** NOTE: NOT SUITABLE FOR DATA "PUT" TO HDFS, ONLY FOR DATA PROCESSED BY AT LEAST ONE JOB
vector<string>* dispatchLocality(const char* inDir, int numSlaves) //remember to "delete[] assignment" after used
{ //considers locality
    vector<string>* assignment = new vector<string>[numSlaves];
//    hdfsFS fs = getHdfsFS();
    int numFiles;
    vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);
    if (fileinfo.empty()) {
        fprintf(stderr, "Failed to list directory %s!\n", inDir);
        exit(-1);
    }
    tOffset* assigned = new tOffset[numSlaves];
    for (int i = 0; i < numSlaves; i++)
        assigned[i] = 0;
    //sort files by size
    vector<sizedFName> sizedfile;
    int avg = 0;
    for (int i = 0; i < numFiles; i++) {
        if (fileinfo[i].mKind == kObjectKindFile) {
            sizedFName cur = { fileinfo[i].mName, fileinfo[i].mSize };
            sizedfile.push_back(cur);
            avg += fileinfo[i].mSize;
        }
    }
    avg /= numSlaves;
    sort(sizedfile.begin(), sizedfile.end());
    //allocate files to slaves
    vector<sizedFName>::iterator it;
    vector<sizedFName> recycler;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        istringstream ss(rfind(it->fname, '/'));
        string cur;
        getline(ss, cur, '_');
        getline(ss, cur, '_');
        int slaveOfFile = atoi(cur.c_str());
        if (assigned[slaveOfFile] + it->size <= avg) {
            assignment[slaveOfFile].push_back(it->fname);
            assigned[slaveOfFile] += it->size;
        } else
            recycler.push_back(*it);
    }
    for (it = recycler.begin(); it != recycler.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < numSlaves; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
    }
    delete[] assigned;
//    hdfsFreeFileInfo(fileinfo, numFiles);
    return assignment;
}

vector<vector<string> >* dispatchRan(const char* inDir) //remember to delete assignment after used
{ //locality is not considered for simplicity
    vector<vector<string> >* assignmentPointer = new vector<vector<string> >(_num_workers);
    vector<vector<string> >& assignment = *assignmentPointer;
//    hdfsFS fs = getHdfsFS();
    int numFiles;
    vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);


    if (fileinfo.empty()) {
        fprintf(stderr, "Failed to list directory %s!\n", inDir);
        exit(-1);
    }

    tOffset* assigned = new tOffset[_num_workers];
    for (int i = 0; i < _num_workers; i++)
        assigned[i] = 0;
    //sort files by size
    vector<sizedFName> sizedfile;

    for (int i = 0; i < numFiles; i++) {
        if (fileinfo[i].mKind == kObjectKindFile) {
            sizedFName cur = { fileinfo[i].mName, fileinfo[i].mSize };
            sizedfile.push_back(cur);
        }
    }

    sort(sizedfile.begin(), sizedfile.end());



    //allocate files to slaves
    vector<sizedFName>::iterator it;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < _num_workers; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
#ifdef _DEBUG
        cout << "assignment " << min << " " << it->fname << endl;
#endif
    }
    delete[] assigned;
//    hdfsFreeFileInfo(fileinfo, numFiles);
    return assignmentPointer;
}

vector<vector<string> >* dispatchRan(vector<string> inDirs) //remember to delete assignment after used
{ //locality is not considered for simplicity
    vector<vector<string> >* assignmentPointer = new vector<vector<string> >(_num_workers);
    vector<vector<string> >& assignment = *assignmentPointer;
//    hdfsFS fs = getHdfsFS();
    vector<sizedFString> sizedfile;
    for (int pos = 0; pos < inDirs.size(); pos++) {
        const char* inDir = inDirs[pos].c_str();
        int numFiles;
        vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);
        if (fileinfo.empty()) {
            fprintf(stderr, "Failed to list directory %s!\n", inDir);
            exit(-1);
        }
        for (int i = 0; i < numFiles; i++) {
            if (fileinfo[i].mKind == kObjectKindFile) {
                sizedFString cur = { fileinfo[i].mName, fileinfo[i].mSize };
                sizedfile.push_back(cur);
            }
        }
//        hdfsFreeFileInfo(fileinfo, numFiles);
    }
    //sort files by size
    sort(sizedfile.begin(), sizedfile.end());
    tOffset* assigned = new tOffset[_num_workers];
    for (int i = 0; i < _num_workers; i++)
        assigned[i] = 0;
    //allocate files to slaves
    vector<sizedFString>::iterator it;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < _num_workers; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
    }
    delete[] assigned;
    return assignmentPointer;
}

//considers locality
//1. compute avg size, define it as quota
//2. sort files by size
//3. for each file, if its slave has quota, assign it to the slave
//4. for the rest, run the greedy assignment
//(libhdfs do not have location info, but we can check slaveID from fileName)
//*** NOTE: NOT SUITABLE FOR DATA "PUT" TO HDFS, ONLY FOR DATA PROCESSED BY AT LEAST ONE JOB
vector<vector<string> >* dispatchLocality(const char* inDir) //remember to delete assignment after used
{ //considers locality
    vector<vector<string> >* assignmentPointer = new vector<vector<string> >(_num_workers);
    vector<vector<string> >& assignment = *assignmentPointer;
//    hdfsFS fs = getHdfsFS();
    int numFiles;
    vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);
    if (fileinfo.empty()) {
        fprintf(stderr, "Failed to list directory %s!\n", inDir);
        exit(-1);
    }
    tOffset* assigned = new tOffset[_num_workers];
    for (int i = 0; i < _num_workers; i++)
        assigned[i] = 0;
    //sort files by size
    vector<sizedFName> sizedfile;
    int avg = 0;
    for (int i = 0; i < numFiles; i++) {
        if (fileinfo[i].mKind == kObjectKindFile) {
            sizedFName cur = { fileinfo[i].mName, fileinfo[i].mSize };
            sizedfile.push_back(cur);
            avg += fileinfo[i].mSize;
        }
    }
    avg /= _num_workers;
    sort(sizedfile.begin(), sizedfile.end());
    //allocate files to slaves
    vector<sizedFName>::iterator it;
    vector<sizedFName> recycler;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        istringstream ss(rfind(it->fname, '/'));
        string cur;
        getline(ss, cur, '_');
        getline(ss, cur, '_');
        int slaveOfFile = atoi(cur.c_str());
        if (assigned[slaveOfFile] + it->size <= avg) {
            assignment[slaveOfFile].push_back(it->fname);
            assigned[slaveOfFile] += it->size;
        } else
            recycler.push_back(*it);
    }
    for (it = recycler.begin(); it != recycler.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < _num_workers; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
    }
    delete[] assigned;
//    hdfsFreeFileInfo(fileinfo, numFiles);
    return assignmentPointer;
}

vector<vector<string> >* dispatchLocality(vector<string> inDirs) //remember to delete assignment after used
{ //considers locality
    vector<vector<string> >* assignmentPointer = new vector<vector<string> >(_num_workers);
    vector<vector<string> >& assignment = *assignmentPointer;
//    hdfsFS fs = getHdfsFS();
    vector<sizedFString> sizedfile;
    int avg = 0;
    for (int pos = 0; pos < inDirs.size(); pos++) {
        const char* inDir = inDirs[pos].c_str();
        int numFiles;
        vector<fileInfo> fileinfo = fsListDirectory(inDir, &numFiles);
        if (fileinfo.empty()) {
            fprintf(stderr, "Failed to list directory %s!\n", inDir);
            exit(-1);
        }
        for (int i = 0; i < numFiles; i++) {
            if (fileinfo[i].mKind == kObjectKindFile) {
                sizedFString cur = { fileinfo[i].mName, fileinfo[i].mSize };
                sizedfile.push_back(cur);
                avg += fileinfo[i].mSize;
            }
        }
//        hdfsFreeFileInfo(fileinfo, numFiles);
    }
    tOffset* assigned = new tOffset[_num_workers];
    for (int i = 0; i < _num_workers; i++)
        assigned[i] = 0;
    //sort files by size
    avg /= _num_workers;
    sort(sizedfile.begin(), sizedfile.end());
    //allocate files to slaves
    vector<sizedFString>::iterator it;
    vector<sizedFString> recycler;
    for (it = sizedfile.begin(); it != sizedfile.end(); ++it) {
        istringstream ss(rfind(it->fname.c_str(), '/'));
        string cur;
        getline(ss, cur, '_');
        getline(ss, cur, '_');
        int slaveOfFile = atoi(cur.c_str());
        if (assigned[slaveOfFile] + it->size <= avg) {
            assignment[slaveOfFile].push_back(it->fname);
            assigned[slaveOfFile] += it->size;
        } else
            recycler.push_back(*it);
    }
    for (it = recycler.begin(); it != recycler.end(); ++it) {
        int min = 0;
        tOffset minSize = assigned[0];
        for (int j = 1; j < _num_workers; j++) {
            if (minSize > assigned[j]) {
                min = j;
                minSize = assigned[j];
            }
        }
        assignment[min].push_back(it->fname);
        assigned[min] += it->size;
    }
    delete[] assigned;
    return assignmentPointer;
}

void reportAssignment(vector<string>* assignment, int numSlaves)
{
    for (int i = 0; i < numSlaves; i++) {
        cout << "====== Rank " << i << " ======" << endl;
        vector<string>::iterator it;
        for (it = assignment[i].begin(); it != assignment[i].end(); ++it) {
            cout << *it << endl;
        }
    }
}

void reportAssignment(vector<vector<string> >* assignment)
{
    for (int i = 0; i < _num_workers; i++) {
        cout << "====== Rank " << i << " ======" << endl;
        vector<string>::iterator it;
        for (it = (*assignment)[i].begin(); it != (*assignment)[i].end(); ++it) {
            cout << *it << endl;
        }
    }
}

#endif
