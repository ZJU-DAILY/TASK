//
// Created by 罗程阳 on 2022/6/8.
//

#ifndef INSTANTQUERY_IOPB_H
#define INSTANTQUERY_IOPB_H

#include "ioUint.h"
#include "PB.h"

struct outBufPB {
    FILE *outFile;
    int numPB;  // the total num of PB
    outBufPB() {}
    outBufPB(char *fileName);
    ~outBufPB();
    void insertInt(int num);
    void insertBackwardList(vector<edgeS> &backwardList);
    void insertNodeSet(set<int> &nodeSet);
    void insertBitMaps(vector<vector<BitMap>> &bitMapSet);
    void insertBitMapID(vector<int> &bitMapID);
    void insertNodeBitMapId(unordered_map<int, vector<int>> &nodeBitMapId);
};

struct inBufPB {
    FILE *inFile;
    inBufPB() {}
    inBufPB(char *fileName);
    void nextInt(int &num);
    void getBackwardList(vector<edgeS> &backwardList);
    void getNodeSet(set<int> &nodeSet);
    void getBitMaps(vector<vector<BitMap>> &bitMaps);
    void getNodeBitMapId(unordered_map<int, vector<int>> &nodeBitMapId);
    ~inBufPB();
};
#endif //INSTANTQUERY_IOPB_H
