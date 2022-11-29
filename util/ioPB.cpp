//
// Created by 罗程阳 on 2022/6/8.
//
//#pragma GCC optimize(2)

#include "ioPB.h"

outBufPB::outBufPB(char *fileName) {
    outFile = fopen(fileName, "wb");
}

void outBufPB::insertInt(int num) {
    fwrite(&num, sizeof(int), 1, outFile);
}

void outBufPB::insertBackwardList(vector<edgeS> &backwardList) {
    int count = backwardList.size();
    edgeS *buf;
    buf = (edgeS *) malloc(sizeof(edgeS) * count);
    for (int i = 0; i < count; i++) {
        buf[i].x = backwardList[i].x;
        buf[i].w = backwardList[i].w;
    }
    fwrite(buf, sizeof(edgeS), count, outFile);
    free(buf);
}

void outBufPB::insertNodeSet(set<int> &nodeSet) {
    int count = nodeSet.size();
    int *buf;
    buf = (int *) malloc(sizeof(int) * count);
    int i = 0;
    for (auto it = nodeSet.begin(); it != nodeSet.end(); it++) {
        buf[i++] = *it;
    }
    fwrite(buf, sizeof(int), count, outFile);
    free(buf);
}

void outBufPB::insertBitMapID(vector<int> &bitMapID) {
    int count = bitMapID.size();
    insertInt(count);
    for (int i = 0; i < count; i++) {
        insertInt(bitMapID[i]);
    }
}

void outBufPB::insertBitMaps(vector<vector<BitMap>> &bitMapSet) {
    int count = bitMapSet.size();
    insertInt(count);
    for (auto &item: bitMapSet) {
        count = item.size();
        insertInt(count);
        for (auto &item1 : item) {
            insertInt(item1.bitNum);
            int count1 = item1.size;
            insertInt(count1);
            unsigned int *buf;
            buf = (unsigned int *) malloc(sizeof(unsigned int) * count1);
            for (int k = 0; k < count1; k++) {
                buf[k] = item1.bits[k];
            }
            fwrite(buf, sizeof(unsigned int), count1, outFile);
            free(buf);
        }
    }
//    for (auto &item : bitMapSet) {
//        insertInt(item.bitNum);
//        int count1 = item.size;
//        insertInt(count1);
//        unsigned int *buf;
//        buf = (unsigned int *) malloc(sizeof(unsigned int) * count1);
//        for (int k = 0; k < count1; k++) {
//            buf[k] = item.bits[k];
//        }
//        fwrite(buf, sizeof(unsigned int), count1, outFile);
//        free(buf);
//    }
}

void outBufPB::insertNodeBitMapId(unordered_map<int, vector<int>> &nodeBitMapId) {
    int count = nodeBitMapId.size();
    insertInt(count);
    for (auto i = nodeBitMapId.begin(); i != nodeBitMapId.end(); i++) {
        insertInt((int &) i->first);
        int count2 = i->second.size();
        insertInt(count2);
        int *buf;
        buf = (int *) malloc(sizeof(int) * count2);
        for (int j = 0; j < count2; j++) {
            buf[j] = i->second[j];
        }
        fwrite(buf, sizeof(int), count2, outFile);
        free(buf);
    }
}

outBufPB::~outBufPB() {
    fflush(outFile);
    if (outFile != NULL) fclose(outFile);
}

inBufPB::inBufPB(char *fileName) {
    inFile = fopen(fileName, "rb");
}

void inBufPB::nextInt(int &num) {
    fread(&num, sizeof(int), 1, inFile);
}

void inBufPB::getBackwardList(vector<edgeS> &backwardList) {
    int size;
    nextInt(size);
    edgeS* buf;
    buf = (edgeS *) malloc(sizeof(edgeS) * size);
    fread(buf, sizeof(edgeS), size, inFile);
    for (int i = 0; i < size; i++) {
        backwardList.push_back(buf[i]);
    }
//    vector<edgeS>(backwardList).swap(backwardList);
    free(buf);
}

void inBufPB::getNodeSet(set<int> &nodeSet) {
    int size;
    nextInt(size);
    int *buf;
    buf = (int *) malloc(sizeof(int) * size);
    fread(buf, sizeof(int), size, inFile);
    for (int i = 0; i < size; i++) {
        nodeSet.insert(buf[i]);
    }
    free(buf);
}

void inBufPB::getBitMaps(vector<vector<BitMap>> &bitMaps) {
    int size;
    nextInt(size);
    bitMaps.resize(size);
    for (int i = 0; i < size; i++) {
        int size1;
        nextInt(size1);
        for (int j = 0; j < size1; j++) {
            int bitNum, bitSize;
            nextInt(bitNum);
            nextInt(bitSize);  // read in the size of unsigned int array
            BitMap tempBitMap(bitNum);
            unsigned int *buf;
            buf = (unsigned int *) malloc(sizeof(unsigned int) * bitSize);
            fread(buf, sizeof(unsigned int), bitSize, inFile);
            for (int k = 0; k < bitSize; k++) {
                tempBitMap.bits[k] = buf[k];
            }
            bitMaps[i].push_back(tempBitMap);
            free(buf);
        }
    }
}

void inBufPB::getNodeBitMapId(unordered_map<int, vector<int>> &nodeBitMapId) {
    int size;
    nextInt(size);
    nodeBitMapId.reserve(size);
    for (int i = 0; i < size; i++) {
        int nodeId, size2;
        vector<int> bitMapIds;
        nextInt(nodeId);
        nextInt(size2);
        bitMapIds.reserve(size2);
        int *buf;
        buf = (int *) malloc(sizeof(int) * size2);
        fread(buf, sizeof(int), size2, inFile);
        for (int j = 0; j < size2; j++) {
            bitMapIds.push_back(buf[j]);
        }
        nodeBitMapId[nodeId] = bitMapIds;
        // free
        free(buf);
        bitMapIds.clear();
        vector<int>().swap(bitMapIds);
    }
}

inBufPB::~inBufPB() {
    if (inFile != NULL)
        fclose(inFile);
}