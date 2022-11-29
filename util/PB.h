//
// Created by 罗程阳 on 2022/5/27.
//

#ifndef INSTANTQUERY_PB_H
#define INSTANTQUERY_PB_H
#include <assert.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <string>
#include <iostream>
#include <climits>
#include "ioUint.h"

using namespace std;

class BitMap {
public:
    vector<unsigned int> bits;
    int size;  // the number of unsigned int to store the bit map
    int bitNum;  // the total number of bit
public:
    BitMap();
    BitMap(int range);
    void setOne(int num);  // set bit num to 1
    void setZero(int num);  // set bit num to 0
    void setAll();  // set all bits to 1
    bool hasZero();  // whether a bitmap contains 0
    bool allZero();  // whether all bits in a bit map are 0
    int bitValueAt(int num);  // get the value at num position
    int findOne(int num);
    int nextOne(int num);
    void nextOne(int &num, int &labelPos, int &countOne);
    void insertOne(int pos);
    void insertZero(int pos);
    void deleteBit(int pos);
    void updateTag(pair<int, int> &tag, int &pos);  // update the tag when search from top to bottom
    int countOne();  // get the total
    int countOne(int pos);
    void compress(BitMap &bitMap);
    void change(BitMap &bitMap);
    void performOr(BitMap &bitMap);
    void performAnd(BitMap &bitMap);
    void uncompress(BitMap &bitMap);
    string show();  // show the concrete content of the bit map
//    bool isEqual(BitMap &bitMap);
//    void assignValue(BitMap &bitMap);
    BitMap& operator= (const BitMap &bm);
    bool operator==(const BitMap &bm);
    bool operator<(const BitMap &bm) const;
    ~BitMap();
};

class PB {
public:
    int vertexId;
    vector<edgeS> backwardList;
    set<int> nodeSet;
//    vector<int> bitMapID;
    vector<vector<BitMap>> bitMaps;
//    vector<BitMap> bitMaps;
    unordered_map<int, vector<int>> nodeBitMapId;  // the bit map from root to a node
public:
    PB() {}
    PB(int id) {
        vertexId = id;
    }
};


#endif //INSTANTQUERY_PB_H
