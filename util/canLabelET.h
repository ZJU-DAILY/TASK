//
// Created by 罗程阳 on 2022/6/27.
//
#ifndef INSTANTQUERY_CANLABELET_H
#define INSTANTQUERY_CANLABELET_H

#include "ioUint.h"
#include <map>
#include <unordered_map>
#include <set>

class labelScore {
public:
    int x;  // vertex id
    double SRD;
    int PED;
    double score;
    vector<int> source;
public:
    labelScore(){}
    labelScore(int x, double SRD, int PED, double score);
    labelScore(int x);
    bool operator<(const labelScore &ls) const;
};

class labelPos {
public:
    int pos;
    int node;
public:
    labelPos(){}
    labelPos(int pos, int node);
    bool operator<(const labelPos &lp) const;
};

class sourceInfo {
public:
    int labelNum;
    vector<pair<int, int>> tag;
    int matchDepth;
    int matchPBNode;
    int qid;
public:
    sourceInfo();
    sourceInfo(int labelNum);
};

class canLabelET {
public:
    set<labelPos> labelPosList;
    map<int, sourceInfo> sourceNode;  // node -> source information
    set<labelScore> labelScores;
public:
    canLabelET(){};
    void insertLabelPos(int pos, int node);
public:
};

#endif //INSTANTQUERY_CANLABELET_H
