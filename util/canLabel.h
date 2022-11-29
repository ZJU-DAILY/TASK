//
// Created by 罗程阳 on 2022/6/17.
//

#ifndef INSTANTQUERY_CANLABEL_H
#define INSTANTQUERY_CANLABEL_H

#include "ioUint.h"
#include <map>
#include <unordered_map>
#include <set>
class canLabel {
public:
    edgeS labelForward;
    int cid;
    set<label> labels;
    set<labelS> labelScore;
    map<int, vector<edgeS>> sourceNode;
public:
    canLabel(){};
    canLabel(edgeS &labelForward);
};

#endif //INSTANTQUERY_CANLABEL_H
