//
// Created by 罗程阳 on 2022/6/27.
//

#include "canLabelET.h"

labelScore::labelScore(int x, double SRD, int PED, double score) {
    this->x = x;
    this->SRD = SRD;
    this->PED = PED;
    this->score = score;
}

labelScore::labelScore(int x) {
    this->x = x;
}

bool labelScore::operator<(const labelScore &ls) const {
    if (this->score != ls.score)
        return this->score < ls.score;
    else if (this->SRD != ls.SRD)
        return this->SRD < ls.SRD;
    else
        return this->x < ls.x;
}

sourceInfo::sourceInfo() {
    labelNum = 0;
    matchDepth = 0;
    qid = 0;
}

sourceInfo::sourceInfo(int labelNum) {
    this->labelNum = labelNum;
}

labelPos::labelPos(int pos, int node) {
    this->pos = pos;
    this->node = node;
}

bool labelPos::operator<(const labelPos &lp) const {
    if (pos != lp.pos)
        return pos < lp.pos;
    else
        return node < lp.node;
}
