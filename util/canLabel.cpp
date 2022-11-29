//
// Created by 罗程阳 on 2022/6/17.
//

#include "canLabel.h"

canLabel::canLabel(edgeS &labelForward) {
    cid = 0;
    this->labelForward.x = labelForward.x;
    this->labelForward.w = labelForward.w;
}