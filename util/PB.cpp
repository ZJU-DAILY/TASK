//
// Created by 罗程阳 on 2022/5/27.
//

//#pragma GCC optimize(2)

#include "PB.h"

BitMap::BitMap() {
    size = 0;
    bitNum = 0;
}

BitMap::BitMap(int range) {
    assert(range >= 0);
    bitNum = range;
    if (range % 32 == 0)
        size = range / 32;
    else
        size = range / 32 + 1;
    bits.resize(size, 0);
}

void BitMap::setOne(int num) {
    assert(num < bitNum);
    int index = num / 32;  // the array position of the bit
    int bitIndex = num % 32;  // the bit position of in the array
    bits[index] |= 1 << bitIndex;
}

void BitMap::setZero(int num) {
    assert(num < bitNum);
    int index = num / 32;  // the array position of the bit
    int bitIndex = num % 32;  // the bit position of in the array
    bits[index] &= ~(1 << bitIndex);
}

void BitMap::setAll() {
    for (int i = 0; i < size; i++) {
        bits[i] = UINT_MAX;
    }
    bits[size - 1] = bits[size - 1] << (32 - bitNum % 32);
    bits[size - 1] = bits[size - 1] >> (32 - bitNum % 32);
}

bool BitMap::hasZero() {
    int wholeSize = bitNum / 32;
    for (int i = 0; i < wholeSize; i++) {
        if (bits[i] < UINT_MAX)
            return true;
    }
    if (size - wholeSize != 0) {
        unsigned int maxLeft = 1 << (bitNum % 32);
        maxLeft -= 1;
        if (bits[size - 1] < maxLeft)
            return true;
    }
    return false;
}

bool BitMap::allZero() {
    for (int i = 0; i < size; i++) {
        if (bits[i] > 0)
            return false;
    }
    return true;
}

int BitMap::bitValueAt(int num) {
    assert(num < bitNum);
//    int index = num / 32;
//    int bitIndex = num % 32;
//    return (bits[index] << (31 - bitIndex)) >> 31;
    if (bits[num/32] & (0x01 << (num % 32)))
        return 1;
    else
        return 0;
}

int BitMap::findOne(int num) {
    // find 1 bit start from num
    if (num >= bitNum)
        return -1;
    int index = num / 32;
    int bitIndex = num % 32;
    // search the int now
    unsigned int x = bits[index];
    x = x >> bitIndex;
    while (x) {
        if (x & 0x01) {
            return num;
        }
        num++;
        x = x >> 1;
    }

    // search the int next
    ++index;
    while (index < size) {
        x = bits[index];
        bitIndex = 0;
        while (x) {
            if (x & 0x01) {
                return 32*index + bitIndex;
            }
            bitIndex++;
            x = x >> 1;
        }
        ++index;
    }
    return -1;
}

int BitMap::nextOne(int num) {
    // find 1 bit after num
    int index = num / 32;
    int bitIndex = num % 32;
    // search the int now
    unsigned int x = bits[index];
    x = x >> bitIndex;
    while (x) {
        if (x & 0x01) {
            return num + 1;
        }
        num++;
        x = x >> 1;
    }

    // search the int next
    ++index;
    while (index < size) {
        x = bits[index];
        bitIndex = 0;
        while (x) {
            if (x & 0x01) {
                return 32*index + bitIndex + 1;
            }
            bitIndex++;
            x = x >> 1;
        }
        ++index;
    }
    return bitNum;
}

void BitMap::nextOne(int &num, int &labelPos, int &countOne) {
    if (labelPos == -1)
        return;
    if (num >= bitNum) {
        labelPos = -1;
        return;
    }
    if (countOne - 1 == labelPos) {
        labelPos = num - 1;
        return;
    }
    int index = num / 32;
    int bitIndex = num % 32;
    // search the int now
    unsigned int *p = bits.data();
    unsigned int x = *(p + index);
    x = x >> bitIndex;
    while (x) {
        if (x & 0x01) {
            ++countOne;
            if (countOne - 1 == labelPos) {
                labelPos = num;
                ++num;
                return;
            }
        }
        ++num;
        x = x >> 1;
    }

    // search the int next
    ++index;
    while (index < size) {
        x = *(p + index);
        bitIndex = 0;
        while (x) {
            if (x & 0x01) {
                ++countOne;
                if (countOne - 1 == labelPos) {
                    labelPos =  32 * index + bitIndex;
                    num = labelPos + 1;
                    return;
                }
            }
            ++bitIndex;
            x = x >> 1;
        }
        ++index;
    }
    num =  bitNum;
    labelPos = -1;
    return;
}

void BitMap::insertOne(int pos) {
    int sizeIndex = pos / 32;
    int bitIndex = pos % 32;

    unsigned int highestValue = 0;
    int x = 0x01 << 31;
    // process the unsigned int that pos in it
    unsigned int low = bits[sizeIndex];
    if (bitIndex == 0){
        low = 0;
    } else {
        low = (bits[sizeIndex] << (32 - bitIndex)) >> (32 - bitIndex);
    }
    unsigned int high = bits[sizeIndex] >> bitIndex << (bitIndex + 1);
    unsigned int insert = 0x01 << bitIndex;
    if (x & bits[sizeIndex])
        highestValue = 1;
    else
        highestValue = 0;
    bits[sizeIndex] = low + high + insert;

    // process the next unsigned int
    sizeIndex++;
    while (sizeIndex < size) {
        unsigned int preHighestValue = highestValue;
        if (x & bits[sizeIndex])
            highestValue = 1;
        else
            highestValue = 0;
        bits[sizeIndex] = bits[sizeIndex] << 1;
        bits[sizeIndex] += preHighestValue;
        sizeIndex++;
    }

    // check whether size needs to be expanded
    if (bitNum % 32 == 0) {
        ++size;
        bits.push_back(highestValue);
    }
    bitNum++;
}

void BitMap::insertZero(int pos) {
    int sizeIndex = pos / 32;
    int bitIndex = pos % 32;

    // special intersection at the end
    if (sizeIndex >= this->size) {
        ++size;
        bits.push_back(0);
        ++bitNum;
        return;
    }

    unsigned int highestValue = 0;
    unsigned int x = 0x01 << 31;
    // process the unsigned int that pos in it
    unsigned int low;
    if (bitIndex == 0){
        low = 0;
    } else {
        low = (bits[sizeIndex] << (32 - bitIndex)) >> (32 - bitIndex);
    }
//    unsigned int high = bits[sizeIndex] >> bitIndex << (bitIndex + 1);
    unsigned int high = bits[sizeIndex] >> bitIndex;
    if ((bitIndex + 1) != 32)
        high = high << (bitIndex + 1);
    else
        high = 0;
    if (x & bits[sizeIndex])
        highestValue = 1;
    else
        highestValue = 0;
    bits[sizeIndex] = low + high;

    // process the next unsigned int
    sizeIndex++;
    while (sizeIndex < size) {
        unsigned int preHighestValue = highestValue;
        if (x & bits[sizeIndex])
            highestValue = 1;
        else
            highestValue = 0;
        bits[sizeIndex] = bits[sizeIndex] << 1;
        bits[sizeIndex] += preHighestValue;
        sizeIndex++;
    }

    // check whether size needs to be expanded
    if (bitNum % 32 == 0) {
        ++size;
        bits.push_back(highestValue);
    }
    bitNum++;
}

void BitMap::deleteBit(int pos) {

    int sizeIndex = pos / 32;
    int bitIndex = pos % 32;

    unsigned int highestValue = 0;

    // process the unsigned int that pos in it
    unsigned int low;
    if (bitIndex == 0){
        low = 0;
    } else {
        low = (bits[sizeIndex] << (32 - bitIndex)) >> (32 - bitIndex);
    }
    unsigned int high = bits[sizeIndex] >> (bitIndex+1) << (bitIndex);
    if (sizeIndex + 1 < size) {
        highestValue = (0x01 & bits[sizeIndex + 1]) << 31;
    }
    bits[sizeIndex] = low + high + highestValue;

    // process the next unsigned int
    sizeIndex++;
    while (sizeIndex < size) {
        highestValue = 0;
        if (sizeIndex + 1 < size) {
            highestValue = (0x01 & bits[sizeIndex + 1]) << 31;
        }
        bits[sizeIndex] = bits[sizeIndex] >> 1;
        bits[sizeIndex] += highestValue;
        sizeIndex++;
    }

    // check whether size needs to be reduced
    if (bitNum % 32 == 1) {
        size--;
        bits.pop_back();
    }
    bitNum--;
}

void BitMap::updateTag(pair<int, int> &tag, int &pos) {
    if (tag.first - 1 == pos) {
        pos = tag.second;
        return;
    }
    int startIndex = tag.first / 32;
    int startBitIndex = tag.first % 32;
    int endIndex = pos / 32;
    int endBitIndex = pos % 32;
    unsigned int x = bits[startIndex];
    x = x >> startBitIndex;
    if (startIndex == endIndex) {
        x = x << (31 + startBitIndex - endBitIndex);
        while (x) {
            tag.second++;
            x = x & (x - 1);
        }
        tag.first = pos + 1;
        pos = tag.second - 1;
        return;
    }

    // search the first unsigned int
    while (x) {
        tag.second++;
        x = x & (x - 1);
    }

    // search the unsigned int before last
    for (int i = startIndex + 1; i <= endIndex - 1; i++) {
        x = bits[i];
        while (x) {
            tag.second++;
            x = x & (x - 1);
        }
    }

    // search the last unsigned int
    x = bits[endIndex];
    x = x << (31 - endBitIndex);
    while (x) {
        tag.second++;
        x = x & (x - 1);
    }
    tag.first = pos + 1;
    pos = tag.second - 1;
    return;
}

int BitMap::countOne() {
    int count = 0;
    for (int i = 0; i < size; i++) {
        unsigned int num = bits[i];
        while (num) {
            count++;
            num = num & (num - 1);
        }
    }

    // move right
//    for (int i = 0; i < size; i++) {
//        unsigned int num = bits[i];
//        while (num) {
//            if (num & 1)
//                count++;
//            num = num >> 1;
//        }
//    }

    // move left
//    for (int i = 0; i < size; i++) {
//        unsigned int flag = 1;
//        while (flag) {
//            if (bits[i] & flag)
//                count++;
//            flag = flag << 1;
//        }
//    }
    return count;
}

int BitMap::countOne(int pos) {
    // count one number before pos
    int count = 0;
    int index = pos / 32;
    int bitIndex = pos % 32;
    for (int i = 0; i <= index - 1; i++) {
        unsigned int x = bits[i];
        while (x) {
            if (x & 0x01) {
                count++;
            }
            x = x >> 1;
        }
    }
    unsigned int x = bits[index];
    for (int i = 0; i < bitIndex; i++) {
        if (x & 0x01) {
            count++;
        }
        x = x >> 1;
    }
    return count;
}


void BitMap::compress(BitMap &bitMap) {
    assert(bitNum == bitMap.bitNum);

    // init the new bit map info
    int newBitNum = bitMap.countOne();
    int newSize;
    if (newBitNum % 32 == 0)
        newSize = newBitNum / 32;
    else
        newSize = newBitNum / 32 + 1;
    vector<unsigned int> newBits(newSize, 0);

    // calculate each bit for new bit map
    int newIndex = 0;
    for (int i = 0; i < bitNum; i++) {
        if (bitMap.bitValueAt(i)) {
            if (bitValueAt(i))
                newBits[newIndex / 32] |= 1 << (newIndex % 32);  // set newIndex bit to 1
            else
                newBits[newIndex / 32] &= ~(1 << (newIndex % 32));  // set newIndex bit to 0
            newIndex++;
        }
    }
    bits = newBits;
    bitNum = newBitNum;
    size = newSize;
    vector<unsigned int>().swap(newBits);
}

void BitMap::change(BitMap &bitMap) {
    assert(this->bitNum >= bitMap.bitNum);

    int pos = 0;  // the index of &bitMap
    for (int i = 0; i < this->bitNum; i++) {
        if (this->bitValueAt(i) == 1) {
            if (bitMap.bitValueAt(pos) == 0)
                this->setZero(i);
            pos++;
        }
    }
}

void BitMap::performOr(BitMap &bitMap) {
    assert(this->bitNum == bitMap.bitNum);

    for (int i = 0; i < size; i++) {
        this->bits[i] |= bitMap.bits[i];
    }
}

void BitMap::performAnd(BitMap &bitMap) {
    assert(this->bitNum == bitMap.bitNum);

    for (int i = 0; i < size; i++) {
        this->bits[i] &= bitMap.bits[i];
    }
}

void BitMap::uncompress(BitMap &bitMap) {
    int index1 = 0;
    for (int i = 0; i < bitNum; i++) {
        if (this->bitValueAt(i) == 1) {
            if (bitMap.bitValueAt(index1) == 0) {
                this->setZero(i);
            }
            index1++;
        }
    }
}
string BitMap::show() {
    string result = "";
    for (int i = 0; i < bitNum; i++) {
        result += to_string(bitValueAt(i));
    }
//    cout << result << endl;
    return result;
}

//bool BitMap::isEqual(BitMap &bitMap) {
//    if (size != bitMap.size)
//        return false;
//    if (bitNum != bitMap.bitNum)
//        return false;
//    for (int i = 0; i < size; i++) {
//        if (bits[i] != bitMap.bits[i])
//            return false;
//    }
//    return true;
//}
//
//void BitMap::assignValue(BitMap &bitMap) {
//    size = bitMap.size;
//    bitNum = bitMap.bitNum;
//    bits = bitMap.bits;
//}

BitMap &BitMap::operator=(const BitMap &bm) {
    size = bm.size;
    bitNum = bm.bitNum;
    bits = bm.bits;
    return *this;
}

bool BitMap::operator==(const BitMap &bm) {
    if (size != bm.size)
        return false;
    if (bitNum != bm.bitNum)
        return false;
    for (int i = 0; i < size; i++) {
        if (bits[i] != bm.bits[i])
            return false;
    }
    return true;
}

bool BitMap::operator<(const BitMap &bm) const {
    if (this->size != bm.size)
        return this->size < bm.size;
    else if (this->bitNum != bm.bitNum)
        return this->bitNum < bm.bitNum;
    else {
        for (int i = 0; i < size; i++) {
            if (bits[i] != bm.bits[i]) {
                return bits[i] < bm.bits[i];
            }
        }
        return false;
    }
}

BitMap::~BitMap() {
    bits.clear();
    vector<unsigned int>().swap(bits);
}