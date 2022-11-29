//
// Created by 罗程阳 on 2022/7/22.
//

#include "util/ioPB.h"
#include "util/ioUint.h"
#include "util/trie.h"
#include "util/util.h"
#include <climits>
#include <time.h>
#include <malloc.h>
#include <algorithm>
#include "util/head.h"

#define MAXIMUM 0x7fffffff

// PLL
Graph g, g1, g2;
string stringFileName;
string destiGraph;
string orderfile;
string PLLfile;
string PLLpointfile;

char *txtName, *labelForwardName, *PBName;
Trie *trie;
vector<string> recs;
vector<string> queries;
int vertexNum;
vector<string> strings;
unordered_map<int, vector<int>> keyVertex;  // key id -> vertex id
unordered_map<int, set<int>> vertexNode;  // vertex id -> node id
unordered_map<int, vector<edgeS>> labelForward;
unordered_map<int, PB> vertexPB;  // vertex id and its corresponding Prefix Backward (PB)
unordered_map<int, string> nodeString;

double insertKeyTime = 0;
double insertLabelTime = 0;
double deleteKeyTime = 0;
double deleteLabelTime = 0;
double EdgeDecreaseTime = 0;
double EdgeIncreaseTime = 0;

int labelCount = 0;
int EdgeDecreaseCount = 0;
int EdgeIncreaseCount = 0;
void LoadLabelForward() {
    labelForwardName = (char *) malloc(50);
    sprintf(labelForwardName, "%s.labelForward", txtName);
    inBufS outLabel(labelForwardName);
    outLabel.start();
    edgeS tempLabel;
    int vertexNow = -1;
    while (1) {
        outLabel.nextEdge(tempLabel);
        if (tempLabel.x == -1 && tempLabel.w == -1)
            break;
        if (tempLabel.w == 0)
            vertexNow = tempLabel.x;
        labelForward[vertexNow].push_back(tempLabel);
    }
}

void LoadPB() {
    PBName = (char *) malloc(50);
    sprintf(PBName, "%s.PB", txtName);
    inBufPB inPB(PBName);
    int numPB;
    inPB.nextInt(numPB);
    for (int i = 0; i < numPB; i++) {
        if (i % 10000 == 0) {
            printf("%d/%d\n", i, numPB);
        }
        int vertexId;
        inPB.nextInt(vertexId);
        vertexPB[vertexId].vertexId = vertexId;
        inPB.getBackwardList(vertexPB[vertexId].backwardList);
        inPB.getNodeSet(vertexPB[vertexId].nodeSet);
        inPB.getBitMaps(vertexPB[vertexId].bitMaps);
        inPB.getNodeBitMapId(vertexPB[vertexId].nodeBitMapId);
    }
}

void LoadVertexNode() {
    // load string
    ifstream srcFile(stringFileName, ios::in);
    string string1;
    while (srcFile >> string1) {
        if (string1.size() > 0)
            strings.push_back(string1);
    }

    // read in the keyVertex file, construct the keyVertex map
    ifstream inputKeyVertex(string(txtName) + ".keyVertex", ios::in);
    int keyNum;
    inputKeyVertex >> keyNum;
    int tempKeyId, tempVertexNum;
    inputKeyVertex >> tempKeyId;  // skip the number of keyword occurrence
    while (keyNum--) {
        inputKeyVertex >> tempKeyId >> tempVertexNum;
        vector<int> tempVertexList;
        for(int i = 0; i < tempVertexNum; i++) {
            int tempVertexId;
            inputKeyVertex >> tempVertexId;
            tempVertexList.push_back(tempVertexId);
        }
        keyVertex[tempKeyId] = tempVertexList;
    }
    inputKeyVertex.close();

    // construct the vertexNode
    for (int i = 0; i < trie->ids.size(); i++) {
        pair<int, int> tempId = trie->ids[i];
        int nodeId = tempId.first;
        int keyId = tempId.second;
        vector<int> vertexList = keyVertex[keyId];
        for (int j = 0; j < vertexList.size(); j++) {
            int tempVertexId = vertexList[j];
            vertexNode[tempVertexId].insert(nodeId);
        }
    }
}

bool comp_id(const TrieNode* n, const int id) {
    return n->id < id;
}

bool comp_node(const TrieNode* t1, const TrieNode* t2) {
    if (t1->id < t2->id)
        return true;
    else if (t1->id > t2->id)
        return false;
    else
        return t1->last > t2->last;
}

void InsertNewKey(int vertex, string keyword) {
    vector<TrieNode*> path;
    int nodePos = trie->appendNewKey(keyword, path);  // the position of the first new node in the trie path
    int offset = keyword.size() - nodePos + 1;
    trie->buildIdx();

    bool flag = false;
    // check if vertex exist in the backward list of PB
    if (vertexNode[vertex].size() == 0) {
        flag = true;
    }
    vertexNode[vertex].insert(path.back()->id);

    // update the corresponding PB
    auto insertLabels = &labelForward[vertex];
    for (auto &item : *insertLabels) {
        auto insertPB = &vertexPB[item.x];

        // update the node id in PB
        set<int> newNodeSet;
        for (auto &item1 : insertPB->nodeSet) {
            if (item1 >= path[nodePos]->id)
                newNodeSet.insert(item1 + offset);
            else
                newNodeSet.insert(item1);
        }
        insertPB->nodeSet = newNodeSet;
//        insertPB->nodeSet.insert(path.back()->id);

        unordered_map<int, vector<int>> newNodeBitMapId;
        for (auto &item1 : insertPB->nodeBitMapId) {
            if (item1.first >= path[nodePos]->id)
                newNodeBitMapId[item1.first + offset] = item1.second;
            else
                newNodeBitMapId[item1.first] = item1.second;
        }
        insertPB->nodeBitMapId = newNodeBitMapId;

        int labelPos = 0;
        // get the location of the label
        if (flag) {
            for (auto &item1: insertPB->backwardList) {
                if (item1.w >= item.w) {
                    break;
                }
                labelPos++;
            }
            // process the first bit map
            insertPB->bitMaps[0][0].insertZero(labelPos);
            insertPB->backwardList.insert(insertPB->backwardList.begin() + labelPos, edgeS(vertex, item.w));
        } else {
            for (auto &item1: insertPB->backwardList) {
                if (item1.x == vertex) {
                    break;
                }
                labelPos++;
            }
        }
        // update the bitmap
        for (int i = 0; i < path.size(); i++) {
            TrieNode *node = path[i];
            // update the bit map of this level
            auto matchNode = insertPB->nodeSet.lower_bound(node->id);
            if (matchNode != insertPB->nodeSet.end() && *matchNode <= node->last) {
                if (i < insertPB->nodeBitMapId[*matchNode].size()) {
                    auto matchBitMap = &insertPB->bitMaps[i][insertPB->nodeBitMapId[*matchNode][i]];
                    if (matchBitMap->bitValueAt(labelPos) == 0) {
                        matchBitMap->setOne(labelPos);
                        // update the bit map of next level
                        set<int> bitMapIdModify;
                        while (matchNode != insertPB->nodeSet.end() && *matchNode <= node->last) {
                            if (i + 1 < insertPB->nodeBitMapId[*matchNode].size()) {
                                int bitMapId = insertPB->nodeBitMapId[*matchNode][i + 1];
                                if (bitMapIdModify.count(bitMapId) == 0) {
                                    matchBitMap = &insertPB->bitMaps[i + 1][insertPB->nodeBitMapId[*matchNode][i + 1]];
                                    matchBitMap->insertZero(labelPos);
                                    bitMapIdModify.insert(bitMapId);
                                }
                            }
                            matchNode++;
                        }
                    }
                    labelPos = matchBitMap->countOne(labelPos);
                }
            } else {
                int newNode = path[path.size() - 1]->id;
                insertPB->nodeSet.insert(newNode);
                node = path[i - 1];
                matchNode = insertPB->nodeSet.lower_bound(node->id);
                for (int j = 0; j < i; j++) {
                    insertPB->nodeBitMapId[newNode].push_back(insertPB->nodeBitMapId[*matchNode][j]);
                }
                auto lastBitMap = &insertPB->bitMaps[i-1][insertPB->nodeBitMapId[*matchNode][i-1]];
                BitMap newBitMap(lastBitMap->countOne());
                newBitMap.setOne(labelPos);
                while (insertPB->bitMaps.size() <= i)
                    insertPB->bitMaps.push_back({});
                insertPB->bitMaps[i].push_back(newBitMap);
                insertPB->nodeBitMapId[newNode].push_back(insertPB->bitMaps[i].size() - 1);
                break;
            }
        }
    }

//    // testPB data
//    ofstream testPBData;
//    testPBData.open("datasets/test/testPBData.txt", ios::out|ios::binary);
//    for (int i = 0; i < vertexNum; i++) {
//        testPBData << "id:" << vertexPB[i].vertexId << endl;
//        testPBData << "backward list:";
//        for (int j = 0; j < vertexPB[i].backwardList.size(); j++) {
//            testPBData << "(" << vertexPB[i].backwardList[j].x << "," << vertexPB[i].backwardList[j].w << ")";
//        }
//        testPBData << endl;
//        for (auto j = vertexPB[i].nodeSet.begin(); j != vertexPB[i].nodeSet.end(); j++) {
//            testPBData << *j << ":";
//            for (int k = 0; k < vertexPB[i].nodeBitMapId[*j].size(); k++) {
//                testPBData << vertexPB[i].bitMaps[k][vertexPB[i].nodeBitMapId[*j][k]].show() << " ";
//            }
//            testPBData << endl;
//        }
//        testPBData << endl << "---------" << endl;
//    }
//    testPBData.close();
}

void InsertKey(int vertex, string keyword) {
    timer tm;
    tm.restart();
    // get the path of the keyword
    vector<TrieNode*> path;
    path.clear();
    path.reserve(keyword.size() + 1);
    path.push_back(trie->root);
    TrieNode* activeNode = trie->root;
    for (int i = 0; i < keyword.size(); i++) {
        TrieNode* childNode = activeNode->child;
        while (childNode != NULL) {
            if (childNode->key == keyword[i]) {
                path.push_back(childNode);
                activeNode = childNode;
                break;
            } else {
                childNode = childNode->next;
            }
        }
    }
    bool flag = false;

    // check if vertex exist in the backward list of PB
    if (vertexNode[vertex].size() == 0) {
        flag = true;
    } else if (vertexNode[vertex].find(path.back()->id) != vertexNode[vertex].end()){
        return;
    }
    vertexNode[vertex].insert(path.back()->id);

    // update the corresponding PB
    auto &insertLabels = labelForward[vertex];
    for (auto &item : insertLabels) {
        auto &insertPB = vertexPB[item.x];
        if (insertPB.backwardList.empty()) {
            insertPB.vertexId = item.x;
            insertPB.backwardList.push_back(edgeS(vertex, item.w));
            BitMap tempBitMap(1);
            tempBitMap.setAll();
            insertPB.bitMaps.clear();
            insertPB.bitMaps.resize(1);
            insertPB.bitMaps[0].push_back(tempBitMap);
            insertPB.nodeSet.clear();
            insertPB.nodeSet.insert(path.back()->id);
            insertPB.nodeBitMapId[path.back()->id].push_back(0);
            insertKeyTime += tm.getTimeMs();
            continue;
        }

        int labelPos = 0;
        // get the location of the label
        if (flag) {
            for (auto &item1 : insertPB.backwardList) {
                if (item1.w >= item.w) {
                    break;
                }
                labelPos++;
            }
            // process the first bit map
            insertPB.bitMaps[0][0].insertZero(labelPos);
            insertPB.backwardList.insert(insertPB.backwardList.begin() + labelPos, edgeS(vertex, item.w));
        } else {
            for (auto &item1 : insertPB.backwardList) {
                if (item1.x == vertex) {
                    break;
                }
                labelPos++;
            }
        }

        set<int> bitMapIdModify;
        // update the bitmap
        for (int i = 0; i < path.size(); i++) {
            TrieNode *node = path[i];
            // update the bit map of this level
            auto matchNode = insertPB.nodeSet.lower_bound(node->id);
            int nextLabelPos;
            if (matchNode != insertPB.nodeSet.end() && *matchNode <= node->last) {
                if (i < insertPB.nodeBitMapId[*matchNode].size()) {
                    auto &matchBitMap = insertPB.bitMaps[i][insertPB.nodeBitMapId[*matchNode][i]];
                    if (matchBitMap.bitValueAt(labelPos) == 0) {
                        matchBitMap.setOne(labelPos);
                        // update the bit map of next level
                        bitMapIdModify.clear();
                        nextLabelPos = matchBitMap.countOne(labelPos);
                        while (matchNode != insertPB.nodeSet.end() && *matchNode <= node->last) {
                            if (i + 1 < insertPB.nodeBitMapId[*matchNode].size()) {
                                int bitMapId = insertPB.nodeBitMapId[*matchNode][i + 1];
                                if (bitMapIdModify.find(bitMapId) == bitMapIdModify.end()) {
                                    auto &bitMap = insertPB.bitMaps[i + 1][insertPB.nodeBitMapId[*matchNode][i + 1]];
                                    bitMap.insertZero(nextLabelPos);
                                    bitMapIdModify.insert(bitMapId);
                                }
                            } else {
                                if (node->id != *matchNode) {
                                    while (insertPB.bitMaps.size() <= i + 1)
                                        insertPB.bitMaps.push_back({});
                                    BitMap tempBM(matchBitMap.countOne());
                                    tempBM.setAll();
                                    tempBM.setZero(nextLabelPos);
                                    insertPB.bitMaps[i + 1].push_back(tempBM);
                                    insertPB.nodeBitMapId[*matchNode].push_back(insertPB.bitMaps[i + 1].size() - 1);
                                    bitMapIdModify.insert(insertPB.bitMaps[i + 1].size() - 1);
                                }
                            }
                            matchNode++;
                        }
                    }
                    labelPos = matchBitMap.countOne(labelPos);
                }
            } else {
                int newNode = path[path.size() - 1]->id;
                node = path[i - 1];
                matchNode = insertPB.nodeSet.lower_bound(node->id);
                insertPB.nodeSet.insert(newNode);
                int size = min(i, (int)insertPB.nodeBitMapId[*matchNode].size());
                for (int j = 0; j < size; j++) {
                    insertPB.nodeBitMapId[newNode].push_back(insertPB.nodeBitMapId[*matchNode][j]);
                }
                while (insertPB.bitMaps.size() <= i)
                    insertPB.bitMaps.push_back({});
                auto &lastBitMap = insertPB.bitMaps[size-1][insertPB.nodeBitMapId[*matchNode][size-1]];
                for (int j = size; j < i; j++) {
                    insertPB.bitMaps[j].push_back(lastBitMap);
                    insertPB.nodeBitMapId[newNode].push_back(insertPB.bitMaps[j].size() - 1);
                }
                BitMap newBitMap(lastBitMap.countOne());
                newBitMap.setOne(labelPos);
                insertPB.bitMaps[i].push_back(newBitMap);
                insertPB.nodeBitMapId[newNode].push_back(insertPB.bitMaps[i].size() - 1);
                break;
            }
        }
    }

    insertKeyTime += tm.getTimeMs();

//    // testPB data
//    ofstream testPBData;
//    testPBData.open("datasets/test/testPBData.txt", ios::out|ios::binary);
//    for (int i = 0; i < vertexNum; i++) {
//        testPBData << "id:" << vertexPB[i].vertexId << "\n";
//        testPBData << "backward list:";
//        for (int j = 0; j < vertexPB[i].backwardList.size(); j++) {
//            testPBData << "(" << vertexPB[i].backwardList[j].x << "," << vertexPB[i].backwardList[j].w << ")";
//        }
//        testPBData << "\n";
//        for (auto j = vertexPB[i].nodeSet.begin(); j != vertexPB[i].nodeSet.end(); j++) {
//            testPBData << *j << ":";
//            for (int k = 0; k < vertexPB[i].nodeBitMapId[*j].size(); k++) {
//                testPBData << vertexPB[i].bitMaps[k][vertexPB[i].nodeBitMapId[*j][k]].show() << " ";
//            }
//            testPBData << "\n";
//        }
//        testPBData << "\n" << "---------" << "\n";
//    }
//    testPBData.close();
};

void DeleteKey(int vertex, string keyword) {
    timer tm;
    tm.restart();
    vector<TrieNode*> path;
    // get the path of the keyword
    path.clear();
    path.reserve(keyword.size() + 1);
    path.push_back(trie->root);
    TrieNode* activeNode = trie->root;
    for (int i = 0; i < keyword.size(); i++) {
        TrieNode* childNode = activeNode->child;
        while (childNode != NULL) {
            if (childNode->key == keyword[i]) {
                path.push_back(childNode);
                activeNode = childNode;
                break;
            } else {
                childNode = childNode->next;
            }
        }
    }

    bool flag = false;  // whether deleting the keyword has no impact
    int pathId;
    for (pathId = 1; pathId < path.size(); pathId++) {
        auto matchNode = vertexNode[vertex].lower_bound(path[pathId]->id);
        if (matchNode != vertexNode[vertex].end() && *matchNode <= path[pathId]->last) {
            ++matchNode;
            if (matchNode != vertexNode[vertex].end() && *matchNode > path[pathId]->last) {
                flag = true;
                break;
            } else if (matchNode == vertexNode[vertex].end()) {
                flag = true;
                break;
            }
        }
    }
    vertexNode[vertex].erase(vertexNode[vertex].find(path.back()->id));

    if (!flag) {
        deleteKeyTime += tm.getTimeMs();
        return;
    }

    auto deleteLabels = &labelForward[vertex];
    for (auto &item : *deleteLabels) {
        auto deletePB = &vertexPB[item.x];
        int deleteVertex = item.x;

        // get the labelPos in the bit map of root level
        int labelPos = 0;
        for (auto &item1 : deletePB->backwardList) {
            if (item1.x == vertex) {
                break;
            }
            labelPos++;
        }

        bool flag = true;
        // get the labelPos of the path id
        for (int i = 1; i < pathId; i++) {
            auto matchNode = deletePB->nodeSet.lower_bound(path[i]->id);
            if (matchNode != deletePB->nodeSet.end() && *matchNode <= path[i]->last) {
                if (i < deletePB->nodeBitMapId[*matchNode].size()) {
                    flag = false;
                    break;
                }
                labelPos = deletePB->bitMaps[i][deletePB->nodeBitMapId[*matchNode][i]].countOne(labelPos);
            }
        }
        if (!flag)
            continue;

        // change the bit map according to the path
        for (int i = pathId; i < path.size(); i++) {
            auto matchNode = deletePB->nodeSet.lower_bound(path[i]->id);
            if (i < deletePB->nodeBitMapId[*matchNode].size()) {
                // check whether the node should be deleted
                auto &matchBitMap = deletePB->bitMaps[i][deletePB->nodeBitMapId[*matchNode][i]];
                if (i == pathId)
                    matchBitMap.setZero(labelPos);
                if (matchBitMap.allZero()) {
                    deletePB->nodeSet.erase(*matchNode);
                    deletePB->nodeBitMapId.erase(*matchNode);
                    break;
                }
                labelPos = matchBitMap.countOne(labelPos);

                set<int> deleteBitMap;
                while (matchNode != deletePB->nodeSet.end() && *matchNode <= path[i]->last) {
                    if ((i + 1) < deletePB->nodeBitMapId[*matchNode].size()) {
                        int bitMapId = deletePB->nodeBitMapId[*matchNode][i + 1];
                        if (deleteBitMap.find(bitMapId) == deleteBitMap.end()) {
                            deletePB->bitMaps[i + 1][bitMapId].deleteBit(labelPos);
                            deleteBitMap.insert(bitMapId);
                        }
                    }
                    matchNode++;
                }
            }
        }
    }
    deleteKeyTime += tm.getTimeMs();

//    // testPB data
//    ofstream testPBData;
//    testPBData.open("datasets/test/testPBData.txt", ios::out|ios::binary);
//    for (int i = 0; i < vertexNum; i++) {
//        testPBData << "id:" << vertexPB[i].vertexId << endl;
//        testPBData << "backward list:";
//        for (int j = 0; j < vertexPB[i].backwardList.size(); j++) {
//            testPBData << "(" << vertexPB[i].backwardList[j].x << "," << vertexPB[i].backwardList[j].w << ")";
//        }
//        testPBData << endl;
//        for (auto j = vertexPB[i].nodeSet.begin(); j != vertexPB[i].nodeSet.end(); j++) {
//            testPBData << *j << ":";
//            for (int k = 0; k < vertexPB[i].nodeBitMapId[*j].size(); k++) {
//                testPBData << vertexPB[i].bitMaps[k][vertexPB[i].nodeBitMapId[*j][k]].show() << " ";
//            }
//            testPBData << endl;
//        }
//        testPBData << endl << "---------" << endl;
//    }
//    testPBData.close();
}

void DeleteLabel(int vertex1, int vertex2) {  // vertex1 is the vertex id of the PB
    timer tm;
    tm.restart();
    auto deletePB = &vertexPB[vertex1];
    int labelPosRoot = 0;
    int labelPos = 0;
    bool flag = false;
    auto it = deletePB->backwardList.begin();
    for (const auto &item : deletePB->backwardList) {
        if (item.x == vertex2) {
            flag = true;
            break;
        }
        labelPos++;
        it++;
    }
    labelPosRoot = labelPos;
    if (!flag)
        return;
    deletePB->backwardList.erase(it);
    auto deleteNode = vertexNode[vertex2];
    map<int, set<int>> deleteBitMap;  // depth -> bit map ids
    vector<TrieNode*> path;
    for (auto &item : deleteNode) {
        if (deleteBitMap[0].find(0) == deleteBitMap[0].end())
            labelPos = labelPosRoot;
        else
            labelPos = labelPosRoot - 1;
        // get the path of each keyword of vertex
        TrieNode *node = trie->root;
        path.clear();
        while (node != NULL) {
            path.push_back(node);
            node = node->child;
            while (node != NULL) {
                if (node->id <= item && node->last >= item) {
                    break;
                } else {
                    node = node->next;
                }
            }
        }

        for (int i = 0; i < path.size(); i++) {
            auto matchNode = deletePB->nodeSet.lower_bound(path[i]->id);
            if (matchNode != deletePB->nodeSet.end() && i < deletePB->nodeBitMapId[*matchNode].size()) {
                int bitMapId = deletePB->nodeBitMapId[*matchNode][i];
                if (deleteBitMap[i].find(bitMapId) == deleteBitMap[i].end()) {
                    deletePB->bitMaps[i][bitMapId].deleteBit(labelPos);
                    deleteBitMap[i].insert(bitMapId);
                }
                if (deletePB->bitMaps[i][deletePB->nodeBitMapId[*matchNode][i]].allZero()) {
                    deletePB->nodeSet.erase(*matchNode);
                    deletePB->nodeBitMapId.erase(*matchNode);
                    break;
                }
                labelPos = deletePB->bitMaps[i][deletePB->nodeBitMapId[*matchNode][i]].countOne(labelPos);
            }
            while (matchNode != deletePB->nodeSet.end() && *matchNode <= path[i]->last) {
                if ((i + 1) < deletePB->nodeBitMapId[*matchNode].size()) {
                    int bitMapId = deletePB->nodeBitMapId[*matchNode][i + 1];
                    if (deleteBitMap[i + 1].find(bitMapId) == deleteBitMap[i + 1].end()) {
                        deletePB->bitMaps[i + 1][bitMapId].deleteBit(labelPos);
                        deleteBitMap[i + 1].insert(bitMapId);
                    }
                }
                matchNode++;
            }
        }
    }

    deleteLabelTime += tm.getTimeMs();
//    // testPB data
//    ofstream testPBData;
//    testPBData.open("datasets/test/testPBData.txt", ios::out|ios::binary);
//    for (int i = 0; i < vertexNum; i++) {
//        testPBData << "id:" << vertexPB[i].vertexId << endl;
//        testPBData << "backward list:";
//        for (int j = 0; j < vertexPB[i].backwardList.size(); j++) {
//            testPBData << "(" << vertexPB[i].backwardList[j].x << "," << vertexPB[i].backwardList[j].w << ")";
//        }
//        testPBData << endl;
//        for (auto j = vertexPB[i].nodeSet.begin(); j != vertexPB[i].nodeSet.end(); j++) {
//            testPBData << *j << ":";
//            for (int k = 0; k < vertexPB[i].nodeBitMapId[*j].size(); k++) {
//                testPBData << vertexPB[i].bitMaps[k][vertexPB[i].nodeBitMapId[*j][k]].show() << " ";
//            }
//            testPBData << endl;
//        }
//        testPBData << endl << "---------" << endl;
//    }
//    testPBData.close();
}

void InsertLabel(int vertex, edgeS label) {
    timer tm;
    tm.restart();
    auto &insertPB = vertexPB[vertex];
    auto &insertNode = vertexNode[label.x];
    if (insertNode.empty())
        return;
    int labelPos = 0;
    // check if the backward list is empty
    if (insertPB.backwardList.empty()) {
        insertPB.vertexId = vertex;
        insertPB.backwardList.push_back(label);
        BitMap tempBitMap(1);
        tempBitMap.setAll();
        insertPB.bitMaps.clear();
        insertPB.bitMaps.resize(1);
        insertPB.bitMaps[0].push_back(tempBitMap);
        insertPB.nodeSet.clear();
        for (auto &item : insertNode) {
            insertPB.nodeSet.insert(item);
            insertPB.nodeBitMapId[item].push_back(0);
        }
        insertLabelTime += tm.getTimeMs();
        labelCount++;
        return;
    }

    // check if vertex exist in the backward list of PB
    for (auto &item : insertPB.backwardList) {
        if (item.x == label.x)
            return;
        if (item.w >= label.w)
            break;
        labelPos++;
    }
    labelCount++;
    insertPB.backwardList.insert(insertPB.backwardList.begin() + labelPos, label);
    map<int, set<int>> insertBitMap;
    vector<TrieNode*> path;
    for (auto &item : insertNode) {
        // get the path of a keyword
        path.clear();
        TrieNode* node = trie->root;
        while (node != NULL) {
            path.push_back(node);
            node = node->child;
            while (node != NULL) {
                if (node->id <= item && node->last >= item) {
                    break;
                } else {
                    node = node->next;
                }
            }
        }

        // update the bitmap in the node of the path
        for (int i = 0; i < path.size(); i++) {
            node = path[i];
            auto matchNode = insertPB.nodeSet.lower_bound(node->id);
            if (matchNode != insertPB.nodeSet.end() && *matchNode <= node->last) {
                if (i < insertPB.nodeBitMapId[*matchNode].size()) {
                    int bitMapId = insertPB.nodeBitMapId[*matchNode][i];
                    auto &matchBitMap = insertPB.bitMaps[i][bitMapId];
                    if (insertBitMap[i].find(bitMapId) == insertBitMap[i].end()) {
                        insertPB.bitMaps[i][bitMapId].insertOne(labelPos);
                        insertBitMap[i].insert(bitMapId);
                    } else {
                        insertPB.bitMaps[i][bitMapId].setOne(labelPos);
                    }
                    labelPos = insertPB.bitMaps[i][bitMapId].countOne(labelPos);

                    // update the bitmap in the child node
                    while (matchNode != insertPB.nodeSet.end() && *matchNode <= node->last) {
                        if (i + 1 < insertPB.nodeBitMapId[*matchNode].size()) {
                            bitMapId = insertPB.nodeBitMapId[*matchNode][i + 1];
                            if (insertBitMap[i + 1].find(bitMapId) == insertBitMap[i + 1].end()) {
                                auto &matchBitMap2 = insertPB.bitMaps[i + 1][bitMapId];
                                matchBitMap2.insertZero(labelPos);
                                insertBitMap[i + 1].insert(bitMapId);
                            }
                        } else {
                            if (node->id != *matchNode) {
                                while (insertPB.bitMaps.size() <= i + 1)
                                    insertPB.bitMaps.push_back({});
                                BitMap tempBitMap(matchBitMap.countOne());
                                tempBitMap.setZero(labelPos);
                                insertPB.bitMaps[i + 1].push_back(tempBitMap);
                                insertPB.nodeBitMapId[*matchNode].push_back(insertPB.bitMaps[i + 1].size() - 1);
                                insertBitMap[i + 1].insert(insertPB.bitMaps[i + 1].size() - 1);
                            }
                        }
                        matchNode++;
                    }
                }
            } else {
                int newNode = path.back()->id;
                node = path[i - 1];
                matchNode = insertPB.nodeSet.lower_bound(node->id);
                insertPB.nodeSet.insert(newNode);
                int size = min(i, (int)insertPB.nodeBitMapId[*matchNode].size());
                for (int j = 0; j < size; j++) {
                        insertPB.nodeBitMapId[newNode].push_back(insertPB.nodeBitMapId[*matchNode][j]);
                }
                while (insertPB.bitMaps.size() <= i)
                    insertPB.bitMaps.push_back({});
                auto &lastBitMap = insertPB.bitMaps[size-1][insertPB.nodeBitMapId[*matchNode][size-1]];
                for (int j = size; j < i; j++) {
                    insertPB.bitMaps[j].push_back(lastBitMap);
                    insertPB.nodeBitMapId[newNode].push_back(insertPB.bitMaps[j].size() - 1);
                }
                BitMap newBitMap(lastBitMap.countOne());
                newBitMap.setOne(labelPos);
                insertPB.bitMaps[i].push_back(newBitMap);
                insertPB.nodeBitMapId[newNode].push_back(insertPB.bitMaps[i].size() - 1);
                break;
            }
        }
    }
    insertLabelTime += tm.getTimeMs();

//    // testPB data
//    ofstream testPBData;
//    testPBData.open("datasets/test/testPBData.txt", ios::out|ios::binary);
//    for (int i = 0; i < vertexNum; i++) {
//        testPBData << "id:" << vertexPB[i].vertexId << endl;
//        testPBData << "backward list:";
//        for (int j = 0; j < vertexPB[i].backwardList.size(); j++) {
//            testPBData << "(" << vertexPB[i].backwardList[j].x << "," << vertexPB[i].backwardList[j].w << ")";
//        }
//        testPBData << endl;
//        for (auto j = vertexPB[i].nodeSet.begin(); j != vertexPB[i].nodeSet.end(); j++) {
//            testPBData << *j << ":";
//            for (int k = 0; k < vertexPB[i].nodeBitMapId[*j].size(); k++) {
//                testPBData << vertexPB[i].bitMaps[k][vertexPB[i].nodeBitMapId[*j][k]].show() << " ";
//            }
//            testPBData << endl;
//        }
//        testPBData << endl << "---------" << endl;
//    }
//    testPBData.close();
}

void InsertKeyTest() {
    srand((unsigned)time(NULL));
    int testNum = 100;
    for (int i = 0; i < testNum; i++) {
        int stringId = rand() % strings.size();
        int vertexId = rand() % vertexNum;
        string keyword = strings[stringId];
        cout << vertexId << " " << keyword << endl;
        InsertKey(vertexId, keyword);
    }
    cout << insertKeyTime / (double) testNum;
}

void InsertLabelTest() {
    int testNum = 1000;
    srand((unsigned)time(NULL));
    for (int i = 0; i < testNum; i++) {
        int v1 = rand() % vertexNum;
        int v2 = rand() % vertexNum;
        int dis = 1 + rand() % 10;
        InsertLabel(v1, edgeS(v2 , dis));
    }
    cout << labelCount << endl;
    cout << insertLabelTime / (double)labelCount << endl;
}

void DeleteKeyTest() {
    srand((unsigned)time(NULL));
    int testNum = 1000;
    for (int i = 0; i < testNum; i++) {
       int vertexId = rand() % vertexNum;
       while (vertexNode[vertexId].size() == 0) {
           vertexId = rand() % vertexNum;
       }
       vector<int> nodes;
       nodes.assign(vertexNode[vertexId].begin(), vertexNode[vertexId].end());
       int nodeId = nodes[rand() % nodes.size()];
       string string1 = "";
       TrieNode* node = trie->root->child;
       while (node != NULL) {
           while (node->id > nodeId || node->last < nodeId)
               node = node->next;
           string1 += node->key;
           if (node->id == nodeId)
               break;
           node = node->child;
       }
       cout << vertexId << " " << string1 << " " << nodeId << endl;
       DeleteKey(vertexId, string1);
   }
   cout << deleteKeyTime / (double)testNum << endl;
}

void DeleteLabelTest() {
    srand((unsigned)time(NULL));
    int testNum = 100;
    for (int i = 0; i < testNum; i++) {
        int v1 = rand() % vertexNum;
        while (vertexPB[v1].backwardList.size() == 0)
            v1 = rand() % vertexNum;
        int v2id = rand() % vertexPB[v1].backwardList.size();
        int v2 = vertexPB[v1].backwardList[v2id].x;
        DeleteLabel(v1, v2);
        cout << v1 << " " << v2 << endl;
    }
    cout << deleteLabelTime / (double)testNum << endl;
}

void InsertNewKeyTest() {

}

double totalTime;
int totalCount;
double changeWeight;

void EdgeDecrease(int v1, int v2, int w) {
    vector<pair<int, pair<int, int>>> allChange = g.PSLdec(v1, v2, w, w - w * changeWeight);
    printf("%d\n", allChange.size());
    if (allChange.size() > changeWeight * 10000)
        return;
    EdgeDecreaseCount++;
    int vertex1, vertex2, weight;
    timer tm;
    for (int i = 0; i < allChange.size(); i++) {
        vertex1 = allChange[i].first;
        vertex2 = allChange[i].second.first;
        weight = allChange[i].second.second;
//        cout << vertex1 << " " << vertex2 << " " << weight << endl;
//        printf("%d %d\n", vertex2, vertex1);
        tm.restart();
        DeleteLabel(vertex2, vertex1);

        InsertLabel(vertex2, edgeS(vertex1, weight));
        totalTime += tm.getTimeMs();
        totalCount++;
    }
}

void EdgeIncrease(int v1, int v2, int w) {
    vector<pair<int, pair<int, int>>> allChange = g.PSLinc(v1, v2, w, w + w * changeWeight);
    printf("%d\n", allChange.size());
    if (allChange.size() > changeWeight * 10000)
        return;
    EdgeIncreaseCount++;
    int vertex1, vertex2, weight;
    timer tm;
    for (int i = 0; i < allChange.size(); i++) {
        vertex1 = allChange[i].first;
        vertex2 = allChange[i].second.first;
        weight = allChange[i].second.second;
//        cout << vertex1 << " " << vertex2 << " " << weight << endl;
//        printf("%d %d\n", vertex2, vertex1);
        tm.restart();
        DeleteLabel(vertex2, vertex1);

        InsertLabel(vertex2, edgeS(vertex1, weight));
        totalTime += tm.getTimeMs();
        totalCount++;
    }
}

void EdgeDecreaseTest() {
    srand((unsigned)time(NULL));
    int testNum = 10000;
    for (int i = 0; i < testNum; i++) {
        int v1 = rand() % vertexNum;
        auto &neighbor = g.Neighbor[v1];
        int neiId = rand() % neighbor.size();
        int v2 = neighbor[neiId].first;
        int w = neighbor[neiId].second;
        EdgeDecrease(v1, v2, w);
    }
    cout << EdgeDecreaseCount << endl;
    cout << EdgeDecreaseTime / (double)EdgeDecreaseCount << endl;
}

void EdgeIncreaseTest() {
    srand((unsigned)time(NULL));
    int testNum = 10;
    for (int i = 0; i < testNum; i++) {
        int v1 = rand() % vertexNum;
        auto &neighbor = g.Neighbor[v1];
        int neiId = rand() % neighbor.size();
        int v2 = neighbor[neiId].first;
        int w = neighbor[neiId].second;
        EdgeIncrease(v1, v2, w);
    }
    cout << EdgeIncreaseCount << endl;
    cout << EdgeIncreaseTime / (double)EdgeIncreaseCount << endl;
}

int main(int argc, char **argv) {
    string filePath = "datasets/" + string(argv[1]) + "/" + string(argv[1]);
    txtName = (char*) malloc(1+filePath.length() + 4);
    strcpy(txtName,filePath.c_str());

    // PLL
    string DesFile="datasets/";
    DesFile += argv[1];
    DesFile += "/";
    destiGraph=DesFile+argv[1];
    orderfile=DesFile+"Order";
    PLLfile=DesFile+"PLL";
    PLLpointfile=DesFile+"PLLPoint";

    g.ReadGraph(destiGraph);
    vertexNum = g.nodenum;
    g.readPLL(PLLfile, PLLpointfile);

    // rebuild trie
    stringFileName = "datasets/" + string(argv[1]) + "/" + string(argv[1]) + ".string";
    readData(stringFileName, recs);

    // construct trie index
    printf("Starting to build trie index...\n");
    trie = new Trie();
    for (auto i = 0; i < recs.size(); i++)
        trie->append(recs[i].c_str(), i);
    trie->buildIdx();
    trie->bfsIndex();
    printf("Trie index building complete.\n");

    // Load index
    printf("Loading index...\n");
    LoadLabelForward();
    LoadPB();
    // init
    LoadVertexNode();
    malloc_trim(0);

    vector<double> varList = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    vector<double> timeResult;
    timeResult.resize(varList.size());
    for (int i = 0; i < varList.size(); i++) {
        totalTime = 0;
        totalCount = 0;
        changeWeight = varList[i];
        EdgeDecreaseTest();
//        EdgeIncreaseTest();
        timeResult[i] = totalTime / (double)totalCount;
    }

    cout << "-----------" << endl;
    for (int i = 0; i < varList.size(); i++) {
//        cout << "Var" << varList[i] << " ";
//        cout << "Time" << timeResult[i] << endl;
        cout << timeResult[i] << endl;
    }
    return 0;
}