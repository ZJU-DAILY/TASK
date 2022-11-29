//
// Created by 罗程阳 on 2022/9/26.
//

#include "util/ioUint.h"
#include "util/util.h"
#include "util/trie.h"
#include "util/PB.h"
#include "util/ioPB.h"
#include <stdio.h>
#include <climits>
#include <omp.h>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <set>
#include <queue>
#include <malloc.h>

using namespace std;

char *txtName, *PBName;
string outName;
char *labelForwardName;  // the name of file storing forward search label
int *forwardLabelID; // this array stores the vertex id list of label forward file
int vertexNum, keyNum;
long long labelBackwardNum = 0; // the number of backward labels
vector<edge> labelBackward;

Trie *trie;
vector<string> recs;
unordered_map<int, vector<int>> keyVertex;  // key id -> vertex id
unordered_map<int, vector<int>> vertexNode;  // vertex id -> node id
vector<PB> vertexPB;  // vertex id and its corresponding Prefix Backward (PB)
vector<int> orderVertex;


bool CompareLabelForward(edgeS &label1, edgeS &label2) {
    return label1.w < label2.w;
}

bool Compare(edge &label1, edge &label2) {
    return label1.x < label2.x || label1.x == label2.x && label1.w < label2.w;
}

struct CompareQueue {
    bool operator () (const labelQ &label1, const labelQ &label2){
        if (label1.w == label2.w)
            return label1.x > label2.x;
        else
            return label1.w > label2.w;
    }
};

void LoadIndex() {
    // read in the keyVertex file, construct the keyVertex map
    ifstream inputKeyVertex(string(txtName) + ".keyVertex", ios::in);
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
            if (vertexNode.count(tempVertexId) == 0) {
                vector<int> tempVertexNode;
                tempVertexNode.push_back(nodeId);
                vertexNode[tempVertexId] = tempVertexNode;
            } else {
                vertexNode[tempVertexId].push_back(nodeId);
            }
        }
    }

    // read in the 2-hop label and vertex number
    ifstream outFile(outName, ios::binary | ios::out);
    outFile.read((char *)(&vertexNum), sizeof(vertexNum));

    // init the PB
    vertexPB.resize(vertexNum);

    edgeS tempLabel;
    int vid, dis;
    int labelNum;
    for (int i = 0; i < vertexNum; i++) {
        bool hasKey = (vertexNode.count(i) == 0)?false:true;
        outFile.read((char *)(&labelNum), sizeof(labelNum));
        for (int j = 0; j < labelNum; j++) {
                outFile.read((char *) (&vid), sizeof(vid));
                outFile.read((char *) (&dis), sizeof(dis));
                tempLabel.x = vid;
                tempLabel.w = dis;
                if (hasKey) {  // whether vertex has keywords(trie node)
                    vertexPB[tempLabel.x].backwardList.push_back(edgeS(i, dis));
                }
        }
    }
    outFile.close();
}

bool ComparePB(const edgeS &l1, const edgeS &l2) {
    if (l1.w < l2.w)
        return true;
    else if (l1.w == l2.w && l1.x < l2.x)
        return true;
    else
        return false;
}

void BuildPB() {
    // init
    for (int i = 0; i < vertexNum; i++) {
        vertexPB[i].vertexId = i;
        vector<edgeS>(vertexPB[i].backwardList).swap(vertexPB[i].backwardList);
        sort(vertexPB[i].backwardList.begin(), vertexPB[i].backwardList.end(), ComparePB);
    }

    unordered_map<int, vector<int>> nodePath;  // the trie node id from root node to a certain node
    unordered_map<int, int> nodeUp; // the upper limit of a node
    for (int i = 0; i < vertexNum; i++) {
        if (vertexNode.count(i) != 0) {
            vector<int> tempNodeIds = vertexNode[i];
            for (int j = 0; j < tempNodeIds.size(); j++) {
                int tempNodeId = tempNodeIds[j];
                if (nodePath.count(tempNodeId) == 0) {
                    TrieNode *node = trie->root;
                    while (node->id != tempNodeId) {
                        nodePath[tempNodeId].push_back(node->id);
                        TrieNode *childNode = node->child;
                        while (childNode != NULL) {
                            if (childNode->id <= tempNodeId && childNode->last >= tempNodeId) {
                                break;
                            } else {
                                childNode = childNode->next;
                            }
                        }
                        node = childNode;
                    }
                    nodePath[tempNodeId].push_back(node->id);
                    nodeUp[tempNodeId] = node->last;
                }
            }
        }
    }

    time_t time1, time2;
    time1 = time(NULL);

    clock_t start,finish;
    start = clock();

    // store the bitMap before compression
    vector<int> indexSize;
    indexSize.resize(vertexNum);
    int threadNum = 4;
    omp_set_num_threads(threadNum);
#pragma omp parallel for
    // calculate bitmaps for each PB
    for (int i = 0; i < vertexNum; i++) {
        auto PBNow = &vertexPB[i];
        if (i % 100000 == 0) {
            malloc_trim(0);
            printf("%d/%d\n", i, vertexNum);
        }
        unordered_map<int, BitMap> nodeBitMap;  // trie node -> BitMap
        auto tempBackwardList = PBNow->backwardList.data();
        int tempLabelNum = PBNow->backwardList.size();  // the number label in backward list of certain vertex PB
        set<int> nodeSet;
        for (int j = 0; j < tempLabelNum; j++) {
            int tempVertex = tempBackwardList[j].x;
            auto tempNodeIds = vertexNode[tempVertex].data();
            int tempNodeIdsSize = vertexNode[tempVertex].size();
            for (int k = 0; k < tempNodeIdsSize; k++) {  // traverse the trie and generate the bit map
                int tempNodeId = tempNodeIds[k];
                for (int l = 0; l < nodePath[tempNodeId].size(); l++) {
                    if (nodeBitMap.count(nodePath[tempNodeId][l]) == 0) {
                        BitMap bitMap(tempLabelNum);
                        nodeBitMap[nodePath[tempNodeId][l]] = bitMap;
                    }
//                    s1 = clock();
                    nodeBitMap[nodePath[tempNodeId][l]].setOne(j);
//                    f1 = clock();
//                    time1 += ((double)(f1-s1) / CLOCKS_PER_SEC);
                }
                nodeSet.insert(tempNodeId);
            }
        }

        // calculate the size before compression
        indexSize[i] += 4;  //
        indexSize[i] += PBNow->backwardList.size() * 8;
        indexSize[i] += 4;
        indexSize[i] += nodeSet.size() * 4;
        indexSize[i] += 4;
        indexSize[i] += nodeBitMap.size() * 4;
        for (auto &item : nodeBitMap) {
            indexSize[i] += item.second.size * 4;
            indexSize[i] += 4; // bitNum
            indexSize[i] += 4; // size
        }

        // delete the unnecessary node
        if (nodeSet.size() != 1 && nodeSet.size() != 0) {
            set<int> tempNodeSet;
            auto it1 = nodeSet.begin();
            auto it2 = nodeSet.begin();
            it2++;
            for (it1, it2; it2 != nodeSet.end(); it1++, it2++) {
                if (nodeUp[*it1] < *it2) {
                    tempNodeSet.insert(*it1);
                }
            }
            tempNodeSet.insert(*it1);
            nodeSet = tempNodeSet;
        }


        // compress the bitmap
        set<int> accessedNodes;
        PBNow->bitMaps.resize(1);
        PBNow->bitMaps[0].push_back(nodeBitMap[1]);
        for (auto &node : nodeSet)  {
            PBNow->nodeSet.insert(node);
            int count = 0;
            for (int k = nodePath[node].size() - 1; k > 0; k--) {
                int nodeNow = nodePath[node][k];
                int nodePre = nodePath[node][k-1];
                auto bitMapNow = nodeBitMap[nodeNow];
                auto bitMapPre = &nodeBitMap[nodePre];
                if (count == 0 && bitMapNow == *bitMapPre && k > 1)  // redundant bitmaps at the tail should be deleted
                    continue;
                count++;
                if (PBNow->nodeBitMapId[node].size() == 0)
                    PBNow->nodeBitMapId[node].resize(k + 1);
                if (accessedNodes.find(nodeNow) == accessedNodes.end()) {
                    bitMapNow.compress(*bitMapPre);
                    while (PBNow->bitMaps.size() <= k)
                        PBNow->bitMaps.push_back({});
                    PBNow->bitMaps[k].push_back(bitMapNow);
                    PBNow->nodeBitMapId[node][k] = PBNow->bitMaps[k].size() - 1;
                    accessedNodes.insert(nodeNow);
                } else {
                    PBNow->nodeBitMapId[node][k] = PBNow->bitMaps[k].size() - 1;
                }
            }
            if (PBNow->nodeBitMapId[node].empty())
                PBNow->nodeBitMapId[node].push_back({});
            PBNow->nodeBitMapId[node][0] = 0;
        }
    }
    finish = clock();
    time2 = time(NULL);
    cout << "Thread Number: " << threadNum << endl;
    cout << "CPU Time: " << ((double)(finish-start) / CLOCKS_PER_SEC) << " (s) "<< endl;
    cout << "Real Time: " << time2 - time1 << " (s)" << endl;

    long long totalIndexSize = 0;
    for (int i = 0; i < vertexNum; i++) {
        totalIndexSize += indexSize[i];
    }
    cout << "Size Before Compression:" << totalIndexSize / 1024.0 / 1024.0 / 1024.0 << "GB" << endl;
////     testPB data
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

void SavePB() {
    outBufPB outPB(PBName);
    outPB.numPB = vertexPB.size();
    outPB.insertInt(outPB.numPB);
    long long comSize = 0;
    comSize += 4;
    for (auto &item : vertexPB) {
        outPB.insertInt(item.vertexId);  // write the vertex id of PB
        int backwardListSize = item.backwardList.size();
        outPB.insertInt(backwardListSize);  // write the backward list size
        outPB.insertBackwardList(item.backwardList);  // write the backward list
        int nodeSetSize = item.nodeSet.size();
        outPB.insertInt(nodeSetSize);  // write the node set size;
        outPB.insertNodeSet(item.nodeSet);  // write the node set;
        outPB.insertBitMaps(item.bitMaps);
        outPB.insertNodeBitMapId(item.nodeBitMapId);

        comSize += 4;
        comSize += 4;
        comSize += 8 * item.backwardList.size();
        comSize += 4;
        comSize += 4 * nodeSetSize;
        comSize += 4;
        for (auto &item1 : item.bitMaps) {
            comSize += 4;
            for (int i = 0; i < item1.size(); i++) {
                comSize += item1[i].size * 4;
                comSize += 4; // bitNum
                comSize += 4; // size
            }
        }
        comSize += 4;
        for (auto &item1 : item.nodeBitMapId) {
            comSize += 4;
            comSize += 4;
            comSize += 4 * item1.second.size();
        }
    }
//    cout << "Size After Compression:" << (double)comSize << "B" << endl;
    cout << "Size After Compression:" << (double)comSize / 1024.0 / 1024.0 / 1024.0 << "GB" << endl;
}

int main(int argc, char **argv) {
    string filePath = "datasets/" + string(argv[1]) + "/" + string(argv[1]);
    txtName = (char*) malloc(1+filePath.length() + 4);
    strcpy(txtName,filePath.c_str());

    outName = "datasets/" + string(argv[1]) + "/" + string(argv[1]) + ".out";
    PBName = (char *) malloc(1 + strlen(txtName) + 50);
    sprintf(PBName, "%s.PB", txtName);

    // construct trie index
    printf("Starting to construct trie index...\n");
    string stringFileName = "datasets/" + string(argv[1]) + "/" + string(argv[1]) + ".string";
    readData(stringFileName, recs);
    trie = new Trie();
    for (auto i = 0; i < recs.size(); i++)
        trie->append(recs[i].c_str(), i);
    trie->buildIdx();
    trie->bfsIndex();
    printf("Trie index constructing complete.\n\n");

    // load 2-hop label index
    printf("Starting loading the 2-hop label index...\n");
    LoadIndex();
    printf("2-hop loading complete.\n\n");


    printf("Starting to construct PB index...\n");
    BuildPB();
    printf("PB constructing complete.\n\n");

    printf("Saving PB.\n\n");
    SavePB();

    free(labelForwardName);
    free(PBName);
    return 0;
}