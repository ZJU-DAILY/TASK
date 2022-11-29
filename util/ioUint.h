#ifndef _IOUINT_H_
#define _IOUINT_H_

#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
#include <sys/time.h>
#include <iostream>

using namespace std;
typedef int wtype;


struct timer {
    timeval stime, etime;

    timer();
    void restart();
    double getTime();
    double getTimeMs();
};

struct edgeS {
    int x;
    wtype w;

    edgeS (int x, wtype w);
    edgeS () {}

    bool operator<(const edgeS &edgeTmp) const;
}__attribute__((packed));

struct edge {
    int x, y;
    wtype w;

    edge() {}
    edge(int xx, int yy, wtype ww);
    bool operator<(const edge &edgeTmp) const;

}__attribute__((packed));

struct labelQ {
    int x;
    wtype w;
    vector<int> source;  // which line the label comes from

    labelQ() {}
    labelQ(int x, wtype w);
    labelQ(int x, wtype w, vector<int> source);
    bool operator<(const labelQ &labelQ1) const;
};

// label in labels of canLabel
struct label {
    int x;
    wtype w;

    label() {}
    label(int x, wtype w);
    bool operator<(const label &label1) const;
    bool operator== (const label &label1) const;
};

// label in label score of canLabel
struct labelS {
    int x;
    double score;

    labelS() {}
    labelS(int x, double score);
    bool operator<(const labelS &labelS1) const;
};

struct labelSQ {
    int x;
    double score;
    vector<int> source;  // which line the label comes from

    labelSQ() {}
    labelSQ(int x, double score);
    labelSQ(int x, double score, vector<int> source);
    bool operator<(const labelSQ &labelSQ1) const;
};

struct inBuf {
    edge *buf;
    int cnt, bufLen;
    FILE *inFile;
    bool isEnd;
    int perBuf;

    inBuf();
    inBuf(char *fileName);
    void init(char *fileName, bool isMem);
    void init(char *fileName);
    void fseek(long long x);
    void start();
    ~inBuf();
    void nextEdge();
    void nextEdge(edge &tmpEdge);
    void nextEdge(edgeS &tmpEdge);
};

struct outBuf {
    edge *buf;
    int cnt;
    FILE *outFile;
    int perBuf;
    long long *sCnt;

    outBuf();
    outBuf(char *fileName);
    outBuf(char *fileName, long long &sCnt_);
    ~outBuf();
    void insert(edge &x);
    void insert(int x, int y, wtype w);
    void flush();
};

struct edgeL {
    long long x, y, w;

    bool operator<(const edgeL &edgeTmp) const {
        return x < edgeTmp.x || x == edgeTmp.x && y < edgeTmp.y || x == edgeTmp.x && y == edgeTmp.y && w < edgeTmp.w;
    }
};

struct inBufL {
    edgeL *buf;
    int cnt, bufLen, perBuf;
    FILE *inFile;
    bool isEnd;

    inBufL() {}
    inBufL(char *fileName);
    void start();
    void init(char *fileName);
    ~inBufL();
    void nextOne(edgeL &tmpEdge);
    void fseek(long long x);
    void nextEdge(edgeL &tmpEdge);
};

struct outBufL {
    edgeL *buf;
    int cnt, perBuf;
    FILE *outFile;

    outBufL() {}
    outBufL(char *fileName);
    ~outBufL();
    void insert(long long x, long long y, long long w);
    void insert(edgeL &x);
    void flush();
};

struct inBufS {
    edgeS *buf;
    int cnt, bufLen;
    FILE *inFile;
    bool isEnd;
    int perBuf;

    inBufS() {}
    inBufS(char *fileName);
    void init(char *fileName);
    void start();
    ~inBufS();
    void fseek(long long x);
    void nextEdge(edgeS &tmpEdge);
};

struct outBufS {
    edgeS *buf;
    int cnt, perBuf;
    FILE *outFile;

    outBufS() {}
    outBufS(char *fileName);
    ~outBufS();
    void insert(int x, wtype w);
    void insert(edgeS &x);
    void flush();
};


struct BTType {
    int w, v;
    long long pos;

    BTType() {}

    BTType(int w_, int v_, long long pos_) { w = w_, v = v_, pos = pos_; }
}__attribute__((packed));

struct BCType {
    int v;
    long long pos;
}__attribute__((packed));

long long checkSize(char *fileName);

void xSort(char *_sName, long long &m, bool deDu);

void xSortL(char *_sName, long long &m);

void toTxtL(char *sName);

long long checkByte(char *fileName);

int rand32();

#endif