//
// Created by 罗程阳 on 2022/5/31.
//
//#pragma GCC optimize(2)

#include "ioUint.h"
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <vector>
#include <map>
#include <sys/time.h>
#include <iostream>

const long long memSize = 4LL * 1024 * 1024 * 1024;
const int bufSize = 1024 * 3 * 256 * 32;//2 * 12 ;
//const int intPerMem = memSize / 2 / sizeof(int);
const int llPerBuf = bufSize / sizeof(long long);
const int fileNameLen = 128;

long long iCnt = 0, oCnt = 0;
double tRead = 0, tWrite = 0;


int iter = 0, allCnt, nxCnt, nyCnt, ilCnt, ihCnt, olCnt, ohCnt;
long long pCnt = 0, c1 = 0, c2 = 0, c3 = 0, c4 = 0, c22 = 0, lenCnt = 0, largeDisk = 0, sumCnt = 0;


double tBi, tRe, tInit, tGen, tRun, tPrune, tMerge1, tMerge2, tGenS = 0, tPruneS = 0, tMerge1S = 0, tMerge2S = 0, tSortS = 0, tSort1 = 0, tSort2 = 0, tQuery, tPruneSort = 0, tPruneCPU = 0, tSortMem = 0;
double tBP = 0;

int ixBuf, oxBuf, iyBuf, oyBuf, ixMerge, oxMerge, iyMerge, oyMerge,
        iSort = 0, oSort = 0, iSort1 = 0, oSort1 = 0, iSort2 = 0, oSort2 = 0;


bool fTest = 0, fBP = 0, fxBuf = 0, fyBuf = 0, fmerge = 0, fmerge2 = 0, fdeb = 1, fiter = 0, fnon = 0, findex = 0, fnew = 0, ftime_ = 1, fio = 0, fprune = 0, fEq = 0, fCheck = 0;

const int BSize = 4096;
const int BCSize = 4096 / sizeof(BCType);

//const int lenSize = 5, lowBit = 65535;
//const int lenPrune[lenSize + 1] = {0, 2, 3, 4, 6, 9};

char *p1, *p2;


// timer
timer::timer() {
    gettimeofday(&stime, NULL);
}

double timer::getTime() {
    gettimeofday(&etime, NULL);
    return (double) etime.tv_sec - stime.tv_sec +
           ((double) etime.tv_usec - stime.tv_usec) / 1000000.0;
}

double timer::getTimeMs() {
    gettimeofday(&etime, NULL);
    return ((double) etime.tv_sec - stime.tv_sec +
           ((double) etime.tv_usec - stime.tv_usec) / 1000000.0) * 1000.0;
}

void timer::restart() {
    gettimeofday(&stime, NULL);
}


// edgeS
edgeS::edgeS(int x, wtype w) {
    this->x = x;
    this->w = w;
}

bool edgeS::operator<(const edgeS &edgeTmp) const {
    return x < edgeTmp.x || x == edgeTmp.x && abs(w) < abs(edgeTmp.w) ||
           x == edgeTmp.x && abs(w) == abs(edgeTmp.w) && w > edgeTmp.w;
}


// edge
edge::edge(int xx, int yy, wtype ww) {
    x = xx;
    y = yy;
    w = ww;
}

bool edge::operator<(const edge &edgeTmp) const {
    return x < edgeTmp.x || x == edgeTmp.x && y < edgeTmp.y ||
           x == edgeTmp.x && y == edgeTmp.y && abs(w) < abs(edgeTmp.w);
}

// labelQ
labelQ::labelQ(int x, wtype w) {
    this->x = x;
    this->w = w;
}

labelQ::labelQ(int x, wtype w, vector<int> source) {
    this->x = x;
    this->w = w;
    this->source = source;
}

bool labelQ::operator<(const labelQ &labelQ1) const {
    if (this->w != labelQ1.w)
        return this->w < labelQ1.w;
    else
        return this->x < labelQ1.x;
}

// label
label::label(int x, wtype w) {
    this->x = x;
    this->w = w;
}

bool label::operator<(const label &label1) const {
    if (this->w != label1.w)
        return this->w < label1.w;
    else
        return this->x < label1.x;
}

bool label::operator==(const label &label1) const {
    if (this->x != label1.x)
        return false;
    if (this->w != label1.w)
        return false;
    return true;
}


// labelS
labelS::labelS(int x, double score) {
    this->x = x;
    this->score = score;
}

bool labelS::operator<(const labelS &labelS1) const {
    if (this->score != labelS1.score)
        return this->score < labelS1.score;
    else
        return this->x < labelS1.x;
}

// labelSQ
labelSQ::labelSQ(int x, double score) {
    this->x = x;
    this->score = score;
}

labelSQ::labelSQ(int x, double score, vector<int> source) {
    this->x = x;
    this->score = score;
    this->source = source;
}

bool labelSQ::operator<(const labelSQ &labelSQ1) const {
    if (this->score != labelSQ1.score)
        return this->score < labelSQ1.score;
    else
        return this->x < labelSQ1.x;
}

// inBuf

inBuf::inBuf() {
    perBuf = bufSize / sizeof(edge);
    inFile = NULL;
}

inBuf::inBuf(char *fileName) {
    perBuf = bufSize / sizeof(edge);
    init(fileName, false);
}

void inBuf::init(char *fileName, bool isMem) {
    if (isMem) perBuf = memSize / 5 / sizeof(edge);
    buf = (edge *) malloc(sizeof(edge) * perBuf);
    inFile = fopen(fileName, "rb");
}

void inBuf::init(char *fileName) {
    buf = (edge *) malloc(sizeof(edge) * perBuf);
    inFile = fopen(fileName, "rb");
}

void inBuf::fseek(long long x) {
    fseeko(inFile, x * sizeof(edge), SEEK_SET);
}

void inBuf::start() {
    timer tm;
    cnt = 0;
    iCnt++;
    bufLen = fread(buf, sizeof(edge), perBuf, inFile);
    isEnd = (bufLen == 0);
    tRead += tm.getTime();
}

void inBuf::nextEdge() {
    if (cnt < bufLen - 1) cnt++;
    else {
        timer tm;
        iCnt++;
        bufLen = fread(buf, sizeof(edge), perBuf, inFile);
        cnt = 0;
        isEnd = (bufLen == 0);
        tRead += tm.getTime();
    }
}

void inBuf::nextEdge(edge &tmpEdge) {
    if (cnt < bufLen - 1) {
        tmpEdge = buf[cnt++];
    } else {
        timer tm;
        iCnt++;
        tmpEdge = buf[cnt];
        bufLen = fread(buf, sizeof(edge), perBuf, inFile);
        cnt = 0;
        isEnd = (bufLen == 0);
        tRead += tm.getTime();
    }
}

void inBuf::nextEdge(edgeS &tmpEdge) {
    if (cnt < bufLen - 1) {
        tmpEdge.x = buf[cnt].y;
        tmpEdge.w = buf[cnt++].w;
    } else {
        timer tm;
        iCnt++;
        tmpEdge.x = buf[cnt].y;
        tmpEdge.w = buf[cnt].w;
        bufLen = fread(buf, sizeof(edge), perBuf, inFile);
        cnt = 0;
        isEnd = (bufLen == 0);
        tRead += tm.getTime();
    }
}

inBuf::~inBuf() {
    if (inFile != NULL) fclose(inFile);
    free(buf);
}


// outBuf

outBuf::outBuf()  {
    perBuf = bufSize / sizeof(edge);
    buf = (edge *) malloc(perBuf * sizeof(edge));
    cnt = 0;
    outFile = NULL;
}

outBuf::outBuf(char *fileName) {
    perBuf = bufSize / sizeof(edge);
    buf = (edge *) malloc(perBuf * sizeof(edge));
    cnt = 0;
    outFile = fopen(fileName, "wb");
}

outBuf::outBuf(char *fileName, long long &sCnt_) {
    sCnt = &sCnt_;
    perBuf = memSize / 5 / sizeof(edge);
    buf = (edge *) malloc(perBuf * sizeof(edge));
    cnt = 0;
    outFile = fopen(fileName, "wb");
}

outBuf::~outBuf() {
    flush();
    fflush(outFile);
    if (outFile != NULL) fclose(outFile);
    free(buf);
}

void outBuf::insert(edge &x) {
    buf[cnt++] = x;
    if (cnt == perBuf) flush();
}

void outBuf::insert(int x, int y, wtype w) {
    buf[cnt].x = x;
    buf[cnt].y = y;
    buf[cnt++].w = w;
    if (cnt == perBuf) flush();
}

void outBuf::flush() {
    if (perBuf != bufSize / sizeof(edge) && cnt > 0) {
        sort(buf, buf + cnt);
        int newCnt = 1;
        for (int px = buf[0].x, py = buf[0].y, i = 1; i < cnt; i++) {
            if (px == buf[i].x && py == buf[i].y) continue;
            buf[newCnt++] = buf[i];
            px = buf[i].x;
            py = buf[i].y;
        }
        *sCnt -= cnt - newCnt;
        cnt = newCnt;
    }

    timer tm;
    oCnt++;
    fwrite(buf, sizeof(edge), cnt, outFile);
    cnt = 0;
    tWrite += tm.getTime();
}


// inBufL
inBufL::inBufL(char *fileName) {
    perBuf = bufSize / sizeof(edgeL);
    buf = (edgeL *) malloc(sizeof(edgeL) * perBuf);
    inFile = fopen(fileName, "rb");
}

void inBufL::start() {
    timer tm;
    cnt = 0;
    iCnt++;
    bufLen = fread(buf, sizeof(edgeL), perBuf, inFile);
    isEnd = (bufLen == 0);
    tRead += tm.getTime();
}

void inBufL::init(char *fileName) {
    inFile = fopen(fileName, "rb");
}

inBufL::~inBufL() {
    if (inFile != NULL) fclose(inFile);
    free(buf);
}

void inBufL::nextOne(edgeL &tmpEdge) {
    fread(&tmpEdge, sizeof(edgeL), 1, inFile);
}

void inBufL::fseek(long long x) {
    fseeko(inFile, x * sizeof(edgeL), SEEK_SET);
}

void inBufL::nextEdge(edgeL &tmpEdge) {
    if (cnt < bufLen - 1) {
        tmpEdge = buf[cnt++];
    } else {
        timer tm;
        iCnt++;
        tmpEdge = buf[cnt];
        bufLen = fread(buf, sizeof(edgeL), perBuf, inFile);
        cnt = 0;
        isEnd = (bufLen == 0);
        tRead += tm.getTime();
    }
}

outBufL::outBufL(char *fileName) {
    perBuf = bufSize / sizeof(edgeL);
    buf = (edgeL *) malloc(sizeof(edgeL) * perBuf);
    cnt = 0;
    outFile = fopen(fileName, "wb");
}

outBufL::~outBufL() noexcept {
    flush();
    fflush(outFile);
    if (outFile != NULL) fclose(outFile);
    free(buf);
}

void outBufL::insert(long long x, long long y, long long w) {
    buf[cnt].x = x;
    buf[cnt].y = y;
    buf[cnt++].w = w;
    if (cnt == perBuf) flush();
}

void outBufL::insert(edgeL &x) {
    buf[cnt++] = x;
    if (cnt == perBuf) flush();
}

void outBufL::flush() {
    timer tm;
    oCnt++;
    fwrite(buf, sizeof(edgeL), cnt, outFile);
    cnt = 0;
    tWrite += tm.getTime();
}

inBufS::inBufS(char *fileName) {
    perBuf = bufSize / sizeof(edgeS);
    buf = (edgeS *) malloc(sizeof(edgeS) * perBuf);
    inFile = fopen(fileName, "rb");
}

void inBufS::init(char *fileName) {
    buf = (edgeS *) malloc(sizeof(edgeS) * perBuf);
    inFile = fopen(fileName, "rb");
}

void inBufS::start() {
    timer tm;
    cnt = 0;
    iCnt++;
    bufLen = fread(buf, sizeof(edgeS), perBuf, inFile);  //将file的数据读进buf中
    isEnd = (bufLen == 0);
    tRead += tm.getTime();
}

inBufS::~inBufS() {
    if (inFile != NULL) fclose(inFile);
    free(buf);
}

void inBufS::fseek(long long x) {
    fseeko(inFile, x * sizeof(edgeS), SEEK_SET);
}

void inBufS::nextEdge(edgeS &tmpEdge) {
    if (cnt < bufLen - 1) {
        tmpEdge = buf[cnt++];
    } else {
        timer tm;
        iCnt++;
        tmpEdge = buf[cnt];
        bufLen = fread(buf, sizeof(edgeS), perBuf, inFile);
        cnt = 0;
        isEnd = (bufLen == 0);
        tRead += tm.getTime();
    }
}

outBufS::outBufS(char *fileName) {
    perBuf = bufSize / sizeof(edgeS);
    buf = (edgeS *) malloc(sizeof(edgeS) * perBuf);
    cnt = 0;
    outFile = fopen(fileName, "wb");
}

outBufS::~outBufS() {
    flush();
    fflush(outFile);
    if (outFile != NULL) fclose(outFile);
    free(buf);
}

void outBufS::insert(int x, wtype w) {
    buf[cnt].x = x;
    buf[cnt++].w = w;
    if (cnt == perBuf) flush();
}

void outBufS::insert(edgeS &x) {
    buf[cnt++] = x;
    if (cnt == perBuf) flush();
}

void outBufS::flush() {
    timer tm;
    oCnt++;
    fwrite(buf, sizeof(edgeS), cnt, outFile);
    cnt = 0;
    tWrite += tm.getTime();
}

long long checkSize(char *fileName) {
    FILE *pFile = fopen(fileName, "rb");
    fseeko(pFile, 0, SEEK_END);
    long long ans = ftello(pFile) / sizeof(edge);
    fclose(pFile);
    return ans;
}

void xSort(char *_sName, long long &m, bool deDu) {

    if (m == 0) return;
    timer tm, tm1;
    tm.restart();
    char *sName = (char *) malloc(1 + strlen(_sName));
    strcpy(sName, _sName);
    char *tmpName = (char *) malloc(1 + strlen(sName) + 4);
    strcpy(tmpName, sName);
    strcat(tmpName, ".tmp");

    int cnt, edgePerMem = memSize / sizeof(edge), edgePerBuf = bufSize / sizeof(edge);

    edge *memEdge = (edge *) malloc(edgePerMem * (long long) sizeof(edge));

    double tin = 0;
    if (deDu) {
        long long tmpM = m;
        outBuf tmpBuf(tmpName);
        inBuf usEdge(sName);
        usEdge.start();
        m = 0;
        while (!usEdge.isEnd) {
            cnt = 0;
            memEdge[0].x = 0;
            while (cnt < edgePerMem && !usEdge.isEnd)
                usEdge.nextEdge(memEdge[cnt++]);
            timer tTmp;
            sort(memEdge, memEdge + cnt);
            tin += tTmp.getTime();
            int px = -1, py = -1;
            for (int i = 0; i < cnt; i++) {
                if (memEdge[i].x == px && memEdge[i].y == py) continue;
                m++;
                tmpBuf.insert(memEdge[i]);
                px = memEdge[i].x;
                py = memEdge[i].y;
            }
        }

        swap(tmpName, sName);


    }


    {
        outBuf tmpBuf(tmpName);
        inBuf usEdge(sName);
        usEdge.start();
        while (!usEdge.isEnd) {
            cnt = 0;
            while (cnt < edgePerMem && !usEdge.isEnd)
                usEdge.nextEdge(memEdge[cnt++]);

            timer tTmp;
            sort(memEdge, memEdge + cnt);
            tin += tTmp.getTime();
            for (int i = 0; i < cnt; i++)
                tmpBuf.insert(memEdge[i]);
        }
    }

    tSort1 += tm1.getTime();

    int iCntTmp, oCntTmp;
    iSort1 += (iCntTmp = iCnt);
    oSort1 += (oCntTmp = oCnt);


    if (m <= edgePerMem) {
        if (1 + strlen(tmpName) > 1 + strlen(sName)) {
            remove(sName);
            rename(tmpName, sName);
        } else {
            remove(sName);
        }
        free(sName);
        free(tmpName);
        free(memEdge);

        tSortS += tm.getTime();

        return;
    }


    timer tm2;


    int h[edgePerMem / edgePerBuf][2], hh, cntBuf[edgePerMem / edgePerBuf];

    long long posFile[edgePerMem / edgePerBuf];
    int ii = 0;
    for (long long len = edgePerMem; len < m; len *= edgePerMem / edgePerBuf) {
        swap(sName, tmpName);
        inBuf sBuf(sName);
        outBuf tmpBuf(tmpName);

//show(len);

        for (long long lastPos = 0; lastPos < m; lastPos += len * (edgePerMem / edgePerBuf)) {
            cnt = 0;
            int maxi;
            for (int i = 0; i < edgePerMem / edgePerBuf && lastPos + i * len < m; i++) {
                maxi = i;

//printf("%d\n", (lastPos+i*len));

                fseeko(sBuf.inFile, (lastPos + i * len) * sizeof(edge), SEEK_SET);
                posFile[i] = 0;
                sBuf.start();
                h[hh = i + 1][0] = cnt;
                h[i + 1][1] = i;
                while (!sBuf.isEnd && cnt - h[hh][0] < edgePerBuf)
                    sBuf.nextEdge(memEdge[cnt++]);

                cntBuf[i] = cnt;
                for (int j = hh, jj; j > 1; j = jj) {
                    jj = (j >> 1);
                    if (memEdge[h[j][0]] < memEdge[h[jj][0]]) {
                        swap(h[j][0], h[jj][0]);
                        swap(h[j][1], h[jj][1]);
                    } else break;
                }
                cntBuf[i] = cnt;
            }
            maxi++;
            while (hh > 0) {
                int minEdge = h[1][0], pi = h[1][1];
                tmpBuf.insert(memEdge[minEdge]);

                if (minEdge < cntBuf[pi] - 1 || posFile[pi] < (len - 1) / edgePerBuf) {
                    if (minEdge < cntBuf[pi] - 1) h[1][0] = minEdge + 1;
                    else {
                        posFile[pi]++;
                        fseeko(sBuf.inFile, (lastPos + pi * len + posFile[pi] * edgePerBuf) * sizeof(edge), SEEK_SET);
                        sBuf.start();
                        cntBuf[pi] = h[1][0] = pi * edgePerBuf;
                        while (!sBuf.isEnd && cntBuf[pi] - h[1][0] < edgePerBuf
                               && edgePerBuf * posFile[pi] + cntBuf[pi] - h[1][0] < len)
                            sBuf.nextEdge(memEdge[cntBuf[pi]++]);

                        if (sBuf.isEnd || edgePerBuf * posFile[pi] + cntBuf[pi] - h[1][0] >= len) posFile[pi] =
                                                                                                          (len - 1) /
                                                                                                          edgePerBuf;

                    }
                } else {
                    h[1][0] = h[hh][0];
                    h[1][1] = h[hh--][1];
                }

                for (int j = 1, jj; (j << 1) <= hh; j = jj) {
                    jj = (j << 1);
                    if (jj < hh && memEdge[h[jj + 1][0]] < memEdge[h[jj][0]]) jj++;
                    if (memEdge[h[jj][0]] < memEdge[h[j][0]]) {
                        swap(h[j][0], h[jj][0]);
                        swap(h[j][1], h[jj][1]);
                    } else break;
                }
            }
        }
//		break;
    }
//	tSort += tm.getTime();

//	printf("%s\n%s\n", tmpName, sName);
//	return;

    if (1 + strlen(tmpName) > 1 + strlen(sName)) {
        remove(sName);
        rename(tmpName, sName);
    } else {
        remove(sName);
    }
    free(sName);
    free(tmpName);
    free(memEdge);

    tSort2 += tm2.getTime();
    tSortS += tm.getTime();
    iSort2 += iCnt - iCntTmp;
    oSort2 += oCnt - oCntTmp;
    iSort += iCnt;
    oSort += oCnt;


//	if (iter > 1)	{	deDu(_sName, m);	printf("\nm edgePerMem time %d %d %lf\n\n", m, edgePerMem, tm.getTime());	}
}

void xSortL(char *_sName, long long &m) {

    if (m == 0) return;
    timer tm, tm1;
    tm.restart();
    iCnt = oCnt = 0;
    char *sName = (char *) malloc(1 + strlen(_sName));
    strcpy(sName, _sName);
    char *tmpName = (char *) malloc(1 + strlen(sName) + 4);
    strcpy(tmpName, sName);
    strcat(tmpName, ".tmp");

    int cnt, edgePerMem = memSize / sizeof(edgeL), edgePerBuf = bufSize / sizeof(edgeL);


    edgeL *memEdge = (edgeL *) malloc(edgePerMem * (long long) sizeof(edgeL));

    double tin = 0;

    {
        outBufL tmpBuf(tmpName);
        inBufL usEdge(sName);
        usEdge.start();
        while (!usEdge.isEnd) {
            cnt = 0;
            while (cnt < edgePerMem && !usEdge.isEnd)
                usEdge.nextEdge(memEdge[cnt++]);

            timer tTmp;
            sort(memEdge, memEdge + cnt);
            tin += tTmp.getTime();
            for (int i = 0; i < cnt; i++)
                tmpBuf.insert(memEdge[i]);
        }
    }

    tSort1 += tm1.getTime();

    int iCntTmp, oCntTmp;
    iSort1 += (iCntTmp = iCnt);
    oSort1 += (oCntTmp = oCnt);


    if (m <= edgePerMem) {
        if (1 + strlen(tmpName) > 1 + strlen(sName)) {
            remove(sName);
            rename(tmpName, sName);
        } else {
            remove(sName);
        }
        free(sName);
        free(tmpName);
        free(memEdge);

        tSortS += tm.getTime();

        return;
    }


    timer tm2;


    int h[edgePerMem / edgePerBuf][2], hh, cntBuf[edgePerMem / edgePerBuf];

    long long posFile[edgePerMem / edgePerBuf];
    int ii = 0;
    for (long long len = edgePerMem; len < m; len *= edgePerMem / edgePerBuf) {
        swap(sName, tmpName);
        inBufL sBuf(sName);
        outBufL tmpBuf(tmpName);

//show(len);

        for (long long lastPos = 0; lastPos < m; lastPos += len * (edgePerMem / edgePerBuf)) {
            cnt = 0;
            int maxi;
            for (int i = 0; i < edgePerMem / edgePerBuf && lastPos + i * len < m; i++) {
                maxi = i;

//printf("%d\n", (lastPos+i*len));

                fseeko(sBuf.inFile, (lastPos + i * len) * sizeof(edgeL), SEEK_SET);
                posFile[i] = 0;
                sBuf.start();
                h[hh = i + 1][0] = cnt;
                h[i + 1][1] = i;
                while (!sBuf.isEnd && cnt - h[hh][0] < edgePerBuf)
                    sBuf.nextEdge(memEdge[cnt++]);

                cntBuf[i] = cnt;
                for (int j = hh, jj; j > 1; j = jj) {
                    jj = (j >> 1);
                    if (memEdge[h[j][0]] < memEdge[h[jj][0]]) {
                        swap(h[j][0], h[jj][0]);
                        swap(h[j][1], h[jj][1]);
                    } else break;
                }
                cntBuf[i] = cnt;
            }
            maxi++;
            while (hh > 0) {
                int minEdge = h[1][0], pi = h[1][1];
                tmpBuf.insert(memEdge[minEdge]);

                if (minEdge < cntBuf[pi] - 1 || posFile[pi] < (len - 1) / edgePerBuf) {
                    if (minEdge < cntBuf[pi] - 1) h[1][0] = minEdge + 1;
                    else {
                        posFile[pi]++;
                        fseeko(sBuf.inFile, (lastPos + pi * len + posFile[pi] * edgePerBuf) * sizeof(edgeL), SEEK_SET);
                        sBuf.start();
                        cntBuf[pi] = h[1][0] = pi * edgePerBuf;
                        while (!sBuf.isEnd && cntBuf[pi] - h[1][0] < edgePerBuf
                               && edgePerBuf * posFile[pi] + cntBuf[pi] - h[1][0] < len)
                            sBuf.nextEdge(memEdge[cntBuf[pi]++]);

                        if (sBuf.isEnd || edgePerBuf * posFile[pi] + cntBuf[pi] - h[1][0] >= len) posFile[pi] =
                                                                                                          (len - 1) /
                                                                                                          edgePerBuf;

                    }
                } else {
                    h[1][0] = h[hh][0];
                    h[1][1] = h[hh--][1];
                }

                for (int j = 1, jj; (j << 1) <= hh; j = jj) {
                    jj = (j << 1);
                    if (jj < hh && memEdge[h[jj + 1][0]] < memEdge[h[jj][0]]) jj++;
                    if (memEdge[h[jj][0]] < memEdge[h[j][0]]) {
                        swap(h[j][0], h[jj][0]);
                        swap(h[j][1], h[jj][1]);
                    } else break;
                }
            }
        }
//		break;
    }
//	tSort += tm.getTime();

//	printf("%s\n%s\n", tmpName, sName);
//	return;

    if (1 + strlen(tmpName) > 1 + strlen(sName)) {
        remove(sName);
        rename(tmpName, sName);
    } else {
        remove(sName);
    }
    free(sName);
    free(tmpName);
    free(memEdge);

    tSort2 += tm2.getTime();
    tSortS += tm.getTime();
    iSort2 += iCnt - iCntTmp;
    oSort2 += oCnt - oCntTmp;
    iSort += iCnt;
    oSort += oCnt;


//	if (iter > 1)	{	deDu(_sName, m);	printf("\nm edgePerMem time %d %d %lf\n\n", m, edgePerMem, tm.getTime());	}
}

void toTxtL(char *sName) {

    char *dName = (char *) malloc(1 + strlen(sName) + 15);
    sprintf(dName, "%s.to.txt", sName);

    inBufL sBuf(sName);
    sBuf.start();
    FILE *dFile = fopen(dName, "w");
    for (edgeL e; !sBuf.isEnd;) {
        sBuf.nextEdge(e);
        fprintf(dFile, "%lld %lld %lld\n", e.x, e.y, e.w);
    }
    fclose(dFile);

}

long long checkByte(char *fileName) {
    FILE *pFile = fopen(fileName, "rb");
    int temp = ftello(pFile);
    fseeko(pFile, 0, SEEK_END);  //move pointer to the end of the file
    long long ans = ftello(pFile) + 0.0;
    fclose(pFile);
    return ans;
}

int rand32() {
    return rand() * 32767LL + rand();
}

