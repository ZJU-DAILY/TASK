#pragma once
#ifndef GRAPH_H
#define GRAPH_H

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <omp.h>
#include <float.h>
#include <cstdio>
#include <iostream>

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

#define Int int
#define Long long long int
#define Double long double
#define Char char

struct edge
{
	Int u; double p;
	Int re; //The index of the reverse edge 
};
#define Pairs struct edge

class Graph
{
public:
	Int n, max_degree, * deg;	//初始图
	Long m, * offs;
	Pairs* data, ** adj;
	vector<pair<pair<int, int>, double> > unselected;	//放缩后的图
	//Int* degn;

	Graph();
	~Graph();

	void read_bin(const string& infile);	//读取图结构，生成deg[n], offs[n], data[m], adj[m]
	//提供两个更新图结构的操作：1）随机选择顶点；若该点被选择，则判断其所有邻居是否被选择，重新构造图结构--统计顶点数；2）随机选择边；若改变被选择，则重新构造图结构--统计边数
	void read_bin(const string& infile, int voe, float scale);
	vector<int> edge_selected_bin(const string& infile, float scale);	//返回未删除边时 顶点的初始度数

	void creat_bin(const string& infile);

	Int get_nm() { return n; }
	Long get_em() { return m; }
	Int get_uem() { return unselected.size(); }
	Int* get_deg() { return deg; }
	Long* get_offs() { return offs; }
	Pairs** get_adj() { return adj; }
};

#endif // !GRAPH_H#pragma once
