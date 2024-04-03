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
	Int n, max_degree, * deg;	//��ʼͼ
	Long m, * offs;
	Pairs* data, ** adj;
	vector<pair<pair<int, int>, double> > unselected;	//�������ͼ
	//Int* degn;

	Graph();
	~Graph();

	void read_bin(const string& infile);	//��ȡͼ�ṹ������deg[n], offs[n], data[m], adj[m]
	//�ṩ��������ͼ�ṹ�Ĳ�����1�����ѡ�񶥵㣻���õ㱻ѡ�����ж��������ھ��Ƿ�ѡ�����¹���ͼ�ṹ--ͳ�ƶ�������2�����ѡ��ߣ����ı䱻ѡ�������¹���ͼ�ṹ--ͳ�Ʊ���
	void read_bin(const string& infile, int voe, float scale);
	vector<int> edge_selected_bin(const string& infile, float scale);	//����δɾ����ʱ ����ĳ�ʼ����

	void creat_bin(const string& infile);

	Int get_nm() { return n; }
	Long get_em() { return m; }
	Int get_uem() { return unselected.size(); }
	Int* get_deg() { return deg; }
	Long* get_offs() { return offs; }
	Pairs** get_adj() { return adj; }
};

#endif // !GRAPH_H#pragma once
