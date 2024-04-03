#pragma once
#include "graph.h"
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <functional>
#include <sys/stat.h>
#include <map>
#include <set>
#include <unordered_set>
#include <queue>

#define vec_i vector<int>
#define vec_b vector<bool>
#define vec_d vector<Double>
#define EPSILON 1e-9
#define EPSILON1 1e-6


//����ڵ�ṹ������ID�Ͷ�Ӧ��doubleֵ
struct Node {
	int id;
	double value;

	Node() : id(0), value(0.0) {}
	Node(int _id, double _value) : id(_id), value(_value) {}
};


// ���رȽϺ���������ָ���ıȽ�
struct NodePtrComparator {
	bool operator()(const Node* a, const Node* b) const {
		if (a->value == b->value) {
			return a->id < b->id;
		}
		return a->value < b->value;
	}
};


class Uncertain_Core : public Graph
{
public:
	Uncertain_Core();
	~Uncertain_Core();
	void get_core();
	void get_core_sorted_adj();

	//ȷ��ͼά��core�����Ƚ����������Ƿ���ͬ
	int certain_graph_inset_edge(int root, int k, vec_b& color);
	void insert_core_compare();
	bool areVectorsEqual(const std::vector<int>& v1, const std::vector<int>& v2);


	//kprob�Զ����¼��㡢�������ж����visited������㡢�����ھӵ�visited������㡢thres�ĳ�ʼ����
	double kprob_comp_topdown(bool flag, int v, vec_d& oldp, int k);
	double kprob_comp(int v, vec_b& visited, int k);
	double kprob_comp_scale(int v, vec_b& visited, int k);
	double kprob_comp_search_candi_increase(int v, int k, int num, vec_d& oldp, double& p);
	void Initial_threshold_compute(vector<vector<double> >& thres);
	void Initial_threshold_compute_map(vector<vector<double> >& thres);	//map��set-�����



	//����һ���ߺ�ȷ��root��range
	void insert_threshold_singlepoint_compute(vector<double>& kthres, double& low, double& upper, int& root, int u, int v, int k);
	//����root��range��Ѱ�Һ�ѡ��
	int insert_search_candi(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	//���º�ѡ����thres
	void insert_update_candidate_thres(vector<double>& kthres, int& k, int& count, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);

	//����mink+1������ά��core����Ǳ仯�ڵ� -- 1��mink=maxk����ʼͼΪcore����Ľڵ㣻2��root�ڵ��thres��ȻΪ0��upper=0
	int insert_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum);
	//���� thres>low || color=true ȷ����ʼͼ
	int insert_search_candi_core_change(vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	//��ȷ��ͼά��core�����Ƚ����������Ƿ���ͬ
	vector<vector<double> > insert_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	void insert_threshold_recompute(vector<vector<double> >& thres, double& time);
	void insert_threshold_compare(vector<vector<double> >& thres);


	//insert�Ż�������0����ʼ������1����̬����range��2��Ѱ�����Ƶ㣻3������ȷ����ѡ���thres��4���Ż������Ľ��
	vector<vector<double> > insert_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int i_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	

	//Ѱ�����Ƶ�
	vector<vector<double> > insert_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	void i_singlepoint_compute_restriction_point(vector<double>& kthres, double& low, double& upper, int& root, int u, int v, int k, vec_b& restRecord);
	

	//����ȷ����ѡ���thres
	vector<vector<double> > insert_threshold_compute_batchUP(vector<vector<double> >& thres, int u, int v, double& time, int& ct);


	//�Ż��������
	vector<vector<double> > insert_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);


	//�ȽϽ���Ƿ���ͬ
	void compareArrays(const std::vector<vector<double> >& array1, const std::vector<vector<double> >& array2, const vector<double>& comp1, const vector<double>& comp2, int mincore);	//�����Ѱ��low��upper���жϱ仯��Χ�Ƿ�Խ��
	bool compareCoreArrays(const std::vector<int>& array1, const std::vector<int>& array2);
	bool compareArraysTrueOrFalse(const std::vector<std::vector<double> >& array1, const std::vector<std::vector<double> >& array2);
	bool compareKsizeArraysTrueOrFalse(const std::vector<std::vector<double> >& array1, const std::vector<std::vector<double> >& array2, int k);	////�����Ƚ�ǰk��


	//��ɾ��
	void selectRandomEdges(double scale);
	void delete_threshold_compare(vector<vector<double> >& thres, double scale);
	double delete_compute_low(vector<double>& kthres, int k, int root);
	//����root[2]��range��Ѱ�Һ�ѡ��
	int delete_search_candi(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	void delete_update_candidate_thres(vector<double>& kthres, int& k, int& count, double low, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	vector<vector<double> > delete_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int delete_one_compute_candidate(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int& k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	int delete_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum);
	int delete_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	//vector<vector<double> > delete_threshold_recompute(double& time);	//����������Ĺ�����ͬ


	//delete�Ż�������0)��ʼ������1����̬����range��2�����ݷ������½磻3������ȷ����ѡ���thres��4���Ż������Ľ��
	vector<vector<double> > delete_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int d_one_compute_update_range(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	int d_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	int d_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);

	//���ݷ������½�
	void delete_compare_range(vector<vector<double> >& thres, double scale);
	vector<vector<double> > delete_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int d_one_compute_shrink_low(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	double d_shrink_lower_bound(vector<double>& kthres, int k, int root);

	//�������º�ѡ��
	vector<vector<double> > delete_threshold_compute_batchUp(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int d_one_compute_batch(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	void d_batch_update_candidate(vector<double>& kthres, int& k, int& count, double low, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);

	//�Ż��������
	vector<vector<double> > delete_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	int d_one_compute_opt(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable);
	Int get_kmax() { return kmax; }

	//private:
	int kmax;
	//int kmax_node_num;	//��kmax��core���нڵ�����
	Int* core;
	pair<int, int>* adj_core_nbrs, ** AdjC;
	//pair<Double, int>* adj_nbrs, ** Adj;
	//vec_i Del, is_Del;

};


#pragma once
