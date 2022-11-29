#ifndef HEAD_H_
#define HEAD_H_

#include <stdio.h>
#include <math.h>
#include <vector>
#include <map>
#include <set>
#include<boost/algorithm/string/split.hpp>
#include<boost/algorithm/string/classification.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include<iostream>
#include<fstream>
#include<math.h>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <boost/thread/thread.hpp>
#include "Semaphore.h"

#include <functional>
#include <bits/stdc++.h>
#include <utility>
//to hash any given pair
struct hash_pair {
   template <class T1, class T2>
   size_t operator()(const pair<T1, T2>& p) const{
      auto hash1 = hash<T1>{}(p.first);
      auto hash2 = hash<T2>{}(p.second);
      return hash1 ^ hash2;
   }
};

#define INF 999999999
using namespace std;
using namespace boost;

struct Nei{
	int nid;
	int w;
	int c;
};

struct tri{
	int u;
	int v;//v is the contracted node
	int w;
};

struct Node{//tree node
	vector<pair<int,pair<int,int>>> vert;//neighID/weight/count
	vector<pair<int,Nei>> neighInf;//posID,neighID,weight,count(for shortcut information maintenance)
	vector<int> pos, pos2;
	vector<int> dis, cnt;//the distance value and corresponding count number
	set<int> changedPos;
	vector<bool> FN;//another succint way of FromNode
	set<int> DisRe;
	vector<int> ch;
	int height, hdepth;//hdepty is the deepest node that a vertex still exists
	int pa;//parent
	int uniqueVertex;
	vector<int> piv;//pivot vetex, used in path retrieval
	Node(){
		vert.clear();
		neighInf.clear();
		pos.clear();
		dis.clear();
		cnt.clear();
		ch.clear();
		pa = -1;
		uniqueVertex = -1;
		height = 0;
		hdepth = 0;
		changedPos.clear();
		FN.clear();
		DisRe.clear();
		piv.clear();
	}
};

class Graph{
public:
	int nodenum;
	int edgenum;

	int threadnum=40;

	int eNum;//(a,b) is a edge with a<b
	vector<pair<int,int>> Edge;//edgeID->(a,b) with a<b
	unordered_map<pair<int,int>, int, hash_pair> EdgeRe;//edge (a,b)->edgeID
	vector<set<int>> E1;

	//bidirectional graph
	vector<vector<pair<int,int>>> Neighbor;
	vector<pair<double,double>> GraphLocation;//location used in A* algorithm

	//vertex order & its inverted list
	vector<int> NodeOrder;
	vector<int> vNodeOrder;
	//in Lattice and Regular network, vertex order decided by Minimum Degree Elimination method
	void CHconsMTOrderGenerate();
	void insertEMTOrderGenerate(int u,int v,int w);
	void NeighborComOrderGenerate(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	//in small-world and scale-free network, vertex order decided by vertex degree
	void NoRegularGraphOrder();

	//ReadGraph
	//all node IDs start from 0; files in the form of (ID1,ID2,weight)
	void ReadGraph(string filename);
	void ReadGraphForCHPinc(string filename);
	//Read update test data
	void ReadIns(string filename, vector<pair<pair<int,int>,int>>& data);

	//correctness check
	void CorCheckDij();
	void CorCheckCHPorder();
	void CorCheckCH();
	void CorCheckH2H();
	void CorCheckPLL();

	//efficiency calculation
	void EffiCheckDij(string filename);
	void EffiCheckCHPorder(string filename);
	void EffiCheckCH(string filename);
	void EffiCheckH2H(string filename);
	void EffiCheckPLL(string filename);

	//Query processing
	int	QueryCHPorder(int ID1, int ID2);
	int	QueryCH(int ID1, int ID2);
	int QueryH2H(int ID1,int ID2);
	int ShortestDisQuery(int ID1,int ID2);

	//BiDijkstra & A*
	int Dij(int ID1, int ID2);
	int BiDij(int ID1, int ID2);
	int Astar(int ID1, int ID2);
	int EuclideanDis(int s, int t);

	//Index size computation
	long long H2Hindexsize();
	long long CHindexsize();
	long long SCconNodesize();
	long long CHPindexsize();
	long long CHPWitPathsize();
	long long PLLindexsize();
	long long PrunePointsize();

	//CHP (Contraction Hierarchy with Pruning)
	//index construction
	/****for decrease***/
	int StartCHPOrderMT(string indexfile, string orderfile);
	void CHPconsorderMT();
	int CHcontractionorderMT(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);
	/****for increase (witness path WP)***/
	int StartCHPOrderMTWP(string indexfile, string orderfile);
	void CHPconsorderMTWP();
	int CHcontractionorderMTWP(int k, int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW, vector<vector<pair<int, int>> >& vvpResult);
	//index write & read
	int writeShortCutorder(string filename);
	int ReadShortCut(string filename);
	void writeCHPIncrease(string indexfile);
	void readCHPIncrease(string indexfile);
	//index update
	void CHPdec(int a, int b, int oldW, int newW);
	void CHPinc(int a, int b, int oldW, int newW);
	/***intermediate variable****/
	Semaphore* smLock = new Semaphore(1);
	vector<map<int,int>> OutEdgesM;
	vector<vector<pair<int, int>>> vvNode;//the adjacent and shortcut during contraction
	vector<vector<pair<int, int>>> vvpShortCut;
	vector<vector<pair<int,pair<int,int>>>> OutEdgesV;//adjacent id / edge length / middle id
	vector<int>	vnContractedNeighbors;//the number of contracted neighbors
	vector<int> vnContractedNeighborsOld;
	vector<vector<pair<int, int>>> AdjaShort;//final adja+shorcut for query answering
	vector<map<int, vector<int>>> SupportNodes;//all possible supportive nodes of shortcut (u,w), uID<wID
	vector<unordered_map<pair<int,int>, int, hash_pair>> PathInfor; //u, (v,w) end points triple, length of path u<w
	vector<vector<tri>> EdgeOnPath;//(a,b) end points of edge, path (u,v,w)
	vector<vector<pair<int, int>>> AdjaShortR;//used for 1-N re-contraction, only store the lower adjacent nodes
	vector<Semaphore*> vSm;
	Semaphore* sm = new Semaphore(threadnum);
	int CHcontractInc(int ID1, int ID2, int dUV, vector<pair<int, int> >& vW, vector<pair<int,int>>& addedShortcut);
	int NewSCweight(int s, int t);
	vector<set<pair<int,int>>> InvalidWP;//record those invalid witness paths

	//CHW (Contraction Hierarchy without pruning)
	//index construction
	int StartCHOrderMT(string indexfile, string orderfile);
	void CHconsorderMT(string orderfile);//easier one
	void NeighborComorder(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	void insertEMTorder(int u,int v,int w);
	//index write & read
	int writeShortCutCH(string filename);
	int ReadShortCutCH(string filename);
	//index update
	void CHdecStr(int a,int b, int oldW, int newW);//UE decrease
	void CHincStrMT(int a,int b, int oldW, int newW);//UE increase
	void CHincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW increase
	void CHdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//CHW decrease
	/****intermediate variables and functions****/
	vector<int> DD,DD2;
	vector<map<int,pair<int,int>>> E;
	vector<vector<pair<int,pair<int,int>>>> NeighborCon;
	vector<map<int, vector<int>>> SCconNodesMT;
	void deleteEorder(int u,int v);
	void insertEorder(int u,int v,int w);
	void deleteE(int u,int v);


	//H2H
	//index construction
	void H2HconOrderMT(string orderfile);
	//index update
	void H2HdecBat(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//decrease
	void H2HincBatMT(vector<pair<pair<int,int>,pair<int,int>>>& wBatch);//increase
	/****intermediate variables and functions****/
	vector<int> rank;
	vector<Node> Tree;
	int heightMax;
	void makeTree();
	int match(int x,vector<pair<int,pair<int,int>>> &vert);
	vector<vector<int>> VidtoTNid;//one vertex exist in those tree nodes (nodeID--->tree node rank)
	vector<int> EulerSeq;//prepare for the LCA calculation
	vector<int> toRMQ;
	vector<vector<int>> RMQIndex;
	void makeRMQDFS(int p, int height);
	void makeRMQ();
	int LCAQuery(int _p, int _q);
	void makeIndex();//make index
	void makeIndexDFS(int p, vector<int> &list);
	void EachNodeProBDis5(int child,vector<int>& line,set<int>& vertexIDChL, map<int,int>& checkedDis);
	void eachNodeProcessIncrease1(int children, vector<int>& line, int& changelabel);

	//PLL
	//index construction
	int StartPLL(string filename, string filenameP, string orderfile);
	void PLLcon(string orderfile);
	void DijksPrune1(int nodeID, vector<pair<int,int>>& vp);
	//index write & read
	void readPLL(string filename, string filenameP);
	void writePLL(string filename, string filenameP);
	//index update
	void PLLdec(int a,int b, int oldW, int newW);//PLL-S decrease
	void PLLinc(int a,int b, int oldW, int newW);//PLL-S increase
    vector<pair<int, pair<int, int>>> PSLdec(int a,int b, int oldW, int newW);//PLL-P decrease
    vector<pair<int, pair<int, int>>> PSLinc(int a,int b, int oldW, int newW);//PLL-P increase
	/****intermediate variables and functions****/
	vector<unordered_map<int,int>> Label;
	vector<unordered_map<int,vector<int>>> PruningPointNew;//v {c,{u}}
	set<pair<int,int>> NoSupportedPair;
	int DisQueryPeak(int ID1, int ID2);
	int DisQueryVally(int ID1, int ID2);
	int DisQueryLower1(int ID1, int ID2);
	int ShortestDisQuery1(int ID1,int ID2,vector<int>& SupNode, int& d);
	void DijkstraPrune1(int NodeID, int ID, int d);
	int PrefixalDisQuery(int ID1, int ID2);
	int DijkstraList(int ID1, vector<int>& distance);
	void PLLaffect(int u, int oldw, vector<int>& PAu, vector<int>& RAu, vector<int> disu, vector<int> disv, vector<int> disvNew);
	void PLLinvalid(int u, vector<int>& ILu, vector<int>& MLu, vector<int> RAu, vector<int> RAv, vector<int> PAv);
	void PLLremove(int u, vector<int> RAu, vector<int> ILu);
	void PLLPrune(int r, set<int> Affect);

	pair<bool, pair<int,int>> VertexPair(int ID, int hopnum);

	//path retrieval
	/*Direct Search*/
	vector<int> BiDijPath(int ID1, int ID2);
	vector<int> AstarPath(int ID1, int ID2);
	/*CHW*/
	vector<int> PathRetriCH(int ID1, int ID2);
	void CHconsPath(string orderfile);
	void NeighborComPath(vector<pair<int,pair<int,int>>> &Neighvec, pair<int,int> p, int x);
	bool insertEPath(int u,int v,int w);
	bool insertEMTPath(int u,int v,int w);
	/*H2H*/
	vector<int> PathRetriH2H(int s, int t);//path retrieve function
	vector<int> PathRetriH2Hpartial(int ID1, int ID2);
	void H2HconPath(string orderfile);
	void makeIndexPath();
	void makeIndexDFSPath(int p, vector<int>& list);
	/*PLL*/
	vector<set<int>> AdjacentNodes;
	vector<vector<pair<int,int>>> label;
	vector<vector<pair<int,pair<int,int>>>> labelPre;//with precedent
	vector<unordered_map<int,int>> labelPos;
	void PLLconsPath(string orderfile, string indexfile);
	void DijkPLLPath(int nodeID, vector<pair<int,pair<int,int>>> &vp);
	void PLLPathIndexWrite(string file);
	void PLLPathIndexRead(string file);
	vector<int> PathRetriPLL(int s, int t);//path retrieve function
	int DistanceCompute(int s, int t);
	vector<int> PointCompute(int ID1, int ID2);
	/*CHP*/
	int StartCHPPath(string orderfile);
	int CHcontractionorderPath(int ID1, int ID2, vector<bool>& vbVisited, int dUV, vector<pair<int, int> >& vW);
	vector<int>	PathRetriCHP(int ID1, int ID2);//path retrieve function

};

namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

template<int log_k, typename k_t, typename id_t>
class heap {

public:

	// Expose types.
	typedef k_t key_t;
	typedef id_t node_t;

	// Some constants regarding the elements.
	//static const node_t NULLINDEX = 0xFFFFFFFF;
	static const node_t k = 1 << log_k;

	// A struct defining a heap element.
	struct element_t {
		key_t key;
		node_t element;
		element_t() : key(0), element(0) {}
		element_t(const key_t k, const node_t e) : key(k), element(e) {}
	};


public:

	// Constructor of the heap.
	heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {
	}

	heap() {

	}

	// Size of the heap.
	inline node_t size() const {
		return n;
	}

	// Heap empty?
	inline bool empty() const {
		return size() == 0;
	}

	// Extract min element.
	inline void extract_min(node_t &element, key_t &key) {
		assert(!empty());

		element_t &front = elements[0];

		// Assign element and key.
		element = front.element;
		key = front.key;

		// Replace elements[0] by last element.
		position[element] = NULLINDEX;
		--n;
		if (!empty()) {
			front = elements[n];
			position[front.element] = 0;
			sift_down(0);
		}
	}

	inline key_t top() {
		assert(!empty());

		element_t &front = elements[0];

		return front.key;

	}

	inline node_t top_value() {

		assert(!empty());

		element_t &front = elements[0];

		return front.element;
	}

	// Update an element of the heap.
	inline void update(const node_t element, const key_t key) {
		if (position[element] == NULLINDEX) {
			element_t &back = elements[n];
			back.key = key;
			back.element = element;
			position[element] = n;
			sift_up(n++);
		} else {
			node_t el_pos = position[element];
			element_t &el = elements[el_pos];
			if (key > el.key) {
				el.key = key;
				sift_down(el_pos);
			} else {
				el.key = key;
				sift_up(el_pos);
			}
		}
	}


	// Clear the heap.
	inline void clear() {
		for (node_t i = 0; i < n; ++i) {
			position[elements[i].element] = NULLINDEX;
		}
		n = 0;
	}

	// Cheaper clear.
	inline void clear(node_t v) {
		position[v] = NULLINDEX;
	}

	inline void clear_n() {
		n = 0;
	}


	// Test whether an element is contained in the heap.
	inline bool contains(const node_t element) const {
		return position[element] != NULLINDEX;
	}


protected:

	// Sift up an element.
	inline void sift_up(node_t i) {
		assert(i < n);
		node_t cur_i = i;
		while (cur_i > 0) {
			node_t parent_i = (cur_i-1) >> log_k;
			if (elements[parent_i].key > elements[cur_i].key)
				swap(cur_i, parent_i);
			else
				break;
			cur_i = parent_i;
		}
	}

	// Sift down an element.
	inline void sift_down(node_t i) {
		assert(i < n);

		while (true) {
			node_t min_ind = i;
			key_t min_key = elements[i].key;

			node_t child_ind_l = (i << log_k) + 1;
			node_t child_ind_u = std::min(child_ind_l + k, n);

			for (node_t j = child_ind_l; j < child_ind_u; ++j) {
				if (elements[j].key < min_key) {
					min_ind = j;
					min_key = elements[j].key;
				}
			}

			// Exchange?
			if (min_ind != i) {
				swap(i, min_ind);
				i = min_ind;
			} else {
				break;
			}
		}
	}

	// Swap two elements in the heap.
	inline void swap(const node_t i, const node_t j) {
		element_t &el_i = elements[i];
		element_t &el_j = elements[j];

		// Exchange positions
		position[el_i.element] = j;
		position[el_j.element] = i;

		// Exchange elements
		element_t temp = el_i;
		el_i = el_j;
		el_j = temp;
	}



private:

	// Number of elements in the heap.
	node_t n;

	// Number of maximal elements.
	node_t max_n;

	// Array of length heap_elements.
	vector<element_t> elements;

	// An array of positions for all elements.
	vector<node_t> position;
};
}

#endif /* HEAD_H_ */
