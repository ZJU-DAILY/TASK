#include "uncertain-core.h"

Uncertain_Core::Uncertain_Core()
{
	kmax = 0;
	core = NULL;
	adj_core_nbrs = NULL;
	AdjC = NULL;
	/*adj_nbrs = NULL;
	Adj = NULL;*/
}

Uncertain_Core::~Uncertain_Core()
{
	if (core != NULL) {
		delete[] core;
	}
	if (adj_core_nbrs != NULL) {
		delete adj_core_nbrs;
	}
	if (AdjC != NULL) {
		delete AdjC;
	}
	/*if (adj_nbrs != NULL) delete[] adj_nbrs;
	if (Adj != NULL) delete[] Adj;*/
}

void Uncertain_Core::get_core() {
	//int corenum;
	//vec_i corenum_record(max_degree + 1, 0);
	int maxcore = 0;
	core = new Int[n];
	memset(core, 0, sizeof(Int) * n);
	vec_i d(n);
	//std::copy(deg, deg + n, d);
	vector<vec_i> D;
	D.resize(max_degree + 1);
	for (int i = 0; i < n; i++) {
		d[i] = deg[i];
		D[deg[i]].push_back(i);
	}
	for (int i = 0; i <= max_degree; i++) {
		//corenum = 0;
		while (!D[i].empty())
		{
			int v = D[i].back();
			D[i].pop_back();
			core[v] = i;
			//corenum++;

			int dd = deg[v];
			for (int k = 0; k < dd; k++) {
				int w = adj[v][k].u;
				if (d[w] > i) {
					auto it = remove(D[d[w]].begin(), D[d[w]].end(), w);	//��ȡw��λ��
					D[d[w]].erase(it, D[d[w]].end());	//ɾ��w
					d[w]--;
					D[d[w]].push_back(w);
				}
			}
			maxcore = i > maxcore ? i : maxcore;
		}
		//corenum_record[i] = corenum;
	}
	kmax = maxcore;
	//kmax_node_num = corenum_record[kmax];
}


void Uncertain_Core::get_core_sorted_adj() {
	adj_core_nbrs = new pair<int, int>[m];
	AdjC = new pair<int, int>*[n];
	AdjC[0] = adj_core_nbrs;
	for (int i = 0; i < n - 1; i++) {
		AdjC[i + 1] = AdjC[i] + deg[i];
	}
	for (int i = 0; i < n; i++) {
		int d = deg[i];
		for (int j = 0; j < d; j++) {
			int w = adj[i][j].u;
			AdjC[i][j].first = core[w];
			AdjC[i][j].second = w;
		}
	}
	for (int i = 0; i < n; i++) {
		sort(AdjC[i], AdjC[i] + deg[i], greater<pair<int, int> >());
	}
}


//����ȷ��ͼcore�㷨Ѱ���Ӽ� -- ͬʱ�޸�core
int Uncertain_Core::certain_graph_inset_edge(int root, int k, vec_b& color) {
	int count = 0;  //core�仯�Ľڵ���
	vec_i cd(n, 0);  //ÿ���ڵ��cd  
	vec_i candiNode(n, -1);
	std::fill(color.begin(), color.end(), false);   //���core����Ľڵ�

	//Ѱ�Һ�ѡ���㼯
	std::queue<int> que;
	vec_b visited(n, false);    //�ýڵ��Ƿ񱻷���
	que.push(root);
	visited[root] = true;
	while (!que.empty()) {
		int r = que.front();
		que.pop();
		//cout << "r:" << r << endl;

		int d = deg[r];
		for (int i = 0; i < d; i++) {
			int w = adj[r][i].u;
			if (core[w] > k) {
				cd[r]++;
			}
			if (core[w] == k) {
				if (!visited[w]) {
					cd[r]++;
				}
				else if (cd[w] > k || cd[w] == 0) { //visited[w]:�Ѿ����ʼ���||�������δ����
					cd[r]++;
				}
			}
		}
		//cout << "cd[r]:" << cd[r] << endl;
		if (cd[r] > k) {
			candiNode[count] = r;
			count++;
			color[r] = true;
			for (int i = 0; i < d; i++) {
				int w = adj[r][i].u;
				if (core[w] == k && !visited[w]) {
					que.push(w);
					visited[w] = true;
				}
			}
		}
		else
		{
			for (int i = 0; i < d; i++) {
				int w = adj[r][i].u;
				if (cd[w] > k) {
					cd[w]--;
				}
			}
		}
	}

	//ȷ�Ϻ�ѡ�����Ƿ�color
	int flag = 1;
	int countNum = count;
	while (flag) {
		int recordCount = count;
		for (int i = 0; i < countNum; i++) {
			int v = candiNode[i];
			if (color[v] && cd[v] <= k) {
				color[v] = false;
				count--;
				for (int l = 0; l < deg[v]; l++) {
					int w = adj[v][l].u;
					if (cd[w] > k) {
						cd[w]--;
					}
				}
				break;  //һ���޸�count���������һ��whileѭ��
			}
		}
		if (recordCount == count) {
			flag = 0;
		}
	}

	//����core
	if (count) {
		for (int i = 0; i < countNum; i++) {
			int v = candiNode[i];
			if (color[v]) {
				core[v]++;
				//cout << "v:" << v << " core:" << core[v] << endl;
			}
		}
	}
	return count;
}

//��ʼͼ�����thres��ά��ͼ�ṹ+Ѱ�Һ�ѡ������+ֱ������
void Uncertain_Core::insert_core_compare() {
	double candidate_tm = 0.0;
	double recompute_tm = 0.0;
	int root;
	int mink;
	vec_b color(n);
	vec_i result_1(n, 0);
	vec_i result_2(n, 0);

	int um = unselected.size();
	cout << "um:" << um << endl;
	for (int i = 0; i < um; i++) {
		int u = unselected[i].first.first;
		int v = unselected[i].first.second;

		//ά��ͼ�ṹ
		int pos = deg[u];
		adj[u][pos].u = v;
		adj[u][pos].p = unselected[i].second;
		adj[u][pos].re = deg[v];

		pos = deg[v];
		adj[v][pos].u = u;
		adj[v][pos].p = unselected[i].second;
		adj[v][pos].re = deg[u];
		deg[u]++;
		deg[v]++;
		//cout << "u:" << u << " v:" << v << endl;

		double tm_1 = omp_get_wtime();
		root = u;
		mink = core[u];
		if (core[v] < core[u]) {
			root = v;
			mink = core[v];
		}
		//cout << "root:" << root << " mink:" << mink << endl;
		int count = certain_graph_inset_edge(root, mink, color);
		tm_1 = omp_get_wtime() - tm_1;
		for (int j = 0; j < n; j++) {
			result_1[j] = core[j];
		}
		//��¼insert�㷨��������core number
		candidate_tm += tm_1;
		cout << "count:" << count << endl;
		//cout << "insert_compute_candidate:" << tm_1 << endl;

		delete[] core;	//��֮ǰ��core��̬�ͷ�
		double tm_2 = omp_get_wtime();
		get_core();
		tm_2 = omp_get_wtime() - tm_2;
		recompute_tm += tm_2;
		for (int j = 0; j < n; j++) {
			result_2[j] = core[j];
		}
		//cout << "insert_threshold_recompute:" << tm_2 << endl;


		//�Ƚϲ���ߺ�ı仯�Ƚ�		
		bool compare = areVectorsEqual(result_1, result_2);
		if (!compare) {
			cout << "insert�������Ƿ���ͬ��" << compare << endl;
		}
	}
	cout << "insert_compute_candidate_time:" << candidate_tm << endl;
	cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
	unselected.clear();
}


//�Ƚ�����vector<int>�Ƿ���ͬ
bool Uncertain_Core::areVectorsEqual(const std::vector<int>& v1, const std::vector<int>& v2) {
	if (v1.size() != v2.size()) {
		return false;  // ������С��ͬ��ֱ�ӷ��� false
	}

	for (size_t i = 0; i < v1.size(); i++) {
		if (v1[i] != v2[i]) {
			return false;  // ���ֲ���ȵ�Ԫ�أ�ֱ�ӷ��� false
		}
	}

	return true;  // ������ȫ��ͬ
}


