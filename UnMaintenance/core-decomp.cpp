#include "uncertain-core.h"


//初始计算--上一轮的概率，顶点v的度数
double Uncertain_Core::kprob_comp_topdown(bool flag, int v, vec_d& oldp, int k) {
    //printf("flag：%d\n", flag);
    //printf("v：%d\n", v);
    vec_d newp(k + 1, 0.0);
    vec_d kprob(k + 1, 0.0);
    kprob[0] = 1;
    int d = deg[v];
    int cp = k + 1;
    if (flag) {
        int count = 0;
        int end = 0;
        for (int i = 0; i < d; i++) {
            int u = adj[v][i].u;
            if (core[u] < k) {
                continue;
            }
            //printf("u：%d\n", u);
            count++;
            end = min(k, count + 1);
            double p = adj[v][i].p;
            for (int j = 1; j <= end; j++) {
                newp[j] = (1 - p) * oldp[j] + p * oldp[j - 1];
            }
            //printf("ending \n");
            std::copy_n(newp.begin(), cp, oldp.begin());
            std::fill(newp.begin(), newp.end(), 0.0);
        }
    }
    else
    {
        for (int i = 0; i < d; i++) {
            int u = adj[v][i].u;
            if (core[u] != k) {
                continue;
            }
            double p = adj[v][i].p;
            for (int j = 1; j <= k; j++) {
                newp[j] = (1 - p) * oldp[j] + p * oldp[j - 1];
            }
            std::copy_n(newp.begin(), cp, oldp.begin());
            std::fill(newp.begin(), newp.end(), 0.0);
        }
    }
    for (int i = 1; i <= k; i++) {
        kprob[i] = kprob[i - 1] - oldp[i];
        //printf("kprob: %lf\n", kprob[i]);
    }
    return kprob[k];
}


//更新计算--直接计算:visited表示所有的顶点
double Uncertain_Core::kprob_comp(int v, vec_b& visited, int k) {
    vec_d kprob(k + 1, 0.0);
    kprob[0] = 1.0;
    vec_d oldp(k + 1, 0.0);
    vec_d newp(k + 1, 0.0);
    oldp[1] = 1.0;

    int edges = 0;
    int end;
    int d = deg[v];
    int cp = k + 1;

    for (int i = 0; i < d; i++) {
        int u = adj[v][i].u;
        if (!visited[u]) {
            continue;
        }
        edges++;
        end = min(k, edges + 1);
        double p = adj[v][i].p;
        for (int j = 1; j <= end; j++) {
            newp[j] = (1 - p) * oldp[j] + p * oldp[j - 1];
        }
        std::copy_n(newp.begin(), cp, oldp.begin());
        std::fill(newp.begin(), newp.end(), 0.0);
    }
    for (int i = 1; i <= k; i++) {
        kprob[i] = kprob[i - 1] - oldp[i];
    }
    return kprob[k];
}




//更新计算--直接计算:visited表示其邻居
double Uncertain_Core::kprob_comp_scale(int v, vec_b& visited, int k) {
    vec_d kprob(k + 1, 0.0);
    kprob[0] = 1.0;
    vec_d oldp(k + 1, 0.0);
    vec_d newp(k + 1, 0.0);
    oldp[1] = 1.0;

    int edges = 0;
    int end;
    int d = deg[v];
    int cp = k + 1;

    for (int i = 0; i < d; i++) {
        if (!visited[i]) {
            continue;
        }
        int u = adj[v][i].u;
        edges++;
        end = min(k, edges + 1);
        double p = adj[v][i].p;
        for (int j = 1; j <= end; j++) {
            newp[j] = (1 - p) * oldp[j] + p * oldp[j - 1];
        }
        std::copy_n(newp.begin(), cp, oldp.begin());
        std::fill(newp.begin(), newp.end(), 0.0);
    }
    for (int i = 1; i <= k; i++) {
        kprob[i] = kprob[i - 1] - oldp[i];
    }
    return kprob[k];
}



double Uncertain_Core::kprob_comp_search_candi_increase(int v, int k, int num, vec_d& oldp, double& p) {
    vec_d newp(k + 1, 0.0);
    int end = min(k, num + 1);
    for (int j = 1; j <= end; j++) {
        newp[j] = (1 - p) * oldp[j] + p * oldp[j - 1];
    }
    std::copy_n(newp.begin(), k + 1, oldp.begin());

    vec_d kprob(k + 1, 0.0);
    kprob[0] = 1;
    for (int i = 1; i <= end; i++) {
        kprob[i] = kprob[i - 1] - oldp[i];
    }
    return kprob[end];
}


void Uncertain_Core::Initial_threshold_compute(vector<vector<double> >& thres) {
    //get_core_sorted_adj();

    vector<vec_d> kprobEqualK(n, vec_d(kmax + 1, 0.0));
    for (int i = 0; i < n; i++) {
        kprobEqualK[i][1] = 1.0;
    }
    vec_b initial(n, true);     //全局数组

    vec_d kprob(n);
    vec_b visited(n, true);		//顶点i的deg(i)邻边是否需要访问        
    vec_i nbr_visited_num(n);  //顶点i需要访问的邻边

    for (int k = kmax; k > 0; k--) {
        printf("k：%d\n", k);
        int count = 0;	//图中剩余顶点数

        std::fill(kprob.begin(), kprob.end(), 0.0);
        std::fill(visited.begin(), visited.end(), true);
        std::fill(nbr_visited_num.begin(), nbr_visited_num.end(), 0);

        for (int i = 0; i < n; i++) {
            if (core[i] < k) {
                visited[i] = false;
            }
            else
            {
                count++;
                nbr_visited_num[i] = deg[i];
                for (int j = 0; j < deg[i]; j++) {
                    int u = adj[i][j].u;
                    if (core[u] < k) {
                        nbr_visited_num[i]--;
                    }
                }
                //计算初始kprob
                kprob[i] = kprob_comp_topdown(initial[i], i, kprobEqualK[i], k);
                initial[i] = false;
                //printf(" kprob：%lf\n", kprob[i]);
            }
        }

        double curThres = 0.0;
        while (count)
        {
            int v = -1;
            double p = 2.0;
            for (int i = 0; i < n; i++) {
                if (visited[i]) {
                    //cout << "need visited i:" << i << endl;
                    if (p > kprob[i]) {
                        p = kprob[i];
                        v = i;
                    }
                }
            }

            curThres = max(curThres, p);
            if (curThres > 0) {
                thres[k - 1][v] = curThres;
                //printf(" thres：%lf\n", curThres);
            }
            visited[v] = false;
            count--;
            //更新相邻点的kprob
            for (int l = 0; l < deg[v]; l++) {
                int u = adj[v][l].u;
                if (visited[u]) {
                    nbr_visited_num[u]--;
                    if (nbr_visited_num[u] < k) {
                        kprob[u] = 0;
                    }
                    else
                    {
                        kprob[u] = kprob_comp(u, visited, k);
                    }
                    //cout << "nei:" << u << " kprob:" << kprob[u] << endl;
                    //更新到队列                     
                    //pq.swap(emptyPQ);                    
                }
            }

        }
    }
}



void Uncertain_Core::Initial_threshold_compute_map(vector<vector<double> >& thres) {
     // 定义平衡二叉搜索树和哈希表
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    double kprob;

    vector<vec_d> kprobEqualK(n, vec_d(kmax + 1, 0.0));
    for (int i = 0; i < n; i++) {
        kprobEqualK[i][1] = 1.0;
    }
    vec_b initial(n, true);     //全局数组
    vec_b visited(n, true);		//顶点i的deg(i)邻边是否需要访问        
    vec_i nbr_visited_num(n);  //顶点i需要访问的邻边

    for (int k = kmax; k > 0; k--) {
        //printf("k：%d\n", k);
        int count = 0;
        std::fill(visited.begin(), visited.end(), true);
        std::fill(nbr_visited_num.begin(), nbr_visited_num.end(), 0);

        for (int i = 0; i < n; i++) {
            if (core[i] < k) {
                visited[i] = false;
            }
            else
            {
                count++;
                nbr_visited_num[i] = deg[i];
                for (int j = 0; j < deg[i]; j++) {
                    int u = adj[i][j].u;
                    if (core[u] < k) {
                        nbr_visited_num[i]--;
                    }
                }
                //计算初始kprob
                kprob = kprob_comp_topdown(initial[i], i, kprobEqualK[i], k);
                Node* initialNode = new Node(i, kprob);
                mySet.insert(initialNode);
                hashTable[i] = initialNode;
                initial[i] = false;
                //printf(" kprob：%lf\n", kprob);
            }
        }

        double curThres = 0.0;
        while (count)
        {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int v = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            if (curThres > 0) {
                thres[k - 1][v] = curThres;
                //printf("v：%d", v);
                //printf("thres：%lf\n", curThres);
            }
            count--;
            mySet.erase(it);
            hashTable.erase(v);
            delete minValueNode;
            //cout << "count:" << mySet.size() << endl;
            visited[v] = false;

            //更新相邻点的kprob
            for (int l = 0; l < deg[v]; l++) {
                int u = adj[v][l].u;
                if (visited[u]) {
                    Node* neiNode = hashTable[u];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    nbr_visited_num[u]--;
                    if (nbr_visited_num[u] < k) {
                        neiNode->value = 0.0;
                    }
                    else
                    {
                        kprob = kprob_comp(u, visited, k);   
                        neiNode->value = kprob;
                    }                    
                    mySet.insert(neiNode);
                    hashTable[u] = neiNode;
                    //printf("%d\n",mySet.find(neiNode)->id);
                }
            }
        }
        //hashTable.clear();
    }

    // 释放剩余节点的内存
    for (auto node : mySet) {
        delete node;
    }
    mySet.clear();

    // 清空哈希表
    for (auto pair : hashTable) {
        delete pair.second;
    }
    hashTable.clear();
}


