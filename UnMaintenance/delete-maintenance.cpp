#include "uncertain-core.h"

//在m条边中标记待删除的边
void Uncertain_Core::selectRandomEdges(double scale) {
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);
    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //确保每条边只处理一次
                if (rand() > scale_size) {  //待删除
                    //std::cout << "unselected-vertex:" << i << " nei:" << u << endl;
                    std::pair<int, int> p1 = std::make_pair(i, u);
                    double p = adj[i][j].p;
                    unselected.push_back(make_pair(p1, p));
                }
            }
        }
    }
    std::cout << "unselected edges:" << unselected.size() << endl;
}



//依次删除ct条边 then新增ct条边
void Uncertain_Core::delete_threshold_compare(vector<vector<double> >& thres, double scale) {
    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //每个方法寻找候选对象个数的累计和

    selectRandomEdges(scale);
    vector<vector<double> > copy_thres = thres;

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "i:" << i << " u:" << u << " v:" << v << endl;

        //维护图结构
        bool flag = false;
        int d = deg[u];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[u][j].u;
                if (w == v) {
                    flag = true;
                }
            }
            else
            {
                adj[u][j - 1].u = adj[u][j].u;
                adj[u][j - 1].p = adj[u][j].p;
            }
        }

        flag = false;
        d = deg[v];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[v][j].u;
                if (w == u) {
                    flag = true;
                }
            }
            else
            {
                adj[v][j - 1].u = adj[v][j].u;
                adj[v][j - 1].p = adj[v][j].p;
            }
        }
        deg[u]--;
        deg[v]--;


        vector<vector<double> > range = delete_threshold_compute_candidate(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_update_range(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_shrink_low(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_batchUp(thres, u, v, tm_1, candi_count);
        //vector<vector<double> > range = delete_threshold_compute_opt(thres, u, v, tm_1, candi_count);   //the best
        candidate_tm += tm_1;
        candi_sum += candi_count;


        /*if ((i + 1) % 10 == 0) {
            std::cout << "delete_compute_candidate:" << candidate_tm << " candi_count:" << candi_sum << endl;
        }*/

        delete[] core;	//将之前的core动态释放   
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);
        recompute_tm += tm_2;
        //cout << "insert_threshold_recompute:" << tm_2 << endl;

        //比较插入边后的变化阈值
        bool compare = compareArraysTrueOrFalse(thres, thres_2);
        if (!compare) {
            cout << "数组是否相同----------------------------------------------：" << compare << endl;
            //compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
        }


        //将newArray复制给thres
        if (thres.size() == thres_2.size()) {
            thres = thres_2;
        }
        else
        {
            thres.resize(thres_2.size());
            for (size_t i = 0; i < kmax; i++) {
                thres[i] = thres_2[i];
            }
        } 
    }
}





//删除边计算low
double Uncertain_Core::delete_compute_low(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //直接邻居数判断
    if (d < k) {
        low = 0;
    }
    else
    {
        int count = 0;  //记录kprob不小于root的邻居数量
        vec_b visited(d, false);
        for (int t = 0; t < d; t++) {
            int w = adj[root][t].u;
            double th = kthres[w];
            if ((th - upper) >= -EPSILON) {
                visited[t] = true;
                count++;
            }
        }

        if (count >= k) {
            low = kprob_comp_scale(root, visited, k);
        }
        else
        {
            low = 0;
        }

    }
    //cout << "root:" << root << " low:" << low << " upper:" << upper << endl;
    return low;
}


int Uncertain_Core::delete_search_candi(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //初始化需要访问的节点
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) > EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;
    vec_b visitedCandi(n, false);
    que.push(root[0]);
    visitedCandi[root[0]] = true;
    if (root[1] != -1) {    //放入第二个root
        que.push(root[1]);
        visitedCandi[root[1]] = true;
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        kprob = kprob_comp(node, visited, k + 1);
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;     
        candiRIndex[node] = count;
        //cout << "candiNode:" << node << " prob:" << candiProb[count] << endl;
        count++;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - low) > EPSILON) {  //kthres[w] > low
                neiNum[node]++;
                if (!visitedCandi[w] && (kthres[w] - upper) < EPSILON) {    //kthres[w] <= upper
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "nei:" << w << endl;
                }
            }
        }
    }

    return count;
}


//更新候选集的thres，且初始值设置为low
void Uncertain_Core::delete_update_candidate_thres(vector<double>& kthres, int& k, int& count, double low, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int index;
    double curThres = low;
    int flag = count;
    double kprob;
    while (flag) {
        auto it = mySet.begin();
        Node* minValueNode = *(it);
        int minV = minValueNode->id;
        double p = minValueNode->value;

        curThres = max(curThres, p);
        if (curThres > 0) {
            kthres[minV] = curThres;
        }
        flag--;
        visited[minV] = false;
        candiRIndex[minV] = -1;
        mySet.erase(it);
        hashTable.erase(minV);
        delete minValueNode;

        //更新其邻居
        for (int t = 0; t < deg[minV]; t++) {
            int w = adj[minV][t].u;
            index = candiRIndex[w];

            //w在candiSet中才需要更新
            if (index != -1) {
                Node* neiNode = hashTable[w];
                it = mySet.find(neiNode);
                mySet.erase(it);
                neiNum[w]--;
                if (neiNum[w] < k) {
                    neiNode->value = 0.0;
                }
                else
                {
                    kprob = kprob_comp(w, visited, k + 1);
                    neiNode->value = kprob;
                }
                mySet.insert(neiNode);
                hashTable[w] = neiNode;
                //cout << "nei:" << w << " kprob:" << candiProb[index] << endl;
            }
        }
    }

    for (auto node : mySet) {
        delete node;
    }
    mySet.clear();

    for (auto pair : hashTable) {
        delete pair.second;
    }
    hashTable.clear();
}



//删除边后 对每个k单独操作 -- core不改变的情况下
int Uncertain_Core::delete_one_compute_candidate(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int& k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //一般只有一个root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //确定thres变化的范围
    if (std::fabs(thres1 - thres2) < EPSILON) {
        root[0] = u;
        up = thres1;
        double kprob1 = delete_compute_low(kthres, k + 1, u);

        root[1] = v;
        double kprob2 = delete_compute_low(kthres, k + 1, v);

        low = min(kprob1, kprob2);
    }
    else if (thres1 > thres2) {
        root[0] = v;
        up = thres2;
        low = delete_compute_low(kthres, k + 1, v);
    }
    else
    {
        root[0] = u;
        up = thres1;
        low = delete_compute_low(kthres, k + 1, u);
    }
    range[0][k] = low;
    range[1][k] = up;

    if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
        return count;
    }

    //寻找候选集
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //更新候选集的thres：在k = kk-1时，维护coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


//删除边 寻找core减小的点 -- 将core减小的点依次放入que，然后寻找thres可能改变的候选点
int Uncertain_Core::delete_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum) {
    int count = 0;
    int xdeg = 0;   //root的core是否减少
    vec_i cd(n, 0);
    std::fill(color.begin(), color.end(), false);   //标记core增大的节点
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);

    std::queue<int> que;
    vec_b visited(n, false);
    int d = deg[root];
    for (int i = 0; i < d; i++) {
        int w = adj[root][i].u;
        if (core[w] >= k) {
            xdeg++;
        }
    }
    if (xdeg < k) { //root的coreness确定减小
        que.push(root);
        visited[root] = true;
    }
    while (!que.empty()) {
        int r = que.front();
        que.pop();
        candiNode[count] = r;
        candiRIndex[r] = count;
        count++;    //与root联通&core=k的子图节点数量
        color[r] = true;

        d = deg[r];
        for (int i = 0; i < d; i++) {
            int w = adj[r][i].u;
            if (core[w] >= k) {
                cd[r]++;
            }
            if (core[w] == k && !visited[w]) {
                que.push(w);
                visited[w] = true;
            }
        }
    }

    //确定候选顶点是否被color
    int flag = 1;
    countNum = count;
    count = 0;
    while (flag) {
        int recordCount = count;
        for (int i = 0; i < countNum; i++) {
            int v = candiNode[i];
            if (color[v] && (cd[v] < k)) {
                color[v] = false;
                count++;
                for (int l = 0; l < deg[v]; l++) {
                    int w = adj[v][l].u;
                    if (color[w] && (cd[w] >= k)) { //在候选集中&&再次判断
                        cd[w]--;
                    }
                }
            }
        }
        if (recordCount == count) {
            flag = 0;
        }
    }

    //更新core
    if (count) {
        for (int i = 0; i < countNum; i++) {
            int v = candiNode[i];
            if (!color[v]) {
                core[v]--;
            }
        }
    }

    return count;
}


//core变化后 维护thres=0 寻找可能改变的候选集
int Uncertain_Core::delete_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double up = 0;  //que中thres最大值
    double kprob;
    std::queue<int> que;
    vec_b visitedCandi(n, false);
    //将core改变影响的点放入que，并维护其thres -- 初始化que：一级扩散
    for (int i = 0; i < countNum; i++) {
        int r = candiNode[i];
        if (!visited[r]) {
            int d = deg[r];
            up = max(up, kthres[r]);
            //double th = kthres[r];
            for (int t = 0; t < d; t++) {
                int w = adj[r][t].u;
                if (!visitedCandi[w] && (core[w] > k)) {
                    if ((kthres[w] - up) < EPSILON) {   //kthres[w] <= th
                        que.push(w);
                        visitedCandi[w] = true;
                    }
                }
            }
            //维护thres
            kthres[r] = 0;
        }
    }

    //初始化候选记录
    int count = 0;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    for (int i = 0; i < n; i++) {
        if (core[i] > k) {
            visited[i] = true;
        }
    }

    if ((root != -1) && (core[root] > k)) {
        if (!visitedCandi[root]) {
            que.push(root); //考虑两个分离的连通分量
            visitedCandi[root] = true;
            up = max(up, kthres[root]);
        }
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        kprob = kprob_comp(node, visited, k + 1);    //初始邻居数可能小于k，计算结果=0   
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;
        candiRIndex[node] = count;
        //cout << "candiNode:" << node << " prob:" << candiProb[count] << endl;
        count++;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if (core[w] > k) {  //kthres[w] > low
                neiNum[node]++;
                if (!visitedCandi[w] && (kthres[w] - up) < EPSILON) {    //kthres[w] <= upper
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "nei:" << w << endl;
                }
            }
        }
    }

    return count;
}


//删除边后 根据候选集进行更新
vector<vector<double> > Uncertain_Core::delete_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "删除边 & 根据原始图计算上下界——————" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //候选对象总数量
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk - 1; k++) {
        //std::cout << "k:" << k << endl;
        count = delete_one_compute_candidate(range, thres[k], u, v, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
    }

    int curK = kk - 1;
    int countNum;   //候选点数量
    //k = kk-1：确定图core维护算法 + 寻找候选集 + 更新thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //确定图core维护 -- kmax必然不变
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量: core变化 + v
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = delete_one_compute_candidate(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else if (core[u] > core[v]) {
        //std::cout << "core[u] > core[v]" << endl;
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = delete_one_compute_candidate(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else
    {
        //std::cout << "core[u] = core[v]" << endl;
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);
        if (count) {    //core[u]变小 
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]、core[v]都变
                    for (int i = 0; i < countNumV; i++) {
                        int node = candiNode_2[i];
                        if (!visited_2[node]) {
                            candiNode[countNum] = node;
                            node_visited[node] = false;
                            countNum++;
                        }
                    }
                }
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都变 -- 没有root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited已考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //维护kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax查询顶点数：" << t << endl;
                            break;
                        }
                    }
                    if (eraseflag) {
                        thres.erase(thres.begin() + kmax - 1);
                        kmax--;
                        std::cout << "kk == kmax && kmax--" << endl;
                    }
                }
            }
            else
            {
                //core[u]改变，core[v]不变
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]不变
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都不变 -- 执行上述循环
                count = delete_one_compute_candidate(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]不变，core[v]减小
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
    }
    ct += count;

    time = omp_get_wtime() - tm;
    return range;
}



//delete优化方案：0)初始方案；1）动态更新range；2）阶梯法收缩下界；3）批量确定候选点的thres；4）优化方案的结合
vector<vector<double> > Uncertain_Core::delete_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "candidate--kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));    //将所求range作为返回值

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk - 1; k++) {
        //std::cout << "k:" << k << endl;
        count = d_one_compute_update_range(range, thres[k], u, v, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
    }

    int curK = kk - 1;
    int countNum;   //候选点数量
    //k = kk-1：确定图core维护算法 + 寻找候选集 + 更新thres
    if (core[u] < core[v]) {
        //确定图core维护 -- kmax必然不变
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量: core变化 + v
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_update_range(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else if (core[u] > core[v]) {
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_update_range(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else
    {
        //core[u] = core[v]
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);
        if (count) {    //core[u]变小 
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]、core[v]都变
                    for (int i = 0; i < countNumV; i++) {
                        int node = candiNode_2[i];
                        if (!visited_2[node]) {
                            candiNode[countNum] = node;
                            node_visited[node] = false;
                            countNum++;
                        }
                    }
                }
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都变 -- 没有root
                count = d_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited已考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //维护kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax查询顶点数：" << t << endl;
                            break;
                        }
                    }
                    if (eraseflag) {
                        thres.erase(thres.begin() + kmax - 1);
                        kmax--;
                        std::cout << "kk == kmax && kmax--" << endl;
                    }
                }
            }
            else
            {
                //core[u]改变，core[v]不变
                count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]不变
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都不变 -- 执行上述循环
                count = d_one_compute_update_range(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]不变，core[v]减小
                count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
    }
    ct += count;

    time = omp_get_wtime() - tm;
    return range;
}

//单次计算--实时更新range上下界
int Uncertain_Core::d_one_compute_update_range(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //一般只有一个root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //确定thres变化的范围
    if (std::fabs(thres1 - thres2) < EPSILON) {
        root[0] = u;
        up = thres1;
        double kprob1 = delete_compute_low(kthres, k + 1, u);

        root[1] = v;
        double kprob2 = delete_compute_low(kthres, k + 1, v);

        low = min(kprob1, kprob2);
    }
    else if (thres1 > thres2) {
        root[0] = v;
        up = thres2;
        low = delete_compute_low(kthres, k + 1, v);
    }
    else
    {
        root[0] = u;
        up = thres1;
        low = delete_compute_low(kthres, k + 1, u);
    }
    range[0][k] = low;
    range[1][k] = up;

    if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
        return count;
    }

    //寻找候选集
    count = d_candi_dynamic_update_range(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //更新候选集的thres：在k = kk-1时，维护coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


//寻找候选集时 动态收缩边界
int Uncertain_Core::d_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double tempUP, kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //初始化需要访问的节点
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) > EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;
    vec_b visitedCandi(n, false);
    que.push(root[0]);
    visitedCandi[root[0]] = true;
    if (root[1] != -1) {
        que.push(root[1]);
        visitedCandi[root[1]] = true;
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        kprob = kprob_comp(node, visited, k + 1);
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;
        candiRIndex[node] = count;
        //cout << "candiNode:" << node << " prob:" << candiProb[count] << endl;
        count++;

        if (kthres[node] < upper) { //upper严格变小
            tempUP = kthres[node];
        }
        else
        {
            tempUP = upper;
        }

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - low) > EPSILON) {
                neiNum[node]++;
                //if (!visitedCandi[w] && (kthres[w] - kthres[node]) < EPSILON) {     //迭代更新upper
                //    que.push(w);
                //    visitedCandi[w] = true;
                //}

                if (!visitedCandi[w] && (kthres[w] - tempUP) < EPSILON) {     //迭代更新upper
                    que.push(w);
                    visitedCandi[w] = true;
                }
            }
        }
    }

    return count;
}

//kk考虑两个联通结构 -- 实时更改range
int Uncertain_Core::d_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double up = 0;
    double kprob;
    std::queue<int> que;
    vec_b visitedCandi(n, false);
    //将core改变影响的点放入que，并维护其thres -- 初始化que：一级扩散
    for (int i = 0; i < countNum; i++) {
        int r = candiNode[i];
        if (!visited[r]) {
            int d = deg[r];
            up = max(up, kthres[r]);
            //double th = kthres[r];
            for (int t = 0; t < d; t++) {
                int w = adj[r][t].u;
                if (!visitedCandi[w] && (core[w] > k)) {
                    if ((kthres[w] - kthres[r]) < EPSILON) {   //kthres[w] <= th
                        que.push(w);
                        visitedCandi[w] = true;
                    }
                }
            }
            //维护thres
            kthres[r] = 0;
        }
    }

    //初始化候选记录
    int count = 0;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    for (int i = 0; i < n; i++) {
        //if (kthres[i] > EPSILON) {    //low=0
        //    visited[i] = true;
        //}
        if (core[i] > k) {
            visited[i] = true;
        }
    }

    if ((root != -1) && (core[root] > k)) {
        if (!visitedCandi[root] && (kthres[root] - up) < EPSILON) {
            que.push(root); //考虑两个分离的连通分量
            visitedCandi[root] = true;
        }
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        kprob = kprob_comp(node, visited, k + 1);
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;
        candiRIndex[node] = count;
        //cout << "candiNode:" << node << " prob:" << candiProb[count] << endl;
        count++;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if (core[w] > k) {  //kthres[w] > low
                neiNum[node]++;
                if (!visitedCandi[w] && (kthres[w] - kthres[node]) < EPSILON) {    //kthres[w] <= upper
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "nei:" << w << endl;
                }
            }
        }
    }

    return count;
}

//删除边 比较两个方法计算的low
void Uncertain_Core::delete_compare_range(vector<vector<double> >& thres, double scale) {
    double tm_1;
    double tm_2;
    int count1, count2;
    selectRandomEdges(scale);

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "第i条边：" << i << " u:" << u << " v:" << v << endl;

        //维护图结构
        bool flag = false;
        int d = deg[u];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[u][j].u;
                if (w == v) {
                    flag = true;
                }
            }
            else
            {
                adj[u][j - 1].u = adj[u][j].u;
                adj[u][j - 1].p = adj[u][j].p;
            }
        }

        flag = false;
        d = deg[v];
        for (int j = 0; j < d; j++) {
            if (!flag) {
                int w = adj[v][j].u;
                if (w == u) {
                    flag = true;
                }
            }
            else
            {
                adj[v][j - 1].u = adj[v][j].u;
                adj[v][j - 1].p = adj[v][j].p;
            }
        }
        deg[u]--;
        deg[v]--;


        int mincore = min(core[u], core[v]);
        //初始计算下界
        vector<vector<double> > range1 = delete_threshold_compute_candidate(thres, u, v, tm_1, count1);
        //阶梯法求取下界
        vector<vector<double> > range2 = delete_threshold_compute_shrink_low(thres, u, v, tm_1, count2);

        delete[] core;
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);

        //比较删除边后的变化阈值
        compareArrays(thres, thres_2, range2[0], range2[1], mincore - 1);   //不能比较core减小的点

        for (int i = 0; i < range1[0].size(); i++) {
            if ((range2[0][i] - range1[0][i]) > EPSILON) {
                std::cout << "low收缩 ------- k：" << i;
                std::cout << "  initial：" << range1[0][i] << " opt:" << range2[0][i] << endl;
            }
        }

        //将newArray复制给thres
        if (thres.size() == thres_2.size()) {
            thres = thres_2;
        }
        else
        {
            thres.resize(thres_2.size());
            for (size_t i = 0; i < kmax; i++) {
                thres[i] = thres_2[i];
            }
        }
    }
    unselected.clear();
}


//优化方案 -- 阶梯法求取下界
vector<vector<double> > Uncertain_Core::delete_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "删除边 & 根据原始图计算上下界——————" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //候选对象总数量
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk - 1; k++) {
        //std::cout << "k:" << k << endl;
        count = d_one_compute_shrink_low(range, thres[k], u, v, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
    }

    int curK = kk - 1;
    int countNum;   //候选点数量
    //k = kk-1：确定图core维护算法 + 寻找候选集 + 更新thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //确定图core维护 -- kmax必然不变
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量: core变化 + v
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_shrink_low(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else if (core[u] > core[v]) {
        //std::cout << "core[u] > core[v]" << endl;
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
            if (count) {
                delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_shrink_low(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else
    {
        //std::cout << "core[u] = core[v]" << endl;
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);
        if (count) {    //core[u]变小 
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]、core[v]都变
                    for (int i = 0; i < countNumV; i++) {
                        int node = candiNode_2[i];
                        if (!visited_2[node]) {
                            candiNode[countNum] = node;
                            node_visited[node] = false;
                            countNum++;
                        }
                    }
                }
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都变 -- 没有root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited已考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //维护kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax查询顶点数：" << t << endl;
                            break;
                        }
                    }
                    if (eraseflag) {
                        thres.erase(thres.begin() + kmax - 1);
                        kmax--;
                        std::cout << "kk == kmax && kmax--" << endl;
                    }
                }
            }
            else
            {
                //core[u]改变，core[v]不变
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]不变
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都不变 -- 执行上述循环
                count = d_one_compute_shrink_low(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]不变，core[v]减小
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
    }
    ct += count;

    time = omp_get_wtime() - tm;
    return range;
}

//单次计算 -- 收缩上下界
int Uncertain_Core::d_one_compute_shrink_low(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //一般只有一个root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //确定thres变化的范围
    if (std::fabs(thres1 - thres2) < EPSILON) {
        root[0] = u;
        up = thres1;
        double kprob1 = d_shrink_lower_bound(kthres, k + 1, u);

        root[1] = v;
        double kprob2 = d_shrink_lower_bound(kthres, k + 1, v);

        low = min(kprob1, kprob2);
    }
    else if (thres1 > thres2) {
        root[0] = v;
        up = thres2;
        low = d_shrink_lower_bound(kthres, k + 1, v);
    }
    else
    {
        root[0] = u;
        up = thres1;
        low = d_shrink_lower_bound(kthres, k + 1, u);
    }
    range[0][k] = low;
    range[1][k] = up;

    if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
        return count;
    }

    //寻找候选集
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //更新候选集的thres：在k = kk-1时，维护coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}



//进一步收缩下界 -- 分为opt1和opt2 两者单独实验 优化性能不同
double Uncertain_Core::d_shrink_lower_bound(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //直接邻居数判断
    if (d < k) {
        low = 0;
    }
    else
    {
        vec_d kprobKqualK(k + 1, 0.0);  //初始kprob值
        kprobKqualK[1] = 1;
        vec_d newp(k + 1, 0.0);
        int end;
        double p, kprob;

        int coreD = 0;  //考虑core不小于k的邻居数量
        int count = 0;  //记录kprob不小于root的邻居数量
        vector<pair<double, int> > unvisitedNei; //记录kprob小于root的邻居
        for (int t = 0; t < d; t++) {
            int w = adj[root][t].u;
            double th = kthres[w];
            if ((th - upper) >= -EPSILON) {     //th >= upper                
                coreD++;
                count++;    //现有邻居数量
                end = min(k, count + 1);
                p = adj[root][t].p;
                for (int j = 1; j <= end; j++) {
                    newp[j] = (1 - p) * kprobKqualK[j] + p * kprobKqualK[j - 1];
                }
                std::copy_n(newp.begin(), k + 1, kprobKqualK.begin());
                std::fill(newp.begin(), newp.end(), 0.0);
            }
            else
            {
                if (core[w] >= k) {
                    unvisitedNei.push_back(make_pair(th, t));
                    coreD++;
                }
            }
        }

        //d >= coreD >= count
        if (coreD < k) {
            low = 0;
        }
        else if (count >= k) {  //利用邻居收缩下界
            low = 1;
            for (int j = 1; j <= end; j++) {
                low = low - kprobKqualK[j]; //初始low
            }

            //if (coreD > count) {    //确保unvisitedNei中有邻居 -- opt1
            //    sort(unvisitedNei.begin(), unvisitedNei.end(), [](const pair<double, int>& p1, const pair<double, int>& p2) {
            //        return p1.first > p2.first;
            //        });

            //    double th = unvisitedNei[0].first;
            //    if ((th - low) > EPSILON) {
            //        bool flag = true;
            //        int countNum = count;
            //        while (flag) {
            //            int index = unvisitedNei[countNum - count].second;
            //            int w = adj[root][index].u;
            //            p = adj[root][index].u;
            //            countNum++;
            //            kprob = kprob_comp_search_candi_increase(root, k, countNum, kprobKqualK, p);

            //            if ((kprob - kthres[w]) >= -EPSILON) {    //确定low取值
            //                low = kthres[w];
            //                flag = false;
            //            }
            //            else
            //            {
            //                low = kprob;    //能否一直赋值？？？
            //                if ((countNum == coreD) || ((kprob - unvisitedNei[countNum - count].first) >= -EPSILON)) {      //意味着新增一条边 kprob必然增大
            //                    flag = false;
            //                }
            //            }
            //        }
            //    }
            //}

        }
        else
        {
            //优先选择剩下的邻居 -- opt2
            sort(unvisitedNei.begin(), unvisitedNei.end(), [](const pair<double, int>& p1, const pair<double, int>& p2) {
                return p1.first > p2.first;
                });

            bool flag = true;
            int countNum = count;
            while (flag) {
                int index = unvisitedNei[countNum - count].second;
                //std::cout << "d:" << d << " countNum:" << countNum << " count:" << count << endl;
                int w = adj[root][index].u;
                double p = adj[root][index].p;
                countNum++;
                kprob = kprob_comp_search_candi_increase(root, k, countNum, kprobKqualK, p);  //增量式计算

                if ((kprob - kthres[w]) >= -EPSILON) {    //确定low取值
                    low = kthres[w];
                    flag = false;
                }
                else
                {
                    low = 0;
                    if (countNum == coreD) {    //没有邻居可以测试
                        flag = false;
                    }
                }
            }
        }
    }
    //cout << "root:" << root << " low:" << low << " upper:" << upper << endl;
    return low;
}


//批量更新候选集
vector<vector<double> > Uncertain_Core::delete_threshold_compute_batchUp(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "删除边 & 根据原始图计算上下界——————" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //候选对象总数量
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk - 1; k++) {
        //std::cout << "k:" << k << endl;
        count = d_one_compute_batch(range, thres[k], u, v, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
    }

    int curK = kk - 1;
    int countNum;   //候选点数量
    //k = kk-1：确定图core维护算法 + 寻找候选集 + 更新thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //确定图core维护 -- kmax必然不变
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量: core变化 + v
            if (count) {
                d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_batch(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else if (core[u] > core[v]) {
        //std::cout << "core[u] > core[v]" << endl;
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
            if (count) {
                d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_batch(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else
    {
        //std::cout << "core[u] = core[v]" << endl;
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);
        if (count) {    //core[u]变小 
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]、core[v]都变
                    for (int i = 0; i < countNumV; i++) {
                        int node = candiNode_2[i];
                        if (!visited_2[node]) {
                            candiNode[countNum] = node;
                            node_visited[node] = false;
                            countNum++;
                        }
                    }
                }
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都变 -- 没有root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited已考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //维护kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax查询顶点数：" << t << endl;
                            break;
                        }
                    }
                    if (eraseflag) {
                        thres.erase(thres.begin() + kmax - 1);
                        kmax--;
                        std::cout << "kk == kmax && kmax--" << endl;
                    }
                }
            }
            else
            {
                //core[u]改变，core[v]不变
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]不变
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都不变 -- 执行上述循环
                count = d_one_compute_batch(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]不变，core[v]减小
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
    }
    ct += count;

    time = omp_get_wtime() - tm;
    return range;
}


//单次计算 -- 批量更新
int Uncertain_Core::d_one_compute_batch(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //一般只有一个root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //确定thres变化的范围
    if (std::fabs(thres1 - thres2) < EPSILON) {
        root[0] = u;
        up = thres1;
        double kprob1 = delete_compute_low(kthres, k + 1, u);

        root[1] = v;
        double kprob2 = delete_compute_low(kthres, k + 1, v);

        low = min(kprob1, kprob2);
    }
    else if (thres1 > thres2) {
        root[0] = v;
        up = thres2;
        low = delete_compute_low(kthres, k + 1, v);
    }
    else
    {
        root[0] = u;
        up = thres1;
        low = delete_compute_low(kthres, k + 1, u);
    }
    range[0][k] = low;
    range[1][k] = up;

    if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
        return count;
    }
    //cout << "k:" << k << " low:" << low << " up:" << up << endl;

    //寻找候选集
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);
    //std::cout << "count:" << count << endl;

    //更新候选集的thres：在k = kk-1时，维护coreness + kmax
    if (count != 0) {
        d_batch_update_candidate(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


void Uncertain_Core::d_batch_update_candidate(vector<double>& kthres, int& k, int& count, double low, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    //std::cout << "---------k:" << k << " candi count:" << count << endl;
    std::queue<int> que;    //利用que批量更新 p < curThres 的点
    double curThres = low;
    int flag = count;
    int d, index, minV; //邻居索引、kprob最小点的index、id、que中数量的大小
    double p;
    //std::unordered_set<int> neiSet; //待更新顶点
    vec_i record(count, 0); //记录候选点是否需要更新
    vec_i updateNum(count);
    int rd;
    //double time = 0;

    while (flag) {
        while (!mySet.empty()) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            minV = minValueNode->id;
            p = minValueNode->value;

            if (p - curThres > EPSILON) {
                break;
            }
            else
            {
                que.push(minV);
                mySet.erase(it);
            }
        }


        if (!que.empty()) { //thres即为curThres
            rd = 0;
            //neiSet.clear();
            while (!que.empty()) {
                int node = que.front();
                que.pop();

                kthres[node] = curThres;
                visited[node] = false;
                hashTable.erase(node);
                candiRIndex[node] = -1;
                flag--;

                //需要更新的邻居
                d = deg[node];
                for (int t = 0; t < d; t++) {
                    int w = adj[node][t].u;
                    index = candiRIndex[w];
                    //w在candiSet中才需要更新
                    if (index != -1) {
                        neiNum[w]--;
                        //neiSet.insert(w);
                        if (!record[index]) {
                            record[index] = 1;
                            updateNum[rd] = w;
                            rd++;
                        }
                    }
                }
            }

            for (int i = 0; i < rd; i++) {
                int w = updateNum[i];
                if (visited[w]) { 
                    Node* neiNode = hashTable[w];
                    auto it = mySet.find(neiNode);
                    mySet.erase(it);
                    if (neiNum[w] < k + 1) {
                        //que.push(w);    //更新时维护que
                        neiNode->value = 0;
                    }
                    else
                    {
                        p = kprob_comp(w, visited, k + 1);
                        neiNode->value = p;                        

                        //if ((p - curThres) < EPSILON) {  //更新时维护que
                        //    que.push(w);
                        //}
                        //else
                        //{
                        //    neiNode->value = p;
                        //    mySet.insert(neiNode);
                        //    hashTable[w] = neiNode;
                        //}
                    }
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                    index = candiRIndex[w];                    
                    record[index] = 0;
                }
            }
        }
        else
        {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            minV = minValueNode->id;
            p = minValueNode->value;

            curThres = max(curThres, p);
            kthres[minV] = curThres;
            visited[minV] = false;
            candiRIndex[minV] = -1;
            flag--;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            
            //更新其邻居
            d = deg[minV];
            for (int t = 0; t < d; t++) {
                int w = adj[minV][t].u;
                index = candiRIndex[w];
                //w在candiSet中才需要更新
                if (index != -1) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    neiNum[w]--;
                    if (neiNum[w] < k + 1) {
                        neiNode->value = 0;
                        //que.push(w);    //更新时维护que
                    }
                    else
                    {
                        p = kprob_comp(w, visited, k + 1);
                        neiNode->value = p;

                        //if ((p - curThres) < EPSILON) {  //更新时维护que
                        //    que.push(w);
                        //}
                        //else
                        //{
                        //    neiNode->value = p;
                        //    mySet.insert(neiNode);
                        //    hashTable[w] = neiNode;
                        //}
                    }
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;                    
                    //cout << "nei:" << w << " kprob:" << candiProb[index] << endl;
                }
            }
        }
    }

    for (auto node : mySet) {
        delete node;
    }
    mySet.clear();

    for (auto pair : hashTable) {
        delete pair.second;
    }
    hashTable.clear();

    //std::cout << "find-time:" << time << endl;
}



//三种优化的组合：动态更新range + 收缩low + 批量处理
vector<vector<double> > Uncertain_Core::delete_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "删除边 & 根据原始图计算上下界——————" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //候选对象总数量
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk - 1; k++) {
        //std::cout << "k:" << k << endl;
        count = d_one_compute_opt(range, thres[k], u, v, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
    }

    int curK = kk - 1;
    int countNum;   //候选点数量
    //k = kk-1：确定图core维护算法 + 寻找候选集 + 更新thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //确定图core维护 -- kmax必然不变
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量: core变化 + v
            if (count) {
                d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_opt(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else if (core[u] > core[v]) {
        //std::cout << "core[u] > core[v]" << endl;
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core变小 && thres=0：在candiNode中 且node_visited=false
        if (count) {    //存在点v core减小
            //寻找候选集 -- 考虑所有core减小的点，up=thres[v]实时变化
            count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
            if (count) {
                d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
        }
        else
        {
            count = d_one_compute_opt(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }
    else
    {
        //std::cout << "core[u] = core[v]" << endl;
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);
        if (count) {    //core[u]变小 
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]、core[v]都变
                    for (int i = 0; i < countNumV; i++) {
                        int node = candiNode_2[i];
                        if (!visited_2[node]) {
                            candiNode[countNum] = node;
                            node_visited[node] = false;
                            countNum++;
                        }
                    }
                }
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都变 -- 没有root
                count = d_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited已考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //维护kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax查询顶点数：" << t << endl;
                            break;
                        }
                    }
                    if (eraseflag) {
                        thres.erase(thres.begin() + kmax - 1);
                        kmax--;
                        std::cout << "kk == kmax && kmax--" << endl;
                    }
                }
            }
            else
            {
                //core[u]改变，core[v]不变
                count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]不变
            if (!node_visited[v] && (core[v] == kk)) {  //u、v为两个不同的连通分量
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]、core[v]都不变 -- 执行上述循环
                count = d_one_compute_opt(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]不变，core[v]减小
                count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //同时考虑两个连通分量
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
    }
    ct += count;

    time = omp_get_wtime() - tm;
    return range;
}



//优化组合
int Uncertain_Core::d_one_compute_opt(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //一般只有一个root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //确定thres变化的范围
    if (std::fabs(thres1 - thres2) < EPSILON) {
        root[0] = u;
        up = thres1;
        double kprob1 = d_shrink_lower_bound(kthres, k + 1, u);

        root[1] = v;
        double kprob2 = d_shrink_lower_bound(kthres, k + 1, v);

        low = min(kprob1, kprob2);
    }
    else if (thres1 > thres2) {
        root[0] = v;
        up = thres2;
        low = d_shrink_lower_bound(kthres, k + 1, v);
    }
    else
    {
        root[0] = u;
        up = thres1;
        low = d_shrink_lower_bound(kthres, k + 1, u);
    }
    range[0][k] = low;
    range[1][k] = up;

    if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
        return count;
    }
    //cout << "k: " << k << " low:" << low << " up:" << up << endl;

    //寻找候选集
    count = d_candi_dynamic_update_range(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //更新候选集的thres：在k = kk-1时，维护coreness + kmax
    if (count != 0) {
        d_batch_update_candidate(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}
