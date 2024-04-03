#include "uncertain-core.h"

//已知插入边端点u、v，计算上下界low、upper、根节点root
void Uncertain_Core::insert_threshold_singlepoint_compute(vector<double>& kthres, double& low, double& upper, int& root, int u, int v, int k) {
    double thres1 = kthres[u];
    double thres2 = kthres[v];
    //std::cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

    //确定range
    if (std::fabs(thres1 - thres2) < EPSILON) {
        low = thres1;
        int d1 = deg[u];
        vec_b visited(d1, true);
        for (int t = 0; t < d1; t++) {
            int w = adj[u][t].u;
            double th = kthres[w];
            if ((thres1 - th) > EPSILON) {  //th < thres1
                //cout << "neibor:" << w << " thres:" << th << endl;
                visited[t] = false;
            }
        }
        double kprob1 = kprob_comp_scale(u, visited, k + 1);

        int d2 = deg[v];
        vec_b visited_2(d2, true);
        for (int t = 0; t < d2; t++) {
            int w = adj[v][t].u;
            double th = kthres[w];
            if ((thres1 - th) > EPSILON) {
                //cout << "neibor:" << w << " thres:" << th << endl;
                visited_2[t] = false;
            }
        }
        double kprob2 = kprob_comp_scale(v, visited_2, k + 1);

        if ((kprob2 - kprob1) > EPSILON) {  //省略一次kprob的初始计算
            root = u;
            upper = kprob1;
        }
        else
        {
            root = v;
            upper = kprob2;
        }
    }
    else if (thres1 < thres2) {
        root = u;
        low = thres1;
        int d = deg[u];
        vec_b visited(d, true);     //初始化threshold不小于root的邻居
        for (int t = 0; t < d; t++) {
            int w = adj[u][t].u;
            double th = kthres[w];
            if ((thres1 - th) > EPSILON) {
                //cout << "neibor:" << w << " thres:" << th << endl;
                visited[t] = false;
            }
        }
        upper = kprob_comp_scale(u, visited, k + 1);
    }
    else
    {
        root = v;
        low = thres2;
        int d = deg[v];
        vec_b visited(d, true);
        for (int t = 0; t < d; t++) {
            int w = adj[v][t].u;
            double th = kthres[w];
            if ((thres2 - th) > EPSILON) {
                //cout << "neibor:" << w << " thres:" << th << endl;
                visited[t] = false;
            }
        }
        upper = kprob_comp_scale(v, visited, k + 1);
    }

    //std::cout << "lower:" << low << " upper:" << upper << endl;
}


//inset寻找候选集
int Uncertain_Core::insert_search_candi(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //初始化需要访问的节点
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) >= -EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;    //候选处理que
    vec_b visitedCandi(n, false);
    que.push(root);
    visitedCandi[root] = true;
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        if (node == root) {
            kprob = upper;
        }
        else
        {
            kprob = kprob_comp(node, visited, k + 1);
            //cout << "node:" << node << " kprob:" << candiProb[count] << endl;
        }
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;
        candiRIndex[node] = count;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - low) >= -EPSILON) {
                neiNum[node]++;
                if (!visitedCandi[w] && ((upper - kthres[w]) > EPSILON)) {
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "候选对象：" << w << endl;
                }
            }
        }
        count++;
    }
    return count;
}


//更新候选集的thres
void Uncertain_Core::insert_update_candidate_thres(vector<double>& kthres, int& k, int& count,vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int index;
    double curThres = 0.0;
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


//根据确定图core算法寻找子集 -- 同时修改core
int Uncertain_Core::insert_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum) {
    int count = 0;  //core变化的节点数
    vec_i cd(n, 0);  //每个节点的cd  
    std::fill(color.begin(), color.end(), false);   //标记core增大的节点
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);

    //寻找候选顶点集
    std::queue<int> que;
    vec_b visited(n, false);    //该节点是否被访问
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
                else if (cd[w] > k || cd[w] == 0) { //visited[w]:已经访问计数||放入队列未计数
                    cd[r]++;
                }
            }
        }
        //cout << "cd[r]:" << cd[r] << endl;
        if (cd[r] > k) {
            candiNode[count] = r;
            candiRIndex[r] = count;
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

    //确认候选顶点是否被color
    int flag = 1;
    countNum = count;
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
                break;  //一旦修改count，则进入下一个while循环
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
            if (color[v]) {
                core[v]++;
                //cout << "v:" << v << " core:" << core[v] << endl;
            }
        }
    }
    return count;
}


//当k = mincore + 1，inset寻找候选集，low必然为0 -- core已维护（visited=true改变）
int Uncertain_Core::insert_search_candi_core_change(vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    int root;
    double upper, kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    //初始化需要访问的节点 -- thres>=low || core变化的顶点(初始thres=0)
    for (int i = 0; i < n; i++) {
        if (core[i] > 0) {
            visited[i] = true;
        }
    }

    //确定root、upper，其中low必然为0
    if (std::fabs(kthres[u] - kthres[v]) < EPSILON) {
        double kprob1 = kprob_comp(u, visited, k + 1);
        double kprob2 = kprob_comp(v, visited, k + 1);

        if ((kprob2 - kprob1) > EPSILON) {
            root = u;
            upper = kprob1;
        }
        else
        {
            root = v;
            upper = kprob2;
        }
    }
    else if (kthres[u] < kthres[v]) {
        root = u;
        upper = kprob_comp(u, visited, k + 1);
    }
    else
    {
        root = v;
        upper = kprob_comp(v, visited, k + 1);
    }

    //寻找候选集
    std::queue<int> que;    //候选处理que
    vec_b visitedCandi(n, false);    //visited包含可能变大的数组
    que.push(root);
    visitedCandi[root] = true;
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        if (node == root) {
            kprob = upper;
        }
        else
        {
            kprob = kprob_comp(node, visited, k + 1);
        }
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;
        candiRIndex[node] = 1;

        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if (visited[w]) {
                neiNum[node]++;
                if (!visitedCandi[w] && ((upper - kthres[w]) > EPSILON)) {
                    que.push(w);
                    visitedCandi[w] = true;
                }
            }
        }
        count++;
    }
    return count;
}


//初始thres、插入边端点、运行时间
vector<vector<double> > Uncertain_Core::insert_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up, kprob;
    int root;

    int count;
    ct = 0;
    int countNum;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += mySet.size();
        //cout << "count:" << count << endl;

        if (count != 0) {
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk：确定图core维护算法 -- 是否有点的coreness改变
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////发生改变的visited标注为true
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "确定图core改变的顶点数：" << count << endl;


    //存在节点coreness变化
    if (kk == kmax) {  //总体core增大 -- 初始图为core增大的顶点
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //初始化kprob
        for (int i = 0; i < countNum; i++) {
            int node = candiNode[i];
            if (node_visited[node]) {
                kprob = kprob_comp(node, node_visited, kk + 1);
                Node* initialNode = new Node(node, kprob);
                mySet.insert(initialNode);
                hashTable[node] = initialNode;

                int d = deg[node];
                for (int j = 0; j < d; j++) {
                    int w = adj[node][j].u;
                    if (node_visited[w]) {
                        node_visited_num[node]++;
                    }
                }
            }
        }

        //只计算core改变的点
        double curThres = 0.0;
        while (count) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int minV = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            thres[kk][minV] = curThres;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            //cout << "minV:" << minV << " thres:" << curThres << endl;

            node_visited[minV] = false;
            count--;
            //更新其邻居
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                if (node_visited[w]) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    node_visited_num[w]--;

                    if (node_visited_num[w] < kk + 1) {
                        kprob = 0;
                    }
                    else
                    {
                        kprob = kprob_comp(w, node_visited, kk + 1);
                    }
                    neiNode->value = kprob;
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                }
            }
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
    else
    {
        //可以直接访问thres，low=0，计算up（必然不为0）
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //更新候选集thres
        if (count != 0) {
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//insert + decomposition
void Uncertain_Core::insert_threshold_recompute(vector<vector<double> >& thres, double& time) {
    double tm = omp_get_wtime();
    get_core();
    //cout << "kmax:" << kmax << endl;
    thres.resize(kmax, vector<double>(n));
    //Initial_threshold_compute(thres);
    Initial_threshold_compute_map(thres);
    time = omp_get_wtime() - tm;
}



//初始图计算的thres，维护图结构+寻找候选集计算+直接重算
void Uncertain_Core::insert_threshold_compare(vector<vector<double> >& thres) {
    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;

    int um = unselected.size();
    cout << "um:" << um << endl;
    for (int i = 670; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        cout << "insert edge: " << i << ":  u:" << u << " v:" << v << endl;

        //维护图结构
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

        //vector<vector<double> > compareThres = thres;
        int mincore = min(core[u], core[v]);
        //vector<vector<double> > range = insert_threshold_compute_candidate(thres, u, v, tm_1, candi_count);     //InsertBase
        //vector<vector<double> > range = insert_threshold_compute_update_range(thres, u, v, tm_1, candi_count);    //-CT
        //vector<vector<double> > range = insert_threshold_compute_restriction_point(thres, u, v, tm_1, candi_count);   //-UPD
        //vector<vector<double> > range = insert_threshold_compute_batchUP(thres, u, v, tm_1, candi_count);     //-BAT
        vector<vector<double> > range = insert_threshold_compute_opt(thres, u, v, tm_1, candi_count);     //-OPT
        candidate_tm += tm_1;
        candi_sum += candi_count;
        //cout << "insert_compute_candidate:" << tm_1 << endl;


        delete[] core;	//将之前的core动态释放   
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);
        recompute_tm += tm_2;
        //cout << "insert_threshold_recompute:" << tm_2 << endl;

        //比较插入边后的变化阈值
        bool compare = compareArraysTrueOrFalse(thres, thres_2);
        if (!compare) {
            cout << "数组是否相同----------------------------------------------：" << compare << endl;
            compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
        } 
        //compareArrays(compareThres, thres_2, range[0], range[1], mincore);

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
    cout << "insert_compute_candidate_time:" << candidate_tm << endl;
    cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
    cout << "candi sum:" << candi_sum << endl;
    unselected.clear();
}



//insert优化方案：0）初始方案；1）动态更新range；2）寻找限制点；3）批量确定候选点的thres；4）优化方案的结合
vector<vector<double> > Uncertain_Core::insert_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up, kprob;
    int root;

    int count;
    ct = 0;
    int countNum;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
        count = i_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk：确定图core维护算法 -- 是否有点的coreness改变
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////发生改变的visited标注为true
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "确定图core改变的顶点数：" << count << endl;


    //存在节点coreness变化
    if (kk == kmax) {  //总体core增大 -- 初始图为core增大的顶点
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //初始化kprob
        for (int i = 0; i < countNum; i++) {
            int node = candiNode[i];
            if (node_visited[node]) {
                kprob = kprob_comp(node, node_visited, kk + 1);
                Node* initialNode = new Node(node, kprob);
                mySet.insert(initialNode);
                hashTable[node] = initialNode;

                int d = deg[node];
                for (int j = 0; j < d; j++) {
                    int w = adj[node][j].u;
                    if (node_visited[w]) {
                        node_visited_num[node]++;
                    }
                }
            }
        }

        //只计算core改变的点
        double curThres = 0.0;
        while (count) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int minV = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            thres[kk][minV] = curThres;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            //cout << "minV:" << minV << " thres:" << curThres << endl;

            node_visited[minV] = false;
            count--;
            //更新其邻居
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                if (node_visited[w]) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    node_visited_num[w]--;

                    if (node_visited_num[w] < kk + 1) {
                        kprob = 0;
                    }
                    else
                    {
                        kprob = kprob_comp(w, node_visited, kk + 1);
                    }
                    neiNode->value = kprob;
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                }
            }
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
    else
    {
        //可以直接访问thres，low=0，计算up（必然不为0）
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //更新候选集thres
        if (count != 0) {
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//动态更新range -- visited = true表示候选子图； 不能一边寻找候选集 一边计算kprob
int Uncertain_Core::i_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double kprob;
    vec_i candiNode(n);
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //初始化th > upper的子图
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - upper) >= -EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;    //候选处理que
    vec_b visitedCandi(n, false);
    que.push(root);
    visitedCandi[root] = true;
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        candiNode[count] = node;    //候选节点
        candiRIndex[node] = count;
        visited[node] = true;
        int d = deg[node];
        for (int i = 0; i < d; i++) {
            int w = adj[node][i].u;
            if ((kthres[w] - kthres[node]) >= -EPSILON) {
                //neiNum[node]++;
                if (!visitedCandi[w] && ((upper - kthres[w]) > EPSILON)) {
                    que.push(w);
                    visitedCandi[w] = true;
                    //cout << "候选对象：" << w << endl;
                }
            }
        }
        count++;
    }

    //初始化candi各顶点的kprob
    for (int i = 0; i < count; i++) {
        int node = candiNode[i];
        if (node == root) {
            kprob = upper;
        }
        else
        {
            kprob = kprob_comp(node, visited, k + 1);
        }
        Node* initialNode = new Node(node, kprob);
        mySet.insert(initialNode);
        hashTable[node] = initialNode;

        int d = deg[node];
        for (int t = 0; t < d; t++) {
            int w = adj[node][t].u;
            if (visited[w]) {
                neiNum[node]++;
            }
        }
    }
    return count;
}



//寻找限制点
vector<vector<double> > Uncertain_Core::insert_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up, kprob;
    int root;

    int count;
    ct = 0;
    int countNum;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);     //记录是否为限制点

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;
        //是否需要判断low、up的大小
        if ((low - up) >= -EPSILON) {
            continue;
        }

        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);   //限制点不包含在visited数组中
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);    //初始值设为low
        }
    }

    //k = kk：确定图core维护算法 -- 是否有点的coreness改变
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////发生改变的visited标注为true
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "确定图core改变的顶点数：" << count << endl;


    //存在节点coreness变化
    if (kk == kmax) {  //总体core增大 -- 初始图为core增大的顶点
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //初始化kprob
        for (int i = 0; i < countNum; i++) {
            int node = candiNode[i];
            if (node_visited[node]) {
                kprob = kprob_comp(node, node_visited, kk + 1);
                Node* initialNode = new Node(node, kprob);
                mySet.insert(initialNode);
                hashTable[node] = initialNode;

                int d = deg[node];
                for (int j = 0; j < d; j++) {
                    int w = adj[node][j].u;
                    if (node_visited[w]) {
                        node_visited_num[node]++;
                    }
                }
            }
        }

        //只计算core改变的点
        double curThres = 0.0;
        while (count) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int minV = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            thres[kk][minV] = curThres;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            //cout << "minV:" << minV << " thres:" << curThres << endl;

            node_visited[minV] = false;
            count--;
            //更新其邻居
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                if (node_visited[w]) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    node_visited_num[w]--;

                    if (node_visited_num[w] < kk + 1) {
                        kprob = 0;
                    }
                    else
                    {
                        kprob = kprob_comp(w, node_visited, kk + 1);
                    }
                    neiNode->value = kprob;
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                }
            }
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
    else
    {
        //可以直接访问thres，low=0，计算up（必然不为0）
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //更新候选集thres
        if (count != 0) {
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


void Uncertain_Core::i_singlepoint_compute_restriction_point(vector<double>& kthres, double& low, double& upper, int& root, int u, int v, int k, vec_b& restRecord) {
    double thres1 = kthres[u];
    double thres2 = kthres[v];
    std::fill(restRecord.begin(), restRecord.end(), false);
    //vec_i restNodes(n, false);
    int ct = 0; //记录限制点的个数
    //std::cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

    //确定range
    if (std::fabs(thres1 - thres2) < EPSILON) {
        low = thres1;
        std::unordered_set<int> testSet;
        for (int t = 0; t < deg[u]; t++) {
            int w = adj[u][t].u;
            if (w != v) {
                if (fabs(kthres[w] - low) < EPSILON) {
                    testSet.insert(w);
                }
            }
        }

        for (int t = 0; t < deg[v]; t++) {
            int w = adj[v][t].u;
            if (w != u) {
                if (fabs(kthres[w] - low) < EPSILON) {
                    testSet.insert(w);
                }
            }
        }
        //std::cout << "限制点需要判断的数量：" << testSet.size() << endl;

        vec_b visited(n, false);    //所有th = low的邻居 都使用该visited数组计算kprob
        for (const auto& node : testSet) {
            for (int t = 0; t < deg[node]; t++) {
                int w = adj[node][t].u;
                if ((kthres[w] - low) >= -EPSILON) {
                    visited[w] = true;
                }
            }
            double kprob = kprob_comp(node, visited, k + 1);
            if (fabs(kprob - low) < EPSILON) {
                //restNodes.push_back(node);
                restRecord[node] = true;
                ct++;
                std::cout << "限制点：" << node << endl;    //遍历过程中 将限制点标记为已访问
            }
        }

        //维护u、v的visited
        for (int t = 0; t < deg[u]; t++) {
            int w = adj[u][t].u;
            if ((kthres[w] - low) >= -EPSILON) {
                visited[w] = true;
            }
            if (restRecord[w]) {
                visited[w] = false;
            }
        }
        for (int t = 0; t < deg[v]; t++) {
            int w = adj[v][t].u;
            if ((kthres[w] - low) >= -EPSILON) {
                visited[w] = true;
            }
            if (restRecord[w]) {
                visited[w] = false;
            }
        }

        //修改限制点的visited
        /*for (const auto& node : restNodes) {
            visited[node] = false;
        }*/

        //计算upper时 不考虑限制点        
        double kprob1 = kprob_comp(u, visited, k + 1);
        double kprob2 = kprob_comp(v, visited, k + 1);

        if ((kprob2 - kprob1) > EPSILON) {  //省略一次kprob的初始计算，upper和root必须相对应
            root = u;
            upper = kprob1;
        }
        else
        {
            root = v;
            upper = kprob2;
        }
    }
    else if (thres1 < thres2) {
        root = u;
        low = thres1;
        int d = deg[u];
        for (int t = 0; t < d; t++) {
            int w = adj[u][t].u;
            if (fabs(kthres[w] - low) < EPSILON) {  //判断是否为限制点
                vec_b visited_w(deg[w], false);
                for (int j = 0; j < deg[w]; j++) {
                    int ww = adj[w][j].u;
                    if ((kthres[ww] - low) >= -EPSILON) {
                        visited_w[j] = true;
                    }
                }
                double kprob = kprob_comp_scale(w, visited_w, k + 1);
                if (fabs(kprob - low) < EPSILON) {
                    //restNodes.push_back(w);
                    restRecord[w] = true;
                    ct++;
                }
            }
        }
        std::cout << "u的限制点数量：" << ct << endl;

        vec_b visited(deg[u], false);
        for (int t = 0; t < deg[u]; t++) {
            int w = adj[u][t].u;
            if (!restRecord[w]) {    //不为限制点
                if ((kthres[w] - low) >= -EPSILON) {
                    visited[t] = true;
                }
            }
        }
        /*for (const auto& node : restNodes) {
            visited[node] = false;
        }*/
        upper = kprob_comp_scale(u, visited, k + 1);
    }
    else
    {
        root = v;
        low = thres2;
        int d = deg[v];
        for (int t = 0; t < d; t++) {
            int w = adj[v][t].u;
            if (fabs(kthres[w] - low) < EPSILON) {  //判断是否为限制点
                vec_b visited_w(deg[w], false);
                for (int j = 0; j < deg[w]; j++) {
                    int ww = adj[w][j].u;
                    if ((kthres[ww] - low) >= -EPSILON) {
                        visited_w[j] = true;
                    }
                }
                double kprob = kprob_comp_scale(w, visited_w, k + 1);
                if (fabs(kprob - low) < EPSILON) {
                    //restNodes.push_back(w);
                    restRecord[w] = true;
                    ct++;
                }
            }
        }
        //std::cout << "v的限制点数量：" << ct << endl;

        vec_b visited(deg[v], false);
        for (int t = 0; t < deg[v]; t++) {
            int w = adj[v][t].u;
            if (!restRecord[w]) {
                if ((kthres[w] - low) >= -EPSILON) {
                    visited[t] = true;
                }
            }
        }
        /*for (const auto& node : restNodes) {
            visited[node] = false;
        }*/
        upper = kprob_comp_scale(v, visited, k + 1);
    }
}



//批量确定候选点的thres
vector<vector<double> > Uncertain_Core::insert_threshold_compute_batchUP(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    std::cout << "kk:" << kk << endl;
    double low, up, kprob;
    int root;

    int count;
    ct = 0;
    int countNum;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);  
            //insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk：确定图core维护算法 -- 是否有点的coreness改变
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////发生改变的visited标注为true
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "确定图core改变的顶点数：" << count << endl;


    //存在节点coreness变化
    if (kk == kmax) {  //总体core增大 -- 初始图为core增大的顶点
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //初始化kprob
        for (int i = 0; i < countNum; i++) {
            int node = candiNode[i];
            if (node_visited[node]) {
                kprob = kprob_comp(node, node_visited, kk + 1);
                Node* initialNode = new Node(node, kprob);
                mySet.insert(initialNode);
                hashTable[node] = initialNode;

                int d = deg[node];
                for (int j = 0; j < d; j++) {
                    int w = adj[node][j].u;
                    if (node_visited[w]) {
                        node_visited_num[node]++;
                    }
                }
            }
        }

        //只计算core改变的点
        double curThres = 0.0;
        while (count) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int minV = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            thres[kk][minV] = curThres;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            //cout << "minV:" << minV << " thres:" << curThres << endl;

            node_visited[minV] = false;
            count--;
            //更新其邻居
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                if (node_visited[w]) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    node_visited_num[w]--;

                    if (node_visited_num[w] < kk + 1) {
                        kprob = 0;
                    }
                    else
                    {
                        kprob = kprob_comp(w, node_visited, kk + 1);
                    }
                    neiNode->value = kprob;
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                }
            }
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
    else
    {
        //可以直接访问thres，low=0，计算up（必然不为0）
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //更新候选集thres
        if (count != 0) {
            //d_batch_update_candidate(thres[kk], kk, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}



//新增边：优化方案的结合
vector<vector<double> > Uncertain_Core::insert_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up, kprob;
    int root;

    int count;
    ct = 0;
    int countNum;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);     //记录是否为限制点

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;
        //是否需要判断low、up的大小
        if ((low - up) >= -EPSILON) {
            continue;
        }


        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
        count = i_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk：确定图core维护算法 -- 是否有点的coreness改变
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////发生改变的visited标注为true
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "确定图core改变的顶点数：" << count << endl;


    //存在节点coreness变化
    if (kk == kmax) {  //总体core增大 -- 初始图为core增大的顶点
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //初始化kprob
        for (int i = 0; i < countNum; i++) {
            int node = candiNode[i];
            if (node_visited[node]) {
                kprob = kprob_comp(node, node_visited, kk + 1);
                Node* initialNode = new Node(node, kprob);
                mySet.insert(initialNode);
                hashTable[node] = initialNode;

                int d = deg[node];
                for (int j = 0; j < d; j++) {
                    int w = adj[node][j].u;
                    if (node_visited[w]) {
                        node_visited_num[node]++;
                    }
                }
            }
        }

        //只计算core改变的点
        double curThres = 0.0;
        while (count) {
            auto it = mySet.begin();
            Node* minValueNode = *(it);
            int minV = minValueNode->id;
            double p = minValueNode->value;

            curThres = max(curThres, p);
            thres[kk][minV] = curThres;
            mySet.erase(it);
            hashTable.erase(minV);
            delete minValueNode;
            //cout << "minV:" << minV << " thres:" << curThres << endl;

            node_visited[minV] = false;
            count--;
            //更新其邻居
            for (int t = 0; t < deg[minV]; t++) {
                int w = adj[minV][t].u;
                if (node_visited[w]) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    node_visited_num[w]--;

                    if (node_visited_num[w] < kk + 1) {
                        kprob = 0;
                    }
                    else
                    {
                        kprob = kprob_comp(w, node_visited, kk + 1);
                    }
                    neiNode->value = kprob;
                    mySet.insert(neiNode);
                    hashTable[w] = neiNode;
                }
            }
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
    else
    {
        //可以直接访问thres，low=0，计算up（必然不为0）
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //更新候选集thres
        if (count != 0) {
            //d_batch_update_candidate(thres[kk], kk, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}



//oldarray、newarray -- 比较上下界
void Uncertain_Core::compareArrays(const std::vector<vector<double> >& array1, const std::vector<vector<double> >& array2, const vector<double>& comp1, const vector<double>& comp2, int mincore) {
    double lower;
    double upper;
    for (int i = 0; i < mincore; i++) {
        //cout << "compare-k:" << i << endl;
        lower = 1.0;
        upper = 0.0;
        for (int j = 0; j < n; j++) {
            if (std::fabs(array1[i][j] - array2[i][j]) >= EPSILON) {
                //cout << "k:" << i << " vertex:" << j << " candi:" << array1[i][j] << " recompute:" << array2[i][j] << endl;
                lower = min(lower, array1[i][j]);
                upper = max(upper, array1[i][j]);
            }
        }
        cout << "compare lower:" << lower << " upper:" << upper << endl;

        //double类型存在精度误差
        if ((comp1[i] - lower) > EPSILON) {
            cout << "Lower Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "实际比较low：" << lower << "单点计算low：" << comp1[i];
        }
        if ((upper - comp2[i]) > EPSILON) {
            cout << "Upper Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "实际比较up：" << upper << "单点计算low：" << comp2[i];
        }
    }
    //记录core新增的点
    /*if (array1.size() != array2.size()) {
        cout << "compare-k:" << kk << endl;
        for (int j = 0; j < n; j++) {
            if (array2[kk][j] != 0) {
                cout << "vertex:" << j << " thres:" << array2[kk][j] << endl;
            }
        }
    }*/
}


bool Uncertain_Core::compareArraysTrueOrFalse(const std::vector<std::vector<double> >& array1, const std::vector<std::vector<double> >& array2) {
    // 检查数组的维度是否相同
    /*if (array1.size() != array2.size()) {
        cout << "数组的维度not相同!" << endl;
        return false;
    }*/

    for (size_t i = 0; i < array2.size(); ++i) {
        // 检查每个子数组的长度是否相同
        if (array1[i].size() != array2[i].size()) {
            cout << "每个子数组的长度not相同!" << endl;
            return false;
        }

        // 比较每个子数组中对应位置的元素
        for (size_t j = 0; j < array1[i].size(); ++j) {
            if (std::fabs(array1[i][j] - array2[i][j]) >= EPSILON1) {
                cout << "i:" << i << " j:" << j << " candidata:" << array1[i][j] << "recompute: " << array2[i][j] << endl;
                return false;
            }
        }
    }

    return true;
}

//比较更改的k个prob是否相等
bool Uncertain_Core::compareKsizeArraysTrueOrFalse(const std::vector<std::vector<double> >& array1, const std::vector<std::vector<double> >& array2, int k) {
    // 比较每个子数组中对应位置的元素    
    bool flag = true;
    for (size_t i = 0; i < k; i++) {
        for (size_t j = 0; j < array1[i].size(); ++j) {
            if (std::fabs(array1[i][j] - array2[i][j]) >= EPSILON1) {
                cout << "k:" << i << " vector:" << j << " candi:" << array1[i][j] << " recompute:" << array2[i][j] << endl;
                flag = false;
            }
        }
    }
    return flag;
}


bool Uncertain_Core::compareCoreArrays(const std::vector<int>& array1, const std::vector<int>& array2) {
    for (int i = 0; i < n; i++) {
        if (array1[i] != array2[i]) {
            cout << "vertex:" << i << " candidate-core:" << array1[i] << " recompute-core:" << array2[i] << endl;
            return false;
        }
    }
    return true;
}