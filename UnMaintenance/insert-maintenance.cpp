#include "uncertain-core.h"

//��֪����߶˵�u��v���������½�low��upper�����ڵ�root
void Uncertain_Core::insert_threshold_singlepoint_compute(vector<double>& kthres, double& low, double& upper, int& root, int u, int v, int k) {
    double thres1 = kthres[u];
    double thres2 = kthres[v];
    //std::cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

    //ȷ��range
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

        if ((kprob2 - kprob1) > EPSILON) {  //ʡ��һ��kprob�ĳ�ʼ����
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
        vec_b visited(d, true);     //��ʼ��threshold��С��root���ھ�
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


//insetѰ�Һ�ѡ��
int Uncertain_Core::insert_search_candi(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //��ʼ����Ҫ���ʵĽڵ�
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) >= -EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;    //��ѡ����que
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
                    //cout << "��ѡ����" << w << endl;
                }
            }
        }
        count++;
    }
    return count;
}


//���º�ѡ����thres
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

        //�������ھ�
        for (int t = 0; t < deg[minV]; t++) {
            int w = adj[minV][t].u;
            index = candiRIndex[w];

            //w��candiSet�в���Ҫ����
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


//����ȷ��ͼcore�㷨Ѱ���Ӽ� -- ͬʱ�޸�core
int Uncertain_Core::insert_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum) {
    int count = 0;  //core�仯�Ľڵ���
    vec_i cd(n, 0);  //ÿ���ڵ��cd  
    std::fill(color.begin(), color.end(), false);   //���core����Ľڵ�
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);

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

    //ȷ�Ϻ�ѡ�����Ƿ�color
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


//��k = mincore + 1��insetѰ�Һ�ѡ����low��ȻΪ0 -- core��ά����visited=true�ı䣩
int Uncertain_Core::insert_search_candi_core_change(vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    int root;
    double upper, kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    //��ʼ����Ҫ���ʵĽڵ� -- thres>=low || core�仯�Ķ���(��ʼthres=0)
    for (int i = 0; i < n; i++) {
        if (core[i] > 0) {
            visited[i] = true;
        }
    }

    //ȷ��root��upper������low��ȻΪ0
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

    //Ѱ�Һ�ѡ��
    std::queue<int> que;    //��ѡ����que
    vec_b visitedCandi(n, false);    //visited�������ܱ�������
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


//��ʼthres������߶˵㡢����ʱ��
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
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += mySet.size();
        //cout << "count:" << count << endl;

        if (count != 0) {
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk��ȷ��ͼcoreά���㷨 -- �Ƿ��е��coreness�ı�
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////�����ı��visited��עΪtrue
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "ȷ��ͼcore�ı�Ķ�������" << count << endl;


    //���ڽڵ�coreness�仯
    if (kk == kmax) {  //����core���� -- ��ʼͼΪcore����Ķ���
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //��ʼ��kprob
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

        //ֻ����core�ı�ĵ�
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
            //�������ھ�
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

        // �ͷ�ʣ��ڵ���ڴ�
        for (auto node : mySet) {
            delete node;
        }
        mySet.clear();

        // ��չ�ϣ��
        for (auto pair : hashTable) {
            delete pair.second;
        }
        hashTable.clear();
    }
    else
    {
        //����ֱ�ӷ���thres��low=0������up����Ȼ��Ϊ0��
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //���º�ѡ��thres
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



//��ʼͼ�����thres��ά��ͼ�ṹ+Ѱ�Һ�ѡ������+ֱ������
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


        delete[] core;	//��֮ǰ��core��̬�ͷ�   
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);
        recompute_tm += tm_2;
        //cout << "insert_threshold_recompute:" << tm_2 << endl;

        //�Ƚϲ���ߺ�ı仯��ֵ
        bool compare = compareArraysTrueOrFalse(thres, thres_2);
        if (!compare) {
            cout << "�����Ƿ���ͬ----------------------------------------------��" << compare << endl;
            compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
        } 
        //compareArrays(compareThres, thres_2, range[0], range[1], mincore);

        //��newArray���Ƹ�thres
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



//insert�Ż�������0����ʼ������1����̬����range��2��Ѱ�����Ƶ㣻3������ȷ����ѡ���thres��4���Ż������Ľ��
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
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = i_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk��ȷ��ͼcoreά���㷨 -- �Ƿ��е��coreness�ı�
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////�����ı��visited��עΪtrue
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "ȷ��ͼcore�ı�Ķ�������" << count << endl;


    //���ڽڵ�coreness�仯
    if (kk == kmax) {  //����core���� -- ��ʼͼΪcore����Ķ���
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //��ʼ��kprob
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

        //ֻ����core�ı�ĵ�
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
            //�������ھ�
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

        // �ͷ�ʣ��ڵ���ڴ�
        for (auto node : mySet) {
            delete node;
        }
        mySet.clear();

        // ��չ�ϣ��
        for (auto pair : hashTable) {
            delete pair.second;
        }
        hashTable.clear();
    }
    else
    {
        //����ֱ�ӷ���thres��low=0������up����Ȼ��Ϊ0��
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //���º�ѡ��thres
        if (count != 0) {
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}


//��̬����range -- visited = true��ʾ��ѡ��ͼ�� ����һ��Ѱ�Һ�ѡ�� һ�߼���kprob
int Uncertain_Core::i_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, int root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double kprob;
    vec_i candiNode(n);
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //��ʼ��th > upper����ͼ
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - upper) >= -EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;    //��ѡ����que
    vec_b visitedCandi(n, false);
    que.push(root);
    visitedCandi[root] = true;
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        candiNode[count] = node;    //��ѡ�ڵ�
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
                    //cout << "��ѡ����" << w << endl;
                }
            }
        }
        count++;
    }

    //��ʼ��candi�������kprob
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



//Ѱ�����Ƶ�
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
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);     //��¼�Ƿ�Ϊ���Ƶ�

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;
        //�Ƿ���Ҫ�ж�low��up�Ĵ�С
        if ((low - up) >= -EPSILON) {
            continue;
        }

        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);   //���Ƶ㲻������visited������
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);    //��ʼֵ��Ϊlow
        }
    }

    //k = kk��ȷ��ͼcoreά���㷨 -- �Ƿ��е��coreness�ı�
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////�����ı��visited��עΪtrue
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "ȷ��ͼcore�ı�Ķ�������" << count << endl;


    //���ڽڵ�coreness�仯
    if (kk == kmax) {  //����core���� -- ��ʼͼΪcore����Ķ���
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //��ʼ��kprob
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

        //ֻ����core�ı�ĵ�
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
            //�������ھ�
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

        // �ͷ�ʣ��ڵ���ڴ�
        for (auto node : mySet) {
            delete node;
        }
        mySet.clear();

        // ��չ�ϣ��
        for (auto pair : hashTable) {
            delete pair.second;
        }
        hashTable.clear();
    }
    else
    {
        //����ֱ�ӷ���thres��low=0������up����Ȼ��Ϊ0��
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //���º�ѡ��thres
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
    int ct = 0; //��¼���Ƶ�ĸ���
    //std::cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

    //ȷ��range
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
        //std::cout << "���Ƶ���Ҫ�жϵ�������" << testSet.size() << endl;

        vec_b visited(n, false);    //����th = low���ھ� ��ʹ�ø�visited�������kprob
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
                std::cout << "���Ƶ㣺" << node << endl;    //���������� �����Ƶ���Ϊ�ѷ���
            }
        }

        //ά��u��v��visited
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

        //�޸����Ƶ��visited
        /*for (const auto& node : restNodes) {
            visited[node] = false;
        }*/

        //����upperʱ ���������Ƶ�        
        double kprob1 = kprob_comp(u, visited, k + 1);
        double kprob2 = kprob_comp(v, visited, k + 1);

        if ((kprob2 - kprob1) > EPSILON) {  //ʡ��һ��kprob�ĳ�ʼ���㣬upper��root�������Ӧ
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
            if (fabs(kthres[w] - low) < EPSILON) {  //�ж��Ƿ�Ϊ���Ƶ�
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
        std::cout << "u�����Ƶ�������" << ct << endl;

        vec_b visited(deg[u], false);
        for (int t = 0; t < deg[u]; t++) {
            int w = adj[u][t].u;
            if (!restRecord[w]) {    //��Ϊ���Ƶ�
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
            if (fabs(kthres[w] - low) < EPSILON) {  //�ж��Ƿ�Ϊ���Ƶ�
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
        //std::cout << "v�����Ƶ�������" << ct << endl;

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



//����ȷ����ѡ���thres
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
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        insert_threshold_singlepoint_compute(thres[k], low, up, root, u, v, k);
        range[0][k] = low;
        range[1][k] = up;


        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);  
            //insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk��ȷ��ͼcoreά���㷨 -- �Ƿ��е��coreness�ı�
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////�����ı��visited��עΪtrue
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "ȷ��ͼcore�ı�Ķ�������" << count << endl;


    //���ڽڵ�coreness�仯
    if (kk == kmax) {  //����core���� -- ��ʼͼΪcore����Ķ���
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //��ʼ��kprob
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

        //ֻ����core�ı�ĵ�
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
            //�������ھ�
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

        // �ͷ�ʣ��ڵ���ڴ�
        for (auto node : mySet) {
            delete node;
        }
        mySet.clear();

        // ��չ�ϣ��
        for (auto pair : hashTable) {
            delete pair.second;
        }
        hashTable.clear();

    }
    else
    {
        //����ֱ�ӷ���thres��low=0������up����Ȼ��Ϊ0��
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //���º�ѡ��thres
        if (count != 0) {
            //d_batch_update_candidate(thres[kk], kk, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}



//�����ߣ��Ż������Ľ��
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
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);     //��¼�Ƿ�Ϊ���Ƶ�

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;
        //�Ƿ���Ҫ�ж�low��up�Ĵ�С
        if ((low - up) >= -EPSILON) {
            continue;
        }


        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = i_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }


    //k = kk��ȷ��ͼcoreά���㷨 -- �Ƿ��е��coreness�ı�
    root = u;
    if (core[v] < core[u]) {
        root = v;
    }
    count = insert_find_core_subcore(root, kk, node_visited, candiNode, candiReverseIndex, countNum);   ////�����ı��visited��עΪtrue
    if (count == 0) {
        time = omp_get_wtime() - tm;
        return range;
    }
    //cout << "ȷ��ͼcore�ı�Ķ�������" << count << endl;


    //���ڽڵ�coreness�仯
    if (kk == kmax) {  //����core���� -- ��ʼͼΪcore����Ķ���
        kmax++;
        cout << "kmax++; " << kmax << endl;
        ct += count;
        vector<double> newRow(n, 0.0);
        thres.push_back(newRow);

        std::fill(node_visited_num.begin(), node_visited_num.end(), 0);

        //��ʼ��kprob
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

        //ֻ����core�ı�ĵ�
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
            //�������ھ�
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

        // �ͷ�ʣ��ڵ���ڴ�
        for (auto node : mySet) {
            delete node;
        }
        mySet.clear();

        // ��չ�ϣ��
        for (auto pair : hashTable) {
            delete pair.second;
        }
        hashTable.clear();
    }
    else
    {
        //����ֱ�ӷ���thres��low=0������up����Ȼ��Ϊ0��
        count = insert_search_candi_core_change(thres[kk], u, v, kk, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "kmin+1 - count:" << count << endl;

        //���º�ѡ��thres
        if (count != 0) {
            //d_batch_update_candidate(thres[kk], kk, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            insert_update_candidate_thres(thres[kk], kk, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    return range;
}



//oldarray��newarray -- �Ƚ����½�
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

        //double���ʹ��ھ������
        if ((comp1[i] - lower) > EPSILON) {
            cout << "Lower Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "ʵ�ʱȽ�low��" << lower << "�������low��" << comp1[i];
        }
        if ((upper - comp2[i]) > EPSILON) {
            cout << "Upper Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            cout << "ʵ�ʱȽ�up��" << upper << "�������low��" << comp2[i];
        }
    }
    //��¼core�����ĵ�
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
    // ��������ά���Ƿ���ͬ
    /*if (array1.size() != array2.size()) {
        cout << "�����ά��not��ͬ!" << endl;
        return false;
    }*/

    for (size_t i = 0; i < array2.size(); ++i) {
        // ���ÿ��������ĳ����Ƿ���ͬ
        if (array1[i].size() != array2[i].size()) {
            cout << "ÿ��������ĳ���not��ͬ!" << endl;
            return false;
        }

        // �Ƚ�ÿ���������ж�Ӧλ�õ�Ԫ��
        for (size_t j = 0; j < array1[i].size(); ++j) {
            if (std::fabs(array1[i][j] - array2[i][j]) >= EPSILON1) {
                cout << "i:" << i << " j:" << j << " candidata:" << array1[i][j] << "recompute: " << array2[i][j] << endl;
                return false;
            }
        }
    }

    return true;
}

//�Ƚϸ��ĵ�k��prob�Ƿ����
bool Uncertain_Core::compareKsizeArraysTrueOrFalse(const std::vector<std::vector<double> >& array1, const std::vector<std::vector<double> >& array2, int k) {
    // �Ƚ�ÿ���������ж�Ӧλ�õ�Ԫ��    
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