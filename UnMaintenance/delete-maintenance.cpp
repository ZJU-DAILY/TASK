#include "uncertain-core.h"

//��m�����б�Ǵ�ɾ���ı�
void Uncertain_Core::selectRandomEdges(double scale) {
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);
    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //ȷ��ÿ����ֻ����һ��
                if (rand() > scale_size) {  //��ɾ��
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



//����ɾ��ct���� then����ct����
void Uncertain_Core::delete_threshold_compare(vector<vector<double> >& thres, double scale) {
    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //ÿ������Ѱ�Һ�ѡ����������ۼƺ�

    selectRandomEdges(scale);
    vector<vector<double> > copy_thres = thres;

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "i:" << i << " u:" << u << " v:" << v << endl;

        //ά��ͼ�ṹ
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

        delete[] core;	//��֮ǰ��core��̬�ͷ�   
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);
        recompute_tm += tm_2;
        //cout << "insert_threshold_recompute:" << tm_2 << endl;

        //�Ƚϲ���ߺ�ı仯��ֵ
        bool compare = compareArraysTrueOrFalse(thres, thres_2);
        if (!compare) {
            cout << "�����Ƿ���ͬ----------------------------------------------��" << compare << endl;
            //compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
        }


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
}





//ɾ���߼���low
double Uncertain_Core::delete_compute_low(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //ֱ���ھ����ж�
    if (d < k) {
        low = 0;
    }
    else
    {
        int count = 0;  //��¼kprob��С��root���ھ�����
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
    //��ʼ����Ҫ���ʵĽڵ�
    for (int i = 0; i < n; i++) {
        if ((kthres[i] - low) > EPSILON) {
            visited[i] = true;
        }
    }

    std::queue<int> que;
    vec_b visitedCandi(n, false);
    que.push(root[0]);
    visitedCandi[root[0]] = true;
    if (root[1] != -1) {    //����ڶ���root
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


//���º�ѡ����thres���ҳ�ʼֵ����Ϊlow
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



//ɾ���ߺ� ��ÿ��k�������� -- core���ı�������
int Uncertain_Core::delete_one_compute_candidate(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int& k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //һ��ֻ��һ��root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //ȷ��thres�仯�ķ�Χ
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

    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        return count;
    }

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


//ɾ���� Ѱ��core��С�ĵ� -- ��core��С�ĵ����η���que��Ȼ��Ѱ��thres���ܸı�ĺ�ѡ��
int Uncertain_Core::delete_find_core_subcore(int root, int k, vec_b& color, vec_i& candiNode, vec_i& candiRIndex, int& countNum) {
    int count = 0;
    int xdeg = 0;   //root��core�Ƿ����
    vec_i cd(n, 0);
    std::fill(color.begin(), color.end(), false);   //���core����Ľڵ�
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
    if (xdeg < k) { //root��corenessȷ����С
        que.push(root);
        visited[root] = true;
    }
    while (!que.empty()) {
        int r = que.front();
        que.pop();
        candiNode[count] = r;
        candiRIndex[r] = count;
        count++;    //��root��ͨ&core=k����ͼ�ڵ�����
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

    //ȷ����ѡ�����Ƿ�color
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
                    if (color[w] && (cd[w] >= k)) { //�ں�ѡ����&&�ٴ��ж�
                        cd[w]--;
                    }
                }
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
            if (!color[v]) {
                core[v]--;
            }
        }
    }

    return count;
}


//core�仯�� ά��thres=0 Ѱ�ҿ��ܸı�ĺ�ѡ��
int Uncertain_Core::delete_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double up = 0;  //que��thres���ֵ
    double kprob;
    std::queue<int> que;
    vec_b visitedCandi(n, false);
    //��core�ı�Ӱ��ĵ����que����ά����thres -- ��ʼ��que��һ����ɢ
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
            //ά��thres
            kthres[r] = 0;
        }
    }

    //��ʼ����ѡ��¼
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
            que.push(root); //���������������ͨ����
            visitedCandi[root] = true;
            up = max(up, kthres[root]);
        }
    }
    while (!que.empty()) {
        int node = que.front();
        que.pop();

        kprob = kprob_comp(node, visited, k + 1);    //��ʼ�ھ�������С��k��������=0   
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


//ɾ���ߺ� ���ݺ�ѡ�����и���
vector<vector<double> > Uncertain_Core::delete_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "ɾ���� & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //��ѡ����������
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
    int countNum;   //��ѡ������
    //k = kk-1��ȷ��ͼcoreά���㷨 + Ѱ�Һ�ѡ�� + ����thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //ȷ��ͼcoreά�� -- kmax��Ȼ����
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����: core�仯 + v
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
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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
        if (count) {    //core[u]��С 
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]��core[v]����
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
                //core[u]��core[v]���� -- û��root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited�ѿ���������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //ά��kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax��ѯ��������" << t << endl;
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
                //core[u]�ı䣬core[v]����
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]����
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]��core[v]������ -- ִ������ѭ��
                count = delete_one_compute_candidate(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]���䣬core[v]��С
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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



//delete�Ż�������0)��ʼ������1����̬����range��2�����ݷ������½磻3������ȷ����ѡ���thres��4���Ż������Ľ��
vector<vector<double> > Uncertain_Core::delete_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "candidate--kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));    //������range��Ϊ����ֵ

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
    int countNum;   //��ѡ������
    //k = kk-1��ȷ��ͼcoreά���㷨 + Ѱ�Һ�ѡ�� + ����thres
    if (core[u] < core[v]) {
        //ȷ��ͼcoreά�� -- kmax��Ȼ����
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����: core�仯 + v
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
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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
        if (count) {    //core[u]��С 
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]��core[v]����
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
                //core[u]��core[v]���� -- û��root
                count = d_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited�ѿ���������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //ά��kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax��ѯ��������" << t << endl;
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
                //core[u]�ı䣬core[v]����
                count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]����
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]��core[v]������ -- ִ������ѭ��
                count = d_one_compute_update_range(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]���䣬core[v]��С
                count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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

//���μ���--ʵʱ����range���½�
int Uncertain_Core::d_one_compute_update_range(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //һ��ֻ��һ��root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //ȷ��thres�仯�ķ�Χ
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

    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        return count;
    }

    //Ѱ�Һ�ѡ��
    count = d_candi_dynamic_update_range(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


//Ѱ�Һ�ѡ��ʱ ��̬�����߽�
int Uncertain_Core::d_candi_dynamic_update_range(vector<double>& kthres, double low, double upper, vector<int>& root, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    int count = 0;
    double tempUP, kprob;
    std::fill(candiRIndex.begin(), candiRIndex.end(), -1);
    std::fill(neiNum.begin(), neiNum.end(), 0);
    std::fill(visited.begin(), visited.end(), false);
    //��ʼ����Ҫ���ʵĽڵ�
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

        if (kthres[node] < upper) { //upper�ϸ��С
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
                //if (!visitedCandi[w] && (kthres[w] - kthres[node]) < EPSILON) {     //��������upper
                //    que.push(w);
                //    visitedCandi[w] = true;
                //}

                if (!visitedCandi[w] && (kthres[w] - tempUP) < EPSILON) {     //��������upper
                    que.push(w);
                    visitedCandi[w] = true;
                }
            }
        }
    }

    return count;
}

//kk����������ͨ�ṹ -- ʵʱ����range
int Uncertain_Core::d_search_candi_core_change(vector<double>& kthres, int root, int k, vec_b& visited, vec_i& candiNode, vec_i& candiRIndex, vec_i& neiNum, int countNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double up = 0;
    double kprob;
    std::queue<int> que;
    vec_b visitedCandi(n, false);
    //��core�ı�Ӱ��ĵ����que����ά����thres -- ��ʼ��que��һ����ɢ
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
            //ά��thres
            kthres[r] = 0;
        }
    }

    //��ʼ����ѡ��¼
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
            que.push(root); //���������������ͨ����
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

//ɾ���� �Ƚ��������������low
void Uncertain_Core::delete_compare_range(vector<vector<double> >& thres, double scale) {
    double tm_1;
    double tm_2;
    int count1, count2;
    selectRandomEdges(scale);

    int um = unselected.size();
    for (int i = 0; i < um; i++) {
        int u = unselected[i].first.first;
        int v = unselected[i].first.second;
        std::cout << "��i���ߣ�" << i << " u:" << u << " v:" << v << endl;

        //ά��ͼ�ṹ
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
        //��ʼ�����½�
        vector<vector<double> > range1 = delete_threshold_compute_candidate(thres, u, v, tm_1, count1);
        //���ݷ���ȡ�½�
        vector<vector<double> > range2 = delete_threshold_compute_shrink_low(thres, u, v, tm_1, count2);

        delete[] core;
        vector<vector<double> > thres_2;
        insert_threshold_recompute(thres_2, tm_2);

        //�Ƚ�ɾ���ߺ�ı仯��ֵ
        compareArrays(thres, thres_2, range2[0], range2[1], mincore - 1);   //���ܱȽ�core��С�ĵ�

        for (int i = 0; i < range1[0].size(); i++) {
            if ((range2[0][i] - range1[0][i]) > EPSILON) {
                std::cout << "low���� ------- k��" << i;
                std::cout << "  initial��" << range1[0][i] << " opt:" << range2[0][i] << endl;
            }
        }

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
    unselected.clear();
}


//�Ż����� -- ���ݷ���ȡ�½�
vector<vector<double> > Uncertain_Core::delete_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "ɾ���� & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //��ѡ����������
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
    int countNum;   //��ѡ������
    //k = kk-1��ȷ��ͼcoreά���㷨 + Ѱ�Һ�ѡ�� + ����thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //ȷ��ͼcoreά�� -- kmax��Ȼ����
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����: core�仯 + v
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
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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
        if (count) {    //core[u]��С 
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]��core[v]����
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
                //core[u]��core[v]���� -- û��root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited�ѿ���������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //ά��kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax��ѯ��������" << t << endl;
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
                //core[u]�ı䣬core[v]����
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
                if (count) {
                    delete_update_candidate_thres(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]����
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]��core[v]������ -- ִ������ѭ��
                count = d_one_compute_shrink_low(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]���䣬core[v]��С
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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

//���μ��� -- �������½�
int Uncertain_Core::d_one_compute_shrink_low(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //һ��ֻ��һ��root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //ȷ��thres�仯�ķ�Χ
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

    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        return count;
    }

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
    if (count != 0) {
        delete_update_candidate_thres(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}



//��һ�������½� -- ��Ϊopt1��opt2 ���ߵ���ʵ�� �Ż����ܲ�ͬ
double Uncertain_Core::d_shrink_lower_bound(vector<double>& kthres, int k, int root) {
    double low;
    double upper = kthres[root];
    int d = deg[root];  //ֱ���ھ����ж�
    if (d < k) {
        low = 0;
    }
    else
    {
        vec_d kprobKqualK(k + 1, 0.0);  //��ʼkprobֵ
        kprobKqualK[1] = 1;
        vec_d newp(k + 1, 0.0);
        int end;
        double p, kprob;

        int coreD = 0;  //����core��С��k���ھ�����
        int count = 0;  //��¼kprob��С��root���ھ�����
        vector<pair<double, int> > unvisitedNei; //��¼kprobС��root���ھ�
        for (int t = 0; t < d; t++) {
            int w = adj[root][t].u;
            double th = kthres[w];
            if ((th - upper) >= -EPSILON) {     //th >= upper                
                coreD++;
                count++;    //�����ھ�����
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
        else if (count >= k) {  //�����ھ������½�
            low = 1;
            for (int j = 1; j <= end; j++) {
                low = low - kprobKqualK[j]; //��ʼlow
            }

            //if (coreD > count) {    //ȷ��unvisitedNei�����ھ� -- opt1
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

            //            if ((kprob - kthres[w]) >= -EPSILON) {    //ȷ��lowȡֵ
            //                low = kthres[w];
            //                flag = false;
            //            }
            //            else
            //            {
            //                low = kprob;    //�ܷ�һֱ��ֵ������
            //                if ((countNum == coreD) || ((kprob - unvisitedNei[countNum - count].first) >= -EPSILON)) {      //��ζ������һ���� kprob��Ȼ����
            //                    flag = false;
            //                }
            //            }
            //        }
            //    }
            //}

        }
        else
        {
            //����ѡ��ʣ�µ��ھ� -- opt2
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
                kprob = kprob_comp_search_candi_increase(root, k, countNum, kprobKqualK, p);  //����ʽ����

                if ((kprob - kthres[w]) >= -EPSILON) {    //ȷ��lowȡֵ
                    low = kthres[w];
                    flag = false;
                }
                else
                {
                    low = 0;
                    if (countNum == coreD) {    //û���ھӿ��Բ���
                        flag = false;
                    }
                }
            }
        }
    }
    //cout << "root:" << root << " low:" << low << " upper:" << upper << endl;
    return low;
}


//�������º�ѡ��
vector<vector<double> > Uncertain_Core::delete_threshold_compute_batchUp(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "ɾ���� & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //��ѡ����������
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
    int countNum;   //��ѡ������
    //k = kk-1��ȷ��ͼcoreά���㷨 + Ѱ�Һ�ѡ�� + ����thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //ȷ��ͼcoreά�� -- kmax��Ȼ����
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����: core�仯 + v
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
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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
        if (count) {    //core[u]��С 
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]��core[v]����
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
                //core[u]��core[v]���� -- û��root
                count = delete_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited�ѿ���������ͨ����
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //ά��kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax��ѯ��������" << t << endl;
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
                //core[u]�ı䣬core[v]����
                count = delete_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]����
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]��core[v]������ -- ִ������ѭ��
                count = d_one_compute_batch(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]���䣬core[v]��С
                count = delete_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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


//���μ��� -- ��������
int Uncertain_Core::d_one_compute_batch(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //һ��ֻ��һ��root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //ȷ��thres�仯�ķ�Χ
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

    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        return count;
    }
    //cout << "k:" << k << " low:" << low << " up:" << up << endl;

    //Ѱ�Һ�ѡ��
    count = delete_search_candi(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);
    //std::cout << "count:" << count << endl;

    //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
    if (count != 0) {
        d_batch_update_candidate(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}


void Uncertain_Core::d_batch_update_candidate(vector<double>& kthres, int& k, int& count, double low, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    //std::cout << "---------k:" << k << " candi count:" << count << endl;
    std::queue<int> que;    //����que�������� p < curThres �ĵ�
    double curThres = low;
    int flag = count;
    int d, index, minV; //�ھ�������kprob��С���index��id��que�������Ĵ�С
    double p;
    //std::unordered_set<int> neiSet; //�����¶���
    vec_i record(count, 0); //��¼��ѡ���Ƿ���Ҫ����
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


        if (!que.empty()) { //thres��ΪcurThres
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

                //��Ҫ���µ��ھ�
                d = deg[node];
                for (int t = 0; t < d; t++) {
                    int w = adj[node][t].u;
                    index = candiRIndex[w];
                    //w��candiSet�в���Ҫ����
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
                        //que.push(w);    //����ʱά��que
                        neiNode->value = 0;
                    }
                    else
                    {
                        p = kprob_comp(w, visited, k + 1);
                        neiNode->value = p;                        

                        //if ((p - curThres) < EPSILON) {  //����ʱά��que
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
            
            //�������ھ�
            d = deg[minV];
            for (int t = 0; t < d; t++) {
                int w = adj[minV][t].u;
                index = candiRIndex[w];
                //w��candiSet�в���Ҫ����
                if (index != -1) {
                    Node* neiNode = hashTable[w];
                    it = mySet.find(neiNode);
                    mySet.erase(it);
                    neiNum[w]--;
                    if (neiNum[w] < k + 1) {
                        neiNode->value = 0;
                        //que.push(w);    //����ʱά��que
                    }
                    else
                    {
                        p = kprob_comp(w, visited, k + 1);
                        neiNode->value = p;

                        //if ((p - curThres) < EPSILON) {  //����ʱά��que
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



//�����Ż�����ϣ���̬����range + ����low + ��������
vector<vector<double> > Uncertain_Core::delete_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "ɾ���� & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    //std::cout << "delete-kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));

    int count;
    ct = 0;     //��ѡ����������
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
    int countNum;   //��ѡ������
    //k = kk-1��ȷ��ͼcoreά���㷨 + Ѱ�Һ�ѡ�� + ����thres
    if (core[u] < core[v]) {
        //std::cout << "core[u] < core[v]" << endl;
        //ȷ��ͼcoreά�� -- kmax��Ȼ����
        count = delete_find_core_subcore(u, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����: core�仯 + v
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
        count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);  //core��С && thres=0����candiNode�� ��node_visited=false
        if (count) {    //���ڵ�v core��С
            //Ѱ�Һ�ѡ�� -- ��������core��С�ĵ㣬up=thres[v]ʵʱ�仯
            count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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
        if (count) {    //core[u]��С 
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                int countV, countNumV;
                vec_b visited_2(n);
                vec_i candiNode_2(n);
                countV = delete_find_core_subcore(v, kk, visited_2, candiNode_2, candiReverseIndex, countNumV);
                if (countV) {
                    //core[u]��core[v]����
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
                //core[u]��core[v]���� -- û��root
                count = d_search_candi_core_change(thres[curK], -1, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //node_visited�ѿ���������ͨ����
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }

                //ά��kmax
                if (kk == kmax) {
                    bool eraseflag = true;
                    for (int t = 0; t < n; t++) {
                        if (core[t] == kmax) {
                            eraseflag = false;
                            std::cout << "kmax��ѯ��������" << t << endl;
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
                //core[u]�ı䣬core[v]����
                count = d_search_candi_core_change(thres[curK], v, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
                if (count) {
                    d_batch_update_candidate(thres[curK], curK, count, 0, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
                }
            }
        }
        else
        {
            //core[u]����
            if (!node_visited[v] && (core[v] == kk)) {  //u��vΪ������ͬ����ͨ����
                count = delete_find_core_subcore(v, kk, node_visited, candiNode, candiReverseIndex, countNum);
            }

            if (core[u] == core[v]) {
                //core[u]��core[v]������ -- ִ������ѭ��
                count = d_one_compute_opt(range, thres[curK], u, v, curK, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
            }
            else
            {
                //core[u]���䣬core[v]��С
                count = d_search_candi_core_change(thres[curK], u, curK, node_visited, candiNode, candiReverseIndex, node_visited_num, countNum, mySet, hashTable);   //ͬʱ����������ͨ����
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



//�Ż����
int Uncertain_Core::d_one_compute_opt(vector<vector<double> >& range, vector<double>& kthres, int u, int v, int k, vec_b& visited, vec_i& candiRIndex, vec_i& neiNum, std::set<Node*, NodePtrComparator>& mySet, std::map<int, Node*>& hashTable) {
    double low = 2;
    double up = 0;
    vector<int> root(2);
    root[1] = -1;   //һ��ֻ��һ��root

    int count = 0;
    double thres1 = kthres[u];
    double thres2 = kthres[v];

    //ȷ��thres�仯�ķ�Χ
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

    if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
        return count;
    }
    //cout << "k: " << k << " low:" << low << " up:" << up << endl;

    //Ѱ�Һ�ѡ��
    count = d_candi_dynamic_update_range(kthres, low, up, root, k, visited, candiRIndex, neiNum, mySet, hashTable);

    //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
    if (count != 0) {
        d_batch_update_candidate(kthres, k, count, low, visited, candiRIndex, neiNum, mySet, hashTable);
    }
    return count;
}
