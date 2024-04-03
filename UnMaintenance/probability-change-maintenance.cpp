#include "probability-change-core.h"


//��Ե��������
void Probability_Core::increase_threshold_compare(vector<vector<double> >& thres, double scale, int count) {
    int ct = 0;
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);

    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //ÿ������Ѱ�Һ�ѡ����������ۼƺ�

    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //ȷ��ÿ����ֻ����һ��
                if (rand() > scale_size) {  //������
                    std::cout << "ct:" << ct << " selected-edge:" << i << " nei:" << u << endl;
                    double p = adj[i][j].p;
                    if (p > 0.9) {
                        continue;
                    }
                    p += 0.05;
                    std::cout << "����p��" << p << endl;
                    //ά��ͼ�ṹ -- ��Ե��������
                    adj[i][j].p = p;

                    int d = deg[u];
                    for (int t = 0; t < d; t++) {
                        int w = adj[u][t].u;
                        if (w == i) {
                            adj[u][t].p = p;
                            break;
                        }
                    }
                    ct++;

                    int mincore = min(core[u], core[i]);
                    //increase_threshold_compute_candidate(thres, u, i, tm_1, candi_count);  //����Ϊu��i
                    //increase_threshold_compute_update_range(thres, u, i, tm_1, candi_count);
                    increase_threshold_compute_restriction_point(thres, u, i, tm_1, candi_count);
                    //increase_threshold_compute_update_range(thres, u, i, tm_1, candi_count);
                    candidate_tm += tm_1;
                    candi_sum += candi_count;

                    //delete[] core;
                    vector<vector<double> > thres_2;
                    insert_threshold_recompute(thres_2, tm_2);
                    recompute_tm += tm_2;

                    //�Ƚϲ���ߺ�ı仯��ֵ
                    bool compare = compareArraysTrueOrFalse(thres, thres_2);
                    if (!compare) {
                        std::cout << "�����Ƿ���ͬ----------------------------------------------��" << compare << endl;
                        compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
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
        }
        if (ct >= count) {  //countΪ��Ե���������ı���
            break;
        }
    }
    std::cout << "insert_compute_candidate_time:" << candidate_tm << endl;
    std::cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
    std::cout << "candi sum:" << candi_sum << endl;
}

//��Ե�������� -- Ѱ�Һ�ѡ�� ������
void Probability_Core::increase_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "��Ե�������� -- ����ԭʼͼ�������½�--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    //return range;
}


//update
void Probability_Core::increase_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "��Ե�������� -- ����ԭʼͼ�������½�--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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

    time = omp_get_wtime() - tm;
}


void Probability_Core::increase_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "��Ե�������� -- ����ԭʼͼ�������½�--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
    vec_i candiReverseIndex(n); //Ϊcandi index��node index��������
    vec_i node_visited_num(n);  //ÿ����ѡ�ڵ�Ķ��� -- ��СΪn
    vec_b node_visited(n);      //��nӳ��ýڵ��Ƿ�Ҫ����
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //����Ӱ�����½磬�����root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;

        //Ѱ�Һ�ѡ�� -- node_visited��ʾ������ѡ��ĳ�ʼ��ͼ������kk��Ѱ��ʱ�ж�deg�Ƿ�����Ҫ��
        count = insert_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;
        //cout << "count:" << count << endl;

        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}


void Probability_Core::increase_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "��Ե�������� -- ����ԭʼͼ�������½�--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
        }
    }

    time = omp_get_wtime() - tm;
}



void Probability_Core::increase_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "��Ե�������� -- ����ԭʼͼ�������½�--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
        }
    }

    time = omp_get_wtime() - tm;
}




//��Ե���ʼ���
void Probability_Core::decrease_threshold_compare(vector<vector<double> >& thres, double scale, int count) {
    int ct = 0;
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);

    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //ÿ������Ѱ�Һ�ѡ����������ۼƺ�

    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //ȷ��ÿ����ֻ����һ��
                if (rand() > scale_size) {  //������
                    std::cout << "selected-edge:" << i << " nei:" << u << endl;
                    double p = adj[i][j].p;
                    if (p < 0.1) {
                        continue;
                    }
                    p -= 0.05;

                    //ά��ͼ�ṹ -- ��Ե��������
                    adj[i][j].p = p;

                    int d = deg[u];
                    for (int t = 0; t < d; t++) {
                        int w = adj[u][t].u;
                        if (w == i) {
                            adj[u][t].p = p;
                            break;
                        }
                    }
                    ct++;

                    vector<vector<double> > copy_thres = thres;
                    int mincore = min(core[u], core[i]);
                    //decrease_threshold_compute_candidate(thres, u, i, tm_1, candi_count);  //����Ϊu��i
                    //decrease_threshold_compute_update_range(thres, u, i, tm_1, candi_count);
                    //decrease_threshold_compute_shrink_low(thres, u, i, tm_1, candi_count);
                    decrease_threshold_compute_batchup(thres, u, i, tm_1, candi_count);
                    //decrease_threshold_compute_opt(thres, u, i, tm_1, candi_count);
                    candidate_tm += tm_1;
                    candi_sum += candi_count;

                    //delete[] core;
                    vector<vector<double> > thres_2;
                    insert_threshold_recompute(thres_2, tm_2);
                    recompute_tm += tm_2;

                    //�Ƚϲ���ߺ�ı仯��ֵ
                    //compareArrays(copy_thres, thres_2, range[0], range[1], mincore);    //���ܱȽ�core��С�ĵ�
                    bool compare = compareArraysTrueOrFalse(thres, thres_2);
                    if (!compare) {
                        cout << "�����Ƿ���ͬ----------------------------------------------��" << compare << endl;
                        compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
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
        }
        if (ct >= count) {  //countΪ��Ե���������ı���
            break;
        }
    }
    std::cout << "insert_compute_candidate_time:" << candidate_tm << endl;
    std::cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
    std::cout << "candi sum:" << candi_sum << endl;
}


//��Ե���ʼ�С -- Ѱ�Һ�ѡ�� ������
void Probability_Core::decrease_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "��Ե���ʼ�С & ����ԭʼͼ�������½硪����������" << endl;

    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    //return range;
}


void Probability_Core::decrease_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}


void Probability_Core::decrease_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = d_shrink_lower_bound(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = d_shrink_lower_bound(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = d_shrink_lower_bound(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = d_shrink_lower_bound(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            delete_update_candidate_thres(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}



void Probability_Core::decrease_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = delete_compute_low(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = delete_compute_low(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = delete_compute_low(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = delete_compute_low(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}


void Probability_Core::decrease_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();
    int kk = min(core[u], core[v]);
    std::cout << "kk:" << kk << endl;
    vector<vector<double> > range(2, vector<double>(kk));
    double low = 2;
    double up = 0;
    vector<int> root(2);

    int count;
    ct = 0;
    vec_i candiReverseIndex(n);
    vec_i node_visited_num(n);
    vec_b node_visited(n);
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;
        double thres1 = thres[k][u];
        double thres2 = thres[k][v];
        root[1] = -1;   //һ��ֻ��һ��root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //ȷ��thres�仯�ķ�Χ
        if (std::fabs(thres1 - thres2) < EPSILON) {
            root[0] = u;
            up = thres1;
            double kprob1 = d_shrink_lower_bound(thres[k], k + 1, u);

            root[1] = v;
            double kprob2 = d_shrink_lower_bound(thres[k], k + 1, v);

            low = min(kprob1, kprob2);
        }
        else if (thres1 > thres2) {
            root[0] = v;
            up = thres2;
            low = d_shrink_lower_bound(thres[k], k + 1, v);
        }
        else
        {
            root[0] = u;
            up = thres1;
            low = d_shrink_lower_bound(thres[k], k + 1, u);
        }
        range[0][k] = low;
        range[1][k] = up;

        if (low >= up || fabs(low - up) < EPSILON) {    //����������ѯrange����û�е���Ҫ�ı�
            continue;
        }

        //Ѱ�Һ�ѡ��
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //���º�ѡ����thres����k = kk-1ʱ��ά��coreness + kmax
        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}