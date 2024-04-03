#include "probability-change-core.h"


//边缘概率增加
void Probability_Core::increase_threshold_compare(vector<vector<double> >& thres, double scale, int count) {
    int ct = 0;
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);

    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //每个方法寻找候选对象个数的累计和

    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //确保每条边只处理一次
                if (rand() > scale_size) {  //待增加
                    std::cout << "ct:" << ct << " selected-edge:" << i << " nei:" << u << endl;
                    double p = adj[i][j].p;
                    if (p > 0.9) {
                        continue;
                    }
                    p += 0.05;
                    std::cout << "新增p：" << p << endl;
                    //维护图结构 -- 边缘概率增加
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
                    //increase_threshold_compute_candidate(thres, u, i, tm_1, candi_count);  //顶点为u、i
                    //increase_threshold_compute_update_range(thres, u, i, tm_1, candi_count);
                    increase_threshold_compute_restriction_point(thres, u, i, tm_1, candi_count);
                    //increase_threshold_compute_update_range(thres, u, i, tm_1, candi_count);
                    candidate_tm += tm_1;
                    candi_sum += candi_count;

                    //delete[] core;
                    vector<vector<double> > thres_2;
                    insert_threshold_recompute(thres_2, tm_2);
                    recompute_tm += tm_2;

                    //比较插入边后的变化阈值
                    bool compare = compareArraysTrueOrFalse(thres, thres_2);
                    if (!compare) {
                        std::cout << "数组是否相同----------------------------------------------：" << compare << endl;
                        compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
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
        }
        if (ct >= count) {  //count为边缘概率新增的边数
            break;
        }
    }
    std::cout << "insert_compute_candidate_time:" << candidate_tm << endl;
    std::cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
    std::cout << "candi sum:" << candi_sum << endl;
}

//边缘概率增大 -- 寻找候选集 并更新
void Probability_Core::increase_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "边缘概率增大 -- 根据原始图计算上下界--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
            insert_update_candidate_thres(thres[k], k, count, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
    //return range;
}


//update
void Probability_Core::increase_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "边缘概率增大 -- 根据原始图计算上下界--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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

    time = omp_get_wtime() - tm;
}


void Probability_Core::increase_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "边缘概率增大 -- 根据原始图计算上下界--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
    vec_i candiReverseIndex(n); //为candi index与node index建立关联
    vec_i node_visited_num(n);  //每个候选节点的度数 -- 大小为n
    vec_b node_visited(n);      //用n映射该节点是否要访问
    std::set<Node*, NodePtrComparator> mySet;
    std::map<int, Node*> hashTable;
    vec_i candiNode(n);
    vec_b rest_node(n);

    for (int k = 0; k < kk; k++) {
        //std::cout << "k:" << k << endl;  
        //计算影响上下界，并求得root
        i_singlepoint_compute_restriction_point(thres[k], low, up, root, u, v, k, rest_node);
        range[0][k] = low;
        range[1][k] = up;

        //寻找候选集 -- node_visited表示包含候选点的初始子图，对于kk，寻找时判断deg是否满足要求
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

    //std::cout << "边缘概率增大 -- 根据原始图计算上下界--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
        }
    }

    time = omp_get_wtime() - tm;
}



void Probability_Core::increase_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    double tm = omp_get_wtime();

    //std::cout << "边缘概率增大 -- 根据原始图计算上下界--------" << endl;
    int kk = min(core[u], core[v]);
    vector<vector<double> > range(2, vector<double>(kk));
    //std::cout << "kk:" << kk << endl;
    double low, up;
    int root;

    int count;
    ct = 0;
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
        }
    }

    time = omp_get_wtime() - tm;
}




//边缘概率减少
void Probability_Core::decrease_threshold_compare(vector<vector<double> >& thres, double scale, int count) {
    int ct = 0;
    size_t scale_size = RAND_MAX * scale;
    std::srand(0);

    double candidate_tm = 0.0;
    double recompute_tm = 0.0;
    double tm_1;
    double tm_2;
    int candi_count;
    long long int candi_sum = 0;    //每个方法寻找候选对象个数的累计和

    for (int i = 0; i < n; i++) {
        int d = deg[i];
        for (int j = 0; j < d; j++) {
            int u = adj[i][j].u;
            if (i < u) {    //确保每条边只处理一次
                if (rand() > scale_size) {  //待增加
                    std::cout << "selected-edge:" << i << " nei:" << u << endl;
                    double p = adj[i][j].p;
                    if (p < 0.1) {
                        continue;
                    }
                    p -= 0.05;

                    //维护图结构 -- 边缘概率增加
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
                    //decrease_threshold_compute_candidate(thres, u, i, tm_1, candi_count);  //顶点为u、i
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

                    //比较插入边后的变化阈值
                    //compareArrays(copy_thres, thres_2, range[0], range[1], mincore);    //不能比较core减小的点
                    bool compare = compareArraysTrueOrFalse(thres, thres_2);
                    if (!compare) {
                        cout << "数组是否相同----------------------------------------------：" << compare << endl;
                        compareKsizeArraysTrueOrFalse(thres, thres_2, mincore);
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
        }
        if (ct >= count) {  //count为边缘概率新增的边数
            break;
        }
    }
    std::cout << "insert_compute_candidate_time:" << candidate_tm << endl;
    std::cout << "insert_threshold_recompute_time:" << recompute_tm << endl;
    std::cout << "candi sum:" << candi_sum << endl;
}


//边缘概率减小 -- 寻找候选集 并更新
void Probability_Core::decrease_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct) {
    //std::cout << "边缘概率减小 & 根据原始图计算上下界——————" << endl;

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
        root[1] = -1;   //一般只有一个root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //确定thres变化的范围
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

        if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
            continue;
        }

        //寻找候选集
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //更新候选集的thres：在k = kk-1时，维护coreness + kmax
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
        root[1] = -1;   //一般只有一个root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //确定thres变化的范围
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

        if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
            continue;
        }

        //寻找候选集
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //更新候选集的thres：在k = kk-1时，维护coreness + kmax
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
        root[1] = -1;   //一般只有一个root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //确定thres变化的范围
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

        if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
            continue;
        }

        //寻找候选集
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //更新候选集的thres：在k = kk-1时，维护coreness + kmax
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
        root[1] = -1;   //一般只有一个root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //确定thres变化的范围
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

        if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
            continue;
        }

        //寻找候选集
        count = delete_search_candi(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //更新候选集的thres：在k = kk-1时，维护coreness + kmax
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
        root[1] = -1;   //一般只有一个root
        //cout << "thres1:" << thres1 << " thres2:" << thres2 << endl;

        //确定thres变化的范围
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

        if (low >= up || fabs(low - up) < EPSILON) {    //如果不满足查询range，则没有点需要改变
            continue;
        }

        //寻找候选集
        count = d_candi_dynamic_update_range(thres[k], low, up, root, k, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        ct += count;

        //更新候选集的thres：在k = kk-1时，维护coreness + kmax
        if (count != 0) {
            d_batch_update_candidate(thres[k], k, count, low, node_visited, candiReverseIndex, node_visited_num, mySet, hashTable);
        }
    }

    time = omp_get_wtime() - tm;
}