#pragma once
#include "uncertain-core.h"

class Probability_Core : public Uncertain_Core
{
public:
	//删除+增加边 -- 保持初始图不变
	//void in_de_compare(vector<vector<double> >& thres, double scale, int count);
	//void in_de_prob_compare(vector<vector<double> >& thres, double scale, int count, double diff);


	//边缘概率增加
	void increase_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	void increase_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//寻找candi动态更新range
	void increase_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//限制点
	void increase_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//批量更新
	void increase_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//三种优化的结合
	void increase_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);


	//边缘概率减少
	void decrease_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	void decrease_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//寻找candi 动态更新range
	void decrease_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//阶梯法收缩下界
	void decrease_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//批量更新
	void decrease_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//三种优化的结合
	void decrease_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);



private:

};
