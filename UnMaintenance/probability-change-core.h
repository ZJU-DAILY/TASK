#pragma once
#include "uncertain-core.h"

class Probability_Core : public Uncertain_Core
{
public:
	//ɾ��+���ӱ� -- ���ֳ�ʼͼ����
	//void in_de_compare(vector<vector<double> >& thres, double scale, int count);
	//void in_de_prob_compare(vector<vector<double> >& thres, double scale, int count, double diff);


	//��Ե��������
	void increase_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	void increase_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//Ѱ��candi��̬����range
	void increase_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//���Ƶ�
	void increase_threshold_compute_restriction_point(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//��������
	void increase_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//�����Ż��Ľ��
	void increase_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);


	//��Ե���ʼ���
	void decrease_threshold_compare(vector<vector<double> >& thres, double scale, int count);
	void decrease_threshold_compute_candidate(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//Ѱ��candi ��̬����range
	void decrease_threshold_compute_update_range(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//���ݷ������½�
	void decrease_threshold_compute_shrink_low(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//��������
	void decrease_threshold_compute_batchup(vector<vector<double> >& thres, int u, int v, double& time, int& ct);
	//�����Ż��Ľ��
	void decrease_threshold_compute_opt(vector<vector<double> >& thres, int u, int v, double& time, int& ct);



private:

};
