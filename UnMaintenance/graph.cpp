#include "graph.h"



Graph::Graph()
{
	n = 0;  m = 0; max_degree = 0;
	offs = NULL; deg = NULL;
	data = NULL; adj = NULL;
}

Graph::~Graph()
{
	if (offs != NULL) delete[] offs;
	if (data != NULL) delete[] data;
	if (adj != NULL) delete[] adj;
	if (deg != NULL) delete[] deg;
}

void Graph::read_bin(const string& infile)
{
	FILE* in = fopen(infile.c_str(), "rb");
	if (in == NULL) {
		printf("No such the input file:%s\n", infile.c_str());
		abort();
	}
	printf("Reading graph ...\nGraph name: %s\n", infile.c_str());
	double tm = omp_get_wtime();
	size_t FRead = 0;
	FRead = fread(&n, sizeof(int), 1, in);
	FRead = fread(&m, sizeof(long long int), 1, in);

	deg = new Int[n]; memset(deg, 0, sizeof(int) * n);
	offs = new Long[n]; memset(offs, 0, sizeof(long long int) * n);
	//Int *edges = new Int[m]; memset(edges, 0, sizeof(Int) * m);
	//double *probability = new double[m]; memset(probability, 0, sizeof(double) * m);

	data = new Pairs[m]; memset(data, 0, sizeof(Pairs) * m);
	adj = new Pairs * [n];

	FRead = fread(deg, sizeof(Int), n, in);
	//fread(edges, sizeof(Int), m, in);
	//fread(probability, sizeof(double), m, in);
	for (Long i = 0; i < m; ++i)
		FRead = fread(&data[i].u, sizeof(Int), 1, in);
	for (Long i = 0; i < m; ++i) {
		FRead = fread(&data[i].p, sizeof(double), 1, in);
		// if (data[i].p + 1e-5 > 1.0)
		// 	data[i].p = 1.0;
	}
	fclose(in);
	for (Int i = 0, s = 0; i < n; ++i) { // Construct offs of all vertices
		offs[i] = s;
		adj[i] = data + s;		//指向了 data 数组中的第 s 个元素
		s += deg[i];
		max_degree = deg[i] > max_degree ? deg[i] : max_degree;
	}
	printf("n = %d, m = %lld, max_deg = %d\n", n, m / 2, max_degree);
	printf("Reading time: %lf s\n", omp_get_wtime() - tm);

	//unselected
	//unselected = { {{916,2634},0.204} };
	//unselected = { {{916,2634},0.204},{{917,1273},0.227}, {{918,1803},0.221},{{918,3531},0.202}, {{920,3535},0.461}, {{922,3026},0.172} };
	return;

	Int* deg_reverse = new Int[n]();
	for (Int i = 0; i < n; ++i) { // The index of the reverse edge
		Int d = deg[i];
		for (Int j = 0; j < d; ++j) {
			Int u = adj[i][j].u;
			assert(u < n&& u >= 0);
			adj[i][j].re = deg_reverse[u]++;
		}
	}
	delete[] deg_reverse;
	//test
	for (Int i = 0; i < n; ++i) {
		int d = deg[i];
		assert(d <= max_degree && d >= 0);
		for (Int j = 0; j < d; ++j) {
			Int u = adj[i][j].u;
			Int re = adj[u][adj[i][j].re].u;
			Long re_pos = adj[i][j].re;
			Double p = adj[i][j].p;
			assert(p - 1e-16 <= 1.0 && p + 1e-16 >= 0.0);
			if (re != i) {
				printf("v=%d, u=%d, re=%d, re_pos=%lld, p=%lf\n", i, u, re, re_pos, adj[i][j].p);
				exit(1);
			}
			if (deg[i] > n) {
				printf("v=%d, u=%d, deg[i]=%d\n", i, u, deg[i]);
				exit(1);
			}
			if (re_pos > deg[u]) {
				printf("v=%d, u=%d, re=%d, re_pos=%lld, deg=%d\n", i, u, re, re_pos, deg[u]);
				exit(1);
			}
			//printf("v = %d, u = %d, p = %f\n", i, u, p);
			//printf("u = %d, v = %d\n", u, re);
		}
		/*if (i > 1)
		break;*/
	}

}

void Graph::read_bin(const string& infile, int voe, float scale)	//scale表示缩放因子
{
	if (voe != 0 && voe != 1) {
		printf("error: voe = 0 or 1 ?\n");
		exit(1);
	}
	if (scale > 1.0) {
		printf("error: scale <= 1 ?\n");
		exit(1);
	}
	FILE* in = fopen(infile.c_str(), "rb");
	if (in == NULL) {
		printf("No such the input file:%s\n", infile.c_str());
		abort();
	}
	printf("Reading graph ...\nGraph name: %s\n", infile.c_str());
	double tm = omp_get_wtime();
	size_t FRead = 0;		//记录成功读取的字节数
	FRead = fread(&n, sizeof(int), 1, in);
	FRead = fread(&m, sizeof(long long int), 1, in);

	deg = new Int[n]; memset(deg, 0, sizeof(int) * n);
	offs = new Long[n]; memset(offs, 0, sizeof(long long int) * n);
	//Int *edges = new Int[m]; memset(edges, 0, sizeof(Int) * m);
	//double *probability = new double[m]; memset(probability, 0, sizeof(double) * m);

	data = new Pairs[m]; memset(data, 0, sizeof(Pairs) * m);
	adj = new Pairs * [n];

	FRead = fread(deg, sizeof(Int), n, in);
	//fread(edges, sizeof(Int), m, in);
	//fread(probability, sizeof(double), m, in);
	for (Long i = 0; i < m; ++i)
		FRead = fread(&data[i].u, sizeof(Int), 1, in);
	for (Long i = 0; i < m; ++i) {
		FRead = fread(&data[i].p, sizeof(double), 1, in);
		// if (data[i].p + 1e-5 > 1.0)
		// 	data[i].p = 1.0;
	}
	fclose(in);
	for (Int i = 0, s = 0; i < n; ++i) { // Construct offs of all vertices
		offs[i] = s;		//顶点i在邻接数组中的起始位置
		adj[i] = data + s;	//顶点的邻接列表
		s += deg[i];
		max_degree = max(deg[i], max_degree);
	}
	printf("n = %d, m = %lld, max_deg = %d\n", n, m / 2, max_degree);
	if (voe == 0)	//按照一定比例随机选择一部分顶点，更新图结构，并记录选中的顶点数和更新后的最大度数
	{
		int* selected = new int[n];
		memset(selected, -1, sizeof(int) * n);
		int cnts = 0;
		size_t scale_size = RAND_MAX * scale;
		//srand(time(NULL));
		srand(0);	//保证每次程序运行时都产生相同的随机数序列
		for (int i = 0; i < n; ++i) {
			if (rand() <= scale_size)		//选中该顶点
				selected[i] = cnts++;
		}
		printf("random nodes, cnts=%d\n", cnts);		//选中的顶点数
		max_degree = 0;
		for (int i = 0; i < n; ++i) {
			int d = deg[i];
			int cnt = 0;
			int vn = selected[i];
			if (vn == -1) {
				deg[i] = 0;
				continue;
			}
			for (int j = 0; j < d; ++j) {
				int u = adj[i][j].u;
				int un = selected[u];	//邻居的选中状态
				if (un >= 0) {
					adj[i][cnt].u = u;
					adj[i][cnt].p = adj[i][j].p;
					cnt++;
				}
			}
			deg[i] = cnt;
			max_degree = max(max_degree, cnt);
		}
		//n = cnts;
		delete[] selected;
	}
	else if (voe == 1) {		//更新图结构，包括顶点度数、邻接数组、选中边数
		size_t scale_size = RAND_MAX * scale;
		int* degn = new int[n]();	//存储每个顶点的度数
		Long cnts = 0;
		//srand(time(NULL));
		srand(0);
		for (int i = 0; i < n; ++i) {
			int d = deg[i];
			int cnt = 0;
			for (int j = 0; j < d; ++j) {
				int u = adj[i][j].u;
				if (i < u) {		//确保每条边只被处理一次
					if (rand() <= scale_size) {
						cnt = degn[i]++;
						adj[i][cnt].u = u;
						adj[i][cnt].p = adj[i][j].p;
						cnt = degn[u]++;
						adj[u][cnt].u = i;
						adj[u][cnt].p = adj[i][j].p;
						cnts += 2;
					}
				}
			}
		}
		for (int i = 0; i < n; ++i)
			deg[i] = degn[i];

		printf("random edges, cnts=%Ld\n", cnts / 2);
		delete[] degn;
	}
	//printf("n = %d, m = %lld, max_deg = %d\n", n, m / 2, max_degree);
	printf("Reading time: %lf s\n", omp_get_wtime() - tm);
	return;

	Int* deg_reverse = new Int[n]();
	for (Int i = 0; i < n; ++i) { // The index of the reverse edge
		Int d = deg[i];
		for (Int j = 0; j < d; ++j) {
			Int u = adj[i][j].u;
			adj[i][j].re = deg_reverse[u]++;
		}
	}
	delete[] deg_reverse;
	//test
	/* for (Int i = 0; i < n; ++i) {
		for (Int j = 0; j < deg[i]; ++j) {
			Int u = adj[i][j].u;
			Int re = adj[u][adj[i][j].re].u;
			Long re_pos = adj[i][j].re;
			Double p = adj[i][j].p;
			if (re != i) {
				printf("v=%d, u=%d, re=%d, re_pos=%lld, p=%lf\n", i, u, re, re_pos, adj[i][j].p);
				exit(1);
			}
			if (deg[i] > n) {
				printf("v=%d, u=%d, deg[i]=%d\n", i, u, deg[i]);
				exit(1);
			}
			if (re_pos > deg[u]) {
				printf("v=%d, u=%d, re=%d, re_pos=%lld, deg=%d\n", i, u, re, re_pos, deg[u]);
				exit(1);
			}
			//printf("v = %d, u = %d, p = %f\n", i, u, p);
			//printf("u = %d, v = %d\n", u, re);
		}
	} */
}


vector<int> Graph::edge_selected_bin(const string& infile, float scale) {
	if (scale > 1.0) {
		printf("error: scale <= 1 ?\n");
		exit(1);
	}
	FILE* in = fopen(infile.c_str(), "rb");
	if (in == NULL) {
		printf("No such the input file:%s\n", infile.c_str());
		abort();
	}
	printf("Reading and scaling graph ...\nGraph name: %s\n", infile.c_str());
	double tm = omp_get_wtime();
	size_t FRead = 0;		//记录成功读取的字节数
	FRead = fread(&n, sizeof(int), 1, in);
	FRead = fread(&m, sizeof(long long int), 1, in);

	deg = new Int[n]; memset(deg, 0, sizeof(int) * n);
	offs = new Long[n]; memset(offs, 0, sizeof(long long int) * n);
	data = new Pairs[m]; memset(data, 0, sizeof(Pairs) * m);
	adj = new Pairs * [n];

	FRead = fread(deg, sizeof(Int), n, in);
	vector<int> deg_record(n);
	for (int i = 0; i < n; i++) {
		deg_record[i] = deg[i];
	}

	//fread(edges, sizeof(Int), m, in);
	//fread(probability, sizeof(double), m, in);
	for (Long i = 0; i < m; ++i)
		FRead = fread(&data[i].u, sizeof(Int), 1, in);
	for (Long i = 0; i < m; ++i) {
		FRead = fread(&data[i].p, sizeof(double), 1, in);
		// if (data[i].p + 1e-5 > 1.0)
		// 	data[i].p = 1.0;
	}
	fclose(in);
	for (Int i = 0, s = 0; i < n; ++i) { // Construct offs of all vertices
		offs[i] = s;		//顶点i在邻接数组中的起始位置
		adj[i] = data + s;	//顶点的邻接列表
		s += deg[i];
		max_degree = max(deg[i], max_degree);
	}
	printf("n = %d, m = %lld, max_deg = %d\n", n, m / 2, max_degree);

	//更新图结构，包括顶点度数、邻接数组、未选中边数
	size_t scale_size = RAND_MAX * scale;
	int* degn = new int[n]();	//存储每个顶点的度数
	Long cnts = 0;
	//srand(time(NULL));
	srand(0);
	for (int i = 0; i < n; ++i) {
		int d = deg[i];
		int cnt = 0;
		for (int j = 0; j < d; ++j) {
			int u = adj[i][j].u;
			if (i < u) {		//确保每条边只被处理一次
				if (rand() <= scale_size) {
					cnt = degn[i]++;
					adj[i][cnt].u = u;
					adj[i][cnt].p = adj[i][j].p;
					//cout << "selected-vertex:" << i << " nei:" << u << endl;

					cnt = degn[u]++;
					adj[u][cnt].u = i;
					adj[u][cnt].p = adj[i][j].p;
					cnts += 2;
				}
				else
				{
					//cout << "unselected-vertex:" << i << " nei:" << u << endl;
					std::pair<int, int> p1 = std::make_pair(i, u);
					double p = adj[i][j].p;
					unselected.push_back(make_pair(p1, p));
				}
			}
		}
	}
	for (int i = 0; i < n; ++i)
		deg[i] = degn[i];
	delete[] degn;

	cout << "unselected edges:" << unselected.size() << endl;
	cout << "random edges count:" << cnts / 2 << endl;
	printf("Reading time: %lf s\n", omp_get_wtime() - tm);
	//return;

	Int* deg_reverse = new Int[n]();
	for (Int i = 0; i < n; ++i) { // The index of the reverse edge
		Int d = deg[i];
		for (Int j = 0; j < d; ++j) {
			Int u = adj[i][j].u;
			adj[i][j].re = deg_reverse[u]++;
		}
	}
	delete[] deg_reverse;
	return deg_record;
}




void Graph::creat_bin(const string& infile)
{
	FILE* in = fopen(infile.c_str(), "r");
	if (in == NULL) {
		printf("No such the input file:%s\n", infile.c_str());
		abort();
	}

	vector< pair< pair<Int, Int>, double> > datas;

	Char buff[1024];
	Int x, y, max_n = 0;
	double p; Long e_len = 1;
	printf("Reading txt...\n");
	//fgets(buff, 1024, in);
	while (fgets(buff, 1024, in) != NULL)
	{
		sscanf(buff, "%d, %d ,%lf", &x, &y, &p);
		//printf("%d %d, %f\n", x, y, p);
		datas.emplace_back(make_pair(x, y), p);
		datas.emplace_back(make_pair(y, x), p);

		max_n = x > max_n ? x : max_n;		//确保为x、y的最大值
		max_n = y > max_n ? y : max_n;
	}
	fclose(in);
	max_n += 1;
	// re-number ids of vertices
	/*if (true) {
	Int vertices = 0;
	vector<Int> mapping(max_n, -1);
	for (size_t i = 0; i < datas.size(); ++i) {
	x = datas[i].first.first;
	y = datas[i].first.second;
	if (mapping[x] < 0)
	mapping[x] = vertices++;
	if (mapping[y] < 0)
	mapping[y] = vertices++;
	datas[i].first.first = mapping[x];
	datas[i].first.second = mapping[y];
	}
	max_n = vertices;
	}*/
	sort(datas.begin(), datas.end());
	for (size_t i = 1; i < datas.size(); ++i) {
		x = datas[i].first.first;
		y = datas[i].first.second;
		if (x == datas[i - 1].first.first && y == datas[i - 1].first.second)
			continue;
		datas[e_len].first.first = x;
		datas[e_len].first.second = y;
		datas[e_len++].second = datas[i].second;
	}
	datas.resize(e_len);
	n = max_n; m = e_len;
	printf("n=%d, m=%lld\n", n, m);
	string outfile;
	size_t pos = infile.find(".");		//找到输入文件名中的第一个(.)，将之间的部分作为文件名
	if (pos > 0)
		outfile = infile.substr(0, pos);
	outfile += ".bin";
	printf("outfile = %s\nSaving to bin...\n", outfile.c_str());
	deg = new Int[n]; memset(deg, 0, sizeof(Int) * n);
	Int* edges = new Int[m]; memset(edges, 0, sizeof(Int) * m);
	//Int *degn = new Int[n]; memset(degn, 0, sizeof(Int) * n);
	double* probability = new double[m]; memset(probability, 0, sizeof(double) * m);
	for (size_t i = 0; i < datas.size(); ++i) {
		x = datas[i].first.first;
		y = datas[i].first.second;
		p = datas[i].second;

		deg[x]++;
		edges[i] = y;
		probability[i] = p;
		//degn[y]++;
	}
	/*for (Int i = 0; i < n; ++i) {
	if (deg[i] != degn[i]) {
	printf("error\n");
	break;
	}
	}*/
	//adj.first = edges; adj.second = probability;

	FILE* out = fopen(outfile.c_str(), "wb");
	fwrite(&n, sizeof(Int), 1, out);
	fwrite(&m, sizeof(Long), 1, out);
	fwrite(deg, sizeof(Int), n, out);
	fwrite(edges, sizeof(Int), m, out);
	fwrite(probability, sizeof(double), m, out);
	fclose(out);
	printf("Saving done\n");

	/*for (size_t i = 0; i < m; ++i) {
	x = datas[i].first.first;
	y = datas[i].first.second;
	p = datas[i].second;
	printf("%d %d %lf\n",x, y, p);
	}*/
	exit(1);
}

