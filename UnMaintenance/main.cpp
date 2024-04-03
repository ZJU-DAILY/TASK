#include "uncertain-core.h"
#include "probability-change-core.h"

std::string getParentDirectory(const std::string& filePath) {
	size_t found = filePath.find_last_of("/\\");
	return filePath.substr(0, found);
}

void load_cores(const char* str, vec_i& cores)
{
	const char* txt = strstr(str, "txt");
	const char* bin = strstr(str, "bin");
	if (txt == NULL && bin == NULL) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}

	int n = 0;
	FILE* in = NULL;
	if (txt != NULL) in = fopen(str, "r");	//以文本读取模式打开文件
	else in = fopen(str, "rb");		//以二进制模式打开文件

	if (in == NULL) {
		printf("No such file: %s\n", str);
		exit(1);
	}
	if (txt != NULL) {
		int x = fscanf(in, "%d", &n);	//x表示成功读取的项数
		printf("file=%s, n=%d\n", str, n);
		cores.resize(n);
		for (int i = 0; i < n; ++i)
			x = fscanf(in, "%d", &cores[i]);
	}
	else {
		size_t x = fread(&n, sizeof(int), 1, in);
		printf("file=%s, n=%d\n", str, n);
		cores.resize(n);
		for (int i = 0; i < n; ++i)
			x = fread(&cores[i], sizeof(int), 1, in);
	}
	fclose(in);
}

void prinf_core(const char* str, vec_i& cores)	//首先存储顶点个数n，然后依次存储c(i)
{
	const char* txt = strstr(str, "txt");
	const char* bin = strstr(str, "bin");
	if (txt == NULL && bin == NULL) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}
	int n = cores.size();
	FILE* in = NULL;
	if (bin != NULL) in = fopen(str, "wb");	//以二进制写入模式打开文件
	else in = fopen(str, "w");		//以文本写入模式打开文件

	if (in == NULL) {
		printf("No such file: %s\n", str);
		exit(1);
	}
	if (bin != NULL) {
		fwrite(&n, sizeof(int), 1, in);
		for (int i = 0; i < n; ++i)
			fwrite(&cores[i], sizeof(int), 1, in);
	}
	else {
		//printf("n=%d\n", n);
		fprintf(in, "%d\n", n);		//以文本形式写入文件
		for (int i = 0; i < n; ++i)
			fprintf(in, "%d\n", cores[i]);
	}
	fclose(in);
}

//大图存储在bin文件中--读取速度更快
void printf_thres(const std::string& filePath, const std::vector<std::vector<double> >& thres) {
	size_t txtPos = filePath.find(".txt");
	size_t binPos = filePath.find(".bin");
	if (txtPos == std::string::npos && binPos == std::string::npos) {
		printf("Usage: file format \"*.txt\" or \"*.bin\"\n");
		exit(1);
	}
	int k = thres.size();
	int n = thres[0].size();
	FILE* in = NULL;
	if (binPos != std::string::npos) in = fopen(filePath.c_str(), "wb");    //以二进制写入模式打开文件
	else in = fopen(filePath.c_str(), "w");        //以文本写入模式打开文件

	if (in == NULL) {
		printf("No such file: %s\n", filePath.c_str());
		exit(1);
	}
	if (binPos != std::string::npos) {
		fwrite(&k, sizeof(int), 1, in);
		//fwrite(&n, sizeof(int), 1, in);
		for (int i = 0; i < k; i++) {
			//fwrite(&i, sizeof(int), 1, in);
			for (int j = 0; j < n; j++) {
				fwrite(&thres[i][j], sizeof(double), 1, in);
			}
		}
	}
	else {
		for (int i = 0; i < k; i++) {	//txt文件 - 可读性更强
			fprintf(in, "%d :\n", i);
			for (int j = 0; j < n; j++) {
				fprintf(in, "%d :", j);
				fprintf(in, " %lf\n", thres[i][j]);
			}
		}
	}
	fclose(in);
}


//使用bin文件读取thres
void load_thres(const std::string& filePath, std::vector<std::vector<double> >& thres) {
	double tm = omp_get_wtime();
	size_t binPos = filePath.find(".bin");
	if (binPos == std::string::npos) {
		printf("Usage: file format \"*.bin\"\n");
		exit(1);
	}

	int kk = 0;
	FILE* in = NULL;
	in = fopen(filePath.c_str(), "rb");
	if (in == NULL) {
		printf("No such file: %s\n", filePath.c_str());
		exit(1);
	}

	size_t x = fread(&kk, sizeof(int), 1, in);
	cout << "read kmax:" << kk << endl;
	for (int i = 0; i < thres.size(); i++) {
		for (int j = 0; j < thres[0].size(); j++) {
			x = fread(&thres[i][j], sizeof(double), 1, in);
		}
	}
	tm = omp_get_wtime() - tm;
	std::cout << "读取初始thres文件---- time：" << tm << endl;
}



//插入边的初始更新算法
vector<vector<double> > insertEdges(string infile, string parentPath, double scale) {
	Uncertain_Core uc;

	//从原数据文件 随机选择 初始边 构建初始图
	uc.edge_selected_bin(infile, scale);	//scale - unselected

	int n = uc.get_nm();
	uc.get_core();
	int kmax = uc.get_kmax();
	cout << "kmax:" << kmax << endl;
	vector<vector<double> > thres(kmax, vector<double>(n));
	uc.Initial_threshold_compute_map(thres);

	//std::string thresInFile = parentPath + "/datas/decomposition/d-initial-Flickr.bin";	//读取分解文件
	//load_thres(thresInFile, thres);	

	//确定图插入边
	//uc.insert_core_compare();
	uc.insert_threshold_compare(thres);
	return thres;
}


//删除边的初始更新算法
void deleteEdges(string readfile, string parentPath, double scale) {
	Uncertain_Core uc;
	uc.read_bin(readfile);
	int n = uc.get_nm();
	uc.get_core();
	int kmax = uc.get_kmax();
	vector<vector<double> > thres(kmax, vector<double>(n));	//初始计算 —— 需存储，then读取即可

	//初始计算thres	
	uc.Initial_threshold_compute_map(thres);
	//std::string thresInFile = parentPath + "/maintence datas/d-initial-try.bin";
	//printf_thres(thresInFile, thres);

	//读取thres	
	//load_thres(thresInFile, thres);

	uc.delete_threshold_compare(thres, scale);
	//uc.delete_compare_range(thres, scale);

	/*std::string outfile = parentPath + "/maintence datas/d-read-try.txt";
	printf_thres(outfile, thres);*/
}


int main()
{
	// 获取 main.cpp 的路径、main所在目录，数据文件相对于main.cpp 的相对路径、构建数据文件的绝对路径
	std::string currentPath = __FILE__;
	std::string parentPath = getParentDirectory(currentPath);
	std::string dataFilePath = "datas/Fruit-Fly.bin";	
	std::string infile = parentPath + "/" + dataFilePath;

	double scale = 0.8;		//选择随机边
	insertEdges(infile, parentPath, scale);	
	//deleteEdges(infile, parentPath, scale);
	return 0;
	
	
	Probability_Core pc;
	//pc.creat_bin(infile);	//创建bin文件

	pc.read_bin(infile);
	int n = pc.get_nm();

	double tm = omp_get_wtime();
	pc.get_core();
	int kmax = pc.get_kmax();
	vector<vector<double> > thres(kmax, vector<double>(n));	//初始计算 —— 需存储，then读取即可
	pc.Initial_threshold_compute_map(thres);
	tm = omp_get_wtime() - tm;
	std::cout << "Time:" << tm << endl;


	//pc.increase_threshold_compare(thres, 0.8, 1000);	//边缘概率增加
	//pc.decrease_threshold_compare(thres, 0.8, 1000);	//边缘概率减小


	//save decomposition result
	/*std::string outfile = parentPath + "/datas/d-Flickr.txt";
	printf_thres(outfile, thres);*/


	system("pause");
	return 0;
}