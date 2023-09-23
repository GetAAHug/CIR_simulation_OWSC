#include"initial_head.h"

static omp_lock_t lock;

void main(void)
{   
	double beta_T = 20*Pi/180;//发散角 20
    double sita_T = 60*Pi/180;//发射顶角 仰角=90-顶角 
	double fai_T = -90*Pi/180;//发射方位角
	double beta_R = 90*Pi/180;//视场角 90 
	double sita_R = 30*Pi/180;//接收顶角 仰角=90-顶角
	double fai_R = 90*Pi/180;//接收方位角
	double Sr = 1.77*1e-4;//接收机面积，单位m^2 
	double tbin = 1e-9;//时域信道冲激响应采样间隔

	
	double radius, density;   //大气浮质半径和密度

	double r;            //发送距离


	bool location = true;// 城市，乡村为false

	double N = 1e8;          //仿真光子数

	cout << "input the number of simulated photons:\n";
	cin >> N;

	int n = 1; //散射次数
	double epsilon = 1e-3;     //以下三个参数和仿真精度有关，不用调整
	double delta = 1e-6;
	int M = 50;

	time_t t_beg;
	char tstr_beg[32];
	t_beg = time(NULL);
	ctime_s(tstr_beg, 32, &t_beg);
	cout << "beginning time:\n" << tstr_beg << endl;

	srand((unsigned)time(0));     //粗生成随机数种子


	const char*filetext_IRF_overall = "C:\\Users\\lby\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\IRF_example.txt";

	int index = 0;

	for (int LE_num = 0; LE_num < 2; LE_num++)          //共仿真的信道环境参数个数
	{

		radius = 0.02 + rand() / (RAND_MAX + 0.0) * 0.1;               //随机选取大气浮质半径
		density = pow(10, 3 + rand() / (RAND_MAX + 0.0) * 1.8);        //随机选取大气浮质分布密度
		r = 100 + rand() / (RAND_MAX + 0.0) * 400;                     //随机选取发送距离
		sita_T = (80.0 - 30.0 * rand() / (RAND_MAX + 0.0))*Pi / 180;   //随机选取发送顶角
		sita_R = (80.0 - 30.0 * rand() / (RAND_MAX + 0.0))*Pi / 180;   //随机选取接收顶角

		index = index + 1;
		cout << "Index:  " << index  << endl;
		cout <<  "distance = " << r << ",  " << "sita_T = " << 90 - sita_T * 180 / Pi << ",  " << "sita_R = " << 90 - sita_R * 180 / Pi << ",  " << "radius = " << radius << ",  " << "density = " << density << endl;

		for (int lam = 0; lam < 2; lam++)               //共仿真的波长数
		{
			double lambda;
			lambda = 0.2 + 0.04*lam;

			double Ke, Ka, Ks;
			double Ksray, Ksmie, Kamie, g;
			double* atmos = new double[4];

			Ksray = ray_scat(lambda);
			double* K_mie = new double[3];
			K_mie = mie_scat(lambda, radius, density, location);
			Ksmie = K_mie[0];
			Kamie = K_mie[1];
			g = K_mie[2];
			Ks = Ksray + Ksmie;
			Ka = Kamie;
			Ke = Ksray + Ksmie + Kamie;

			atmos[0] = Ksray;
			atmos[1] = Ksmie;
			atmos[2] = Kamie;
			atmos[3] = g;

			double* geo = new double[7];
			geo[0] = sita_T;
			geo[1] = sita_R;
			geo[2] = beta_T;
			geo[3] = beta_R;
			geo[4] = fai_T;
			geo[5] = fai_R;
			geo[6] = r;

			double* parameters = new double[2];
			parameters[0] = epsilon;
			parameters[1] = delta;

			cout << "lambda = " << lambda << " um;" << endl;

			const int length = ((5e-6) / (1e-9) + 1);    //注意该行定义了时间间隔，且1e-9需要与tbin相同
			double* T_i = new double[length];
			for (int l = 0; l < length; l++)
			{
				T_i[l] = l*tbin;
			}
			double(*h_i_j)[length] = new double[5][length];
			double*h_j = new double[length];

			//赋初值
			for (int col = 0; col < length; col++)
			{
				h_j[col] = 0;
				for (int row = 0; row < n + 1; row++)
				{
					h_i_j[row][col] = 0;
				}
			}

			double PD_i[5] = { 0 };
			double u_T[3] = { cos(fai_T)*sin(sita_T),sin(fai_T)*sin(sita_T),cos(sita_T) };
			double u_R[3] = { cos(fai_R)*sin(sita_R),sin(fai_R)*sin(sita_R),cos(sita_R) };
			double r_R = sqrt(Sr / Pi);
			const double c = 2.997046e8;

			if ((-u_T[1] >= cos(beta_T / 2)) && (u_R[1] >= cos(beta_R / 2)))
			{
				PD_i[0] = exp(-Ke*r)*(1 - r*pow((pow(r, 2) + pow(r_R, 2)), (-1.0 / 2))) / (1 - cos(beta_T / 2));
				h_i_j[0][compute_j(r / c, T_i, length)] = PD_i[0];
			}
			int i, j;
			double location[3], u[3], u_new[3], location_new[3], location_minus_new[3];
			double d_total, p, alpha, psi, temp, d, Pi_R, PD_cur, t;
			bool flag;

			for (int k = 1; k <= N; k++)
			{
				location[0] = 0;
				location[1] = r;
				location[2] = 0;
				u[0] = u_T[0];
				u[1] = u_T[1];
				u[2] = u_T[2];
				d_total = 0;
				p = 1 - PD_i[0];
				i = 0;
				while ((i < n) && (p >= epsilon))
				{
					i = i + 1;
					flag = false;
					while (!flag)
					{
						if (i == 1)
						{
							alpha = acos(1 - rand() / (RAND_MAX + 0.0)*(1 - cos(beta_T / 2)));
							psi = 2 * Pi*rand() / (RAND_MAX + 0.0);
						}
						else
						{
							alpha = scat_angle(M, delta, atmos);
							psi = 2 * Pi*rand() / (RAND_MAX + 0.0);
						}
						temp = sqrt(1 - pow(u[2], 2));
						if (!(temp == 0.0))
						{
							u_new[0] = sin(alpha)*(u[0] * u[2] * cos(psi) - u[1] * sin(psi)) / temp + u[0] * cos(alpha);
							u_new[1] = sin(alpha)*(u[1] * u[2] * cos(psi) + u[0] * sin(psi)) / temp + u[1] * cos(alpha);
							u_new[2] = -sin(alpha)*cos(psi)*temp + u[2] * cos(alpha);
						}
						else
						{
							u_new[0] = sin(alpha)*cos(psi);
							u_new[1] = sin(alpha)*sin(psi);
							u_new[2] = u[2] * cos(alpha);
						}
						d = -log(1 - rand() / (RAND_MAX + 0.0)) / Ke;
						location_new[0] = location[0] + d*u_new[0];
						location_new[1] = location[1] + d*u_new[1];
						location_new[2] = location[2] + d*u_new[2];
						location_minus_new[0] = location[0] - location_new[0];
						location_minus_new[1] = location[1] - location_new[1];
						location_minus_new[2] = location[2] - location_new[2];
						double* det = new double[3];
						if ((((norm(cross(location, location_new, det))) / norm(location_minus_new)) > r_R) || (dot(location, u_new) >= 0) || (d < norm(location)))
							flag = true;
						delete[]det;
					}
					location[0] = location_new[0];
					location[1] = location_new[1];
					location[2] = location_new[2];
					u[0] = u_new[0];
					u[1] = u_new[1];
					u[2] = u_new[2];
					d_total = d_total + d;
					Pi_R = min(1.0, 2 * Pi*(1 - norm(location)*pow((pow(norm(location), 2) + pow(r_R, 2)), -1.0 / 2))
						*p_total(-dot(location, u) / norm(location), atmos))*exp(-Ke*norm(location));
					if ((dot(location, u_R) / norm(location)) >= cos(beta_R / 2))
					{
						PD_cur = p * Ks / Ke*Pi_R;
						PD_i[i] = PD_i[i] + PD_cur;
						t = (d_total + norm(location)) / c;
						j = compute_j(t, T_i, length);
						if (j != 0)              //如果没找到光子，放弃赋值
						{
							h_i_j[i][j] = h_i_j[i][j] + PD_cur;
						}
					}
					p = p * Ks / Ke*(1 - Pi_R);
				}
			}

			
			double PD_overall = PD_i[0], PL_overall;
			for (int m = 1; m < n + 1; m++)
			{
				PD_i[m] = PD_i[m] / N;
				PD_overall = PD_overall + PD_i[m];
			}

			PL_overall = 10 * log10(1.0 / PD_overall);

			/*
			PL_first = 10 * log10(1.0 / PD_i[1]);
			PL_second = 10 * log10(1.0 / PD_i[2]);
			PL_third=10*log10(1.0/PD_i[3]);
			PL_fourth=10*log10(1.0/PD_i[4]);
			PL_higher=10*log10(1.0/(PD_i[2]+PD_i[3]+PD_i[4]));
			*/
			for (int col = 0; col < length; col++)
			{
				h_j[col] = *(h_i_j[0] + col);
			}

			for (int col = 0; col < length; col++)
			{
				for (int row = 1; row < n + 1; row++)
				{
					*(h_i_j[row] + col) = *(h_i_j[row] + col) / N;
					h_j[col] = h_j[col] + (*(h_i_j[row] + col));
				}
			}

			//cout<<"单次散射路径损耗为"<<PL_first<<"dB"<<endl;
			//cout<<"多次路径损耗为"<<PL_higher<<"dB"<<endl;
			/*
			const char*filetxt_PL_overall="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\PassLoss_overall.txt";
			const char*filetxt_PL_first="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\PassLoss_first.txt";
			const char*filetxt_PL_second="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\PassLoss_second.txt";
			const char*filetxt_PL_third="E:\\PassLoss_third.txt";
			const char*filetxt_PL_fourth="E:\\PassLoss_fourth.txt";
			const char*filetxt_PL_higher="E:\\PassLoss_higher.txt";
			*/
			//const char*filetext_IRF_overall = "C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\IRF_overall.txt";
			/*
			const char*filetext_IRF_single="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\IRF_single.txt";
			const char*filetext_IRF_multiple="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\IRF_multiple.txt";
			const char*filetext_IRF_second="C:\\Users\\1004\\Desktop\\lby\\Paper\\宽光谱\\MatlabSim\\散射信道模型\\UV_channel\\IRF_second.txt";
			const char*filetext_IRF_third="E:\\IRF_third.txt";
			const char*filetext_IRF_fourth="E:\\IRF_fourth.txt";
			*/
			//fstream ftxt;
			/*ftxt.open(filetxt_PL_overall,ios::app);
			ftxt<<PL_overall<<endl;
			ftxt.close();
			ftxt.open(filetxt_PL_first,ios::app);
			ftxt<<PL_first<<endl;
			ftxt.close();
			ftxt.open(filetxt_PL_second,ios::app);
			ftxt<<PL_second<<endl;
			ftxt.close();
			ftxt.open(filetxt_PL_third,ios::app);
			ftxt<<PL_third<<endl;
			ftxt.close();
			ftxt.open(filetxt_PL_fourth,ios::app);
			ftxt<<PL_fourth<<endl;
			ftxt.close();
			ftxt.open(filetxt_PL_higher,ios::app);
			ftxt<<PL_higher<<endl;
			ftxt.close();
			*/

			double lambda0 = lambda;
			fstream ftxt;
			ftxt.open(filetext_IRF_overall, ios::app);
			ftxt << LE_num << "\t" << lambda << "\t" << radius << "\t" << density << "\t" << r << "\t" << 90 - sita_T * 180 / Pi << "\t" << 90 - sita_R * 180 / Pi << "\t" << PL_overall << "\t";

			for (int i = 0; i < length; i++)
			{
				ftxt << h_j[i] / Sr / tbin << '\t';
			}
			ftxt << endl;
			ftxt.close();

			/*
				ftxt.open(filetext_IRF_single,ios::app);
				for (int i = 0; i < length; i++)
				{
					ftxt<<(*(h_i_j[1]+i))/ (1e-8) / Sr <<endl;
				}
				ftxt.close();

				ftxt.open(filetext_IRF_multiple,ios::app);
				for (int i = 0; i < length; i++)
				{
					ftxt<<(*(h_i_j[2]+i)+*(h_i_j[3]+i)+*(h_i_j[4]+i))/ (1e-8) / Sr <<endl;
				}
				ftxt.close();

				ftxt.open(filetext_IRF_second,ios::app);
				for (int i = 0; i < length; i++)
				{
					ftxt<<(*(h_i_j[2]+i))/ (1e-8) / Sr <<endl;
				}
				ftxt.close();

				ftxt.open(filetext_IRF_third,ios::app);
				for (int i = 0; i < length; i++)
				{
					ftxt<<(*(h_i_j[3]+i))/ (1e-8) / Sr <<endl;
				}
				ftxt.close();

				ftxt.open(filetext_IRF_fourth,ios::app);
				for (int i = 0; i < length; i++)
				{
					ftxt<<(*(h_i_j[4]+i))/ (1e-8) / Sr <<endl;
				}
				ftxt.close();
				*/
			delete[]T_i;
			delete[]h_j;
			delete[]h_i_j;
			delete[]atmos;
			delete[]K_mie;
		}
	}

	time_t t_end;
	char tstr_end[32];
	t_end=time(NULL);
	ctime_s(tstr_end,32,&t_end);
	cout<<"ending time:\n"<<tstr_end<<endl;
	system("pause");
}




		 
