#include <iostream>

#include <cstdlib>

#include <cmath>

#include <ctime>

#include <cmath>

#include <fstream>

using namespace std;

const double M_PI = 3.14159365368;

double L, t, tw, w, wp, E0 ,m ,I, F, stdev, kse, T; 

const double g = 9.81;

const double P0 = 100.;

int n, j, inc, num, nsim, nt, mesh;

double unif() {
return rand() / double(RAND_MAX);
 }

 double unif(double a, double b) {
 return (b-a)*unif() + a;
 }

long unif(long a) {
 if (a < 0) a = -a;
 if (a==0) return 0;
 return long(unif()*a) + 1;
 }

 double randn() {
 static bool has_saved = false;
 static double saved;

 if (has_saved) {
 has_saved = false;
 return saved;
 }

 double x,y,r;
 do {
 x = unif(-1.,1.);
 y = unif(-1.,1.);
 r = x*x + y*y;
 } while (r >= 1.);

 double mu = sqrt(-2.0 * log(r) / r);

 saved = mu*y;
 has_saved = true;

 return mu*x;
 }

 double SDF(double mr, double st)
 {
	return 250. * st*st * pow(mr, 2.) * exp(-10. * mr);
 }
void main()
{
	const char* output_file_name_1 = "displacement.out";
	const char* output_file_name_2 = "variance.out";
	const char* output_file_name_3 = "skewu.out";
	const char* output_file_name_4 = "kurtu.out";
	const char* output_file_name_5 = "u.out";
	ofstream my_out_1(output_file_name_1);
	ofstream my_out_2(output_file_name_2);
	ofstream my_out_3(output_file_name_3);
	ofstream my_out_4(output_file_name_4);
	ofstream my_out_5(output_file_name_5);
	if (my_out_1.fail() || my_out_2.fail()|| my_out_3.fail()|| my_out_4.fail() || my_out_5.fail()) {
		cerr << "Unable to open output or input file " << endl;
	}


// PARAMETER SPECIFICATION

/*	cout << "DEFINE WAVE LENGTH DISCRETIZATION" << endl;
	cin >>*/ num = 200;

	/*cout << "DEFINE NUMBER OF SIMULATIONS" << endl;
	cin >> */nsim = 50000;

	/*cout << "DEFINE STOCHASTIC FIELD STANDARD DEVIATION" << endl;
	cin >> */stdev = .6;

	/*cout << "DEFINE HORIZONTAL AXIS MESH" << endl;
	cin >>*/ mesh = 150;

	/*cout << "DEFINE ELASTIC MODULUS MEAN" << endl;
	cin >>*/ E0= 1.25e+08;

	/*cout << "DEFINE BEAM MOMENT OF INERTIA" << endl;
	cin >>*/ I = .1;

	/*cout << "DEFINE DAMPING RATIO" << endl;
	cin >>*/ kse = 0.05;

	/*cout << "DEFINE TIME WINDOW" << endl;
	cin >>*/ tw = 10.;

	/*cout << "DEFINE NUMBER OF TIME STEPS" << endl;
	cin >>*/ nt = 1000;

	/*cout << "DEFINE LENGTH OF BEAM IN METERS" << endl;
	cin >> */L = 4.;

	cout << "DEFINE LOAD CYCLICAL FREQUENCY" << endl;
	cin >> wp;

	cout << "DEGINE MODEL EIGEN PERIOD IN SEC" << endl;
	cin >> T;

	m = 3. / 4. / M_PI/M_PI / pow(L,3) * T*T * E0 *I;
	
	// ARRAY MEMORY ALLOCATION

	double u, *phi, *sumu, *sumu2, *sumu3, *sumu4, *mu, *varu, *eskewu, *ekurtu; 

	sumu = new double[nt];
	sumu2 = new double[nt];
	sumu3 = new double[nt];
	sumu4 = new double[nt];
	mu = new double[nt];
	varu = new double[nt];
	eskewu = new double[nt];
	ekurtu = new double[nt];
	phi = new double[num];
	
	double dx = L/float(mesh);
	double dk = 1.22718464/float(num);
	double dt = tw/double(nt);
	double lamda = 1.8751/L;
	double x, k, Ke, sumust, sumust2, beta, sumw, mw;

	// RANDOM SEEDS' GENERATION & TRANSLATION IN DISPLACEMENT VALUES
	sumw = 0.;
	for (int ij=1; ij<=nt; ij++)
	{sumu[ij-1]=0.; sumu2[ij-1]=0.; sumu3[ij-1]=0.; sumu4[ij-1]=0.;
	mu[ij-1]=0.; varu[ij-1]=0.; eskewu[ij-1]=0.; ekurtu[ij-1]=0.;}
	
	sumust = 0.;
	sumust2 = 0.;

	for (j=1; j<=nsim; j++)
	{
		for (int i=1; i<=num; i++)
		{
			phi[i-1] = unif(0., (2.*M_PI));
		}
		
		x= dx / 2.; Ke= 0.; double u_st=0.;
	
		for (int jj=1; jj<=mesh; jj++)
		{
			double f = 0.;
			k = dk/2.;
			
			for (int i=1; i<=num; i++)
			{
				f += sqrt(2.) * sqrt(2. * SDF(k, stdev) * dk) * double(cos(k * x + phi[i-1]));
				k+= dk;
			}

			if (nsim==1)
			{
				f =0.;
			}
			else if ((abs(f) >= .9) && (f <= 0.))
			{
				f = -0.9;
			}
			else if ((abs(f) >= .9) && (f >= 0.))
			{
				f = 0.9;
			}
		
			//my_out_1 << f << endl;
			u_st += (L-x) * (L-x) * (1+f) /E0/I * dx;
			x += dx;
		}
		sumust += u_st;
		sumust2 += u_st*u_st;
		Ke = 1. / u_st;
		w = sqrt(abs(Ke / m));
		sumw += w;
		beta = wp / w;
		double wd = w*sqrt(1-kse*kse);
		t=dt/2.;
		
		for (int ij=1; ij<=nt; ij++)
		{
			double omikron = 1./(pow((1.-beta*beta), 2)+pow((2.*kse*beta),2));
			double X = P0/Ke*omikron*(3.*beta*beta-1.)*kse/sqrt(1.-kse*kse);
			double Y = -P0/Ke*omikron*(1.-beta*beta);
			double C1 = P0 / Ke * omikron *2*kse*beta;
			double C2 = P0 / Ke * omikron *(1-beta*beta);
			double u0 = exp(-kse*w*t)*(X*sin(wd*t)+Y*cos(wd*t));
			double up = C1*sin(wp*t) + C2*cos(wp*t);
			u = u0 + up;
			sumu[ij-1] += u;
			sumu2[ij-1] += u*u;
			sumu3[ij-1] += pow(u, 3);
			sumu4[ij-1] += pow(u, 4);
			t+=dt;
		}
	}
	mw = sumw / float(nsim);
	double mu_st = sumust/float(nsim);
	double varu_st = (sumust2 - mu_st*mu_st * float(nsim))/float(nsim-1);
	my_out_5 << mw << endl;
	my_out_5 << mu_st << endl;
	my_out_5 << varu_st << endl;

	for (int ij=1; ij<=nt; ij++)
	{
		mu[ij-1] = sumu[ij-1] / float (nsim);
		varu[ij-1] = sumu2[ij-1] / double(nsim) - mu[ij-1] * mu[ij-1] ;
		eskewu[ij-1] = sumu3[ij-1] / float(nsim);
		ekurtu[ij-1] = sumu4[ij-1] / float(nsim);
		my_out_1 << mu[ij-1] << endl;
		my_out_2 << varu[ij-1] << endl;
		my_out_3 << eskewu[ij-1] << endl;
		my_out_4 << ekurtu[ij-1] << endl;
	}

	delete[] sumu;
	delete[] sumu2;
	delete[] sumu3;
	delete[] sumu4;
	delete[] mu;
	delete[] varu;
	delete[] eskewu;
	delete[] ekurtu;
	delete[] phi;
}