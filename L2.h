#include<bits/stdc++.h>

using namespace std;

double L2(int Ptype, int N, double R[], double D[]){
	double z;
	for(int i=0; i<N; ++i){
//		if(i>0.25*N && i<0.75*N){
		double e=D[i]-DA(Ptype, R[i]);
		z+=e*e;
//		}
	}
	return sqrt(z);
}
/*
void errorDensidad(int Ptype, int Stype, int Gtype, int N){
	double m[N], R[N], h[N], D[N];
	double Na, Nb;
	bool ERes=false, NPV=false;
	Na=(double)N;
	if(NPV==false){
		Nb=(double)N;
	}else{
	}
  	if(Gtype==0){
  		Masa(N, m);
  		Dis0(Ptype, N, m, R);
	}
 	if(Gtype==1){
  		grid(Ptype, N, m, R);
  	}
  	Suavizado(Ptype, Stype, N, m, R, h);
	Densidad0(ERes, N, m, R, h, D);
	ofstream file("ErrorDensity.dat");
	for(int i=0; i<N; i++){
		cout << DA(Ptype, R[i]) << " " << D[i] << '\n';
		double ED=abs(D[i]-DA(Ptype, R[i]));
		file << R[i] << " " << ED << '\n';
	}
	file.close();
}
*/ 
/*
void error(int Ptype, int Stype, int Gtype, int NError){
	int N=100;
	ofstream file("Error.dat");
	for(int i=0; i<NError; ++i){
		double m[N], R[N], h[N], D[N];
		bool ERes=false, NPV=false;
  		if(Gtype==0){
  			Masa(N, m);
  			Dis0(Ptype, N, m, R);
		}
 		if(Gtype==1){
  			grid(Ptype, N, m, R);
  		}
		Suavizado(Ptype, Stype, N, m, R, h);
		Densidad0(ERes, N, m, R, h, D);
		file << N << " " << L2(Ptype, N, R, D) << '\n';
		cout << N << '\n';
		N+=100;
	}
	file.close();
}
*/
