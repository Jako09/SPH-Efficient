//#include<bits/stdc++.h>
#include"hAdaptative.h" 

using namespace std;

void Masa(int N, double m[]){
	for(int i=0; i<N; i++){
		m[i]=0.0;
		m[i]=1.0/(double)N;
	}
}

void Suavizado(int Ptype, int Stype, int N, double m[], double R[], double D[], double h[], double Zh[] ,double dZh[]){
	if(Stype==1){
	double q0=1.4;
		for(int i=0; i<N; i++){
			h[i]=0.0;
			h[i]=q0*m[i]/DA(2,R[i]);
		}
	}
	if(Stype==2){	
	double q0=1.1;
		for(int i=0; i<N; i++){
			h[i]=0.0;
			h[i]=q0*m[i]/D[i];
		}
	}
	if(Stype==3){
		NRh("hAdaptative", N, m, R, h, D, Zh , dZh);
	}
}
double alcance(double tol, double h){
	double q=sqrt(-log(tol*h*sqrt(PI)));	
	cout << q << '\n';
	return q*h;
}
/*
void hAdaptative(int N, double h[]){
	for(int i=0; i<N; i++){
		h[i]+=-
	}
}
*/
