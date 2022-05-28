#define PI 3.1415926535897932 
#include"Densidad.h"
using namespace std;

void DimF(string Str, int N, double m[], double h[], double R[], double D[], double Zh[]){
	if(Str=="hAdaptative"){
		double eta=1.4; //Parámetro a variar
		for(int i=0; i<N;i++){
			Zh[i]=0.0;
			D[i]=0.0;
			double etap=(eta/h[i]);
			for(int j=0; j<N;j++){
//				D[i]=DA(2,R[i]);
				D[i]+=m[j]*ker(h[i], R[i], R[j]); //esto es para la densidad no analítica
			}
			Zh[i]=m[i]*(etap)-D[i];	
		}
	}
}
void Function(string Str, int N, int i, double m[], double h[], double R[], double D[], double Zh[]){
	if(Str=="hAdaptative"){
		double eta=1.4;
			Zh[i]=0.0;
			D[i]=0.0;
			double etap=(eta/h[i]);
			for(int j=0; j<N;j++){
//				D[i]=DA(2,R[i]);
				D[i]+=m[j]*ker(h[i], R[i], R[j]); //esto es para la densidad no analítica
			}
			Zh[i]=m[i]*(etap)-D[i];	
	}
}
double dhker(double xi, double xj, double hi){
	double R=(xi-xj)/hi;
	return -(1.0-2.0*R*R)*exp(-R*R)/(hi*hi*sqrt(PI));
}
void Omegai(string Str, int N, double m[], double x[], double h[], double D[], double Omega[]){
	if(Str=="hAdaptative"){
		double fact1;
		for(int i=0; i<N; i++){
			Omega[i]=0.0;
			Omega[i]=1.0;
			fact1=h[i]/D[i];
			for(int j=0; j<N; j++){
				Omega[i]+=fact1*m[j]*dhker(x[i],x[j],h[i]);
			}
		}	
	}
}
void Omegaiaux(string Str, int N, int i, double m[], double x[], double h[], double D[], double Omega[]){
	if(Str=="hAdaptative"){
		double fact1;
			Omega[i]=0.0;
			Omega[i]=1.0;
			fact1=h[i]/D[i];
			for(int j=0; j<N; j++){
				Omega[i]+=fact1*m[j]*dhker(x[i],x[j],h[i]);
			}	
	}
}
void NRh(string Str, int N, double m[], double x[], double h[], double D[], double Zh[] ,double Omega[]){
	double EPS=1.0e-3;
	double hnew;
		DimF(Str, N, m, h, x, D, Zh);
		Omegai(Str, N, m, x, h, D, Omega);
//		cout << "pass" << '\n';
	for(int i=0; i<N; i++){
		double e=0.0,b=0.0, hnew=0.0, h0=h[i];
		e=h[i]*Zh[i]/(D[i]*Omega[i]);
		hnew=h[i]+e;
		b=abs(hnew-h[i])/h0;
//		b=abs(e/h[i]);
//		cout << Zh[i] << " " << Omega[i] << " " << abs(b) <<'\n';
		while(b>EPS){
//			hold=h[i];
			Function(Str, N,i, m, h, x,  D, Zh);
			Omegaiaux(Str, N,i, m, x, h, D, Omega);
			e=h[i]*Zh[i]/(D[i]*Omega[i]);
			hnew=h[i]+e;
//			b=abs(e/h[i]);
			b=abs(hnew-h[i])/h0;
			h[i]=hnew;
			cout << Zh[i] << "\t" << i << "\t" << Omega[i] << "\t" << b <<   '\n';
		}
	}
}

