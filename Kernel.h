//#include<bits/stdc++.h>
#define PI 3.1415926535897932
#define h1 1.5
using namespace std;
double ker(double hi, double ri, double rj){
	double alpha=(1.0/(hi*sqrt(PI)));
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R);
}
double dker(double hi, double ri, double rj){
	double alpha=(1.0/(hi*sqrt(PI))), beta=(2.0/hi);
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(double hi, double ri, double rj){
	double alpha=(1.0/(hi*sqrt(PI))), beta=(2.0/hi), delta=2.0/(hi*hi);
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
double dddker(double hi, double ri, double rj){
	double alpha=(1.0/(hi*sqrt(PI))), gamma=4.0/(hi*hi*hi);   
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*gamma*R*(-2.0*R*R+3.0);
	}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////Kernel for 2D	
/*
double ker2D(double hi, double xi, double xj, double Rj){
	double alpha=(1.0/(hi*sqrt(PI)));
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R);
}
double dker2D(double hi, double xi, double xj, double yi, double yj){
	double alpha=(1.0/(hi*sqrt(PI))), beta=(2.0/hi);
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*(-R*beta);
}
double ddker(double hi, double xi, double xj, double yi, double yj){
	double alpha=(1.0/(hi*sqrt(PI))), beta=(2.0/hi), delta=2.0/(hi*hi);
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*(R*R*beta*beta-delta);
}
double dddker(double hi, double xi, double xj, double yi, double yj){
	double alpha=(1.0/(hi*sqrt(PI))), gamma=4.0/(hi*hi*hi);   
	double R=(ri-rj)/hi;
	return alpha*exp(-R*R)*gamma*R*(-2.0*R*R+3.0);
	}
*/
