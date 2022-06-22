#include"Dis0.h"
#include"Aceleracion.h"
using namespace std;


void SPH(int search, int Ptype, double x, int Dtype, int N, int hf, double g, double R[], double m[], double h[], double V[], double D[], double Dx[], double Dxx[], double Dxxx[], double Pxx[], double Aq[], double Agp[], double Av[], double A[], double Zh[], double Omega[]){
	double Na, Nb, Vc;
	double xmin=-x, xmax=x;
	if(Dtype==1){
		if(Ptype==1||Ptype==2||Ptype==3){
			grid(2,N,m,R); //Distribución analitica mediante un Grid
			}
		if(Ptype==4){
			grid(6,N,m,R);
			}
		}
	if(Dtype==2){
			Glasslike(N, xmin, xmax, R, m); //Distribución equidistante
		}
	if(Dtype==3){
			UniformD(N,xmin,xmax,R,m); //Distribución Uniforme Aleatoria
		}
	
	//Suavizado(2,1, N, m, R, D, h, Zh, Omega);
	for(int i=0; i<N;++i){
		//h[i]=1.0;
		//h[i]=200.0/(double)N; only for initial data
		h[i]=(double)hf/(double)N;	//for BEC h=400 for 005/////---> h=620/N for HO-007 and h=650/N for BEC-007//----->for 009 we use h=500/N, 490.0/N for HO	
//		printf("h[i]:%lf \n",h[i]);
		//for initial data we use the standar examples 000 with h=200/N and the discretization
		//h[i]=150.0/(double)N;
		V[i]=0.0;
		}
	int keyS[N]={}, idx[N]={};
		double h_hash=(1.0/sqrt(2.0))*h[0];
		int nclass=int((xmax-xmin)/h_hash);
//		printf("pass h[0] %lf \n", h[0]);
		int idxmin[nclass]={}, idxmax[nclass]={};
        int act[nclass]={};
		int xf=5;
	if(search==1){
		Densidad0(N, m, R, h, D);
		Densidad1( N, m, R, h, D, Dx);
		Densidad2( N, m, R, h, D, Dx, Dxx);
		Pressxx( N, m, R, h, D, Dx, Dxx,Dxxx, Pxx);
		AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq);
		AceGP(N,g, m, R, h, D,Dx, Agp);
		AceV(N,m,h, R, D,Av);
	}
	if(search==2){
		CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
		Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
		Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
		Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
		Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
		AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
		AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
		//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
		//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
		AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
	}

  
  for(int i=0;i<N;++i){
    A[i] = Aq[i] + Agp[i] + Av[i];
//	cout << "Aq " << Aq[i] << "Agp " << Agp[i] << "Av " << Av[i] << '\n';
//    A[i] = Aq_2[i] + Av_2[i];
//    V[i]=A[i]/4.0; 
//	cout << " D " << D[i] << " keyS " << keyS[i] << " idx "<<  idx[i]<< " nclass " << nclass << '\n';

  }  
/*  
//  AceDamp(N,m,h,R,D, 4.0, V, Ad_2);
  for(int i=0;i<N;++i){
//  	cout << Ad_2[i] << "\t\t" << A[i] <<  '\n';
    A[i] = A[i]+Ad_2[i];
//    A[i] = Aq_2[i] + Av_2[i];
  }
  */  
}

double Qenergy(int N, double g, double m[], double R[], double V[], double D[], double Dx[],double &EKin, double &EPot, double &EQn, double &Enl){
	double E=0.0;
	EKin=0.0;
	EPot=0.0; 
	EQn=0.0;
	Enl=0.0;
	for(int i=0; i<N; i++){
	EKin=EKin+0.5*m[i]*V[i]*V[i];
	EPot=EPot+0.5*m[i]*R[i]*R[i]; 
	EQn=EQn+0.5*m[i]*0.25*Dx[i]*Dx[i]/(D[i]*D[i]);
	Enl=Enl+0.5*m[i]*g*D[i];	
		}
	return E=EKin+EPot+EQn+Enl;
	}
	
double ChePotential(int N, double g, double m[], double R[], double V[], double D[], double Dx[]){
	double Mu=0.0;
	for(int i=0; i<N; i++){
		Mu=Mu+0.5*m[i]*(V[i]*V[i]+R[i]*R[i]+0.25*Dx[i]*Dx[i]/(D[i]*D[i])+2.0*g*D[i]);
		}
	return Mu;
	}

void Qenergyi(int N, double m[], double R[], double V[], double D[], double Dx[], double E[]){
	for(int i=0; i<N; i++){
		E[i]=0.0;
		E[i]=0.5*m[i]*(V[i]*V[i]+R[i]*R[i]+0.25*Dx[i]*Dx[i]/(D[i]*D[i]));
		}
	}
