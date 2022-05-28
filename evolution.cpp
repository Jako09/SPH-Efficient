#include <iostream>
#include"SPH.h"
#include"EulerInt.h"
#include"LeapFrog.h"
#include"gnu.h"
#include<chrono>


using namespace std;

int main(){
	//We define the functions to follow the time
	using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
	// ****Initial variables to external modification****
	int Ptype; // Ptype:=Problem type. For this case we have -Harmonic Oscillator->"-HarmOsc", -Dynamics of HO->"HODynamic", -BEC ground states->"GSBEC", -Solitons Collition->"SolCol"
	int Graphtype; // Graphtype:= Graph type. We have two types of graph, Animation and final results.
	int N; //# of particles
	double g, DV; // nonlinear parameter.
	int Dtype;// initial distribution type
	string search;
	cout << "Which search engine you want to use? standar or efficient" << '\n';
	cin >> search;
	cout << "Which problem do you want?" << '\n';
	cout << "Problems:" << '\n';
	cout << "1.-Harmonic Oscillator" << '\n';
	cout << "2.-Harmonic Oscillator Dynamics" << '\n';
	cout << "3.-Groud States of Bose-Einstein Condensate" << '\n';
	cout << "4.-Soliton Collision" << '\n';
	cout << "Write the number" << '\n';
	cin >> Ptype;
	while((Ptype!=1)&&(Ptype!=2)&&(Ptype!=3)&&(Ptype!=4)){
			cout << "Wrong number! Writte the right number"<<'\n';
			cin >> Ptype;
			}
	cout << "You choose the problem " << to_string(Ptype)  << '\n';
	cout << "How many particles do you want to use?" << '\n';
	cin >> N;
	while(N<=0){
		cout << "The number of particles need to be a positive, try again:"<< '\n';
		cin >> N;
		}
	cout << "Enter the value of damping value DV:"<< '\n';
	cin >> DV;
	while(DV<=0.0){
		cout << "The damping value need to be a positive value, write a right value:"<< '\n';
		cin >> DV;
		}
	if(Ptype==1|| Ptype==2){
		cout << "For this problem the nonlinear parameter is cero g=0.0"<< '\n';
		g=0.0;		
		}
	if(Ptype==3){
			cout << "Enter the nonlinear parameter g"<< '\n';
			cin >> g;			
	}
	if(Ptype==4){
			cout << "For this problem the value for nonlinear parameter is g=-1.0" << '\n';
			g=-1.0;
		}
	cout << "Which initial distribution do you want?" << '\n';
	cout << "1.-Trapecious type-Analytic"<< '\n';
	cout << "2.-Equidistant distribution"<< '\n';
	cout << "3.-Random distribution" << '\n';
	cin >> Dtype;
	while(Dtype!=1&&Dtype!=2&&Dtype!=3){
		cout << "The value is wrong, choose the right number" << '\n';
		cin >> Dtype;
		}
	//  ****Now we implement the SPH algorithm****
	if(Ptype==1){
		cout << "Which kind of graph do you want?" << '\n';
		cout << "1.- Relaxation to the ground state" << '\n';
		cout << "2.- Static with different times" << '\n';
		cout << "Write the number" << '\n';
		cin  >> Graphtype;
		while(Graphtype!=1&&Graphtype!=2){
			cout << "The number is incorrect, try again";
			cout << "1.- Relaxation to the ground state" << '\n';
			cout << "2.- Static with different times" << '\n';
			cout << "Write the number" << '\n';
			cin >> Graphtype;
			}
			double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
			// R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, 
			//Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
				
			double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure , 
			//Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
				
			//for time integration leap frog
			double step=4.0e-3; //lenght of step in the time 
				
			//for adaptative smoothing length
			double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function.
			int hf=490;  //hf-->h-factor. 490 es el estadar y funciona con N=640 // for 1280 particles we use hf=790
			
			//for energy
			double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
			// EQn----> energy quantum nature, Mu-------->Chemical Potential 
				
			//for error
//			printf("pass\n");
			double error=0.0;
				double xmin, xmax;
				//start the evolution

			//*****Starting the SPH*****
			SPH(search,Ptype,4.0,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
			int keyS[N]={}, idx[N]={};
			double h_hash=(1.0/sqrt(2.0))*h[0];
			int nclass=int((8.0)/h_hash); //Initial number of classes
			int xf=5;
//
				
			E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
			if(g!=0.0){
				Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
			}
			//Qenergyi(N, m, R, V, D, Dx, Ei);
				 
			//data
			string dataname1;
			string dataname2;
			if(Graphtype==1){
				if(search=="efficient"){
					dataname1="RelaxationGroundStateDEff.xxx";
					dataname2="eRelaxationGroundStateDEff.xxx";
				}
				if(search=="standar"){
					dataname1="RelaxationGroundStateD.xxx";
					dataname2="eRelaxationGroundStateD.xxx";
				}
			}
			if(Graphtype==2){
				if(search=="efficient"){
					dataname1="RelaxationGroundStateEff.xxx";
					dataname2="eRelaxationGroundStateEff.xxx";
				}
				if(search=="standar"){
					dataname1="RelaxationGroundState.xxx";
					dataname2="eRelaxationGroundState.xxx";
				}
			}
			ofstream file(dataname1); //open file to data
			ofstream fileE(dataname2); //open file to data for energy
			//print initial data 
			//file << "\n\n\n"; //print in data file the initial values
			//file1 << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
			for(int i=0; i < N; i++){
				file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
				}
			
			if(Graphtype==1){	
				int itmax=7500; // # of iterations of evolution Leap Froag
				//start the evolution
				auto start = high_resolution_clock::now();
				for(int t=0; t<itmax; t++){
					if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
					for(int i=0; i<N; i++){
							R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
							Aold[i]=A[i];
//							cout << "R[i] " << R[i] << " t " << t << '\n';
						}
					xmin=getMind(N,R)-h[0];
					xmax=getMaxd(N,R)+h[0];
					nclass=int((xmax-xmin)/h_hash);
					int idxmin[nclass]={}, idxmax[nclass]={};
					int act[nclass]={};
//					cout << " xmin "<< xmin <<" xmax "<<xmax<<" nclass "<<nclass<< '\n';
					if(search=="standar"){
						Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
						Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
						Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
						AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}
					if(search=="efficient"){
//						printf("pass if %d \n", t);
						CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
//						printf("pass censo %d \n", t);
						Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
//						printf("pass densidadEff %d \n", t);
						Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
						Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
//						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
						AceDampEff(N, nclass, xf, m, R, h, D, DV, V, Ad, keyS, idx, idxmin, idxmax, act); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}
					for(int i=0; i<N; i++){
							A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
							V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}

					if(t%50==0){//print the data only for multiples to 50
						file << "\n\n\n";
							for(int i=0; i < N; i++){
								file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
								}
							E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
							if(g!=0.0){
								Mu=ChePotential(N, g, m, R, V,  D, Dx);
								}
							fileE.open("eRelaxationGroundStateD.xxx",std::fstream::app);// necesary to write the data in the file
							fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
							fileE.close();
						}
					}
				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<microseconds>(stop - start);
				fileE.close(); // close the datafile for energy
				file.close(); //we close the datafile
				cout<< "The process has been finished and the time it takes is " << duration.count() << '\n';
				return 0;
				}
			if(Graphtype==2){
				int itmax=5000;
				int Ngraphs=5;
				double Nit[Ngraphs];
				//start the evolution
//				printf("pass\n"); //////////////
				for(int i=0; i<Ngraphs; i++){
					Nit[i]=4.0*(double)i;
					}
				//open data for static graph
				ofstream file0("RelaxationGroundState0.xxx"); //open file to data for time 1
				ofstream file1("RelaxationGroundState1.xxx"); //open file to data for time 2
				ofstream file2("RelaxationGroundState2.xxx"); //open file to data for time 3
				ofstream file3("RelaxationGroundState3.xxx"); //open file to data for time 4
				ofstream file4("RelaxationGroundState4.xxx"); //open file to data for time 5
				//start the evolution
				auto start = high_resolution_clock::now();
				for(int t=0; t<itmax; t++){
//					printf("pass %d \n", t);
					if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
					for(int i=0; i<N; i++){
							R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
							Aold[i]=A[i];
//						    cout << "R: " << R[i]  << "V: " << V[i]  << "A: " << A[i] << '\n';
						}
					xmin=getMind(N,R)-h[0];
					xmax=getMaxd(N,R)+h[0];
					nclass=int((xmax-xmin)/h_hash);
					int idxmin[nclass]={}, idxmax[nclass]={};
					int act[nclass]={};
//					cout << " xmin "<< xmin <<" xmax "<<xmax<<" nclass "<<nclass<< '\n';


						// aqui pondremos las funciones que nos ayuden a calcular los valores minimos y maximos de las R para creear las nuevas clases
					if(search=="standar"){
						Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
						Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
						Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
						AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}
					if(search=="efficient"){
//						printf("pass if %d \n", t);
						CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
//						printf("pass censo %d \n", t);
						Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
//						printf("pass densidadEff %d \n", t);
						Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
						Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
						AceDampEff(N, nclass, xf, m, R, h, D, DV, V, Ad, keyS, idx, idxmin, idxmax, act); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}

					for(int i=0; i<N; i++){
							A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
							V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}
						if(t*step==Nit[0]){//print the data only for multiples to 50
							file0 << "\n\n\n";
								for(int i=0; i < N; i++){
									file0 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eRelaxationGroundState.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==Nit[1]){
								file1 << "\n\n\n";
								for(int i=0; i < N; i++){
									file1 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eRelaxationGroundState.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==Nit[2]){
								file2 << "\n\n\n";
								for(int i=0; i < N; i++){
									file2 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eRelaxationGroundState.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==Nit[3]){
								file3 << "\n\n\n";
								for(int i=0; i < N; i++){
									file3 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eRelaxationGroundState.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==Nit[4]){
								file4 << "\n\n\n";
								for(int i=0; i < N; i++){
									file4 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eRelaxationGroundState.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						
					}
				auto stop = high_resolution_clock::now();
				auto duration = duration_cast<microseconds>(stop - start);
				fileE.close(); 
				file.close(); 
				file0.close(); 
				file1.close(); 
				file2.close(); 
				file3.close(); 
				file4.close(); 
				printfilegnu1(Nit);
				cout<< "The process has been finished and the time it takes is " << duration.count() << '\n';
				return 0;
				}
			
		}
	if(Ptype==2){
		cout << "Which kind of graph do you want?" << '\n';
		cout << "1.- Dynamic to the Coherent states for HO" << '\n';
		cout << "2.- Static with different times" << '\n';
		cout << "Write the number" << '\n';
		cin  >> Graphtype;
		while(Graphtype!=1&&Graphtype!=2){
			cout << "The number is incorrect, try again";
			cout << "1.- Dynamic to the Coherent states for HO" << '\n';
			cout << "2.- Static with different times" << '\n';
			cout << "Write the number" << '\n';
			cin >> Graphtype;
			}
		double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
		// R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, 
		//Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
				
		double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure , 
		//Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
				
		//for time integration leap frog
		double step=4.0e-3; //lenght of step in the time 
				
		//for adaptative smoothing length
		double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
		int hf=490;  //hf-->h-factor.
			
		//for energy
		double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
		// EQn----> energy quantum nature, Mu-------->Chemical Potential 
				
		//for error
		double error=0.0;
		//for effsearch
		double x=4.0;
		int keyS[N]={}, idx[N]={};
		int h_hash=2*h[0];
		int nclass=int((2*x)/h_hash);
		int idxmin[nclass]={}, idxmax[nclass]={};
		int act[nclass]={};
		int xf=1;
				
		//*****Starting the SPH*****
		SPH(search,Ptype,x,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
				
		E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
		if(g!=0.0){
			Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
		}
		//Qenergyi(N, m, R, V, D, Dx, Ei);
				 
		//data
		string dataname1;
		string dataname2;
		if(Graphtype==1){
			dataname1="CoherentHOD.xxx";
			dataname2="eCoherentHOD.xxx";
			}
		if(Graphtype==2){
			dataname1="CoherentHO.xxx";
			dataname2="eCoherentHO.xxx";
			}
		ofstream file(dataname1); //open file to data
		ofstream fileE(dataname2); //open file to data for energy
		//print initial data 
		//file << "\n\n\n"; //print in data file the initial values
		//file1 << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
		for(int i=0; i < N; i++){
			file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
			}
		if(Graphtype==1){	
			int itmax=10000; // # of iterations of evolution Leap Froag
			//start the evolution
			for(int t=0; t<itmax; t++){
				if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
				for(int i=0; i<N; i++){
					R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
					Aold[i]=A[i];
					}
				if(search=="standard"){
					Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
					Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
					Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
					Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
					Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
					AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
					AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
					//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
					//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
					AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
					AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
					//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
				}
				for(int i=0; i<N; i++){
						A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
						if(t*step==28.0){
							V[i]=V[i]+0.5*(Aold[i]+A[i])*step+1.0;
							DV=0.0;
							}else{
								V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
								}
						}

				if(t%50==0){//print the data only for multiples to 50
					file << "\n\n\n";
					for(int i=0; i < N; i++){
						file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
						}
					E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
					if(g!=0.0){
						Mu=ChePotential(N, g, m, R, V,  D, Dx);
						}
					fileE.open("eCoherentHOD.xxx",std::fstream::app);// necesary to write the data in the file
					fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
					fileE.close();
					}
				}
				fileE.close(); // close the datafile for energy
				file.close(); //we close the datafile
				cout<< "The process has been finished"<< '\n';
				return 0;
			}
		if(Graphtype==2){
			int itmax=10000;
				int Ngraphs=5;
				double t0=28.0, t1=28.2, t2=28.4, t3=28.6, t4=28.8, t5=29.0; // times for print the density
				
				//open data for static graph
				ofstream file0("CoherentHO0.xxx"); //open file to data for time 1
				ofstream file1("CoherentHO1.xxx"); //open file to data for time 2
				ofstream file2("CoherentHO2.xxx"); //open file to data for time 3
				ofstream file3("CoherentHO3.xxx"); //open file to data for time 4
				ofstream file4("CoherentHO4.xxx"); //open file to data for time 5
				ofstream file5("CoherentHO5.xxx"); //open file to data for time 6
				
				//start the evolution
				for(int t=0; t<itmax; t++){
					if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
					for(int i=0; i<N; i++){
							R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
							Aold[i]=A[i];
						}
					if(search=="standar"){
						Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
						Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
						Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
						AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}
					for(int i=0; i<N; i++){
						A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
						if(t*step==28.0){
							V[i]=V[i]+0.5*(Aold[i]+A[i])*step+1.0;
							DV=0.0;
							}else{
								V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
								}
						}
						if(t*step==t0){//print the data only for multiples to 50
							file0 << "\n\n\n";
								for(int i=0; i < N; i++){
									file0 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==t1){
								file1 << "\n\n\n";
								for(int i=0; i < N; i++){
									file1 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==t2+0.004){
								file2 << "\n\n\n";
								for(int i=0; i < N; i++){
									file2 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==t3){
								file3 << "\n\n\n";
								for(int i=0; i < N; i++){
									file3 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==t4){
								file4 << "\n\n\n";
								for(int i=0; i < N; i++){
									file4 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
						if(t*step==t5){
								file5 << "\n\n\n";
								for(int i=0; i < N; i++){
									file5 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
									}
								E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
								if(g!=0.0){
									Mu=ChePotential(N, g, m, R, V,  D, Dx);
									}
								fileE.open("eCoherentHO.xxx",std::fstream::app);// necesary to write the data in the file
								fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
								fileE.close();
							}
					}
				fileE.close(); 
				file.close(); 
				file0.close(); 
				file1.close(); 
				file2.close(); 
				file3.close(); 
				file4.close(); 
				file5.close();
				printfilegnu2(t0, t1, t2, t3, t4, t5);
				cout<< "The process has been finished"<< '\n';
				return 0;
				}
			}	
		if(Ptype==3){
			double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
			// R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, 
			//Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
				
			double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure , 
			//Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
				
			//for time integration leap frog
			double step=4.0e-3; //lenght of step in the time 
				
			//for adaptative smoothing length
			double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
			int hf=490;  //hf-->h-factor.
			
			//for energy
			double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
			// EQn----> energy quantum nature, Mu-------->Chemical Potential 
				
			//for error
			double error=0.0;
			//for effsearch
			double x=4.0;
			int keyS[N]={}, idx[N]={};
			int h_hash=2*h[0];
			int nclass=int((2*x)/h_hash);
			int idxmin[nclass]={}, idxmax[nclass]={};
			int act[nclass]={};
			int xf=1;
				
			//*****Starting the SPH*****
			SPH(search,Ptype,x,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
				
			E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
			if(g!=0.0){
				Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
			}
			//Qenergyi(N, m, R, V, D, Dx, Ei);
				 
			//data
			string dataname1="BECgroundstate"+to_string(g)+".xxx";
			string dataname2="eBECgroundstate"+to_string(g)+".xxx";
			ofstream file(dataname1); //open file to data
			ofstream fileE(dataname2); //open file to data for energy
			//print initial data 
			//file << "\n\n\n"; //print in data file the initial values
			//fileE << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
			for(int i=0; i < N; i++){
				file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
				}
			int itmax=20000; // # of iterations of evolution Leap Froag
			print_file_gnu_BEC_dynamic(itmax, g);
			//start the evolution
				for(int t=0; t<itmax; t++){
					if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
					for(int i=0; i<N; i++){
							R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
							Aold[i]=A[i];
						}
					if(search=="standar"){
						Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
						Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
						Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
						Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
						Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
						AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
						AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
						AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
						AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					}
					for(int i=0; i<N; i++){
							A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
							V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}

					if(t%50==0){//print the data only for multiples to 50
						file << "\n\n\n";
							for(int i=0; i < N; i++){
								file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
								}
							E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
							if(g!=0.0){
								Mu=ChePotential(N, g, m, R, V,  D, Dx);
								}
							fileE.open("eBECgroundstate"+to_string(g)+".xxx",std::fstream::app);// necesary to write the data in the file
							fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
							fileE.close();
						}
					}
				fileE.close(); // close the datafile for energy
				file.close(); //we close the datafile
				cout<< "The process has been finished"<< '\n';
				return 0;
			}
		if(Ptype==4){
			double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
			// R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, 
			//Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
				
			double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure , 
			//Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
				
			//for time integration leap frog
			double step=4.0e-3; //lenght of step in the time 
				
			//for adaptative smoothing length
			double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
			int hf=640;  //hf-->h-factor.
			
			//for energy
			double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
			// EQn----> energy quantum nature, Mu-------->Chemical Potential 
				
			//for error
			double error=0.0;

			//for effsearch
			double x=4.0;
			int keyS[N]={}, idx[N]={};
			int h_hash=2*h[0];
			int nclass=int((2*x)/h_hash);
			int idxmin[nclass]={}, idxmax[nclass]={};
			int act[nclass]={};
			int xf=1;
				
			//*****Starting the SPH*****
			SPH(search,Ptype,x,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
				
			E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
			if(g!=0.0){
				Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
			}
			//Qenergyi(N, m, R, V, D, Dx, Ei);
				 
			//****times
			double t1=10.0, t2=14.0, t3=15.0, t4=18.0, t5=22.0;
			
			//data
			string dataname1="CollisionSoliton.xxx";
			string dataname2="eCollisionSoliton.xxx";
			string dataname00="CollisionSoliton"+to_string(t1)+".xxx";
			string dataname01="CollisionSoliton"+to_string(t2)+".xxx";
			string dataname02="CollisionSoliton"+to_string(t3)+".xxx";
			string dataname03="CollisionSoliton"+to_string(t4)+".xxx";
			string dataname04="CollisionSoliton"+to_string(t5)+".xxx";
			ofstream file(dataname1); //open file to data
			ofstream fileE(dataname2); //open file to data for energy
			ofstream file0(dataname00);
			ofstream file1(dataname01);
			ofstream file2(dataname02);
			ofstream file3(dataname03);
			ofstream file4(dataname04);
			//print initial data 
			//file << "\n\n\n"; //print in data file the initial values
			//fileE << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
			for(int i=0; i < N; i++){
				file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
				}
				int itmax=5000; // # of iterations of evolution Leap Froag
				//start the evolution
				for(int t=0; t<itmax; t++){
					if(t%100==0){cout << step*t << '\n';}//print only the values that have multiples of 100
					for(int i=0; i<N; i++){
							R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
							Aold[i]=A[i];
						}
					Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
					Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
					Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
					Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
					Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
					AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure 
					AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
					//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive  
					//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
					AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
					AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping 
					//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
					
					for(int i=0; i<N; i++){
						A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
						if(R[i]<0.0){
							if(t*step==10.0){
								V[i]=V[i]+0.5*(Aold[i]+A[i])*step+1.0;
								DV=0.0;
								}else{
									V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
									}
								}
						if(R[i]>0.0){
							if(t*step==10.0){
								V[i]=V[i]+0.5*(Aold[i]+A[i])*step-1.0;
								DV=0.0;
								}else{
									V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
									}
							}
						}

					if(t*step==10.0){
						file0 << "\n\n\n"; //print in data file the initial values
						for(int i=0; i < N; i++){
							file0 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
							}
						E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
						if(g!=0.0){
							Mu=ChePotential(N, g, m, R, V,  D, Dx);
							}
						fileE.open("eCollisionSoliton.xxx",std::fstream::app);
						fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
						fileE.close();
						}
					if(t*step==14.0){
						file1 << "\n\n\n"; //print in data file the initial values
						for(int i=0; i < N; i++){
							file1 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
							}
						E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
						if(g!=0.0){
							Mu=ChePotential(N, g, m, R, V,  D, Dx);
							}
						fileE.open("eCollisionSoliton.xxx",std::fstream::app);
						fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
						fileE.close();
						}
					if(t*step==15.0){
						file2 << "\n\n\n"; //print in data file the initial values
						for(int i=0; i < N; i++){
							file2 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
							}
						E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
						if(g!=0.0){
							Mu=ChePotential(N, g, m, R, V,  D, Dx);
							}
						fileE.open("eCollisionSoliton.xxx",std::fstream::app);
						fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
						fileE.close();
						}
					if(t*step==18.0){
						file3 << "\n\n\n"; //print in data file the initial values
						for(int i=0; i < N; i++){
							file3 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
							}
						E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
						if(g!=0.0){
							Mu=ChePotential(N, g, m, R, V,  D, Dx);
							}
						fileE.open("eCollisionSoliton.xxx",std::fstream::app);
						fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
						fileE.close();
						}
					if(t*step==22.0){
						file4 << "\n\n\n"; //print in data file the initial values
						for(int i=0; i < N; i++){
							file4 << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
							}
						E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
						if(g!=0.0){
							Mu=ChePotential(N, g, m, R, V,  D, Dx);
							}
						fileE.open("eCollisionSoliton.xxx",std::fstream::app);
						fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
						fileE.close();
						}
				}
				printfilegnu3(t1,t2,t3,t4,t5);
				fileE.close(); // close the datafile for energy
				file.close(); //we close the datafile
				cout<< "The process has been finished"<< '\n';
				return 0;
			}
	/*
	if(Dim=="1Dimensional"){
//		int N=300; // # of particles -----> for BEC and Harmonic Oscillator
//			N=704;
//			N=pow(2,5)*20;	// # of particles   -----> for Harmonic Oscillator
		int itmax=12000; // # of iterations of evolution Leap Froag
	   // nonlinear parameter for NLSE ----> BEC
	  //Initial values
	  double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N]; // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, Pxx-> component xx of press tensor, A-> Total Acceleration. 
	  //Different acelerations Aold --> old acceleration
	  double Aq[N],Agp[N],Av[N],Ad[N], DV=4.0; // Aq -> acceleration due  for Quantum potential or quantum pressure , Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
	  //for time integration leap frog
	//	double step=0.1;
	  double step=4.0e-3;
	//  double Vp[N],Vf[N],Vc[N],Xp[N],Xf[N],Xc[N],step, DV=4.0, tol=1.0e-16, alc; // Vp-> velocity t-1/2, Vf-> velocity t+1/2, Vc-> center velocity, Xp-> position t-1/2, Xf-> position t+1/2, Xc-> central postion, step->time step ,DV->Damping parameter.
	  //for adaptative smoothing length
	  double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
	  //for energy
	  double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average 
	  //for error
	  double error=0.0;
	  SPH(Ptype,N,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual partícles with the initial values of  g,R,m,h,V,D,Dx,Dxx,Pxx,A, Zh, Omega
	  //In this case -h because we are goint to a back step.
	  E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
	  if(g!=0.0){
		Mu=ChePotential(N, g, m, R, V,  D, Dx);
	  } 
	//  Qenergyi(N, m, R, V, D, Dx, Ei);
	//	data
	  ofstream file("Solitons.xxx"); //open file to data
		ofstream file1("eSolitons.xxx");
//		file << "\n\n\n"; //print in data file the initial values
//		file1 << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
		  for(int i=0; i < N; i++){
			file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
		  }
		//start the evolution
		for(int t=0; t<itmax; t++){
			if(t%100==0){cout << step*t << '\n';}
			for(int i=0; i<N; i++){
				R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
				Aold[i]=A[i];
			}
		
			Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
			Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
			Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density
			Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
			AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure 
			AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
	//		AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive  
	//    	AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
			AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
			AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping 
	//		Qenergyi(N, m, R, V, D, Dx, Ei);
// this is for a HO evolution o soliton collition
			///////////////////////////
			
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
					if(t*step==28.1){
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step+1.0;
						DV=0.1;
					}else{
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}
			}
			
			////////////////////////////////This is for BEC
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
					V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						
			}
			
			////////////////////////////////
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
				if(R[i]<0.0){
					if(t*step==10.0){
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step+1.0;
						DV=0.0;
					}else{
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}
				}
				if(R[i]>0.0){
					if(t*step==10.0){
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step-1.0;
						DV=0.0;
					}else{
						V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
						}
				}
					
			}
		
			 //////////////////////////////////
			///////////////////////////////////////////////////////////////////
			if(t*step==28.0){//this part correspond to the convergence of the method
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
				fstream file2;
				file2.open("Error.xxx",std::fstream::app); 
					if (!file2.is_open()) {
						cout << "failed to open " << '\n';
					} else {
						double w=0.0,z=0.0;
						for(int i=0; i<N; i++){
								z=D[i]-DA(2, R[i]);
								w+=z*z;
							}	
							error=sqrt(w);
							file2 << N << "\t" << error << '\n';
							file2.close();
					}					
			}
			////////////////////////////////////////////////////////////////////
				////////////////////////////// for dynamic to the coherent states
			if(t*step==28.0 || t*step==28.2 || t*step==28.4 || t*step==28.6 || t*step==28.8){ // we use conditional only for the dynamics of HO , in other case we use t%50==0
				file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
				E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
				if(g!=0.0){
					Mu=ChePotential(N, g, m, R, V,  D, Dx);
				}
				file1.open("eDynamics01.xxx",std::fstream::app);
				file1 << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
				file1.close();
			}
				/////////////////////////////////
			/////////////////
			
			if(t%50==0){
				  file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
				E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
				if(g!=0.0){
					Mu=ChePotential(N, g, m, R, V,  D, Dx);
				}
				file1.open("eSolitons2.xxx",std::fstream::app);
				file1 << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
				file1.close();
			}
			
			////////////////This part use to do a collision to sollitons
			
			if(t*step==10.0 || t*step==12.0 || t*step==14.0 || t*step==15.0 || t*step==16.0 || t*step==18.0 || t*step==22.0){
				  file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
				E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
				if(g!=0.0){
					Mu=ChePotential(N, g, m, R, V,  D, Dx);
				}
				file1.open("eSolitons.xxx",std::fstream::app);
				file1 << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
				file1.close();
			}
			
		}
	//	file1.close(); // close the datafile for energy
		file.close(); //we close the datafile
			cout<< "El proceso ha terminado"<< '\n';
				return 0;	
		
		}
///// 2Dimensional

	if(Dim="2Dimensional"){
				//	int N=100; // # of particles -----> for BEC and Harmonic Oscillator
		int N=pow(2,5)*20;	// # of particles   -----> for Harmonic Oscillator
		int itmax=10000; // # of iterations of evolution Leap Froag
	  double g = 10.0; // nonlinear parameter for NLSE ----> BEC
	  //Initial values
	  double R[N], X[N], Y[N], m[N], V[N], VX[N], VY[N], h[N], D[N], Dx[N],Dy[N], Dxx[N],Dyy[N], Dyx[N], Dxxx[N], Dyyy[N], Dxyy[N],Dxxy[N], Pxx[N], A[N],Ax[N],Ay[N], Aold[N],Aoldx[N],Aoldy[N]; // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density, Pxx-> component xx of press tensor, A-> Total Acceleration. 
	  //Different acelerations Aold --> old acceleration
	  double Aq[N],Agp[N],Av[N],Ad[N], DV=4.0; // Aq -> acceleration due  for Quantum potential or quantum pressure , Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
	  //for time integration leap frog
	//	double step=0.1;
	  double step=4.0e-3;
	//  double Vp[N],Vf[N],Vc[N],Xp[N],Xf[N],Xc[N],step, DV=4.0, tol=1.0e-16, alc; // Vp-> velocity t-1/2, Vf-> velocity t+1/2, Vc-> center velocity, Xp-> position t-1/2, Xf-> position t+1/2, Xc-> central postion, step->time step ,DV->Damping parameter.
	  //for adaptative smoothing length
	  double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function. 
	  //for energy
	  double E, Eff; // E----> energy average 
	  SPH(N,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual partícles with the initial values of  g,R,m,h,V,D,Dx,Dxx,Pxx,A, Zh, Omega
	  //In this case -h because we are goint to a back step.
	  E=Qenergy(N, m, R, V, D, Dx); //Initial Energy 
	//  Qenergyi(N, m, R, V, D, Dx, Ei);
			
		//data
	  ofstream file("data.xxx"); //open file to data
		ofstream file1("energydata.xxx");
		file << "\n\n\n"; //print in data file the initial values
		file1 << 0 << "\t\t" << E << "\t\t" << Eff << '\n';
		  for(int i=0; i < N; i++){
			file << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] <<  "\t\t"<< Ei[i] <<"\n";
		  }
		//start the evolution
		for(int t=0; t<itmax; t++){
			if(t%100==0){cout << step*t << '\n';}
			for(int i=0; i<N; i++){
					R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
					Aold[i]=A[i];
			}
			////part for 
			Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
			Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
			Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density
			Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
			AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure 
			AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
	//		AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive  
	//    	AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
			AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
			AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping 
	//		Qenergyi(N, m, R, V, D, Dx, Ei);
			for(int i=0; i<N; i++){
					A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
					V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
			}

			if(t%10==0){ // 
				file << "\n\n\n"; //print in data file the initial values
				for(int i=0; i < N; i++){
					file << R[i] << "\t\t" << D[i] << "\t\t"<< Dx[i] << "\t\t"<< Dxx[i] << "\t\t"<< Dxxx[i] << "\t\t"  << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] << "\t\t" << Zh[i] << "\t\t" << dZh[i] << "\t\t"<< Omega[i] << "\t\t"<< Ei[i] <<'\n';
				}
			}
			E=Qenergy(N, m, R, V, D, Dx);
			file1 << t*step  << "\t\t" << E  << "\t\t" << Eff << '\n';
		}
			file1.close(); // close the datafile for energy
		file.close(); //we close the datafile
			cout<< "El proceso ha terminado"<< '\n';
				return 0;
		}
		*/ 
}
