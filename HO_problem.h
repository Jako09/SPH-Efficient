#include <iostream>
#include"SPH.h"
#include"EulerInt.h"
#include"LeapFrog.h"
#include"gnu.h"
#include"L2.h"
#include <chrono>
#include <iomanip>

//#include<chrono>

using namespace std;

void RelaxHO(int N, int search);
void L2_time_HO(int N, int search);
void InitialDataHO(int N);



void HarmonicOscillator(int N,int search){
    int graphtype;
    cout << "Which graph you need? " << '\n';
    cout << "1.-Relaxation" << '\n';
    cout << "2.-Norm L2 And Time operation" << '\n';
    cout << "3.-Initial data" << '\n';
    cin >> graphtype;

    if(graphtype==1){
        cout << "You choose the Relaxation, you'll obtain 4 type of graphs (Dynamic, Static, EnergyvsTime, Time of search per iteration)"<< '\n';
        RelaxHO(N,search);
    }

    if(graphtype==2){
        cout << "You'll obtain two data sets, one for the L2 norm, to print the graph, you need to use a respective app"<< '\n';
        L2_time_HO(N,search);
    }

    if(graphtype==3){
        cout << "you choose the 3"<< '\n';
        InitialDataHO(N);
    }


}

void RelaxHO(int N, int search){
    using std::chrono::time_point;
    using std::chrono::system_clock;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds;

    double g=0.0, DV=8.0; // nonlinear parameter.
    //---------------
    double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
    // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density,
    // Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
    //---------------
    double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure ,
    //Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
    //---------------
    double step=4.0e-3; //lenght of step in the time
	//---------------
    //for adaptative smoothing length and h factor
    double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function.
    int hf=1605;  //hf-->h-factor. 490 es el estadar y funciona con N=640 // for 1280 particles we use hf=790
    //---------------
    //for energy
    double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
    // EQn----> energy quantum nature, Mu-------->Chemical Potential
	//---------------
    //for error
    double error=0.0;
    double xmin, xmax;
    //---------------
    //*****Starting the SPH*****
    // The Ptype is HO, it means ew use the number 2
    int Ptype=2, Dtype=2;
    // The distrbution type
    SPH(search,Ptype,4.0,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
    int keyS[N]={}, idx[N]={};
    double h_hash=(1.0/sqrt(2.0))*h[0];
    int nclass=int((8.0)/h_hash); //Initial number of classes
    int xf=5;
    E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
    if(g!=0.0){
        Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
    }
    //data
    string dataname1;
    string dataname2;
    string datanamet[5];
    string dataname3;
    if(search==1){
        dataname1="RelaxHO_D_S.xxx";
        dataname2="E_Relax_HO_S.xxx";
        datanamet[0]="Relax_HO_S_0.xxx";
        datanamet[1]="Relax_HO_S_1.xxx";
        datanamet[2]="Relax_HO_S_2.xxx";
        datanamet[3]="Relax_HO_S_3.xxx";
        datanamet[4]="Relax_HO_S_4.xxx";
        dataname3="Relax_HO_it_S.xxx";
    }
    if(search==2){
        dataname1="RelaxHO_D_E.xxx";
        dataname2="E_Relax_HO_E.xxx";
        datanamet[0]="Relax_HO_E_0.xxx";
        datanamet[1]="Relax_HO_E_1.xxx";
        datanamet[2]="Relax_HO_E_2.xxx";
        datanamet[3]="Relax_HO_E_3.xxx";
        datanamet[4]="Relax_HO_E_4.xxx";
        dataname3="Relax_HO_it_E.xxx";
    }

    ofstream fileD(dataname1); //open file to data
    ofstream fileE(dataname2); //open file to data for energy
    ofstream filet0(datanamet[0]);
    ofstream filet1(datanamet[1]);
    ofstream filet2(datanamet[2]);
    ofstream filet3(datanamet[3]);
    ofstream filet4(datanamet[4]);
    ofstream fileit(dataname3);
    //print initial data
    //file << "\n\n\n"; //print in data file the initial values
    fileE << 0 << "\t\t" << E << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
    for(int i=0; i < N; i++){
        fileD << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
        filet0 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
    }

    int itmax=7500; // # of iterations of evolution Leap Froag
    //start the evolution
    std::chrono::time_point<std::chrono::system_clock> startC, stopC, startL, startD0, stopD0, startD1, stopD1, startD2, stopD2, startAQ, stopAQ, startAV, stopAV, startAD, stopAD, stopL;
/*    double durationD0;
    double durationD1;
    double durationD2;
    double durationAQ;
    double durationAV;
    double durationAD;
    //---To compare
    double durationL0 = double(stopAD - startD0);
    double durationL = double(stopL - startL);
*/
    int Ngraphs=5;
    double Nit[Ngraphs];
    for(int i=0; i<Ngraphs; i++){
        Nit[i]=4.0*(double)i;
    }
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
        if(search==1){
            startL=std::chrono::system_clock::now();
            startD0=std::chrono::system_clock::now();
            Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
            stopD0=std::chrono::system_clock::now();
            startD1=std::chrono::system_clock::now();
            Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
            stopD1=std::chrono::system_clock::now();
            startD2=std::chrono::system_clock::now();
            Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
            stopD2=std::chrono::system_clock::now();
//            Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
            Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
            startAQ=std::chrono::system_clock::now();
            AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
            stopAQ=std::chrono::system_clock::now();
            AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
						//AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
						//AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
            startAV=std::chrono::system_clock::now();
            AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
            stopAV=std::chrono::system_clock::now();
            startAD=std::chrono::system_clock::now();
            AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
						//Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
            stopAD=std::chrono::system_clock::now();
            stopL=std::chrono::system_clock::now();
        }
        if(search==2){
            startL=std::chrono::system_clock::now();
//			printf("pass if %d \n", t);
            startC=std::chrono::system_clock::now();
            CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
            stopC=std::chrono::system_clock::now();
//			printf("pass censo %d \n", t);
            startD0=std::chrono::system_clock::now();
            Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
            stopD0=std::chrono::system_clock::now();
            startD1=std::chrono::system_clock::now();
//			printf("pass densidadEff %d \n", t);
            Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
            stopD1=std::chrono::system_clock::now();
            startD2=std::chrono::system_clock::now();
            Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
            stopD2=std::chrono::system_clock::now();
//			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
            Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
            startAQ=std::chrono::system_clock::now();
            AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
            stopAQ=std::chrono::system_clock::now();
            AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
            //AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
            //AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
            startAV=std::chrono::system_clock::now();
            AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
            stopAV=std::chrono::system_clock::now();
            startAD=std::chrono::system_clock::now();
            AceDampEff(N, nclass, xf, m, R, h, D, DV, V, Ad, keyS, idx, idxmin, idxmax, act); // The AceDamp function gives values to Acceleration due for damping
            stopAD=std::chrono::system_clock::now();
            //Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
            stopL=std::chrono::system_clock::now();
        }

        std::chrono::duration<double> durationC = stopC - startC;
        std::chrono::duration<double> durationD0 = stopD0 - startD0;
        std::chrono::duration<double> durationD1 = stopD1 - startD1;
        std::chrono::duration<double> durationD2 = stopD2 - startD2;
        std::chrono::duration<double> durationAQ = stopAQ - startAQ;
        std::chrono::duration<double> durationAV = stopAV - startAV;
        std::chrono::duration<double> durationAD = stopAD - startAD;
        //---To compare
        std::chrono::duration<double> durationL0 = stopAD - startD0;
        std::chrono::duration<double> durationL = stopL - startL;


        for(int i=0; i<N; i++){
            A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
            V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
        }


        if(t%50==0){//print the data only for multiples to 50
            //For the field values-------
            fileD << "\n\n\n";
            for(int i=0; i < N; i++){
                fileD << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
            //---------------------------
            //For the energy-------------
            E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl);
            if(g!=0.0){
                Mu=ChePotential(N, g, m, R, V,  D, Dx);
            }
            fileE.open(dataname2,std::fstream::app);// necesary to write the data in the file
            fileE << t*step  << "\t\t" << E   << "\t\t" << Enl <<"\t\t" << EKin <<"\t\t" << EPot <<"\t\t" << EQn <<"\t\t" << Mu << '\n';
            fileE.close();
            //---------------------------
            //For the time iterations----
            if(search==1){
                fileit  << t << "\t\t"  << durationD0.count() << "\t\t" << durationD1.count() << "\t\t" << durationD2.count() << "\t\t" << durationAQ.count() << "\t\t" << durationAV.count() << "\t\t" << durationAD.count() << "\t\t" << durationL.count() << "\t\t" << durationL0.count() << '\n';
            }
            if(search==2){
                fileit << t << "\t\t" << durationC.count() << "\t\t" << durationD0.count() << "\t\t" << durationD1.count() << "\t\t" << durationD2.count() << "\t\t" << durationAQ.count() << "\t\t" << durationAV.count() << "\t\t" << durationAD.count() << "\t\t" << durationL.count() << "\t\t" << durationL0.count() <<  '\n';
            }
            //---------------------------
        }

        if(t*step==Nit[0]){//print the data only for multiples to 50
            filet0 << "\n\n\n";
            for(int i=0; i < N; i++){
                filet0 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }

        if(t*step==Nit[1]){
            filet1 << "\n\n\n";
            for(int i=0; i < N; i++){
                filet1 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }

        if(t*step==Nit[2]){
            filet2 << "\n\n\n";
            for(int i=0; i < N; i++){
                filet2 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }

        if(t*step==Nit[3]){
            filet3 << "\n\n\n";
            for(int i=0; i < N; i++){
                filet3 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }

		if(t*step==Nit[4]){
            filet4 << "\n\n\n";
            for(int i=0; i < N; i++){
                filet4 << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }
    }
    printGNU_Energy_HO(search, dataname2, "EPS");
    printGNU_Energy_HO(search, dataname2, "PNG");
    printGNU_iteration_HO(search, dataname3, "EPS");
    printGNU_iteration_HO(search, dataname3, "PNG");
    printGNU_Static_HO(search, datanamet,"EPS",Nit);
    printGNU_Static_HO(search, datanamet,"PNG",Nit);

    filet0.close();
    filet1.close();
    filet2.close();
    filet3.close();
    filet4.close();
    fileit.close();
    fileE.close(); // close the datafile for energy
    fileD.close(); //we close the datafile
    cout<< "The process has been finished " << '\n';
}

void L2_time_HO(int N,int search){
        string nameDynamic, nameConvergence, nameTime, nameL2F, nameL2NonF, pathtime, pathL2, pathConvergence, pathDynamic;
        pathtime="DataTime/";
        pathL2="DataL2/";
        pathConvergence="DataConvergence/";
        pathDynamic="DataDynamicRelax/";
        if(search==1){
            nameDynamic="dynamic_"+to_string(N)+"_std.dat";
            nameConvergence="convergence_"+to_string(N)+"_std.dat";
            nameTime="time_operation_std_HO.dat";
            nameL2F="L2_norm_std_HO.dat";
            nameL2NonF="L2nF_norm_std_HO.dat";
        }
        if(search==2){
            nameDynamic="dynamic_"+to_string(N)+"_eff.dat";
            nameConvergence="convergence_"+to_string(N)+"_eff.dat";
            nameTime="time_operation_eff_HO.dat";
            nameL2F="L2_norm_eff_HO.dat";
            nameL2NonF="L2nF_norm_eff_HO.dat";
        }
        ofstream fileConvergence(pathConvergence+nameConvergence);
        ofstream fileDynamic(pathDynamic+nameDynamic);
        ofstream fileTime, fileL2F, fileL2NonF;
        double g=0.0, DV=8.0; // nonlinear parameter.
        //---------------
        double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
        // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density,
        // Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
        //---------------
        double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure ,
        //Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
        //---------------
        double step=4.0e-3; //lenght of step in the time
        //---------------
        //for adaptative smoothing length and h factor
        double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function.
        int hf=490;  //hf-->h-factor. 490 es el estadar y funciona con N=640 // for 1280 particles we use hf=790
        //---------------
        //for energy
        double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
        // EQn----> energy quantum nature, Mu-------->Chemical Potential
        //---------------
        //for error
        double error=0.0;
        double xmin, xmax;
        //---------------
        //*****Starting the SPH*****
        // The Ptype is HO, it means ew use the number 2
        int Ptype=2, Dtype=2;
        // The distrbution type
        SPH(search,Ptype,4.0,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
        fileDynamic << "\n\n\n";
        for(int i=0; i<N; i++){
            fileDynamic << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
        }
        int keyS[N]={}, idx[N]={};
        double h_hash=(1.0/sqrt(2.0))*h[0];
        cout << h_hash << '\n';
        int nclass=int((8.0)/h_hash); //Initial number of classes
        int xf=5;
        E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
        if(g!=0.0){
            Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
        }
        //data
        //print initial data
        //file << "\n\n\n"; //print in data file the initial values

        int itmax=7500; // # of iterations of evolution Leap Froag
        double  time_process=0.0;
        //start the evolution
        std::chrono::time_point<std::chrono::system_clock> startL, stopL;
        std::chrono::duration<double> durationL;
        for(int t=0; t<itmax; t++){
            if(t%100==0){cout << step*t << " " << N  << '\n';}//print only the values that have multiples of 100
            for(int i=0; i<N; i++){
//                cout << "BR[i] " << R[i] << " t " << t << '\n';
                R[i]=R[i]+V[i]*step+0.5*A[i]*(step*step);
                Aold[i]=A[i];
//                cout << "V[i] " << V[i] << " A[i] " << A[i] << '\n';
            }
            xmin=getMind(N,R)-h[0];
//            cout << "pass getMind " << xmin << '\n';
            xmax=getMaxd(N,R)+h[0];
//            cout << "pass getMaxd " << xmax << '\n';
            nclass=int((xmax-xmin)/h_hash);
//            cout << "pass nclass " << nclass<< '\n';
            int idxmin[nclass]={}, idxmax[nclass]={};
            int act[nclass]={};
//            cout << " xmin "<< xmin <<" xmax "<<xmax<<" nclass "<<nclass<< '\n';
            if(search==1){
                startL=std::chrono::system_clock::now();
                Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
                Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
                Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
    //            Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
                Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
                AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
                AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
                            //AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
                            //AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
                AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
                AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
                            //Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
                stopL=std::chrono::system_clock::now();
            }
            if(search==2){
                startL=std::chrono::system_clock::now();
    //			printf("pass if %d \n", t);
                CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
    //			printf("pass censo %d \n", t);
                Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
    //			printf("pass densidadEff %d \n", t);
                Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
                Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
    //			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
                Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
                AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
                AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
                //AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
                //AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
                AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
                AceDampEff(N, nclass, xf, m, R, h, D, DV, V, Ad, keyS, idx, idxmin, idxmax, act); // The AceDamp function gives values to Acceleration due for damping
                //Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
                stopL=std::chrono::system_clock::now();
            }
            durationL = stopL - startL;
            time_process+=(double)durationL.count();
            for(int i=0; i<N; i++){
                A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
                V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
            }
            fileDynamic << "\n\n\n";
            for(int i=0; i<N; i++){
                fileDynamic << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            }
        }
        fileDynamic.close();
        fileTime.open(pathtime+nameTime,std::fstream::app);
        fileTime <<  N << "\t\t " << time_process << '\n';
        fileTime.close();
        fileL2F.open(pathL2+nameL2F,std::fstream::app);
        fileL2F << N << "\t\t " << L2(2,1, N, R, D) << '\n';
        fileL2F.close();
        fileL2NonF.open(pathL2+nameL2NonF,std::fstream::app);
        fileL2NonF << N << "\t\t " << L2(2,2, N, R, D) << '\n';
        fileL2NonF.close();
        for(int i=0; i<N; i++){
            fileConvergence.open(pathConvergence+nameConvergence,std::fstream::app);
            fileConvergence << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
            fileConvergence.close();
        }
}

void L2_ti(int search){
    using std::chrono::time_point;
    using std::chrono::system_clock;
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    int N, countmax=6;
//    int hf[countmax];
    N=20;
    string nameF[countmax];
    string nameL2F, nameL2NonF, nameT;
    if(search==1){
        nameL2F="HO_Norm_L2_Std_f.xxx";
        nameL2NonF="HO_Norm_L2_Std_nf.xxx";
        nameT="HO_Time_Std.xxx";
    }
    if(search==2){
        nameL2F="HO_Norm_L2_Eff_f.xxx";
        nameL2NonF="HO_Norm_L2_Eff_nf.xxx";
        nameT="HO_Time_Eff.xxx";
    }
    ofstream fileL2F(nameL2F);
    ofstream fileL2NonF(nameL2NonF);
    ofstream fileT(nameT);
    for(int count=0; count<=countmax;count++){
        if(search==1){
            nameF[count]="HO_Field_Std_"+to_string(count)+".xxx";
        }
        if(search==2){
            nameF[count]="HO_Field_Eff_"+to_string(count)+".xxx";
        }
        ofstream filefield;
        double g=0.0, DV=8.0; // nonlinear parameter.
        //---------------
        double R[N], m[N], V[N], h[N], D[N], Dx[N], Dxx[N], Dxxx[N], Pxx[N], A[N], Aold[N], Ei[N];
        // R -> position , m -> mass, V -> velocity, h -> smoothing length, D-> density, Dx-> derivative density, Dxx-> second derivative density,
        // Pxx-> component xx of press tensor, A-> total Acceleration, Aold --> old acceleration.
        //---------------
        double Aq[N],Agp[N],Av[N],Ad[N]; // Aq -> acceleration due  for Quantum potential or quantum pressure ,
        //Agp-> aceleration due for Nonlinear or Gross-Pitaevskii , Av-> Acceleration due to potential, Ad-> damping aceleration.
        //---------------
        double step=4.0e-3; //lenght of step in the time
        //---------------
        //for adaptative smoothing length and h factor
        double Zh[N], Omega[N], dZh[N]; // Zh-> Zeta function for h, Omega-> h-grad, dZh-> derivative zeta function.
        int hf=790;  //hf-->h-factor. 490 es el estadar y funciona con N=640 // for 1280 particles we use hf=790
        //---------------
        //for energy
        double E, Enl, EKin, EPot, EQn, Mu=0.0; // E----> energy average, Enl----> energy nonlinear, Ekin-------> energy kinetics, Epot----->energy potential,
        // EQn----> energy quantum nature, Mu-------->Chemical Potential
        //---------------
        //for error
        double error=0.0;
        double xmin, xmax;
        //---------------
        //*****Starting the SPH*****
        // The Ptype is HO, it means ew use the number 2
        int Ptype=2, Dtype=2;
        // The distrbution type
        SPH(search,Ptype,4.0,Dtype,N,hf,g,R,m,h,V,D,Dx,Dxx, Dxxx, Pxx,Aq,Agp, Av,A, Zh, Omega); // SPH function generate N virtual particles with the initial values.
        int keyS[N]={}, idx[N]={};
        double h_hash=(1.0/sqrt(2.0))*h[0];
        int nclass=int((8.0)/h_hash); //Initial number of classes
        int xf=5;
        E=Qenergy(N, g, m, R, V, D, Dx, EKin, EPot, EQn, Enl); //Initial Energy
        if(g!=0.0){
            Mu=ChePotential(N, g, m, R, V,  D, Dx);//ChemicalPotential
        }
        //data
        //print initial data
        //file << "\n\n\n"; //print in data file the initial values

        int itmax=7500; // # of iterations of evolution Leap Froag
        double  time_process=0.0;
        //start the evolution
        std::chrono::time_point<std::chrono::system_clock> startL, stopL;
        std::chrono::duration<double> durationL;
        for(int t=0; t<itmax; t++){
//            if(t%100==0){cout << step*t << " " << N << " " << count << '\n';}//print only the values that have multiples of 100
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
            if(search==1){
                startL=std::chrono::system_clock::now();
                Densidad0(N, m, R, h, D); //The Densidad0 function gives values to Density
                Densidad1(N, m, R, h, D, Dx); // The Densidad1 function gives values to derivative of Density
                Densidad2(N, m, R, h, D, Dx, Dxx); // The Densidad2 function gives values to second derivative of Density
    //            Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
                Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
                AceQ(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx, Aq); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
                AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
                            //AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
                            //AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
                AceV(N, m, h, R, D, Av); // the AceV function gives values to acceleration due for potential term
                AceDamp(N,m, h, R, D,  DV, V, Ad); // The AceDamp function gives values to Acceleration due for damping
                            //Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
                stopL=std::chrono::system_clock::now();
            }
            if(search==2){
                startL=std::chrono::system_clock::now();
    //			printf("pass if %d \n", t);
                CensoSPH(N, h_hash, xmin, R, keyS, idx, idxmin,  idxmax, act);
    //			printf("pass censo %d \n", t);
                Densidad0Eff( N, nclass, xf, m, R, h, D, keyS, idx, idxmin, idxmax, act);
    //			printf("pass densidadEff %d \n", t);
                Densidad1Eff( N, nclass, xf, m, R, h, D, Dx, keyS, idx, idxmin, idxmax, act);// The Densidad1 function gives values to derivative of Density
                Densidad2Eff( N, nclass, xf, m, R, h, D, Dx, Dxx, keyS, idx, idxmin, idxmax, act); // The Densidad2 function gives values to second derivative of Density
    //			Densidad3(N , m, R, h,  D, Dx,Dxx,Dxxx); // The new density, which is necesary only for the quatum force
                Pressxx(N, m, R, h, D, Dx, Dxx, Dxxx, Pxx); // the Pressxx function gives values to xx component of tensor pressure
                AceQEff( N, nclass, xf, m, R, h, D, Dx, Dxx,Pxx,Aq, keyS, idx, idxmin, idxmax, act); // The AceQ function gives the values to acceleration due for quantum potential or quantum pressure
                AceGP(N,g, m, R, h, D, Dx, Agp); // the AceGP function gives the values to nonlinear term  of acceleration with parameter g
                //AceQAdaptative(N, m, Xc, h, D, Pxx,Omega, Aq); // The AceQAdaptative function gives values to acceleration due for quantum pressure special for h adaptive
                //AceGPAdaptative(N,g, m, Xc, h, D,Omega, Agp); // the AceGPAdaptative function gives values to acceleration due for Non linear term g, special for h adaptive
                AceVEff( N, nclass, xf, m, R, h, D, Av, keyS, idx, idxmin, idxmax, act); // the AceV function gives values to acceleration due for potential term
                AceDampEff(N, nclass, xf, m, R, h, D, DV, V, Ad, keyS, idx, idxmin, idxmax, act); // The AceDamp function gives values to Acceleration due for damping
                //Qenergyi(N, m, R, V, D, Dx, Ei);//Energy for each particle i-esima
                stopL=std::chrono::system_clock::now();
            }
            durationL = stopL - startL;
            time_process+=(double)durationL.count();
            for(int i=0; i<N; i++){
                A[i] = Aq[i]+Agp[i]+Av[i]+Ad[i];
                V[i]=V[i]+0.5*(Aold[i]+A[i])*step;
            }
        }
        time_process=time_process/(double)itmax;
        fileT.open(nameT,std::fstream::app);
        filefield.open(nameF[count],std::fstream::app);
        fileL2F.open(nameL2F,std::fstream::app);
        fileL2NonF.open(nameL2NonF,std::fstream::app);

        fileT << "N=" <<  N << "\t\t time=" << time_process << '\n';
        for(int i=0; i<N; i++){
            filefield << R[i] <<  "\t\t"<< D[i] << "\t\t"<< Dx[i] <<  "\t\t"<< Dxx[i]  << "\t\t" << V[i] << "\t\t"<< Aq[i]<< "\t\t" << Agp[i]<<  "\t\t" <<Av[i]<< "\t\t" << Ad[i]<<"\t\t" << A[i] << "\t\t" << Pxx[i] << "\t\t" << h[i] <<  "\t\t"<< Ei[i] <<"\n";
        }
        //-------We need to obtain the iterations here, for example
        //-------Print three kinds of graphs, one that show us the norm L2, in the stable state and second, the time operations
        fileL2F << "N=" << N << "\t\t L2F=" << L2(2,1, N, R, D) << '\n';
        fileL2NonF << "N=" << N << "\t\t L2nF=" << L2(2,2, N, R, D) << '\n';
        N=20*pow(2.0,count+1);
        filefield.close();
        fileL2F.close();
        fileL2NonF.close();
        fileT.close();
    }
    printGNU_time_operation(search,nameT,"EPS");
    printGNU_time_operation(search,nameT,"PNG");
    printGNU_L2(search, nameL2F, nameL2NonF, "EPS");
    printGNU_L2(search, nameL2F, nameL2NonF, "PNG");
    printGNU_field(countmax, search , nameF,"EPS");
    printGNU_field(countmax, search, nameF,"PNG");
    cout<< "The process has been finished " << '\n';
}
void InitialDataHO(int N){
    //In this function we obtain two graphs, one for the initial data for the problem the particle in a Box, second the use of the searching algorithm

}
