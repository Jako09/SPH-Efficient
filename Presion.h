
using namespace std;

void Pressxx(int N, double m[], double h[], double R[], double D[], double Dx[], double Dxx[], double Dxxx[], double Pxx[]){

	double Q=0.0, P=0.0;
	//paper's discretization of pressure tensor 
	/*
  for(int i=0; i<N; ++i){
    Pxx[i]=0.0;
    for(int j=0; j<N; ++j){
      P=0.25*(Dxx[j]-Dx[j]*Dx[j]/D[j]);
      Pxx[i]=Pxx[i]+(m[j]/D[j])*ker(h[i], R[i], R[j])*P;
    }
  }
// 
*/

//-------------->007, 008, 009
	for(int i=0;i<N; i++){
		Pxx[i]=0.0;
		Pxx[i]=0.25*(Dxx[i]-Dx[i]*Dx[i]/D[i]);
		}
//----------->007, 008, 009

////////////////////////////Bohm Representation
/*	//second discretization for quantum potential for Bohm representation 
	for(int i=0; i<N; ++i){
		Pxx[i]=0.0;
		for(int j=0; j<N; ++j){
			Q=0.25*(Dxx[j]-0.5*Dx[j]*Dx[j]/D[j])/D[j];
			Pxx[i]+=(m[j]/D[j])*ker(h[i], R[i], R[j])*Q;
		}
	}
*/
/*
//------------>004, 005, 006
for(int i=0; i<N; i++){
	Pxx[i]=0.0;
 Pxx[i]=-(0.25/D[i])*(Dxx[i]-0.5*Dx[i]*Dx[i]/D[i]);
}
*/ 
//------------->004, 005, 006
//Verify pressure
/*
//------->Resultado 001
 for(int i=0; i<N; i++){
	 Pxx[i]=0.5*D[i];
	 }
//------->Resultado 001
*/
}
/*
void Pressxx(int N,int nclass, int xf, double m[], double R[], double h[], double D[], double Dx[],double Dxx[],double Pxx[], int keyS[], int idx[], int idxmin[], int idxmax[], int act[]){

}
*/
