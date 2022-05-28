
using namespace std;

void LeapFrog(int N, double s,double A[],double Xc[],double Vp[], double Xf[], double Vf[]){
    for(int i =  0; i < N; ++i){
        Vf[i] = Vp[i] + s*A[i];  //Vf<-(s/2) Vp<-(-s/2)
    }

    for(int i =  0; i < N; ++i){
        Xf[i] = Xc[i] + s*Vf[i];//Vf <- (s/2) Xf<-(s)
    }
}