using namespace std;

void EulerInt(int N,double s, double A[], double V[], double X[],double Vn[], double Xn[]){
    for(int i =  0; i < N; ++i){
        Xn[i] = X[i] + s*V[i];
        Vn[i] = V[i] + s*A[i];
    }
}