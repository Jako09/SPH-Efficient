// We use the radix sort algorithm
// C++ implementation of Radix Sort
#include<iostream>
using namespace std;

// A utility function to get maximum value in arr[]
int getMax(int arr[], int size)
{
    int max = arr[0];
    for (int i = 1; i < size; i++)
        if (arr[i] > max)
            max = arr[i];
    return max;
}
 ///// necesitamos buscar el minimo y maximo de todas las posiciones
double getMaxd(int N, double R[]){
	double max=R[0];
	for(int i=0; i<N;i++){
		if(R[i]>max){
			max=R[i];
		}
	}
	return max;
}
double getMind(int N, double R[]){
	double min=R[0];
	for(int i=0; i<N; i++){
		if(R[i]<min){
			min=R[i];
		}
	}
	return min;
}

void CountingSort(int arr[],int idx[], int size, int div)
{
    int output[size]={};
    int idxout[size]={};
    int count[10] = {0};

    for (int i = 0; i < size; i++)
        count[ (arr[i]/div)%10 ]++;

    for (int i = 1; i < 10; i++)
        count[i] += count[i - 1];

    for (int i = size - 1; i >= 0; i--){
        output[count[ (arr[i]/div)%10] - 1] = arr[i];
        idxout[count[ (arr[i]/div)%10] - 1] = idx[i];
        count[(arr[i]/div)%10]--;
    }

    for (int i = 0; i < size; i++){
        arr[i] = output[i];
        idx[i] = idxout[i];
    }
}


void RadixSort(int arr[],int idx[], int size)
{
    int m = getMax(arr, size);
    for (int div = 1; m/div > 0; div *= 10){
            CountingSort(arr,idx, size, div);
    }
}


void RadixSPH(int N,int key[], int keyS[],int idx[]){
    for(int i=0; i<N; i++){
        keyS[i]=key[i];
    }
	RadixSort(keyS, idx, N);
}

void Discretx(int N, double h_hash, double xmin, double R[], int Rd[]){//In this case we use the h_hash same like h_smoothing lenght
		for(int i=0; i<N;i++){
			Rd[i]=0;
			Rd[i]=int((R[i]-xmin)/h_hash);
//			cout << "Para la particula "<< i <<" la posicion continua es " << R[i] << ", la posicion discreta es " << Rd[i] << "\n";
		}
}

void CensoSPH(int N, double h_hash, double xmin, double R[], int keyS[], int idx[], int idxmin[], int idxmax[], int act[]){
    int Rd[N]={};
	Discretx(N, h_hash, xmin, R, Rd);//With this function, we obtain the position in discret coordinates
	//define the index max for the
	//Define the key value
	int key[N]={}; //Define the key array and the key sort array and the array with the permuted index
	for(int i=0; i<N;i++){
		key[i]=Rd[i]; //Here define the value of key for each particle
		idx[i]=i;
    //        printf("Rd: %d, R: %lf, i: %d, key %d, idx %d\n", Rd[i], R[i], i, key[i], idx[i]);
	}
	//We need 4 arrays:
	// 1.- Index original for particles idx[]
	// 2.- New Index for sort keyS: KeyS[i]: i-> is the New Index
	// 3.-
	RadixSPH(N,key,keyS,idx); //Ordered the values of the Key and give value for the permuted index, with the noral index is the original index, so for

	// the array idx[i] provide the number of the particle (original): idx[i] and i is the permuted index where the key is ordered.
	//The number of classes in this case 1D is the same that the H_x
	// Define the classes
	//POdemos hacer mass eficiente partiendo de un loop sobre todas las particulas e ir revisando cada cuando cambia un escalon
	 // the number of the class is in the index of this arrays
    //	bool act[nclass]={false}; //This array only say us if the cell is ocuped or not
	int numclass;
	numclass=keyS[0];
	act[numclass]=1;
	idxmin[numclass]=0;
	for(int i=1; i<N; i++){
		if(keyS[i]!=numclass){
			idxmax[numclass]=i-1;
			numclass=keyS[i];
			act[numclass]=1;
			idxmin[numclass]=i;
		}
	}
    idxmax[numclass]=N-1;
    for(int i=0; i<N;i++){
//        printf("Rd: %d, R: %lf, i: %d, key %d, idx %d\n", Rd[i], R[i], i, keyS[i], idx[i]);
	}
}



