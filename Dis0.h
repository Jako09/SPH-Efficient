#include<random>
#include<iostream>
#include"m&h.h"

using namespace std;
//Se regresa un valor arbitrario para las posiciones entre un rango [xmin,xmax], tal que la distribución es uniforme.
void UniformD(int N, double xmin, double xmax, double R[], double m[]){
  random_device rd; /*Con esto se obtiene la semilla*/
  mt19937 gen(rd());
  uniform_real_distribution<> dis(xmin,xmax);
  
  for(int n=0;n<N;++n){
    R[n]=dis(gen);
    m[n]=1.0/(double)N;  
  }

}
//Se definen los puntos del grid, donde z es la distancia entre los puntos del grid
double grid0(int i, double min, double z){
	return min+z*(double)(i);
}
//Se obtienen los valores de la densidad en los puntos del Grid y obtenemos el área por definición del paralelogramo
double grid1(int Ptype, int i, double min, double z){
	double B1=DA(Ptype, grid0(i, min, z));
	double B2=DA(Ptype, grid0(i+1, min, z));
	return z*(B1+B2)/2.0;
}
//Se obtienen los valores donde se van a colocar las partículas de acuerdo a la distribución analítica. 
void grid(int Ptype, int N, double m[], double R[]){
	int dots=10000*N;
	double z, l1, l2, min, max, sum=0.0, mparcial;
	if(Ptype==1){
		min=0.0;
		max=5.0;
	}
	if(Ptype==2){
		min=-4.0;
		max=4.0;
	}
	if(Ptype==3){
		min=0.0;
		max=1.0;
	}
	if(Ptype==4){
		min=0.0;
		max=4*PI;
	}
	if(Ptype==5){
		min=-0.85;
		max=0.85;
	}
	if(Ptype==6){
		min=-15.0;
		max=15.0;
		}
	if(Ptype==7){
		min=-15.0;
		max=15.0;
		}
//Se define la separación entre los puntos que definen el grid
	z=(max-min)/((double)dots);
//Se suma el área de cada uno de los paralelogramos dando como resultado el área total aproximada, debajo de la curva y entre el eje coordenado, entre los limites xmax y xmin.	
	for(int i=0; i<dots; ++i){
		sum+=grid1(Ptype, i, min, z);
//		cout << sum << '\n';
	}
//Se obtiene la masa parcial para cada una de las partículas
	mparcial=sum/((double)N);
	cout << mparcial << '\n';
	int ia=0;
//Se suman las masas comenzando del primer paralelogramo, susecivamente hasta llegar a un valor	mparcial, y se detiene el contador ia, que nos dice en que nodo del grid estamos, en donde se logró detener, usamos el valor ai para hallar  la posción en donde se va colocar la partícula.
	for(int i=0; i<N; i++){
		m[i]=0.0;
		double ka=0.0, l1=grid0(ia, min, z);
//		cout << m[i] << "masa" << '\n';
		while(ka<mparcial){
			ka+=grid1(Ptype, ia, min, z);
			ia++;
		}
		ia=ia-1; //
		l2=grid0(ia, min, z);
		R[i]=(l2+l1)/2.0;
		m[i]=mparcial;
	}
}
//Las partículas equidistantes a lo largo de un intervalo con extremos xmax y xmin.
void Glasslike(int N, double xmin, double xmax, double R[], double m[]){
	double z=0.0;
//Se define la separación entre partículas
	z=(xmax-xmin)/(double)(N);
//Se colan las partículas equidistantes	
	for(int i=0; i<N; i++){
		R[i]=xmin+((double)(i)+0.5)*z;
		m[i]=1.0/(double)N;	
	}	
}


