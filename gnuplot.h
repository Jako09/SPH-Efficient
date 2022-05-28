#include<fstream>

using namespace std;

void printGraphtimeEff(string format){
    ofstream file("timeEff.gnu");
    if(format=="EPS"){
        file << "set term postscript eps enhanced color" << '\n';
        file << "set out 'timeEff.eps'" << '\n';
    }if(format=="PNG"){
        file << "set terminal png size 1000,700" << '\n';
        file << "set out 'timeEff.png'" << '\n';
    }
    file << "set grid" << '\n';
    file << "set title 'Tiempo vs Numero de particulas'" << '\n';
    file << "set ylabel 't:= tiempo de computo'"<<'\n';
    file << "set xlabel 'N:= numero de particulas'"<<'\n';
    file << "f(x)=a*x+b" << '\n';
    file << "g(x)=c*x**2+d*x+e" << '\n';
    file << "FIT_LIMIT=1e-6" << '\n';
    file << "fit f(x) 'Efftime.dat' u ($1):($3) via a,b" << '\n';
    file << "fit g(x) 'Efftime.dat' u ($1):($2) via c,d,e" << '\n';
    file << "p 'Efftime.dat' u ($1):($2) t 'Estandar' lt 1, '' u ($1):($3) t 'Eficiente' lt 2, f(x) t 'Fit Eficiente' lt 2, g(x) t 'Fit Estandar' lt 1" << '\n';
    file << "pause -1" << '\n';
    file.close();
}

void printGraphDensity(string format, int iterations, string term){
    ofstream file("densidadEff.gnu");
    file << "set term "+term << '\n';
    file << "set grid" << '\n';
    file << "set title 'Posicion vs Densidad'" << '\n';
    file << "set ylabel '{/Symbol r}(r):=Densidad'" << '\n';
    file << "set xlabel 'r:=Posicion'" << '\n';
    file << "n="+to_string(iterations) << '\n';
    file << "do for [it=0:n]{" << '\n';
    file << "p 'Densidad.dat' i it u ($1):($2) lt 1 t 'Estandar','' i it u ($1):($3) lt 2 t 'Eficiente' "<< '\n';
    file << "pause 1" << '\n';
    file << "}" << '\n';
}
