#include <string.h>
#include <fstream>
#include <iostream>

using namespace std;

void L2graph(int nfiles, string names[], string kind);
void TimeGraph(int nfiles, string names[], string kind);

int main(){
    int nfiles;
    int kind;
    cout << "This program is for help to print all the file for L2 and time process that you have in the respective folders" << '\n';
    cout << "What graph do you want? 1.- time or 2.-L2" << '\n';
    cin >> kind;
    cout << "How many files there are?" << '\n';
    cin >> nfiles;
    string names[nfiles];
    if(kind==1){
        cout << "Put the names of each datafile" << '\n';
        for(int i=0; i<nfiles;i++){
            cin >> names[i];
            cout << '\n';
        }
        TimeGraph(nfiles, names, "EPS");
        TimeGraph(nfiles, names, "PNG");
    }
    if(kind==2){
        cout << "Put the names of each datafile" << '\n';
        for(int i=0; i<nfiles;i++){
            cin >> names[i];
            cout << '\n';
        }
        L2graph(nfiles, names, "EPS");
        L2graph(nfiles, names, "PNG");
    }
    return 0;
}

void L2graph(int nfiles, string names[], string kind){
    ofstream file;
	if(kind=="EPS"){
		file.open("NormL2_HO_eps.gnu");
		file << "set term postscript eps enhanced color" << '\n';
        file << "set out 'NormL2.eps'" << '\n';
	}
	if(kind=="PNG"){
		file.open("NormL2_HO_png.gnu");
		file << "set terminal png size 1000,700" << '\n';
        file << "set out 'NormL2.png'" << '\n';
    }
	file << "set grid" << '\n';
	file << "set title 'Norma L2'" << '\n';
	file << "set ylabel 'L2({/Symbol r}):= Error Norma L2'" << '\n';
	file << "set xlabel 'N:=Numero de particulas'" << '\n';
	file << "  p ";
	for(int i=0; i<nfiles-1; i++){
        file << "\""+names[i]+"\"  u ($1):($2) t '', \\" << '\n';
    }
    file << "\""+names[nfiles-1]+"\"  u ($1):($2) t ''" << '\n';
	file << "pause -1";
	file.close();
}

void TimeGraph(int nfiles, string names[], string kind){
    ofstream file;
	if(kind=="EPS"){
		file.open("Time_Operation_HO_eps.gnu");
		file << "set term postscript eps enhanced color" << '\n';
        file << "set out 'time_operation.eps'" << '\n';
	}
	if(kind=="PNG"){
		file.open("Time_Operation_HO_png.gnu");
		file << "set terminal png size 1000,700" << '\n';
        file << "set out 'time_operation.png'" << '\n';
	}
	file << "set grid" << '\n';
	file << "set xrange [-4.0:4.0]" << '\n';
	file << "set title 'Time operation'" << '\n';
	file << "set ylabel 't:= tiempo'" << '\n';
	file << "set xlabel 'N:=Numero de particulas'" << '\n';
//	if(search==1){
//		file << "	p   \""+names+"\"  u ($1):($2) t 'Estandar'" << '\n';
//	}
//	if(search==2){
//		file << "	p   \""+names+"\"  u ($1):($2) t 'Eficiente'" << '\n';
//	}
	file << "pause -1";
	file.close();
}
