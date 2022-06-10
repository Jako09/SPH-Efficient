//This code solve Nonlinear Schrödinger equation using smoothed particle hydrodynamics. For diferent physics problems, starting with the harmonic oscillator (HO)
#include "HO_problem.h"
#include "CS-HO_problems.h"
#include "Soliton_problems.h"
#include "BEC_problems.h"

using namespace std;

int main(){
    int Ptype; // Ptype:=Problem type. For this case we have -Harmonic Oscillator->"-HarmOsc", -Dynamics of HO->"HODynamic", -BEC ground states->"GSBEC", -Solitons Collition->"SolCol"
	int N; //# of particles
    int search;
    cout << "Which search you want to use? 1.-standar or 2.-efficient. Write the number" << '\n';
	cin >> search;
    while((search!=1) && (search!=2)){
        cout << "Choose a right answer" << '\n';
        cin >> search;
    }
	cout << "Which problem do you want?" << '\n';
	cout << "Problems:" << '\n';
	cout << "1.-Harmonic Oscillator" << '\n'; //Give the relaxation base states for HO, in different times, and animated
	cout << "2.-Coherent States" << '\n';
	cout << "3.-Groud States of Bose-Einstein Condensate" << '\n';
	cout << "4.-Bright Soliton Collision" << '\n';
	cout << "Write the number" << '\n';
	cin >> Ptype;
	while((Ptype!=1)&&(Ptype!=2)&&(Ptype!=3)&&(Ptype!=4)){
			cout << "Wrong number! Writte the right number"<<'\n';
			cin >> Ptype;
			}
	cout << "You choose the problem " << to_string(Ptype)  << '\n';
	cout << "How many particles do you want to use?" << '\n';
	cout << "We recomend you 640, 1280 and 2560, because, we know which smoothing lenght is requiered" << '\n';
	cout << "You can choose another particle number, but previusly you need search the smoothing length factor (hf in the code)" << '\n';
	cin >> N;
	while(N<=0){
		cout << "The number of particles need to be a positive, try again:"<< '\n';
		cin >> N;
		}

    if(Ptype==1){
        HarmonicOscillator(N,search);
    }
    if(Ptype==2){
        CoherentStates(N,search);
    }
    if(Ptype==3){
        SolitonCollision(N,search);
    }
    if(Ptype==4){
        BECGroundStates(N,search);
    }
    return 0;
}

