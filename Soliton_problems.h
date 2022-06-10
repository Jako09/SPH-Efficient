
using namespace std;

void SolitonCollision(int N,int search){
    int Skind;
    cout << "Which kind of soliton do you want?" << '\n';
    cout << "1.-Bright or 2.-Dark? Write the number" << '\n';
    cin >> Skind;
    while(Skind!=1&&Skind!=2){
        cout << "You choose a wrong number, try again" << '\n';
        cin >> Skind;
    }
    cout << "For this problem you obtain two graphs" << '\n';
    cout << "First the evolution of the soliton collision, second the evolution of the energy for this collision" << '\n';
    if(Skind==1){
        cout << "You choose bright soliton"<<'\n';
    }
    if(Skind==2){
        cout << "You choose dark soliton" << '\n';
    }
}
