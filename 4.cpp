#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

vector<double> racun_funkcije(const double&, const double&, const int&, const double&);
double trapez_pravilo(const vector<double>&, const double&);
double Simpson_pravilo(const vector<double>&, const double&);

int main()
{
    double a=1, b=100;
    vector<int> N;
    N.push_back(10);
    N.push_back(20);
    N.push_back(40);
    N.push_back(100);
    N.push_back(1000);
    N.push_back(10000);

    vector<double> trapez_integrali, Simpson_integrali;

    for(size_t i=0; i<N.size(); i++)
    {
        double h = (b-a)/N[i];
        vector<double> podint_funkcija = racun_funkcije(a, b, N[i], h);
        trapez_integrali.push_back(trapez_pravilo(podint_funkcija, h));
        Simpson_integrali.push_back(Simpson_pravilo(podint_funkcija, h));
    }

    // ispis
    cout << setw(10) << "N ";
    for(size_t i=0; i<N.size(); i++){
        cout << setw(11) << N[i] << " ";
    }
    cout << endl << setw(10) << "I_trapez ";
    for(size_t i=0; i<trapez_integrali.size(); i++){
        cout << setw(11) << setprecision(8) << trapez_integrali[i] << " ";
    }
    cout << endl << setw(10) << "I_Simpson ";
    for(size_t i=0; i<Simpson_integrali.size(); i++){
        cout << setw(11) << setprecision(8) << Simpson_integrali[i] << " ";
    }

    return 0;
}

vector<double> racun_funkcije(const double& a, const double& b, const int& N_i, const double& h_i)
{
    vector<double> podintegralna_funkcija;
    for(int i=0; i<=N_i; i++)
    {
        double x = a + i*h_i;
        podintegralna_funkcija.push_back(exp(-x)/x);
    }
    return podintegralna_funkcija;
}

double trapez_pravilo(const vector<double>& funkcija, const double& korak)
{
    double integral = 0;
    integral += ( funkcija[0] + funkcija[funkcija.size()] ) / 2; // pribrajamo srednju vrijednost u početnoj i završnoj točki
    for(size_t i=1; i<funkcija.size(); i++){ // pribrajamo sumu iz analitičkog izraza
        integral += funkcija[i];
    }
    integral *= korak;

    return integral;
}

double Simpson_pravilo(const vector<double>& funkcija, const double& korak)
{
    double integral = 0;
    integral += funkcija[0] + funkcija[funkcija.size()]; // pribrajamo vrijednosti funkcije u početnoj i završnoj točki
    for(size_t i=1; i<funkcija.size(); i++){ // pribrajamo obje sume iz analitičkog izraza
        if(i%2 == 0){
            integral += 2*funkcija[i]; // suma za parne vrijednosti k
        }
        else{
            integral += 4*funkcija[i]; // suma za neparne vrijednosti k
        }
    }
    integral *= korak/3;

    return integral;
}