#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double legendre(const int&, const double&);
double legendre_derivacija(const int&, const double&);
double newton_metoda(const double&, const double&, const double&);

int main()
{
    int N;
    double granica;
    double preciznost = 1.e-10;
    vector<double> nultocke;
    cout << "Unesite broj N: " << endl;
    cin >> N;

    if (int(N)%2 == 0){
        granica = N/2; 
    }
    else {
        granica = (N+1)/2;
    }
    // petlja po svim nultočkama
    for(int i=1; i<=granica; i++)
    {
        double x0 = cos(M_PI*(double(i)-0.25) / (N+0.25));
        nultocke.push_back(newton_metoda(x0, preciznost, N));
    }

    // ispis nultočaka
    sort(nultocke.begin(), nultocke.end());
    for (size_t i=0; i<nultocke.size(); i++){
        cout << nultocke[i] << endl;
    }

    return 0;
}

double legendre(const int& n, const double& x)
{
    vector<double> P;
    P.push_back(1.);
    P.push_back(x);
    for(int i = 2; i<=n; i++){
        P.push_back(1/(double)i * (double)(2*i-1)*x*P[i-1] - (double)(i-1)/(double)i*P[i-2]);
    }
    return P[n];
}
double legendre_derivacija(const int& n, const double& x)
{
    vector<double> P;
    P.push_back(1.);
    P.push_back(x);
    for(int i = 2; i<=n; i++){
        P.push_back(1/(double)i * (double)(2*i-1)*x*P[i-1] - (double)(i-1)/(double)i*P[i-2]);
    }

    return (double)n / (1-x*x) * (P[n-1]-x*P[n]);
}

double newton_metoda(const double& x0, const double& preciznost, const double& n)
{
    int max_iteracija = 1000;
    double tocka = x0;
    for (int i=0; i<max_iteracija; i++)
    {
        if(abs(legendre(n, tocka)) <= preciznost){
            return tocka;
        }
        else {
            tocka -= legendre(n, tocka)/legendre_derivacija(n, tocka); 
        }
    }
    return tocka;
}