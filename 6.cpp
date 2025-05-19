#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

double legendre_polinom(const int&, const double&);
double legendre_derivacija(const int&, const double&);
double newton_metoda(const double&, const double&, const double&);
double podint_funkcija(const double&);
double integral_Gauss(const double&);
double trapez_pravilo(const vector<double>&, const double&);
double Simpson_pravilo(const vector<double>&, const double&);

int main()
{
    vector<double> N;
    N.push_back(10);
    N.push_back(20);
    N.push_back(40);
    N.push_back(100);
    N.push_back(1000);
    cout << setw(10) << "N" << setw(10) << "trapez" << setw(10) << "Simpson" << setw(10) << "Gauss" << endl;

    for(int i=0; i<N.size(); i++){
        cout << setw(10) << N[i];
        double korak = (100-1)/N[i];
        vector<double> podintegralna_funkcija;
        for(int j=0; j<=N[i]; j++)
        {
            double x = 1 + j*korak;
            podintegralna_funkcija.push_back(podint_funkcija(x));
        }
        cout << setw(10) << setprecision(6) << trapez_pravilo(podintegralna_funkcija, korak);
        cout << setw(10) << setprecision(6) << Simpson_pravilo(podintegralna_funkcija, korak);
        cout << setw(10) << setprecision(6) << integral_Gauss(N[i]) << endl;
    }

    return 0;
}

double legendre_polinom(const int& n, const double& x)
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
        if(abs(legendre_polinom(n, tocka)) <= preciznost){
            return tocka;
        }
        else {
            tocka -= legendre_polinom(n, tocka)/legendre_derivacija(n, tocka); 
        }
    }
    return tocka;
}

double podint_funkcija(const double& x)
{
    return exp(-x)/x;
}

double integral_Gauss(const double& red)
{
    double N = red;
    double preciznost = 1.e-6;
    vector<double> nultocke_pozitivne, nultocke;

    // račun nultočki u 2 slučaja (parni i neparni N)
    if (int(N)%2 == 0){
        for(int i=1; i<=N/2; i++){
            double x0 = cos(M_PI*(double(i)-0.25) / (N+0.25));
            nultocke_pozitivne.push_back(newton_metoda(x0, preciznost, N));
        }
        sort(nultocke_pozitivne.begin(), nultocke_pozitivne.end());

        for(int i=1; i<=N/2; i++){
            nultocke.push_back(-nultocke_pozitivne[nultocke_pozitivne.size()-i]);
        }
        for(int i=0; i<N/2; i++){
            nultocke.push_back(nultocke_pozitivne[i]);
        }
    }
    else{
        for(int i=1; i<=(N+1)/2; i++){
            double x0 = cos(M_PI*(double(i)-0.25) / (N+0.25));
            nultocke_pozitivne.push_back(newton_metoda(x0, preciznost, N));
        }
        sort(nultocke_pozitivne.begin(), nultocke_pozitivne.end());

        for(int i=1; i<=(N-1)/2; i++){
            nultocke.push_back(-nultocke_pozitivne[nultocke_pozitivne.size()-i]);
        }
        for(int i=0; i<(N+1)/2; i++){
            nultocke.push_back(nultocke_pozitivne[i]);
        }
    }
    
    // račun težinskih faktora
    vector<double> tezinski_faktori;
    for(int i=0; i<nultocke.size(); i++){
        tezinski_faktori.push_back( 2 / ( (1-pow(nultocke[i],2)) * pow(legendre_derivacija(red, nultocke[i]), 2) ) );
    }

    //račun integrala
    double integral = 0;
    double a=0, b=1;
    for(int i=0; i<nultocke.size(); i++){ // promjena granica integrala iz 0 --> 1 u -1 --> 1
        integral += (b-a)/2 * tezinski_faktori[i] * podint_funkcija( (b-a)/2*nultocke[i] + (b+a)/2 );
    }
    return integral;
}

double trapez_pravilo(const vector<double>& funkcija, const double& korak)
{
    double integral = 0;
    integral += ( funkcija[0] + funkcija[funkcija.size()] ) / 2;
    for(size_t i=1; i<funkcija.size(); i++){
        integral += funkcija[i];
    }
    integral *= korak;

    return integral;
}

double Simpson_pravilo(const vector<double>& funkcija, const double& korak)
{
    double integral = 0;
    integral += funkcija[0] + funkcija[funkcija.size()];
    for(size_t i=1; i<funkcija.size(); i++){
        if(i%2 == 0){
            integral += 2*funkcija[i];
        }
        else{
            integral += 4*funkcija[i];
        }
    }
    integral *= korak/3;

    return integral;
}