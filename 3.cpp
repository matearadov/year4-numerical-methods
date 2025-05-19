#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
using namespace std;

double analiticka_derivacija(const double&);
double num_derivacija(const double&, const double&);
void funkcija(const vector<double>&, const vector<double>&, const vector<double>&, vector<double>&, vector<int>&);

int main()
{
    vector<double> tocke, koraci, rel_greska;
    vector<int> indeksi;
    for(int i=-5; i<=5; i++){
        tocke.push_back( (double) i ); // popunjavamo vektor tocke
        rel_greska.push_back(0);
        indeksi.push_back(0);
    }
    for(int j=0; j<15; ++j){ // popunjavamo vektor koraci
        koraci.push_back( (double) pow(10., (double) (-j-1)));
    }

    vector<double> analiticko;
    for(size_t i=0; i<tocke.size(); i++){ // popunjavamo vektor s analitičkim vrijednostima derivacije za svaki element vektora tocke
        analiticko.push_back( analiticka_derivacija(tocke[i]) );
    }

    funkcija(tocke, koraci, analiticko, rel_greska, indeksi);

    // ispis
    cout << setw(6) << "x ";
    for(size_t i=0; i<tocke.size(); i++){
        cout << setw(12) << setprecision(7) << tocke[i] << " ";
    }
    cout << endl << setw(6) << "f'(x) ";
    for(size_t i=0; i<tocke.size(); i++){
        cout << setw(12) << setprecision(7) << analiticko[i] << " ";
    }
    cout << endl << setw(6) << "r ";
    for(size_t i=0; i<rel_greska.size(); i++){
        cout << setw(12) << setprecision(7) << rel_greska[i] << " ";
    }
    cout << endl << setw(6) << "h ";
    for(size_t i=0; i<indeksi.size(); i++){
        cout << setw(12) << setprecision(7) << koraci[int(indeksi[i])] << " ";
    }

    return 0;
}

double analiticka_derivacija(const double& x)
{
    return ( exp(x)*sin(x) + exp(x)*cos(x) );
}
double num_derivacija(const double& x, const double& h)
{
    return ( exp(x+h)*sin(x+h) - exp(x-h)*sin(x-h) ) / (2.0*h);
}

void funkcija(const vector<double>& tocke, const vector<double>& koraci, const vector<double>& analiticko, vector<double>& rel_greska, vector<int>& indeksi)
{
    for(size_t i=0; i<tocke.size(); i++) // petlja po elementima od x
    {
        double min_greska = 100;
        int min_indeks; // indeks vektora h za koji je relativna greška minimalna
        for(size_t j=0; j<koraci.size(); j++) // petlja po elementima od h
        {
            double num = num_derivacija(tocke[i], koraci[j]);
            double rel_temp; // relativna greška za trenutni par x,h
            rel_temp = abs(num - analiticko[i]) / analiticko[i];
            if (rel_temp < min_greska){
                min_greska = rel_temp;
                min_indeks = (int) j;
            }
        }
        rel_greska[i] = min_greska;
        indeksi[i] = min_indeks;
    }
    return;
}