#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double f_ro(const double&, const double&, const double&, const double&);
double f_m(const double&, const double&);
void Runge_Kutta_4(const double&, const double&, const double&, vector<double>&, vector<double>&, vector<double>& );
double podint_funkcija(const double&, vector<double>&, vector<double>&, bool);
double Simpson_pravilo(const vector<double>&, const double&);

int main()
{
    int N = 4000; // broj koraka
    double h = 0.001; // veličina koraka
    vector<double> ro_central;
    ro_central.push_back(0.1);
    ro_central.push_back(10.);
    ro_central.push_back(1000.);
    ro_central.push_back(1e6);

    vector<double> r(N);
    vector<double> ro(N);
    vector<double> m(N);

    // otvaranje datoteke
    ofstream outputFile("output.txt");
    if (!outputFile.is_open()) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    for(size_t i = 0; i<ro_central.size(); i++){
        double ro_c = ro_central[i];
        Runge_Kutta_4(h, N, ro_c, m, ro, r);

        // upisujemo vektore ro i m dobivene RK metodom u tekstualnu datoteku kako bi rezultate mogli nacrtati u Python-u
        for (size_t i = 0; i < N; i++) {
                outputFile << ro[i] << endl;
                outputFile << m[i] << endl;
        }

        // nakon određenog elementa vektora ro i m, svi elementi su "NaN", tražimo zadnji koji nije i odgovarajući element u vektoru r je tada radijus zvijezde
        size_t zadnji_indeks = ro.size();
        for(size_t i = 0; i<ro.size(); i++){
            if (isnan(ro[i])){
                zadnji_indeks = i;
                break; // izlazak iz petlje nakon što smo našli prvu "NaN" vrijednost
            }
        }
        double radijus = r[zadnji_indeks]; // radijus zvijezde dobiven RK metodom
        cout << "ro_c = " << ro_c << endl;
        cout << "R = " << radijus << endl;
        cout << "M = " << m[zadnji_indeks-1] << endl;

        // integracija ukupne kinetičke i gravitacijske potencijalne energije
        vector<double> podint_funkcija_U(zadnji_indeks-1);
        vector<double> podint_funkcija_W(zadnji_indeks-1);
        for(size_t i=0; i<zadnji_indeks-1; i++){
            podint_funkcija_U.push_back(podint_funkcija(r[i], ro, m, true));
            podint_funkcija_W.push_back(podint_funkcija(r[i], ro, m, false));
        }
        cout << "U = " << Simpson_pravilo(podint_funkcija_U, h) << endl;
        cout << "W = " << Simpson_pravilo(podint_funkcija_W, h) << "\n" << endl;
    }

    return 0;
}

double f_ro(const double& m, const double& r, const double& ro)
{
    double x = pow(ro, 1./3); // u vektoru ro podrazumijevamo bezdimenzionalne vrijednosti, pa je u njemu zapravo sadržan parametar ro0
    double gamma = x*x / (3.*sqrt(1.+x*x));
    return -m/gamma * ro/(r*r);
}
double f_m(const double& r, const double& ro)
{
    return r*r * ro;
}

void Runge_Kutta_4(const double& h, const double& N, const double& ro_c, vector<double>& m, vector<double>&ro, vector<double>&r)
{
    m[0] = 0.0;
    ro[0] = ro_c;
    r[0] = 0.001; // ne možemo početi od 0.0 zbog greške koja se događa pri dijeljenju s nulom u funkciji f_ro

    for(int i=0; i<N; i++){
        double k1_ro = h * f_ro(m[i], r[i], ro[i]);
        double k1_m = h * f_m(r[i], ro[i]);
        double k2_ro = h * f_ro(m[i] + k1_m/2.0, r[i] + h/2.0, ro[i] + k1_ro/2.0);
        double k2_m = h * f_m(r[i] + h/2.0, ro[i] + k1_ro/2.0);
        double k3_ro = h * f_ro(m[i] + k2_m/2.0, r[i] + h/2.0, ro[i] + k2_ro/2.0);
        double k3_m = h * f_m(r[i] + h/2.0, ro[i] + k2_ro/2.0);
        double k4_ro = h * f_ro(m[i] + k3_m, r[i] + h, ro[i] + k3_ro);
        double k4_m = h * f_m(r[i] + h, ro[i] + k3_ro);

        r[i + 1] = r[i] + h;
        ro[i + 1] = ro[i] + 1.0/6.0 * (k1_ro + 2*k2_ro + 2*k3_ro + k4_ro);
        m[i + 1] = m[i] + 1.0/6.0 * (k1_m + 2*k2_m + 2*k3_m + k4_m);

    }
    return;
}

double podint_funkcija(const double& r, vector<double>& ro, vector<double>& m, bool UW)
{
    // računamo indeks r-a koji smo dali funkciji da znamo koji ćemo element od ro i m uzeti
    double h = 0.001; // korak
    int i = int(r/h);
    double x = pow(ro[i], 1./3);
    double epsilon = 3./(8.) * (x*(1.+2.*x*x) * sqrt(1.+x*x) - log(x + sqrt(1.+x*x)));
    if (UW == true){ // bool UW određuje koji integral rješavamo, ako je True funkcija će vraćati podintegralnu funkciju od U, a ako je False od W
        return epsilon*r*r;
    }
    else{
        return m[i]*ro[i]*r;
    }
}
double Simpson_pravilo(const vector<double>& funkcija, const double& korak)
{
    double integral = 0.;
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