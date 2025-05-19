#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace chrono;

double funkcija(const double&);
MatrixXd LU(const int&);
vector<double> Thomas(const int&);
double u_analiticko(const double&);

int main()
{
    int N = 10;

    // ovaj dio je napravljen da bi dobili vremena izvođenja dvaju algoritama koja su kopirana u python gdje je napravljen traženi graf
    vector<int> n;
    vector<double> vremena1, vremena2;
    n.push_back(10);
    n.push_back(100);
    n.push_back(1000);
    n.push_back(2000);
    n.push_back(3000);
    n.push_back(4000);
    n.push_back(5000);
    
    for(int i=0; i<n.size(); i++){
        MatrixXd u(n[i],1);
        auto start1 = high_resolution_clock::now();
        u = LU(n[i]);
        auto stop1 = high_resolution_clock::now();
        auto duration1 = duration_cast<duration<double> >(stop1-start1);
        vremena1.push_back(duration1.count());

        auto start2 = high_resolution_clock::now();
        vector<double> u_Thomas;
        u_Thomas = Thomas(n[i]);
        auto stop2 = high_resolution_clock::now();
        auto duration2 = duration_cast<duration<double> >(stop2-start2);
        vremena2.push_back(duration2.count()); 
    }
 
    // rješenje za N = 10
    MatrixXd u(N,1);
    u = LU(N); // rješenja za u ispisana su u prikladnom formatu i kopirana u python kod gdje su nacrtani odgovarajući grafovi
    vector<double> u_Thomas;
     u_Thomas = Thomas(N);
    
    vector<double> xi;
    for(int i = 0; i<N+2; i++){
        xi.push_back(i/(double(N)+1));
    }

    cout << setw(15) << "x_i" << setw(15) << "u_i" << setw(15) << "u_i (LU)" << setw(15) << "u_i (Thomas)" << endl;
    // rezultate za x_0 i x_n ćemo ispisati posebno jer početni uvjeti nisu uključeni u algoritam, pa te vrijednosti nisu u vektorima koje smo dobili kao rezultat funkcija koje provode te algoritme
    // ali znamo da su = 0 jer je ta pretpostavka nužna da algoritam funkcionira
    cout << setw(15) << setprecision(10) << xi[0];
    cout << setw(15) << setprecision(10) << u_analiticko(xi[0]);
    cout << setw(15) << setprecision(10) << 0.0;
    cout << setw(15) << setprecision(10) << 0.0 << endl;
    for(int i=1; i<N+1; i++){
        cout << setw(15) << setprecision(10) << xi[i];
        cout << setw(15) << setprecision(10) << u_analiticko(xi[i]);
        cout << setw(15) << setprecision(10) << u(i-1,0);
        cout << setw(15) << setprecision(10) << u_Thomas[i-1] << endl;
    }
    cout << setw(15) << setprecision(10) << xi[N+1];
    cout << setw(15) << setprecision(10) << u_analiticko(xi[N+1]);
    cout << setw(15) << setprecision(10) << 0.0;
    cout << setw(15) << setprecision(10) << 0.0 << endl;

    return 0;
}

double funkcija(const double& x)
{
    return (3*x+x*x) * exp(x);
}

MatrixXd LU(const int& N)
{
    MatrixXd A = MatrixXd::Zero(N,N);
    for(int i=0; i<N; i++){
        A(i,i) = 2;
        if(i != N-1){
            A(i, i+1) = -1;
            A(i+1, i) = -1;
        }
    }

    double h = 1/(double(N)+1);
    MatrixXd w(N,1);
    for(int i=0; i<N; i++){
        w(i,0) = h*h*funkcija((i+1)*h);
    }
    MatrixXd u(N,1);
    u = A.partialPivLu().solve(w);
    return u;
}

vector<double> Thomas(const int& N)
{
    MatrixXd A = MatrixXd::Zero(N,N); // matrica A s koeficijentima je ista kao i kod LU metode, svi a_i = c_i = -1, a svi b_i = 2
    for(int i=0; i<N; i++){
        A(i,i) = 2;
        if(i != N-1){
            A(i, i+1) = -1;
            A(i+1, i) = -1;
        }
    }

    vector<double> d, a, b, c, x;
    double h = 1/(double(N)+1);
    for(int i=0; i<N; i++){
        d.push_back( h*h*funkcija((i+1)*h) );
        a.push_back(-1.0);
        c.push_back(-1.0);
        b.push_back(2.0);
        x.push_back(0.0);
    }
    for(int k=1; k<N; k++){
        double m = a[k]/b[k-1];
        b[k] -= m*c[k-1];
        d[k] -= m*d[k-1];
    }
    x[N-1] = d[N-1] / b[N-1];
    for(int k=N-2; k>=0; k--)
    {
        x[k] = (d[k] - c[k]*x[k+1]) / b[k];
    }
    return x;
}

double u_analiticko(const double& x)
{
    return x*(1-x) * exp(x);
}