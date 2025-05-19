# year4-numerical-methods
Skup domaćih zadaća i završni projekt izrađeni u sklopu kolegija Numeričke metode i matematičko modeliranje na četvrtoj godini integriranog studija fizike.

## Projekt - Struktura bijelih patuljaka
Cilj projekta bio je modeliranje bijelog patuljka - iznimno gustih svemirskih objekata koji su završni stadij evolucije zvijezda nedovoljno masivnih da postanu crne rupe ili neutronske zvijezde. Takvi objekti opisani su sustavom vezanih diferencijalnih jednadžbi prvog reda, koje su u ovom projektu u svojem bezdimenzionalnom obliku riješene numerički koristeći Runge-Kutta metodu četvrtog reda. Provjerena je stabilnost rješenja, te su ovim putem dobiveni profili gustoće, mase i radijusi različitih teoretskih zvijezda (različite centralne gustoće). Osim toga, Simpsonovim pravilom riješeni su integrali koji predstavljaju ukupnu kinetičku i gravitacijsku potencijalnu energiju zvijezde. Dobivena rješenja eksportirana su u .txt datoteku koja je ponovno pročitana u Python-u radi jednostavnijeg grafičkog prikaza rezultata.

## Domaće zadaće
### Zadatak 3
Program računa prvu derivaciju funkcije $e^x \ \sin(x)$ koristeći numerički izraz za derivaciju pomoću 3 točke.
Dobivene vrijednosti uspoređene su s poznatim analitičkim vrijednostima derivacije, te je izračunata veličina koraka za koju je relativna pogreška minimalna.

### Zadatak 4
Program računa integral $\int_1^{100} \frac{\exp(-x)}{x}dx$ koristeći trapezno i Simpsonovo pravilo integracije za broj točaka $N=10, 20, 40, 100, 1000, 10000.$

### Zadatak 5
Koristeći Newton-Raphsonovu metodu, program računa i ispisuje nenegativne nultočke Legendreovog polinoma $N$-tog stupnja $P_N(x)$ na intervalu $[-1,1]$.

### Zadatak 6
Program računa integral iz Zadatka 3, no sada Gauss-Legendreovom kvadraturom, te uspoređuje rezultate s onim dobivenim koristeći prethodne metode.

### Zadatak 8
Program rješava Poissonovu jednažbu u 1D: $-u''(x) = (3x+x^2)e^x$, $x\in [0,1]$, s početnim uvjetima $u(0)=u(1)=1$. Koristeći Eigen biblioteku u C++, jednadžba je rješena u obliku $Au=w$, gdje je nad matricom $A$ napravljena LU dekompozicija. Rješenja su dobivena za 2 različite diskretizacije intervala $[0,1]$, te su grafički uspoređena s poznatim analitičkim rješenjem. Osim toga, kod implementira i Thomasov algoritam za rješavanje iste jednadžbe, te uspoređuje vremena izvršavanja ove dvije metode.
