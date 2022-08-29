/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

// DICHIARAZIONE DELLE FUNIONI PER CALCOLO DI INTEGRALE DI rho*f
// - rho densità di probabilità;
// - f funzione integranda;
// - V potenziale;
// - e_plus, e_minus, funzioni per codice più compatto
double V(double x);
double V_OA(double x);
double e_plus(double x, double mu, double sigma);
double e_minus(double x, double mu, double sigma);
double rho(double x, double mu, double sigma);
double f(double x, double mu, double sigma); 
double f_OA(double x, double mu, double sigma); 
double min(double a, double b);

 
int main (int argc, char *argv[]){

	// RACCOLGO VARIABILI INSERITE DA LINEA DI COMANDO

	int POT=0;
	double mu=0;
	double sigma=0;
	double walk=0;
	if(argc!=5){
		cerr << "Usage: ./main.exe POT mu sigma walk" << endl;
		cerr << "POT=0 --> OA, POT=1 --> V=x^4-2.5x^2" << endl;
		return -1;
	}
	else{
		POT = atoi(argv[1]);
		mu = atof(argv[2]);
		sigma = atof(argv[3]);
		walk = atof(argv[4]);
	}
	

	// PREPARAZIONE GENERATORE NUMERI CASUALI	
	
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


  // PREPARAZIONE DATA-BLOCKING E ALGORTIMO METROPOLIS

	 int M = 100000; //total attemps
	 int N = 100; // total blocks
	 int L = M/N; // attemps per block

	 double x_old = 0; // initial point
   double x_new = 0; //final point

	 double p=0; // probabilità di muoversi

	 double dato = 0;
	 double somma_prog = 0;
	 double somma_prog_2 = 0;

	 double accepted=0; // to verify 50/50 rule
	 double attempted = (double) L; // mosse per blocco

	 ofstream output("integrale.txt"); // stampa stima, media progressiva, errore progressivo

	
   // PREPARO PER RECORD PUNTI ESTRATTI

   int nbins= 500;
   double * bin = new double[nbins];
   for(int i=0; i<nbins; i++) bin[i]=0;	
	 int appo =0;

	
	// DATA BLOCKING
	
	 for(int i=0; i<N; i++){ // ciclo esterno su blocchi
     dato = 0;
		 accepted = 0;
		 
		 for(int j=0; j<L; j++){ // ciclo interno: Metropolis nel blocco
			 x_new = x_old + rnd.Rannyu(-walk, +walk);
			 p = min(1, rho(x_new,mu,sigma)/rho(x_old,mu,sigma));
			 if(p>rnd.Rannyu()){
				 x_old = x_new;
				 accepted = accepted + 1.;
			 }

			 appo=(int) floor( (x_old+3.)/6. * nbins);
       if(appo>=0 && appo<nbins) bin[appo] ++;
			 
			 if(POT==0) dato += f_OA(x_old,mu,sigma);	
			 else dato += f(x_old,mu,sigma);
		 }
		 
	   cout << "Block n° " << i+1 << ". % accepted: " << accepted/attempted << endl;
		 cout << endl << "------------------------------" << endl;

		 somma_prog += dato/(double)L;
		 somma_prog_2 += pow(dato/(double)L,2);	 

		 output << setw(15) << dato/(double)L;
		 output << setw(15) << somma_prog/(double)(i+1);
		 if (i==0) output << setw(15) << 0;
		 else output << setw(15) << sqrt((somma_prog_2/(double)(i+1) - pow(somma_prog/(double)(i+1),2))/(double)i);
		 if(i!=N-1) output << endl;
	
	 } 


  // STAMPO VALORI PUNTI ESTRATTI PER ISTOGRAMMA IN (-3,3)
   
	 double normal=0;
   ofstream output2("histo.txt");
   for(int i=0; i<nbins; i++){
	 	normal+=bin[i];
	 }
	
   for(int i=0; i<nbins; i++){
	   bin[i]=bin[i]/normal*nbins/6.;
	   output2 << setw(15) << -3.+i*6./(double)nbins << setw(15) << bin[i] << endl;
}

	 output.close();
	 output2.close();
	 delete[] bin;
   rnd.SaveSeed();
   return 0;
}


// DEFINIZIONE DI FUNZIONI

double V(double x){
	return ( pow(x,4) - 2.5*pow(x,2) );
}
double V_OA(double x){
	return ( pow(x,2)/2. );
}

double e_minus(double x, double mu, double sigma){
	return exp(- pow(x-mu,2) / (2*sigma*sigma) );
}
double e_plus(double x, double mu, double sigma){
	return exp(- pow(x+mu,2) / (2*sigma*sigma) );
}

double rho(double x, double mu, double sigma){
	return  pow( ( e_minus(x,mu,sigma) + e_plus(x,mu,sigma) ) , 2 ); 
}

double f(double x, double mu, double sigma){
	double a = V(x) + 1./(2.*sigma*sigma);
	double k = 1./(2.*pow(sigma,4));
	double b = pow(x-mu,2)*e_minus(x,mu,sigma) + pow(x+mu,2)*e_plus(x,mu,sigma);
	double c = e_minus(x,mu,sigma) + e_plus(x,mu,sigma);
	return a-k*b/c;
}

double f_OA(double x, double mu, double sigma){
	double a = V_OA(x) + 1./(2*sigma*sigma);
	double k = 1./(2*pow(sigma,4));
	double b = pow(x-mu,2)*e_minus(x,mu,sigma) + pow(x+mu,2)*e_plus(x,mu,sigma);
	double c = e_minus(x,mu,sigma) + e_plus(x,mu,sigma);
	return a-k*b/c;
}

double min(double a, double b){
	if(a<b) return a;
	else return b;
}