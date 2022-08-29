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
// - integral fa l'integrale per una data trial function

double V(double x);
double V_OA(double x);
double e_plus(double x, double mu, double sigma);
double e_minus(double x, double mu, double sigma);
double rho(double x, double mu, double sigma);
double f(double x, double mu, double sigma); 
double f_OA(double x, double mu, double sigma); 
double min(double a, double b);

double integral(double mu, double sigma, Random* rnd, int M, int N, double* x_in, double walk, bool move, bool print);

 
int main (){
	
	// PREPARAZIONE GENERATORE NUMERI CASUALI	
	
   Random* rnd = new Random;
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
            rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

	// FILE OUPUT

	 ofstream output1("stime_energia.txt"); // energia per oguno dei (100 times)*(41 steps)
	 ofstream output2("stime_parametri.txt"); // parametri per ognuno dei (100 times)* (41 steps)
	 ofstream output3("parametriSA.txt"); // parametri per ognuno dei 41 steps
	 ofstream output5("best.txt");
	
  // PREPARAZIONE SA

	 int M = 100000; // param per calcolo integrale
	 int N = 100; 
	 double* x_in = new double;
	 *x_in=0;
	 double walk = 2.5;
		 
	 double T = 10; // temp iniz, tasso diminuzione, indici
	 double R = 0.75;
	 int iT = 1;
	 int times = 100; // stati campionati per ogni temp
	 double Tstop = 0.0001;
	 int accepted = 0;
   int attempted = 0;

   double mu_old = 1; // parametri
	 double sigma_old = 1;
	 double en_old = integral(mu_old,sigma_old,rnd,M,N,x_in,walk,1,0);
	 double mu_new, sigma_new, en_new;
	 double delta=0.05;
   double mu_min = mu_old; // record minimo e minimizzanti
   double sigma_min = sigma_old;
	 double minimo = en_old;


	 while(T>Tstop){
		 for(int i=0; i<times; i++){
			 mu_new = rnd-> Rannyu(mu_old-delta, mu_old+delta);
			 sigma_new = sigma_old * rnd-> Rannyu(1-delta, 1+delta);
       en_new = integral(mu_new,sigma_new,rnd,M,N,x_in,walk,1,0);
       if(en_new<minimo){
				 minimo = en_new;
				 mu_min = mu_new;
				 sigma_min = sigma_new;
				 cout << "better" << endl;
			 }
			 if(rnd->Rannyu()<min(1,exp((en_old-en_new)/T))){
				 mu_old = mu_new;
				 sigma_old = sigma_new;
				 en_old = en_new;
				 accepted ++;
			 }
			 attempted ++;
			 output1 << en_old << endl;
			 output2 << setw(15) << mu_old << setw(15) << sigma_old << setw(15)<< *x_in << endl;
		 }
		 cout << "STEP " << iT << ". Temp: " << T << endl;
		 cout << "Acc rate: " << (double)accepted/(double)attempted << endl;
		 cout << "En: " << integral(mu_old,sigma_old,rnd, M, N, x_in, walk, 0, 1) << endl;
		 output3 << setw(15) << iT << setw(15) << T << setw(15) << mu_old << setw(15) << sigma_old << endl;
		 accepted =0;
		 attempted =0;
     T*=R;
		 iT=iT+1;
	 }

	 output5 << setw(15) << minimo << setw(15) << mu_min << setw(15) << sigma_min;

	 output1.close();
	 output2.close();
	 output3.close();
   rnd->SaveSeed();
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
	double a = V(x) + 1./(2*sigma*sigma);
	double k = 1./(2*pow(sigma,4));
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


double integral(double mu, double sigma, Random* rnd, int M, int N, double* x_in, double walk,bool move, bool print){
	
	 int L = M/N; // attemps per block
	 double x_old;
	 if(move) x_old = *x_in;
	 else x_old=0; // initial point
   double x_new = 0; //final point
	 double p=0; // probabilità di muoversi
	 double dato = 0;
	 double somma_prog = 0;
	 double somma2_prog = 0;

	 for(int i=0; i<N; i++){ // ciclo esterno su blocchi
     dato = 0;
		 for(int j=0; j<L; j++){ // ciclo interno: Metropolis nel blocco
			 x_new = x_old + rnd->Rannyu(-walk, +walk);
			 p = min(1, rho(x_new,mu,sigma)/rho(x_old,mu,sigma));
			 if(p>rnd->Rannyu()){
				 x_old = x_new;
			 }
     dato += f(x_old,mu,sigma);	
		 }
		 somma_prog += dato/L;
     if(print) somma2_prog += pow(dato/L,2);
	 }
  	if(print){
		ofstream output("energiaSA.txt", ios::app);
		output << setw(15) << somma_prog/(double)N;
		output << setw(15) << sqrt((somma2_prog/(double)N - pow(somma_prog/(double)N,2))/(double)(N-1)) << endl;
		output.close();
	}

	if(move)*x_in=x_old;
	
	return somma_prog/(double)N;
}