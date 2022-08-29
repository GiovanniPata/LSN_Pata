/****************************************************************
*****************************************************************
		GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include <iomanip>

using namespace std;

double max(double, double);
 
int main (int argc, char *argv[]){

	
	//PREPARAZIONE GENERATORE NUMERI CASUALI
	
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

	
  // ORA DICHIARO LE VARIABILI CHE MI SERVIRANNO NEI DUE CICLI
	
	 int M = 100000; // numero estrazioni MC
	 int N = 100; // numero dati finali
	 int L = M/N; // lunghezza di blocco

	 double T = 1; // delivery time
	 double K = 100; // strike prize
	 double mu = 0.1; // risk-free interest rate
	 double sigma = 0.25; // volatility
	 double S_0 = 100; // asset price

	 double r = 0; // in questa variabile piazzerò i numeri random
	
	 double dato_c = 0;
	 double dato_c_2 = 0; // dato e dato al quadrato, da usare nel ciclo esterno
	 double dato_p = 0;
	 double dato_p_2 = 0;
	
   double sum_c_prog = 0;
	 double sum_c_2_prog = 0; //somma progr dato e dato al quadrato
	 double sum_p_prog = 0;
	 double sum_p_2_prog = 0;
	
	 ofstream output("direct_call.txt"); 
	 ofstream output2("direct_put.txt");
   ofstream output3("discrete_call.txt");
	 ofstream output4("discrete_put.txt");
	 

  // 1) DIRECT SAMPLING

	 for(int i=0; i<N; i++){

		 dato_c = 0;
		 dato_p = 0;

		 for(int j=0; j<L; j++){
			 r = S_0 * exp((mu-0.5*sigma*sigma)*T + sigma*rnd.Gauss(0,T));
			 dato_c += exp(-mu*T) * max(0, r-K);
			 dato_p += exp(-mu*T) * max(0, K-r);	 
		 }

		 dato_c = dato_c/L;
		 dato_c_2 = pow(dato_c, 2);
  	 sum_c_prog += dato_c;
		 sum_c_2_prog += dato_c_2;
 		 
 		 dato_p = dato_p/L;
		 dato_p_2 = pow(dato_p, 2);
  	 sum_p_prog += dato_p;
		 sum_p_2_prog += dato_p_2;


		 output << setw(15) << sum_c_prog/(i+1);
		 if (i==0) output << setw(15) << 0;
		 else output << setw(15) << sqrt((sum_c_2_prog/(i+1) - pow(sum_c_prog/(i+1),2))/i);
		 if(i!=N-1) output << endl;

		 output2 << setw(15) << sum_p_prog/(i+1);
		 if (i==0) output2 << setw(15) << 0;
		 else output2 << setw(15) << sqrt((sum_p_2_prog/(i+1) - pow(sum_p_prog/(i+1),2))/i);
		 if(i!=N-1) output2 << endl;		 
	 
	 }

	sum_c_prog = 0;
	sum_c_2_prog = 0;
	sum_p_prog = 0;
	sum_p_2_prog = 0;
	

  // 2) DISCRETE SAMPLING

	 int S = 100; // steps
	 double t = (double) T/S; // lunghezza step
	
	 for(int i=0; i<N; i++){

		 dato_c = 0;
		 dato_p = 0;

		 for(int j=0; j<L; j++){

			 r = S_0;

			 for(int k=0; k<S; k++)// devo inserire ciclo interno sugli S steps
				 r = r * exp((mu-0.5*sigma*sigma)*t + sigma*sqrt(t)*rnd.Gauss(0,1));
			 
			 dato_c += exp(-mu*T) * max(0, r-K);
			 dato_p += exp(-mu*T) * max(0, K-r);
    
		 }

		 dato_c = dato_c/L;
		 dato_c_2 = pow(dato_c, 2);
  	 sum_c_prog += dato_c;
		 sum_c_2_prog += dato_c_2;
 		 
 		 dato_p = dato_p/L;
		 dato_p_2 = pow(dato_p, 2);
  	 sum_p_prog += dato_p;
		 sum_p_2_prog += dato_p_2;


		 output3 << setw(15) << sum_c_prog/(i+1);
		 if (i==0) output3 << setw(15) << 0;
		 else output3 << setw(15) << sqrt((sum_c_2_prog/(i+1) - pow(sum_c_prog/(i+1),2))/i);
		 if(i!=N-1) output3 << endl;

		 output4<< setw(15) << sum_p_prog/(i+1);
		 if (i==0) output4 << setw(15) << 0;
		 else output4 << setw(15) << sqrt((sum_p_2_prog/(i+1) - pow(sum_p_prog/(i+1),2))/i);
		 if(i!=N-1) output4 << endl;
	 
	 }

	
	 output.close();
	 output2.close();
	 output3.close();
   output4.close();
   rnd.SaveSeed();
	
   return 0;
	
}



double max(double a, double b){
	if(a>b) return a;
	else return b;
}