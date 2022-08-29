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

	 double r = 0; // in questa variabile piazzerò i numeri random

	 double dato = 0;
	 double dato_2 = 0; // dato e dato al quadrato, da usare nel ciclo esterno
   double sum_prog = 0;
	 double sum_2_prog = 0; //somma progr dato e dato al quadrato
	
	 ofstream output("uniforme.txt"); 
	 ofstream output2("sampling.txt"); // accept reject
	 ofstream output3("bad_sampling.txt");
	 

  // 1) DISTRIBUZIONE UNIFORME 

	 for(int i=0; i<N; i++){

		 dato = 0;

		 for(int j=0; j<L; j++){
			 r = rnd.Rannyu();
			 r = M_PI/2*cos(M_PI/2*r);
       dato += r; 	 
		 }

		 dato = dato/L;
		 dato_2 = pow(dato, 2);
  	 sum_prog += dato;
		 sum_2_prog += dato_2;

		 output << setw(15) << sum_prog/(i+1);
		 if (i==0) output << setw(15) << 0;
		 else output << setw(15) << sqrt((sum_2_prog/(i+1) - pow(sum_prog/(i+1),2))/i);
		 if(i!=N-1) output << endl;	 
	 
	 }

	sum_prog = 0;
	sum_2_prog = 0;

	
  // 2) IMPORTANCE SAMPLING (ACCEPT-REJECT)

	 for(int i=0; i<N; i++){

	 dato = 0;

	 for(int j=0; j<L; j++){
		 r = rnd.Parabola();
		 r = (M_PI/2*cos(M_PI/2*r)) / (1.5*(1 - pow(r,2)));
     dato += r; 
	 }

	 dato = dato/L;
	 dato_2 = pow(dato, 2);
 	 sum_prog += dato;
	 sum_2_prog += dato_2;

	 output2 << setw(15) << sum_prog/(i+1);
	 if (i==0) output2 << setw(15) << 0;
	 else output2 << setw(15) << sqrt((sum_2_prog/(i+1) - pow(sum_prog/(i+1),2))/i);
		 if(i!=N-1) output2 << endl;
	 
   }

	 sum_prog = 0;
	 sum_2_prog = 0;


  // 3) BAD SAMPLING (p(x)=3x^2)

	 for(int i=0; i<N; i++){

	 dato = 0;

	 for(int j=0; j<L; j++){
		 r = rnd.Bad_sampling();
		 r = (M_PI/2*cos(M_PI/2*r)) / (3*r*r);
     dato += r; 
	 }

	 dato = dato/L;
	 dato_2 = pow(dato, 2);
 	 sum_prog += dato;
	 sum_2_prog += dato_2;

	 output3 << setw(15) << sum_prog/(i+1);
	 if (i==0) output3 << setw(15) << 0;
	 else output3 << setw(15) << sqrt((sum_2_prog/(i+1) - pow(sum_prog/(i+1),2))/i);
		 if(i!=N-1) output3 << endl;
	 
  }

	
	 output.close();
	 output2.close();
	 output3.close();
   rnd.SaveSeed();
   return 0;
}