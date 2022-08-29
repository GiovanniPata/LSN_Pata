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

	
  // ORA DICHIARO LE VARIABILI CHE MI SERVIRANNO NEI DUE CICLI
	
	 int M = 1000000; // numero estrazioni MC
	 int N = 100; // numero dati finali
	 int L = M/N; // lunghezza di blocco

	 double d = 2; // distanza righe
	 double l = 1; // lunghezza ago
	 int N_hit = 0; // lanci a buon fine

	 double x = 0; // in questa variabile piazzerò i numeri random tra 0 e d/2(=1)
   double	x_p = 0; // x del punto a caso nel cerchio (in 0,1)
	 double y_p = 0; // y del punto a caso nel cerchio (in -1,1)
	 double cos_theta = 0; // tramite x_p, y_p genero il coseno di un angolo estratto casualmente uniformemente
	
	 double dato = 0;
	 double sum_prog = 0;
	 double dato_2 = 0;
   double sum_2_prog = 0;  
	
	 ofstream output("pi.txt");


  // DOPPIO CICLO: 
	// INTERNO: L LANCI RNADOM, PRODUZIONE DI UN DATO PER PI
	// ESTERNO: MEDIO SU TUTTI GLI N DATI DI PI E METTO SU FILE IN MODO PROGRESSIVO

	 for(int i=0; i<N; i++){

		 N_hit = 0;

		 for(int j=0; j<L; j++){
			 x = rnd.Rannyu();
       x_p = rnd.Rannyu();
			 y_p = rnd.Rannyu(-1,1);
			 while(pow(x_p, 2) + pow(y_p, 2) > 1){
				 x_p = rnd.Rannyu();
			   y_p = rnd.Rannyu(-1,1);
			 } 
			 cos_theta = sqrt(1./(1 + pow( y_p/x_p, 2)));

			 if(x <= l/2 * cos_theta) N_hit ++;
		 }
		 dato = (2*l)/d * ((double) L)/ ((double)N_hit);
		 dato_2 = pow(dato, 2);
  	 sum_prog += dato;
		 sum_2_prog += dato_2;

	
		 output << setw(15) << sum_prog/(i+1);
		 if (i==0) output << setw(15) << 0;
		 else output << setw(15) << sqrt((sum_2_prog/(i+1) - pow(sum_prog/(i+1),2))/i);
		 if(i!=N-1) output << endl;
		 
	 }
	
	 output.close();
   rnd.SaveSeed();
   return 0;
}