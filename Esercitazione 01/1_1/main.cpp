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
	 double dato_2 = 0; // media e media al quadrato, prodotte dentro al ogni blocco
   double sum_prog = 0;
	 double sum_2_prog = 0; // somma progr media e media al quadrato

	 double dato_var = 0;
	 double dato_2_var = 0; // varianza e varianza al quadrato, prodotte dentro a ogni blocco
   double sum_prog_var = 0;
	 double sum_2_prog_var = 0; //somma progr varianza e varianza al quadrato
	
	 ofstream output("media.txt"); 
	 ofstream output2("varianza.txt");


  // DOPPIO CICLO: INTERNO GENERO RANDOM PER GENERAZIONE  BLOCCHI
	//               ESTERNO USO DATI E METTO SU FILE

	// OSS: ciò che c'è nel ciclo esterno lo posso riciclare
	//      per sempre. Ci sarà ogni volta da cambiare il 
	//      codice nel ciclo più interno.

	 for(int i=0; i<N; i++){

		 dato = 0;
		 dato_var = 0;

		 for(int j=0; j<L; j++){
			 r = rnd.Rannyu();
       dato += r; 
			 dato_var += pow((r-0.5), 2);
		 }

		 dato = dato/L;
		 dato_2 = pow(dato, 2);
  	 sum_prog += dato;
		 sum_2_prog += dato_2;

		 dato_var = dato_var/L;
		 dato_2_var = pow(dato_var, 2);
  	 sum_prog_var += dato_var;
		 sum_2_prog_var += dato_2_var;


		 output << setw(15) << sum_prog/(i+1);
		 if (i==0) output << setw(15) << 0;
		 else output << setw(15) << sqrt((sum_2_prog/(i+1) - pow(sum_prog/(i+1),2))/i);
		 if(i!=N-1) output << endl;

		 output2 << setw(15) << sum_prog_var/(i+1);
		 if (i==0) output2 << setw(15) << 0;
		 else output2 << setw(15) << sqrt((sum_2_prog_var/(i+1) - pow(sum_prog_var/(i+1),2))/i);
		 if(i!=N-1) output2 << endl;
		 
	 }


	// PARTE FINALE: TEST DEL CHI QUADRO

	 ofstream output3("chiquadro.txt");

	 int N_INT = 100; // numero intervalli in cui dividere (0,1)
	 int n = 10000; // numero estrazioni casuali per calcolo chi quadro
	 int* v = new int[N_INT]; // vettore che conta occorrenze in intervallo
   int RIP = 10000; // ripetizioni calcolo chi quadro
	 double chi=0;

	 for(int i=0; i<RIP; i++){ //ciclo esterno: calcolo RIP=10000 chiquadri

		 for(int j=0; j<N_INT; j++) v[j]=0;

		 for(int k=0; k<n; k++){ // ciclo interno1: conto le occorrenze nei vari bin
			 r = rnd.Rannyu();
			 v[(int)floor(r*N_INT)]++;
		 }

		 for(int j=0; j<N_INT; j++) // ciclo interno2: calcolo chiquadro
		   chi += (pow((v[j]-n/N_INT),2))/(n/N_INT);
		 
		 output3 << chi << endl;
		 chi = 0;
	 }
	
	
   delete [] v;
	 output.close();
	 output2.close();
	 output3.close();
   rnd.SaveSeed();
   return 0;
}
