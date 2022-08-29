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
#include <vector>
#include "punto.h"

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
	 int P = 100; // passi RW

	 Punto O; // punto(0,0,0)

	 vector<Punto> RW(P, O); // vettore di 100 punti (passi RW)
	 vector<double> r(P, 0); // moduli della distanza raggiunta

	 vector<double> dato(P, 0);
	 vector<double> dato_2(P, 0); // dato (varianza) e dato al quadrato, da usare nel ciclo esterno
   vector<double> sum_prog(P, 0);
	 vector<double> sum_2_prog(P, 0); //somma progr dato e dato al quadrato

	 double dev_std = 0;
	 double err_var = 0;
	 double err_std = 0;

	 ofstream output("lattice.txt"); 
	 ofstream output2("latticefinale.txt");
	 ofstream output3("continuo.txt");
	 ofstream output4("continuofinale.txt");
	 const char* filename1("latticeRW.txt");
	 const char* filename2("continuoRW.txt");
	// per cancellare file a ogni esecuzione:
	 ofstream out("latticeRW.txt");
	 out.close();
	 ofstream out2("continuoRW.txt");
	 out.close();
	 

	
  // 1) LATTICE
  // Idea: 3 cicli uno dentro l'altro
	// I. ciclo sui blocchi (i in range(N))
	//   II. ciclo sui throws nel blocco (j in range(L)) 
	//     III. ciclo sui passi del RW (k in range (P))

	 for(int i=0; i<N; i++){ 

		 dato.assign(100, 0); //risetto a zero

		 for(int j=0; j<L; j++){
			 for(int k=0; k<P; k++){ // questo ciclo è un singolo RW			   
         if(k==0) RW[k] = rnd.StepLattice();
			   else RW[k] = rnd.StepLattice() + RW[k-1];
				 r[k] = RW[k].Modulo();
				 dato[k] += (pow(r[k], 2));
				// stampo esempio RW
				 if(i==0 && j==0) RW[k].Print(filename1);
			 }
		 } // fine ciclo su j: ho simulato tutti i RW. Ora analisi dati

		 
		 for(int k=0; k<P; k++){ // ciclo per analisi dati di ogni passo

		   dato[k] = dato[k]/L;
			 dato_2[k] = pow(dato[k],2);
			 sum_prog[k] += dato[k];
			 sum_2_prog[k] += dato_2[k];

			// Stampo dati progressivi (in stile es 1_1) del 100° passo
			 if(k==P-1){
				 dev_std = sqrt(sum_prog[k]/(i+1));
				 output2 << setw(15) << dev_std;
				 if(i==0) output2 << setw(15) << 0;
				 else{
					 err_var = sqrt((sum_2_prog[k]/(i+1)-pow(sum_prog[k]/(i+1),2))/(i)); 
					 err_std = err_var/(2*dev_std); // propagazione errore
					 output2 << setw(15) << err_std;
				 }
				 output2 << endl;
			 }
		 
	   } // fine ciclo su k per analisi dati
	 
	 } // fine ciclo su i (blocchi)

	// Ultimo ciclo su passi RW per dev_std finale a ogni passo (equazione di diffusione)
	 for(int k=0; k<P; k++){
		 dev_std = sqrt(sum_prog[k]/N);
		 output << setw(15) << dev_std;
		 err_var = sqrt((sum_2_prog[k]/N-pow(sum_prog[k]/N,2))/(N-1)); 
		 err_std = err_var/(2*dev_std); // propagazione errore 
		 output << setw(15) << err_std;
		 output << endl;
	 }

	
  // 2) CONTINUO
	// Come prima, ma stavolta i punti vengono
	// estratti con distribuzione angolare continua

	 RW.assign(P, O); // rimetto tutto a zero
	 r.assign(P, 0); 
	 dato.assign(P, 0);
   dato.assign(P, 0);
	 dato_2.assign(P, 0); 
   sum_prog.assign(P, 0);
	 sum_2_prog.assign(P, 0); 


	 for(int i=0; i<N; i++){ 

		 dato.assign(100, 0);

		 for(int j=0; j<L; j++){
			 for(int k=0; k<P; k++){			 
			   if(k==0) RW[k] = rnd.StepContinuum();
			   else RW[k] = rnd.StepContinuum() + RW[k-1];
				 r[k] = RW[k].Modulo();
				 dato[k] += (pow(r[k], 2));
				// stampo esempio RW
				 if(i==0 && j==0) RW[k].Print(filename2);
		   }
		 }

		 for(int k=0; k<P; k++){

		   dato[k] = dato[k]/L;
			 dato_2[k] = pow(dato[k],2);
			 sum_prog[k] += dato[k];
			 sum_2_prog[k] += dato_2[k];

			 if(k==P-1){
				 dev_std = sqrt(sum_prog[k]/(i+1));
				 output4 << setw(15) << dev_std;
				 if(i==0) output4 << setw(15) << 0;
				 else{
					 err_var = sqrt((sum_2_prog[k]/(i+1)-pow(sum_prog[k]/(i+1),2))/(i)); 
					 err_std = err_var/(2*dev_std); // propagazione errore
					 output4 << setw(15) << err_std;
				 }
				 output4 << endl;
			 }
		
		 }
		 
	 }

	 for(int k=0; k<P; k++){
		 dev_std = sqrt(sum_prog[k]/N);
		 output3 << setw(15) << dev_std;
		 err_var = sqrt((sum_2_prog[k]/N-pow(sum_prog[k]/N,2))/(N-1)); 
		 err_std = err_var/(2*dev_std); // propagazione errore 
		 output3 << setw(15) << err_std;
		 output3 << endl;
	 }

	 output.close();
	 output2.close();
	 output3.close();
	 output4.close();
	
   rnd.SaveSeed();
   return 0;
}