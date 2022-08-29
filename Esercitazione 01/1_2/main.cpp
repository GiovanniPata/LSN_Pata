/****************************************************************
*****************************************************************
		GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"

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


	// DEFINISCO LE VARIABILI COMUNI AL PROBLEMA. PER OGNUNA
	// DELLE 4 DISTRIBUZIONI ESTRAGGO PER M VOLTE:
	// I. UN NUMERO
	// II. MEDIA FRA 2 NUMERI
	// III. MEDIA FRA 10 NUMERI
	// IV. MEDIA FRA 100 NUMERI

	ofstream output1("uniform.txt");
	ofstream output2("exp.txt");
	ofstream output3("lorentian.txt");
	ofstream output4("dice.txt");
	// Creerò quattro file (uno per distribuzione). Ognuno conterrà
	// 10000x4 dati

  int M = 10000; // numero ripetizioni

	int * N = new int[4]; // vettore che mi dice tra quanti dati mediare
	N[0] = 1;
	N[1] = 2;
	N[2] = 10;
	N[3] = 100;

	double sum1 = 0; // distr uniforme
	double r1 = 0; 
	
	double sum2 = 0; // distr exp
	double r2 = 0;
	double lambda = 1;
	
	double sum3 = 0; // distr lorentziana
	double r3 = 0;
	double mu = 0;
	double gamma = 1;

	double sum4 = 0; // distr dado
	double r4 = 0;
	int faces = 6;

	// APRO 3 CICLI FOR UNO DENTRO L'ALTRO:
	// 1) SULLE 4 POSSIBILITA' DI NUMERO DI PUNTI NELLA MEDIA
	// 2) SULLE M RIPETIZIONI
	// 3) SULLE N[i] MISURE PER RIPETIZIONE

	for(int i=0; i<4; i++){

		for(int j=0; j<M; j++){

			for(int k=0; k<N[i]; k++){
				r1 = rnd.Rannyu();
				r2 = rnd.Exp(lambda);
				r3 = rnd.Lorentian(gamma, mu);
				r4 = rnd.Dice(faces);
				sum1 += r1;
				sum2 += r2;
				sum3 += r3;
				sum4 += r4;
    	}

			output1 << sum1/N[i] << endl;
			output2 << sum2/N[i] << endl;
			output3 << sum3/N[i] << endl;
			output4 << sum4/N[i] << endl;

		  sum1 = 0;
			sum2 = 0;
			sum3 = 0;
			sum4 = 0;

		}
  }
   
	delete [] N;
  output1.close();
	output2.close();
	output3.close();
	output4.close();
	
  rnd.SaveSeed();
  return 0;
}