/****************************************************************
*****************************************************************
GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "random.h"
#include <cmath>
#include "punto.h"
#include <cassert>
#include <string>

using namespace std;


// DICHIARAZIONE DELLE VARIABILI GLOBALI

 const int G = 10000; //generazioni
 const int M = 5000; //numero cromosomi
 const int n = 34; //cities

 int seed[4]; // per random gen
 Random rnd;


// DICHIARAZIONE DI FUNZIONI

// prepara randomgen genera vettore di città e prima generazione casuale
 void Initialize(bool mode, Punto* city, int** pop); 

// permuta vettore di n interi lasciando invariato il primo elemento
 void Permuta(int* v); 

// dà la lunghezza della strada per un cromosoma
 double Lenght(Punto* city, int* v);

// calcola il fitness di una popolazione e registra su fit e sum_fit
 void Fitness(Punto* city, int** pop, double* fit, double* sum_fit);

// ordina vettore fit in fit_ordered
 void Order(double* fit, double* fit_ordered);

// dà la media dei primi M/2 elementi di fit_ordered
 double Ave (double* fit_ordered);

// seleziona la label di un genitore
 int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered);

// attua crossing over e posiziona in pop_new
 void Crossover(int mum, int dad, int label, int** pop, int** pop_new);

// muta in pop_new
 void Mutation(int label, int** pop_new);

// controlla se la popolazione è ben definita
 void Check(int** pop);

// stampa il resoconto della generazione
void Report(int g, double* fit_ordered);

// dà la label del miglior percorso in una generazione
 int BestPoint(double* fit);



int main (int argc, char *argv[]){

	if(argc!=2){
		cerr << "Usage: ./main.exe mode" << endl;
		cerr << "mode=0: Circonferenza R=1" << endl;
		cerr << "mode=1: Quadrato l=1" << endl;
		return -1;
	}
	
  // ALLOCAZIONE DINAMICA DELLE VARIABILI

	 Punto* city = new Punto[n]; // vettore di città
	 bool mode = atoi(argv[1]);
	
	 int** pop = new int*[M]; // popolazione
	 for(int i=0; i<M; i++) pop[i] = new int[n];
	 double* fit = new double[M]; 
	 double* sum_fit = new double;

	 int** pop_new = new int*[M]; // per registare la nuova pop
	 for(int i=0; i<M; i++) pop_new[i] = new int[n];

	 double* fit_ordered = new double[M];
	 int mum, dad;

	 ofstream output("besties.txt");
	 ofstream output2("path.txt");
	
	 
  // INIZIALIZZO E OPERO SULLA GENERAZIONE ZERO

	 Initialize(mode, city, pop);
	 Check(pop);
	 Fitness(city, pop, fit, sum_fit);
	 Order(fit, fit_ordered);
	 output << setw(15) << "0" <<  setw(15) << fit_ordered[0] << setw(15) << Ave(fit_ordered) << endl;
	 Report(0, fit_ordered);


	// CICLO SU TUTTE LE GENERAZIONI
	
	 for(int g=1; g<G; g++){

		// CICLO SU TUTTI I CROMOSOMI
		 for(int i=0; i<M/2; i++){ 
		   mum = Selection(pop, fit, sum_fit, fit_ordered);
		   dad = Selection(pop, fit, sum_fit, fit_ordered);
       Crossover(mum, dad, i, pop, pop_new);
			 Mutation(i, pop_new);
     }

		// CONTROLLO POPOLAZIONE BUONA
		 Check(pop_new);
		// SOSTITUISCO LA POP NUOVA ALLA VECCHIA
		 for(int i=0; i<M; i++){
			 delete[] pop[i];
			 pop[i] = pop_new[i];
		 }
		 delete [] pop;
		 pop = pop_new;
		 // RI-CREO SPAZIO PER LA POP NUOVA
		 pop_new = new int*[M]; // per registare la nuova pop
	   for(int i=0; i<M; i++) pop_new[i] = new int[n];

	  // NUOVO CALCOLO FITNESS
		 Fitness(city, pop, fit, sum_fit);
	   Order(fit, fit_ordered);
	   output << setw(15) << g << setw(15) << fit_ordered[0] << setw(15) << Ave(fit_ordered) << endl;
		 if(g<10 or g%50==0) Report(g, fit_ordered);
	 }
	
	// REGISTRO IL MIGLIOR PERCORSO PER L'ULTIMA GENERZIONE 
   int best = BestPoint(fit);
	 int label = 0;
   for(int i=0; i<n; i++){
		 label = pop[best][i]-1;
		 output2 << setw(15) << city[label].GetX();
		 output2 << setw(15) << city[label].GetY();
		 output2 << endl;
	 }
	 output2 << setw(15) << city[pop[best][0]-1].GetX();
	 output2 << setw(15) << city[pop[best][0]-1].GetY();
	 output2 << endl;
	
  // PULIZIA FINALE
	 for(int i=0; i<M; i++) delete[] pop[i];
	 delete[] pop;
	 delete[] fit;
	 delete[] fit_ordered;
	 delete[] city;
   output.close();
	 output2.close();

	
   return 0;

}

//_______________________________________________________________
// DEFINIZIONE DELLE FUNZIONI UTILIZZATE
//_______________________________________________________________


void Initialize(bool mode, Punto* city, int** pop){
	// Preparo random gen
	 int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
 	 Primes.close();
	 ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
	// Genero vettore di città
	 double R=1;
	 double theta;
	 for(int i=0; i<n; i++){
	   if(mode){ 
		   city[i].SetPoint(rnd.Rannyu(0,R), rnd.Rannyu(0,R));
	   }
	   else{
		   theta = rnd.Rannyu(0,2*M_PI);
		   city[i].SetPoint(R*cos(theta), R*sin(theta));
	   }
   }

	/*ifstream in;
	in.open("American_capitals.dat", ios::in);
	string s;
	double x;
	double y;
	in>>s>>s>>s>>s;
	for(int i=0; i<50; i++){
		in>>s>>s>>x>>y;
		city[i].SetPoint(x,y);
		cout << i<< endl;
	}
	in.close();
  */
	// Preparo randomly la generazione zero
   for(int i=0; i<M; i++){
		 for(int j=0; j<n; j++) pop[i][j] = j+1;
	   Permuta(pop[i]);
	 }
}


void Permuta(int* v){
	 for(int i=1; i<n; i++){
	 int label = rnd.Dice(n-1); // numero da 1 a n-1
	 int appo = v[i];
	 v[i] = v[label];
	 v[label] = appo;
	}
}


double Lenght(Punto* city, int* v){
	 double l=0;
	 for(int i=0; i<n-1; i++)  l += (city[v[i]-1]-city[v[i+1]-1]).Modulo();
	 l += (city[v[0]-1]-city[v[n-1]-1]).Modulo();
	 return l;
}


void Fitness(Punto* city, int** pop, double* fit, double* sum_fit){
	 *sum_fit = 0;
	 for(int i=0; i<M; i++){
		 fit[i]=Lenght(city, pop[i]);
		 *sum_fit += pow(1./fit[i],3);
	 }
}


void Order(double* fit, double* fit_ordered){
	 for(int i=0; i<M; i++) fit_ordered[i]=fit[i]; //copio fit in fit_ordered
	 int pos=0;
	 double appo=0;
   for(int i=0; i<M-1; i++){
		 pos=i;
		 for(int j=i+1; j<M; j++)
		 	 if(fit_ordered[j]<fit_ordered[pos]) pos=j;
	   appo = fit_ordered[i];
		 fit_ordered[i] = fit_ordered[pos];
		 fit_ordered[pos] = appo;
	 }
}


double Ave (double* fit_ordered){
	 double accu = 0;
	 for(int i=0; i<M/2; i++) accu+=fit_ordered[i];
	 return accu/(M/2);
}


int Selection(int** pop, double* fit, double* sum_fit, double* fit_ordered){
	/* double F = *sum_fit;
	 double accu = (pow(1./fit[0],3))/F;
	 int conta = 0;
	 double r = rnd.Rannyu();
	 while(accu<r){
		 conta++;
		 accu += (pow(1./fit[conta],3))/F;
	 }
	 if(conta<0 or conta>M-1){
		 cout << "Errore: conta= " << conta << " in 'Selection': setto conta=M-1" << endl;
		 conta = M-1;
	 }*/
  double r=rnd.Rannyu(); //IMPLEMENTAZIONE ALTERNATIVA
	int j = floor(M*pow(r,5));
	double appo=fit_ordered[j];
	int conta=0;
	while(appo!=fit[conta]){
		conta++;
	} 
	return conta;	
}



void Crossover(int mum, int dad, int label, int** pop, int** pop_new){
	 if(rnd.Rannyu()<0.65){ // avviene crossover
		 int pos = rnd.Dice(n); //num da 1 a n
		 //cout << "taglio in posizione " << pos << endl;
		 //cout << pos << " numeri buoni a sx" << endl;
		 //cout << n-pos << " numeri mancanti a dx" << endl;
		 int pos_copy = pos;
		 bool check = 0;
		 for(int i=0; i<pos; i++){ //ricopio a sx del taglio
			 pop_new[2*label][i] = pop[mum][i];
			 pop_new[2*label+1][i] = pop[dad][i];
		 }

		 //cout << "Inizio setaccio: " << endl;
		 for(int i=1; i<n; i++){ //riempio a dx il primo figlio
			 check=0;
			 //cout << i << ") ";
			 for(int j=1; j<pos; j++){
				 if(pop[dad][i]==pop_new[2*label][j]){
					 check=1;
					 //cout << "Il  " << pop[dad][i] << " c'è" << endl;
				 }
			 }	 
			 if(!check){
				 pop_new[2*label][pos] = pop[dad][i];
				// cout << "manca il " << pop[dad][i] << " in pos " << pos << endl;
				 pos++;
				 //cout <<"incremento: pos=" << pos << endl;	 
			 }
		 }
		 //cerr << "fin qui OK";
		 //cerr << "Inizio setaccio: " << endl;
		 for(int i=1; i<n; i++){ // riempio a dx il secondo figlio
			 check=0;
			 //cerr << i << ") ";
			 for(int j=1; j<pos_copy; j++){
				 if(pop[mum][i]==pop_new[2*label+1][j]){
					 check=1;
					 //cerr << "Il  " << pop[mum][i] << " c'è" << endl;
				 }
			 }
			 if(!check){
				 pop_new[2*label+1][pos_copy] = pop[mum][i];
				 //cout << "manca il " << pop[mum][i] << " in pos " << pos_copy << endl;
				 pos_copy++;
				 //cout <<"incremento: pos_copy=" << pos_copy << endl;	
			 }
		 }
		 if(pos!=n or pos_copy!=n) cout << "Errore in crossover" << endl;
	 }
	 else{ // non avviene crossover
		 for(int i=0; i<n; i++){
			 pop_new[2*label][i] = pop[2*label][i];
			 pop_new[2*label+1][i] = pop[2*label+1][i];
		 }
	 }
 }

	
void Mutation(int label, int** pop_new){
	//prima mutazione: single swap
	 if(rnd.Rannyu()<0.1){
		 for(int i=0; i<5; i++){
			 int a = rnd.Dice(n-1); //pos 1
		   int b = rnd.Dice(n-1); // pos 2
		   swap( pop_new[2*label][a], pop_new[2*label][b] );
		 }
	 }
	 if(rnd.Rannyu()<0.1){
		 for(int i=0; i<5; i++){
			 int a = rnd.Dice(n-1); //pos 1
		   int b = rnd.Dice(n-1); // pos 2
		   swap( pop_new[2*label+1][a], pop_new[2*label+1][b] );
		 }
	 }
	//seconda mutazione: shift di m celle contigue
	 if(rnd.Rannyu()<0.1){
		 int m = rnd.Dice(n-1); // numero celle contigue
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<n; i++){
			 if(i+m<n) v[i]=pop_new[2*label][i+m];
			 else v[i]=pop_new[2*label][i+m+1-n];
		 }
		 for(int i=0; i<n; i++) pop_new[2*label][i]=v[i];
		 delete[] v;
		 }
	 if(rnd.Rannyu()<0.1){
		 int m = rnd.Dice(n-1); // numero celle contigue
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<n; i++){
			 if(i+m<n) v[i]=pop_new[2*label+1][i+m];
			 else v[i]=pop_new[2*label+1][i+m+1-n];
		 }
		 for(int i=0; i<n; i++) pop_new[2*label+1][i]=v[i];
		 delete[] v;
		 }
		// terza mutazione: permutazione celle contigue con altre
	 if(rnd.Rannyu()<0.1){
		 int start = rnd.Dice(n-11);
		 int m = rnd.Dice(5);
		 int start2 = start + m + rnd.Dice(n-m-start-1-6);
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<start; i++) v[i]=pop_new[2*label][i];
		 for(int i=start; i<start+m; i++) v[i-start+start2]=pop_new[2*label][i];
		 for(int i=start+m; i<start2; i++) v[i]=pop_new[2*label][i];
		 for(int i=start2; i<start2+m; i++) v[i-start2+start]=pop_new[2*label][i];
		 for(int i=start2+m; i<n; i++) v[i]=pop_new[2*label][i];
		 for(int i=0; i<n; i++) pop_new[2*label][i]=v[i];
		 delete[] v;
	 }
   if(rnd.Rannyu()<0.1){
		 int start = rnd.Dice(n-11);
		 int m = rnd.Dice(5);
		 int start2 = start + m + rnd.Dice(n-m-start-1-6);
		 int* v = new int[n];
		 v[0]=1;
		 for(int i=1; i<start; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=start; i<start+m; i++) v[i-start+start2]=pop_new[2*label+1][i];
		 for(int i=start+m; i<start2; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=start2; i<start2+m; i++) v[i-start2+start]=pop_new[2*label+1][i];
		 for(int i=start2+m; i<n; i++) v[i]=pop_new[2*label+1][i];
		 for(int i=0; i<n; i++) pop_new[2*label+1][i]=v[i];
		 delete[] v;
	 }
	// quarta mutazione: inversione
	 if(rnd.Rannyu()<0.1){
		 int m = rnd.Dice(n-1);
		 int start = rnd.Dice(n-1);
		 int* v = new int[m];
		 v[0]=1;
		 for(int i=0; i<m; i++){
			 if(start+m-1-i<n) v[i]=pop_new[2*label][start+m-i-1];
			 else v[i]=pop_new[2*label][start+m-1-i-n+1];
		 }
		 for(int i=0; i<m; i++){
			 if(start+i<n) pop_new[2*label][start+i]=v[i];
			 else pop_new[2*label][start+i-n+1]=v[i];
		 } 
		 delete[]v;
	 }
	if(rnd.Rannyu()<0.1){
		 int start = rnd.Dice(n-3);
		 int m = rnd.Dice(n-start-1);
		 int* v = new int[m];
		 v[0]=1;
		 for(int i=0; i<m; i++){
			 if(start+m-1-i<n) v[i]=pop_new[2*label+1][start+m-i-1];
			 else v[i]=pop_new[2*label+1][start+m-1-i-n+1];
		 }
		 for(int i=0; i<m; i++){
			 if(start+i<n) pop_new[2*label+1][start+i]=v[i];
			 else pop_new[2*label+1][start+i-n+1]=v[i];
		 }
	 }
}

	
void Check(int** pop){
	 bool check=0;
	 for(int i=0; i<M; i++){
		 if(pop[i][0]!=1){
			 cout << "Error in creating new gen: doesn't start with 1" << endl;
			 exit(1);
		 }
		 for(int num=2; num<n+1; num++){
			 check=0;
			 for(int j=1; j<n; j++){
				 if(pop[i][j]==num) check=1;
			 }
			 if(!check){
				 cout << "Error in creating new gen: absent " << num << endl;
				 exit(1);
			 }
		 }
	 }
} 


void Report(int g, double* fit_ordered){
	cout << endl;
	cout << "GENEREATION " << g << " OF "<< G-1 << endl;
	cout << "Best fit: " << fit_ordered[0] << endl;
	cout << "______________________________________" << endl;
}

	
 int BestPoint(double* fit){
	 int best = 0;
	 for(int i=1; i<M; i++){
		 if(fit[i]<fit[best]) best = i;
	 }
	 return best;
 }