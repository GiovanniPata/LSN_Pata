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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

	ReadInput >> restart;
	ReadInput >> initial; //0=random 1=aligned
	
	ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

	if(metro==1) cout << "Algorithm: Metropolis" << endl;
	else cout << "Algorithm: Gibbs" << endl;

	if(restart==1) cout << "Initial config: from previous simulation" << endl;
	else if(initial==0) cout << "Initial config: random (T=inf)" << endl;
	else cout << "Initial config: spin aligned (T=0)" << endl;
	
  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
	ic2 = 4;
	ix2 = 5;
 
  n_props = 6; //Number of observables

//initial configuration
	ifstream ReadConf;
  if(restart==1){
		ReadConf.open("config.final");
		for(int i=0; i<nspin; i++) ReadConf >> s[i];
		ReadConf.close();
	}
	else if(initial==0){
		for (int i=0; i<nspin; ++i){
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
    }
	}
	else{
		if(h<0) for(int i=0; i<nspin; i++) s[i]=-1;
		else for(int i=0; i<nspin; i++) s[i] = 1;
	}
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy (per spin) = " << walker[iu]/(double)nspin << endl;
	cout << "Heat capacity = " << beta*beta*(walker[ic]-pow(walker[iu],2)) << endl;
	cout << "Magnetization = " << beta*walker[im] << endl;
	cout << "Magnetic susceptibility = " << walker[ix] << endl;

	for(int i=0; i<50; i++){
		cout << s[i] << endl;
	}
	 cout << endl << "----------------------------" << endl << endl;
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;
	double delta_e;
	double A;

  for(int i=0; i<50; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
			delta_e = -2*Boltzmann(s[o], o); // variazione enrgia dopo spinflip
			A = Min(1, exp(-beta*delta_e)); // prob accettazione

			if(rnd.Rannyu() < A){
				s[o] = - s[o];
				accepted = accepted + 1.0;
			}

			attempted = attempted + 1.0;
			
    }
    else //Gibbs sampling
    {
      delta_e = -2*Boltzmann(s[o], o);
			A = 1./(1.+exp(beta*delta_e)); // prob scambio
			if(rnd.Rannyu() < A) s[o] = -s[o];
			accepted = accepted + 1.0;
			attempted = attempted + 1.0;
    }
/*
OSS: SEMBRA SBAGLIATO PERCHE' SIMILE AL METROPOLIS NORMALE, MA NON E' COSI': IL FATTO E' CHE LA DISTRIBUZIONE CAMPIONATA
DAL GIBBS E' DISCRETA E TRA SOLI DUE VALORI, QUINDI PER QUESTA
COINCIDENZA SEMBRANO UGUALI. NON LO SONO!!!
DIFFERENZA:
1) METROPOLIS: propongo mossa e vedo se accettare
2) GIBBS: fissato il resto della configurazione, determino
          l'ultimo spin. Di fatto, anche qui, o rimango uguale,
          o spinflippo. Ma se la distribuzione fosse discreta 
          a più uscite o continua apprezzerei la differenza
*/
  }
}

double Min(double a, double b){
	if(a<b) return a;
	else return b;
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = m;
	walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props-2; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props-2; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props-2; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi, Heat2, Chi2;
   const int wd=15;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    if(iblk==40) Ene.open("output.ene.dat",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if(iblk==40) Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    if(iblk==40) Ene.close();

	  if(iblk==40) Heat.open("output.heat.dat",ios::app);
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(stima_u*nspin, 2) )/ (double)nspin; //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    if(iblk==40) Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    if(iblk==40) Heat.close();

	  Mag.open("output.mag.dat",ios::app);
    stima_m = blk_av[im]/blk_norm/ (double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    if(iblk==40) Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    if(iblk==40) Mag.close();

	  if(iblk==40) Chi.open("output.chi.dat",ios::app);
    stima_x = beta*blk_av[ix]/blk_norm/ (double)nspin; //Magnetic susceptibility
    glob_av[ix] += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    if(iblk==40) Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    if(iblk==40) Chi.close();

	  if(iblk==40) Heat2.open("output.2heat.dat",ios::app);
    stima2_c = beta*beta*(blk_av[ic]/blk_norm - pow(nspin*U_teo(J,temp,nspin), 2) )/ (double)nspin; //Heat capacity
    glob_av[ic2]  += stima2_c;
    glob_av2[ic2] += stima2_c*stima2_c;
    err2_c=Error(glob_av[ic2],glob_av2[ic2],iblk);
    if(iblk==40) Heat2 << setw(wd) << iblk <<  setw(wd) << stima2_c << setw(wd) << glob_av[ic2]/(double)iblk << setw(wd) << err2_c << endl;
    if(iblk==40) Heat2.close();

	  if(iblk==40) Chi2.open("output.2chi.dat",ios::app);
    stima2_x = beta*(blk_av[ix]/blk_norm-pow(stima_m*nspin,2))/ (double)nspin; //Magnetic susceptibility
    glob_av[ix2] += stima2_x;
    glob_av2[ix2] += stima2_x*stima2_x;
    err2_x=Error(glob_av[ix2],glob_av2[ix2],iblk);
    if(iblk==40) Chi2 << setw(wd) << iblk <<  setw(wd) << stima2_x << setw(wd) << glob_av[ix2]/(double)iblk << setw(wd) << err2_x << endl;
    if(iblk==40) Chi2.close();
	
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double U_teo(double J, double T, double N){
	double b=1./T;
	double th = tanh(b*J);
	return -J*th*((1+pow(th, N-2))/(1+pow(th, N)));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
