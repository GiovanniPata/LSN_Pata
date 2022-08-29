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
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}


double Random :: Exp(double lambda){
  return -1/lambda * log(1-Rannyu());
}


double Random :: Lorentian(double gamma, double mu){
	double r = Rannyu();
	while(r==0) r =Rannyu();
  return mu + gamma * tan((r-0.5)*M_PI);
}


int Random :: Dice(int faces){
  return (int) (1 + floor(Rannyu()*faces));
}

double Random :: Parabola(void){ // per esercizio
	double x = Rannyu();
	double y = Rannyu(0, 1.5);
	while(y > (1.5*(1-pow(x,2))) ){
	  x = Rannyu();
		y = Rannyu(0, 1.5);
	}
	return x;	
}

Punto Random :: StepLattice(void){
	double r = Dice(6);
	Punto p;
	if(r==1) p.SetX(1);
	else if(r==2) p.SetX(-1);
	else if(r==3) p.SetY(1);
	else if(r==4) p.SetY(-1);
	else if(r==5) p.SetZ(1);
	else p.SetZ(-1);
  return p;
}

Punto Random :: StepContinuum(void){
	double phi = Rannyu(0, 2*M_PI);
  double theta = acos(1-2*Rannyu());
  Punto p( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) );
  return p;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
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
