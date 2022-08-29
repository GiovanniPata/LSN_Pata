/****************************************************************
*****************************************************************
		GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "punto.h"

using namespace std;

//costruttori

Punto :: Punto(){
	x = 0;
	y = 0;
  z = 0;
}

Punto :: Punto(double val_x, double val_y, double val_z){
	x = val_x;
	y = val_y;
	z = val_z;
}

// distruttore

Punto :: ~Punto(){}

// metodi set

void Punto :: SetX(double value){
  x = value;
  return;
}

void Punto :: SetY(double value){
  y = value;
  return;
}

void Punto :: SetZ(double value){
  z = value;
  return;
}

void Punto :: SetPoint(double val_x, double val_y, double val_z) {
	x = val_x;
	y = val_y;
	z = val_z;
	return;
}

// metodi get

double Punto :: GetX(void) const{
	return x;
}

double Punto :: GetY(void) const{
	return y;
}

double Punto :: GetZ(void) const{
	return z;
}

double Punto :: Modulo(void) const{
	return sqrt(x*x+y*y+z*z);
}

void Punto :: Print(const char* a) const{
	ofstream output(a, ios::app);
	output << x << " " << y << " " << z << endl;
	output.close();
	return;
}

// operatori

  Punto Punto :: operator+(const Punto& P) const{
	Punto sum(P);
	sum.x += x;
	sum.y += y;
	sum.z += z;
	return sum;
}
	
Punto & Punto :: operator+=(const Punto& P){
	return (*this = *this + P);
}