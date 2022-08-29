/****************************************************************
*****************************************************************
GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "punto.h"

using namespace std;

//costruttori

Punto :: Punto(){
	x = 0;
	y = 0;
}

Punto :: Punto(double val_x, double val_y){
	x = val_x;
	y = val_y;
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

void Punto :: SetPoint(double val_x, double val_y) {
	x = val_x;
	y = val_y;
	return;
}

// metodi get

double Punto :: GetX(void) const{
	return x;
}

double Punto :: GetY(void) const{
	return y;
}


double Punto :: Modulo(void) const{
	return sqrt(x*x+y*y);
}

double Punto :: Modulo2(void) const{
	return (x*x+y*y);
}

void Punto :: Print(void) const{
	cout << x << " " << y << endl;
	return;
}

// operatori

Punto Punto :: operator+(const Punto& P) const{
	Punto sum(P);
	sum.x += x;
	sum.y += y;
	return sum;
}

Punto Punto :: operator-(const Punto& P) const{
	Punto sum(P);
	sum.x = sum.x-x;
	sum.y = sum.y-y;
	return sum;
}  
	
Punto & Punto :: operator+=(const Punto& P){
	return (*this = *this + P);
}


