/****************************************************************
*****************************************************************
GIOVANNI PATA 944333
*****************************************************************
*****************************************************************/


#ifndef __Punto__
#define __Punto__

class Punto {

protected:
  double x,y;


public:

  // constructors
  Punto();
  Punto(double, double);

  // destructor
  ~Punto();

  // methods
  void SetX(double);
  void SetY(double);
  void SetPoint(double, double);

  double GetX(void) const;
  double GetY(void) const;
  double Modulo(void) const;
  double Modulo2(void) const;
  void Print(void) const;

  Punto operator+(const Punto&) const;
  Punto operator-(const Punto&) const;
	Punto & operator+=(const Punto&);

};



#endif