#ifndef __Punto__
#define __Punto__

class Punto {

protected:
  double x,y,z;


public:

  // constructors
  Punto();
  Punto(double, double, double);

  // destructor
  ~Punto();

  // methods
  void SetX(double);
  void SetY(double);
  void SetZ(double);
  void SetPoint(double, double, double);

  double GetX(void) const;
  double GetY(void) const;
  double GetZ(void) const;
  double Modulo(void) const;
  void Print(const char*) const;

  Punto operator+(const Punto&) const;
	Punto & operator+=(const Punto&);

};



#endif