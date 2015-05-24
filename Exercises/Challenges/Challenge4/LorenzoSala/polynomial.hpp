#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <vector>
#include <utility>
#include <iostream>

class Polynomial{

   private:
	std::vector<double> coeffs;
	unsigned int max_degree;
	inline void P_normalize();

   public:
	// Costruttore
	Polynomial(std::vector<double> _coeffs, unsigned int _max_degree);
	Polynomial(std::vector<double> _coeffs);
	Polynomial();

	
	// Costruttore di copia (default)
	Polynomial(Polynomial const &) = default;
	
	// Distruttore di default
	
	// Metodi per le operazioni
	Polynomial &operator +=(Polynomial const &);
	Polynomial &operator -=(Polynomial const &);
	Polynomial &operator /=(Polynomial const &);
  	Polynomial &operator *=(Polynomial const &);
  	
  	
	// Funzioni per le operazioni
	friend Polynomial operator+(Polynomial const &,Polynomial const &);
  	friend Polynomial operator-(Polynomial const &,Polynomial const &);
  	friend Polynomial operator*(Polynomial const &,Polynomial const &);
  	friend std::pair<Polynomial,Polynomial> operator/(Polynomial const &,Polynomial const &);
	
	// Input/Output streaming operators
	friend std::ostream & operator << (std::ostream &, Polynomial const &);
	friend std::istream & operator >> (std::istream &, Polynomial &);	

	// Valutare un polinomio in un punto
	double operator()(double const &);
	// Calcola la derivata di un polinomio
	friend Polynomial der_pol(Polynomial const & p);

	// Solutore per le radici Newton-Horner
	std::vector<double> solve(double x0, double toll = 1e-8, unsigned const MAX_ITER=1000);
	// x0 punto iniziale, toll tollerenza e MAX_ITER numero massimo di iterazioni
};



#endif
