#ifndef __POLYNOMIAL_HPP__
#define __POLYNOMIAL_HPP__

#include <vector>
#include <utility>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <complex>
#include <cmath>


template <unsigned D=0, typename T = double>
class Polynomial{

   private:
	std::vector<T> coeffs;
	inline void P_normalize();

   public:
	// Costruttore
	Polynomial<D,T>(std::vector<T> _coeffs){coeffs.resize(D+1); coeffs = _coeffs;this->P_normalize();}
	Polynomial<D,T>(){coeffs.resize(D+1);coeffs[0] = 0;}

	
	// Costruttore di copia (default)
	Polynomial<D,T>(Polynomial<D,T> const & p) = default;

	// Distruttore di default
	
	// Metodi per le operazioni
	template <unsigned U>
	Polynomial<D,T> &operator +=(Polynomial<U,T> const & p){*this = *this + p;return *this;}
	template <unsigned U>
	Polynomial<D,T> &operator -=(Polynomial<U,T> const & p){*this = *this - p;return *this;}
	template <unsigned U>
  	Polynomial<D,T> &operator *=(Polynomial<U,T> const & p){*this = *this * p;return *this;}
	template <unsigned U>
	Polynomial<D,T> &operator /=(Polynomial<U,T> const & p){auto result = *this / p;*this = result.first;return *this;}
	// Ricordare che l'operazione di divisione restituisce quoziente e resto, in questo caso verrà restituito solo il quoziente.
  	
  	
	// Funzioni per le operazioni
	template <unsigned S,unsigned U,typename W>
	friend Polynomial<S,W> operator+(Polynomial<S,W> const &,Polynomial<U,W> const &);
  	
	template <unsigned S,unsigned U,typename W>
	friend Polynomial<S,W> operator-(Polynomial<S,W> const &,Polynomial<U,W> const &);

	template <unsigned S,unsigned U,typename W>
  	friend Polynomial<S+U,W> operator*(Polynomial<S,W> const &,Polynomial<U,W> const &);
	
	template <unsigned S,unsigned U, typename W>
  	friend std::pair<Polynomial<S,W>,Polynomial<S,W>> operator/(Polynomial<S,W> const &,Polynomial<U,W> const &);
	
	// Input/Output streaming operators
	template <unsigned int U,typename S>
	friend std::ostream & operator << (std::ostream &, Polynomial<U,S> const &);
	template <unsigned int U,typename S>
	friend std::ostream & operator << (std::ostream &, Polynomial<U,std::complex<S>> const &);

	
	template <unsigned int U,typename S>
	friend std::istream & operator >> (std::istream &, Polynomial<U,S> &);	

	// ALTRE FUNZIONI:
	// Valutare un polinomio in un punto
	T operator()(T const &);

	// Calcola la derivata di un polinomio
	template <unsigned U,typename W>
	friend Polynomial<U-1,W> der_pol(Polynomial<U,W> const & p);

	// Solutore per le radici Newton-Horner
	template <unsigned int U,typename S>
	friend std::vector<S> solve(Polynomial<U,S> const & p, double x0, const double toll, const unsigned MAX_ITER);
	// x0 punto iniziale, toll tollerenza e MAX_ITER numero massimo di iterazioni
	template <unsigned U>
	friend std::vector<double> solve(Polynomial<U,double> const & p, double x0, const double toll, const unsigned MAX_ITER);
	template <unsigned U>
	friend std::vector<std::complex<double>> solve(Polynomial<U,std::complex<double>> const & p, std::complex<double> x0,const double toll, const unsigned MAX_ITER);

};
// NORMALIZE
template <unsigned D,typename T>
void Polynomial<D,T>::P_normalize(){
	while(coeffs.size()-1 > 0){
		// if(fabs(this->coeffs[coeffs.size()-1]) < 1e-100){
		if(this->coeffs[coeffs.size()-1] == T(0)){
			this->coeffs.erase(coeffs.end()-1);
		}else{
			break;
		}	
		
	}
	return;
}

// Operations:

template <unsigned S, unsigned U,typename T>
Polynomial<S,T> operator+(Polynomial<S,T> const & pl,Polynomial<U,T> const & pr){
	Polynomial<S,T> tmp = pl;	
	if(S<U){
		throw std::runtime_error("This sum is not possible, apply the commutative property because the first addend has to be greater than the second");
		return tmp;
	}
	for (unsigned int i=0 ; i < pr.coeffs.size();++i)
		tmp.coeffs[i] += pr.coeffs[i];
	tmp.P_normalize();	
	return tmp;
}	

template <unsigned S,unsigned U,typename T>
Polynomial<S,T> operator-(Polynomial<S,T> const & pl,Polynomial<U,T> const & pr){
	Polynomial<S,T> tmp = pl;
	for (unsigned int i=0 ; i < pr.coeffs.size();++i)
		tmp.coeffs[i] -= pr.coeffs[i];	
	tmp.P_normalize();
	return tmp;
}

template <unsigned S,unsigned U,typename T>
Polynomial<S+U,T> operator*(Polynomial<S,T> const & pl ,Polynomial<U,T> const & pr){
	Polynomial<S+U,T> tmp;
	for (unsigned int i = 0; i < pl.coeffs.size() ; ++i){
		for (unsigned int j = 0; j < pr.coeffs.size(); ++j){
			tmp.coeffs[i+j] += pl.coeffs[i] * pr.coeffs[j];
		}
	}
	tmp.P_normalize();
	return tmp;
}

template <unsigned S,unsigned U, typename T>
std::pair<Polynomial<S,T>,Polynomial<S,T>> operator/(Polynomial<S,T> const & pl,Polynomial<U,T> const & pr){
	Polynomial<S,T> q; // quoziente	
	Polynomial<S,T> r; // resto
	r = pl;

	if(S < U){
		return std::make_pair(q,r);
	}
	while(r.coeffs.size()-1 >= U){
		q.coeffs[r.coeffs.size()-1-U] = r.coeffs[r.coeffs.size()-1]/pr.coeffs[U];
		for (unsigned int i = 0; i <= U; ++i){
			r.coeffs[r.coeffs.size()-1-U+i] -= q.coeffs[r.coeffs.size()-1-U]*pr.coeffs[i];
		}
		r.P_normalize();
	}			
	q.P_normalize();
	return std::make_pair(q,r);
}

// specializzazione per gli interi
template <unsigned S,unsigned U>
std::pair<Polynomial<S,int>,Polynomial<S,int>> operator/(Polynomial<S,int> const & pl,Polynomial<U,int> const & pr){

	throw std::runtime_error("This operation is not possible with integers");

}
template <unsigned S,unsigned U>
std::pair<Polynomial<S,std::complex<int>>,Polynomial<S,std::complex<int>>> operator/(Polynomial<S,std::complex<int>> const & pl,Polynomial<U,std::complex<int>> const & pr){

	throw std::runtime_error("This operation is not possible with integers");

}
// Streaming operators
template <unsigned int U,typename T>
std::ostream & operator << (std::ostream & str, Polynomial<U,T> const & p){
	if (U != 0){
		for(auto i = p.coeffs.size()-1; i>1; --i){
			if(fabs(p.coeffs[i]) > 1e-10) 	
				str << std::showpos << p.coeffs[i] << "*x^" << std::noshowpos << i << "  ";
		}
		if(fabs(p.coeffs[1]) > 1e-10)
			str << std::showpos << p.coeffs[1] << "*x" << "  ";
	}	
 	str << std::showpos << p.coeffs[0];
	return str;
}

template <unsigned int U,typename T>
std::ostream & operator << (std::ostream & str, Polynomial<U,std::complex<T>> const & p){
	if (U != 0){
		for(auto i = p.coeffs.size()-1; i>1; --i){
			if(std::abs(p.coeffs[i]) > 1e-10) 	
				str << std::showpos << p.coeffs[i] << "*x^" << std::noshowpos << i << "  ";
		}
		if(std::abs(p.coeffs[1]) > 1e-10)
			str << std::showpos << p.coeffs[1] << "*x" << "  ";
	}	
 	str << std::showpos << p.coeffs[0];
	return str;
}



template <unsigned int U,typename T>
std::istream & operator >> (std::istream & str, Polynomial<U,T> & p){
	// Vanno scritti dal coefficiente con potenza più alta. ES: 3x^4+2x^1-3x^0 
	// e ricordarsi di mettere uno spazio alla fine.

	std::string tmp;
	T help;
	unsigned int power;
	std::size_t found;
	// Parte per ricavare il primo coefficiente e il grado massimo	
	std::getline(str,tmp,'^');		
	found = tmp.find('x');
	std::istringstream s(tmp.substr(0,found));
	s >> help;
	std::getline(str,tmp,' ');
	std::istringstream t(tmp.substr(0,tmp.size()));
	t >> power;
	
	if (U < power){
		throw std::runtime_error("Dimensions must agree");
		return str;
	}
	p.coeffs[power] = help;
	
	// Ciclo per ricevere tutti gli altri coefficienti
	while(power != 0){
		std::getline(str,tmp,'^');		
		found = tmp.find('x');
		std::istringstream s(tmp.substr(0,found));
		s >> help;
		std::getline(str,tmp,' ');
		std::istringstream t(tmp.substr(0,tmp.size()));
		t >> power;
				
		// std::cout << "Valore: " << help << std::endl;
		// std::cout << "Potenza: " << power << std::endl;
		p.coeffs[power]=help;
	}
	
	p.P_normalize();
	return str;
}

// ALTRE FUNZIONI

template <unsigned int D,typename T>
T Polynomial<D,T>::operator()(T const & x){
	T u = this->coeffs.back();
	for (int i = coeffs.size()-2; i >= 0; --i){
		u = u*x + this->coeffs[i];
	}
	return u;	
}



// Calcola la derivata di un polinomio
template <unsigned D,typename T>
Polynomial<D-1,T> der_pol(Polynomial<D,T> const & p){
	Polynomial<D-1,T> tmp;
		
	for (int i = p.coeffs.size()-2; i >=0  ; --i){
		tmp.coeffs[i] = p.coeffs[i+1]*static_cast<T>(i+1);
	}
	tmp.P_normalize();
	return tmp;
}

// Solutore per le radici Newton-Horner
template <unsigned int U,typename S>
std::vector<S> solve(Polynomial<U,S> const & p, double x0, const double toll=1e-8, const unsigned MAX_ITER=1000){
	throw std::runtime_error("Solve is defined only for type double");
}

template <unsigned D>
std::vector<double> solve(Polynomial<D,double> const & p, double x0, double const toll=1e-5, unsigned const MAX_ITER=100000){
	unsigned N = p.coeffs.size()-1;
	std::vector<double> root(N);
	std::complex<double> complex_x0(x0,1.1);	
	std::vector<std::complex<double>> help(p.coeffs.size());

	for (unsigned i=0; i< help.size();i++){
		help[i] = p.coeffs[i];
	}

	Polynomial<D,std::complex<double>> tmp(help);
	tmp.P_normalize();
	auto root2 = solve(tmp,std::complex<double>(1.1,1.1));
	bool flag = 0;
	for (unsigned i=0;i<root.size();++i){
		if (fabs(std::imag(root2[i])) > 1e-5)
			flag = 1;
		root[i] = std::real(root2[i]);
	}
	if(flag == 1)
		std::cerr << "ATTENZIONE: questo polinomio presenta radici complesse, le radici calcolate potrebbero non essere esatte" << std::endl;	
	return root;
}


template <unsigned D>
std::vector<std::complex<double>> solve(Polynomial<D,std::complex<double>> const & p, std::complex<double> x0, double const toll=1e-5, unsigned const MAX_ITER=1000){
	unsigned N = p.coeffs.size()-1;
	Polynomial<D,std::complex<double>> tmp = p;
	std::vector<std::complex<double>> help(2);
	help[1] = 1;
	unsigned int niter = 0;
	std::complex<double> x,xnew;
	double diff;
	std::complex<double> pz,dpz;
	std::vector<std::complex<double>> root(N);
	for (unsigned k=0; k<N; ++k){
		niter = 0;
		x = x0;
		diff = 1+toll;
		while(niter < MAX_ITER && diff > toll){
			pz = tmp(x);
			dpz = (der_pol(tmp))(x);
			xnew = x - pz/dpz;
			diff = std::abs(xnew-x);						
			++niter;						
			x = xnew;
		}
		if (niter == MAX_ITER)
			std::cerr << "Raggiunto numero massimo di iterazioni" << std::endl;
		// Deflazione
		root[k] = x;
		help[0] = -x;
		tmp /= Polynomial<1,std::complex<double>>(help);		
	}	
	return root;
}



#endif
