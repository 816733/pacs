#include "polynomial.hpp"
#include <string>
#include <sstream>
#include <cmath>

// Costruttori
Polynomial::Polynomial(std::vector<double> _coeffs, unsigned int _max_degree) : coeffs(_coeffs), max_degree(_max_degree){this->P_normalize();}

Polynomial::Polynomial(std::vector<double> _coeffs) : coeffs(_coeffs), max_degree(_coeffs.size() - 1){this->P_normalize();}

Polynomial::Polynomial():max_degree(0){this->coeffs.push_back(0);}



// Normalizzazione
void Polynomial::P_normalize(){
	auto i = this->coeffs.end()-1;
	while(this -> max_degree > 0){
		if(fabs(this->coeffs[max_degree]) < 1e-10){
			this->coeffs.erase(i);
			this->max_degree--;
			--i;
		}else{
			break;
		}	
		
	}
	return;
}

// Metodi operazioni
Polynomial &Polynomial::operator +=(Polynomial const & p){
	if (this->max_degree < p.max_degree){
		this -> coeffs.resize(p.max_degree+1);
		this -> max_degree = p.max_degree;		
	}
	for (unsigned int i=0;i < this->coeffs.size();++i)
			this->coeffs[i] += p.coeffs[i];
	this -> P_normalize();	
	return (*this);
}

Polynomial &Polynomial::operator -=(Polynomial const & p){
	if (this->max_degree < p.max_degree){
		this->coeffs.resize(p.max_degree+1);
		this->max_degree = p.max_degree;		
	}
	for (unsigned int i=0 ; i < this->coeffs.size();++i)
		this->coeffs[i] -= p.coeffs[i];
	this -> P_normalize();
	return (*this);
}


Polynomial & Polynomial::operator /=(Polynomial const & p){
	// Ricordare che l'operazione di divisione restituisce quoziente e resto, in questo caso verrà restituito solo il quoziente.
	auto result = *this / p;
	*this = result.first;
	return *this;
}

Polynomial & Polynomial::operator *=(Polynomial const & p){
	Polynomial tmp;
	tmp.coeffs.resize(max_degree+p.max_degree+1);	
	tmp.max_degree = max_degree+p.max_degree+1;
	for (unsigned int i = 0; i < this->coeffs.size() ; ++i){
		for (unsigned int j = 0; j < p.coeffs.size(); ++j){
			tmp.coeffs[i+j] += this -> coeffs[i] * p.coeffs[j];
		}
	}
	tmp.P_normalize();		
	this -> coeffs.swap(tmp.coeffs);
	this -> max_degree = tmp.max_degree;
	return *this;
}

// Funzioni per le operazioni
Polynomial operator+(Polynomial const & pl ,Polynomial const & pr){
	Polynomial tmp = pl;
	tmp += pr;
	return tmp;
}

Polynomial operator-(Polynomial const & pl,Polynomial const & pr){
	Polynomial tmp = pl;
	tmp -= pr;
	return tmp;
}

Polynomial operator*(Polynomial const & pl ,Polynomial const & pr){
	Polynomial tmp = pl;
	tmp *= pr;
	return tmp;
}

std::pair<Polynomial,Polynomial> operator/(Polynomial const & pl,Polynomial const & pr){
	Polynomial q; // quoziente
	Polynomial r; // resto
	r = pl;
	if(pl.max_degree < pr.max_degree){
		return std::make_pair(q,r);
	}

	q.coeffs.resize(pl.max_degree - pr.max_degree+1);
	q.max_degree = pl.max_degree - pr.max_degree;

	while(r.max_degree >= pr.max_degree){
		q.coeffs[r.max_degree-pr.max_degree] = r.coeffs[r.max_degree]/pr.coeffs[pr.max_degree];
		for (unsigned int i = 0; i < pr.coeffs.size(); ++i){
			r.coeffs[r.max_degree-pr.max_degree+i] -= q.coeffs[r.max_degree-pr.max_degree]*pr.coeffs[i];
		}
		r.P_normalize();
	}			
	q.P_normalize();
	return std::make_pair(q,r);
}


// Streaming operators
std::ostream & operator << (std::ostream & str, Polynomial const & p){
	if (p.max_degree != 0){
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



std::istream & operator >> (std::istream & str, Polynomial & p){
	// Vanno scritti dal coefficiente con potenza più alta. ES: 3x^4+2x^1-3x^0 
	// e ricordarsi di mettere uno spazio alla fine.

	std::string tmp;
	double help;
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
	
	p.coeffs.resize(power+1);
	p.max_degree = power;
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


// Valutare un polinomio
double Polynomial::operator()(double const & x){
	double u = this->coeffs.back();
	for (int i = coeffs.size()-2; i >= 0; --i){
		u = u*x + this->coeffs[i];
	}
	return u;	
}

// Calcola la derivata di un polinomio
Polynomial der_pol(Polynomial const & p){
	Polynomial tmp;
	tmp.coeffs.resize(p.coeffs.size()-1);
	tmp.max_degree = p.max_degree-1;
		
	for (int i = tmp.coeffs.size()-1; i >=0  ; --i){
		tmp.coeffs[i] = p.coeffs[i+1]*(i+1);
	}
	tmp.P_normalize();
	return tmp;
}


// Solutore per le radici Newton-Horner
std::vector<double> Polynomial::solve(double x0, double toll, unsigned const MAX_ITER){
	unsigned N = this->coeffs.size()-1;
	Polynomial tmp = *this;
	std::vector<double> help(2);
	help[1] = 1;
	unsigned int niter = 0;
	double x,xnew, diff;
	double pz,dpz;
	std::vector<double> root(N);
	for (unsigned k=0; k<N; ++k){
		niter = 0;
		x = x0;
		diff = 1+toll;
		while(niter <= MAX_ITER && diff >= toll){
			pz = tmp(x);
			dpz = (der_pol(tmp))(x);
			xnew = x - pz/dpz;
			diff = fabs(xnew-x);			
			++niter;						
			x = xnew;
		}
		// Deflazione
		root[k] = x;
		help[0] = -x;
		tmp /= Polynomial(help);
	}
	
	return root;
}


