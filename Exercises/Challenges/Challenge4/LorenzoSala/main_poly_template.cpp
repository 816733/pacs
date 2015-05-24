#include "polynomial_template.hpp"

int main(){
	using real= std::complex<double>;// std::complex<float>;

	std::vector<real> v;
	
	v.emplace_back(1,1);
	v.emplace_back(2,1);
	v.emplace_back(2,0);
	v.emplace_back(1,-1);
	
	/*
	v.emplace_back(1);
	v.emplace_back(2);
	v.emplace_back(2);
	v.emplace_back(1);
	*/
	Polynomial<3,real> p1(v);

	std::cout << "p1: " << p1 << std::endl;
	
	
	std::cout << "Valuto p1 in 3: " << p1(3) << std::endl;
	
	v[3] = 0;
	v[2] = 0;
	v[1] = real(1,0);
	v[0] = 1;
	Polynomial<1,real> p2(v);
	std::cout << "p2: " << p2 << std::endl;
	
	std::cout << "Moltiplicazione: " << p1*p2 << std::endl;
	std::cout << "Addizione: " << p1+p2 << std::endl;
	std::cout << "Sottrazione: " << p1-p2 << std::endl;	
	
	
	try{
		std::cout << "Divisione: " << std::endl;
		auto result = p1/p2;
		std::cout << "Quoziente: " << result.first << "\t" << "Resto: " << result.second << std::endl; 
	}catch(const std::runtime_error & re){
		std::cerr << "Runtime error: " << re.what() << std::endl;	
	}
	
	try{
		
		std::cout << std::endl << "Inserisci un nuovo polinomio per trovare le sue radici (grado massimo 100)" << std::endl;
		Polynomial<100> p;
		std::cin >> p;
	 	std::cout << std::endl << "p: " << p << std::endl;
		
		/*
		std::vector<double> v2;
		v2.emplace_back(1);
		v2.emplace_back(-2);
		v2.emplace_back(1);
		Polynomial<9> p(v2);
		*/
	
		std::cout << std::endl << "Radici di " << p << " : " << std::endl;
		auto root2 = solve(p,1.1);
		for (auto i:root2)
			std::cout << i << " ";
		std::cout << std::endl;

	} catch (const std::runtime_error & re){
		std::cerr << "Runtime error: " << re.what() << std::endl;	
	}
	
	return 0;	 
}
