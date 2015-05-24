#include "polynomial.hpp"
#include <iostream>

int main(){
	std::vector<double> v;
	v.push_back(1);
	v.push_back(-2);
	v.push_back(1);
	v.push_back(-1);
	
	
	std::cout << "Inserire un polinomio" << std::endl;
	Polynomial p1;
	std::cin >> p1;
	std::cout << "p1: " << p1 << std::endl;
	
	auto root = p1.solve(1.2);
	std::cout << std::endl << "Radici: " << std::endl;
	for (auto i:root)
		std::cout << i << " ";
	std::cout << std::endl;

		
	Polynomial p2(v , 3u);
	std::cout << std::endl << "p2: " <<  p2 << std::endl << std::endl;
	
	
	std::cout << "Addizione: " << std::endl;
	std::cout << p2 + p1 << std::endl;
	std::cout << "Sottrazione: " << std::endl;
	std::cout << p1 - p2 << std::endl;
	std::cout << "Moltiplicazione: " << std::endl;
	std::cout << p1*p2 << std::endl;
	std::cout << "Divisione: " << std::endl;
	auto result = p1/p2;
	std::cout << "Quoziente: " << result.first << "\t" << "Resto: " << result.second << std::endl; 
	
	
		
	return 0;	 
}
