//===========================================================================================
// Name        : NumInt.cpp
// Author      : lindorfer
// Version     :
// Copyright   : Your copyright notice
// Description : Numerical Integration using Gauss-Chebychev Integration and Lambda Functions
//===========================================================================================

#include <iostream>
#include "ublas.hpp"

using namespace std;

template<typename T, typename Vec_Type, typename limit_type>
double Gauss_Chebychev_Quadrature(T& function, Vec_Type& nodes, Vec_Type& weights, const limit_type& lower, const limit_type& upper, const bool& init_grid = false){

	//-----Inititalize a Sample Grid of 1e6 Points if Grid is NOT Pre-Defined, i.e. init_grid != true-----
	if(init_grid == true){
		double n = 100000;
		double pi = 3.1415926535897932384626434;

		nodes.resize(n);
		weights.resize(n);

		for (long double i = 1; i <= n; i++) {
			nodes(i - 1) = cos((2.0 * i - 1.0) / (2.0 * n) * pi);
			weights(i - 1) = pi / n;
		}
	}

	//-----Do the Gauss-Chebychev Integration-----
	double sum = 0.0;
	int n = nodes.size();

	//-----Transform the Upper / Lower Limits to [-1, 1]-----
	auto f = [&function](const double& xpara, const double& a, const double& b) {
		double x = 0;
		x = (b - a) / 2.0 * xpara + (b + a) / 2.0;
		return ((b - a) / 2.0 * function(x));
	};

	//-----Perform the Gauss-Chebychev-Quadrature Summation-----
	for (int i = 0; i < n; i++) {
		sum += weights(i) * f(nodes(i), lower, upper) * sqrt(1 - nodes(i) * nodes(i));
	}

	return sum;
}


int main() {

	//-----Number of Nodes & Pi-----
	double pi = 3.1415926535897932384626434;
	double n = 10000;

	//-----Initialize the Integration Grid-----
    Vektor<long double> nodes(n);
    Vektor<long double> weights(n);

    for(long double i = 1; i <=n; i++){
        nodes(i-1) = cos((2.0 * i - 1.0)/(2.0*n) * pi);
        weights(i-1) = pi / n;
    }

    //-----Test Lambda Function-----
	auto test_func = [](double x) {return exp(- x * x);};

	//-----Integrate the Function in the Region [0,10]-----
	double sum = Gauss_Chebychev_Quadrature(test_func, nodes, weights, 0, 10);
	cout << "Integral Result with Quadrature: " << sum << endl;

	return 0;
}
