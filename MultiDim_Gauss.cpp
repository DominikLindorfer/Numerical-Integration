//============================================================================
// Name        : MultiDim_Gauss.cpp
// Author      : lindorfer
// Version     :
// Copyright   : Your copyright notice
// Description : Numerical Integration in 3D using Gauss-Chebychev Quadrature
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

template<typename T, typename Vec_Type, typename limit_type>
double Gauss_Chebychev_Quadrature(T& function, Vec_Type& nodes, Vec_Type& weights, const limit_type& lower, const limit_type& upper, const bool& init_grid = false){

	//-----Inititalize a Sample Grid of 1e6 Points if Grid is NOT Pre-Defined, i.e. init_grid != true-----
	if(init_grid == true){
		double n = 100000;
		double pi = 3.1415926535897932384626434;

		nodes.resize(n);
		weights.resize(n);

		for (int i = 1; i <= n; i++) {
			nodes[i - 1] = cos((2.0 * i - 1.0) / (2.0 * n) * pi);
			weights[i - 1] = pi / n;
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
		sum += weights[i] * f(nodes[i], lower, upper) * sqrt(1 - nodes[i] * nodes[i]);
	}

	return sum;
}

template<typename T, typename Vec_Type, typename limit_type>
double Gauss_Chebychev_Quadrature_2D(T& function, Vec_Type& nodes_x, Vec_Type& nodes_y, Vec_Type& weights_x, Vec_Type& weights_y, const limit_type& lower_x, const limit_type& upper_x, const limit_type& lower_y, const limit_type& upper_y){

	//-----Do the Gauss-Chebychev Integration-----
	double sum = 0.0;
	int n = nodes_x.size();
	int m = nodes_y.size();

	//-----Transform the Upper / Lower Limits to [-1, 1]-----
	auto f = [&function](const double& xpara, const double& ypara, const double& a, const double& b, const double& c, const double& d) {

		double x, y = 0;
		x = (b - a) / 2.0 * xpara + (b + a) / 2.0;
		y = (d - c) / 2.0 * ypara + (d + c) / 2.0;

		return ((b - a) / 2.0 * (d - c) / 2.0 * function(x, y));
	};

	//-----Perform the Gauss-Chebychev-Quadrature Summation-----
	for(int j = 0; j < m; j++){
		for(int i = 0; i < n; i++){
			sum += weights_x[i] * weights_y[j] * f(nodes_x[i], nodes_y[j], lower_x, upper_x, lower_y, upper_y) * sqrt(1 - nodes_x[i] * nodes_x[i]) * sqrt(1 - nodes_y[j] * nodes_y[j]);
		}
	}

	return sum;
}

template<typename T, typename Vec_Type, typename limit_type>
double Gauss_Chebychev_Quadrature_3D(T& function, Vec_Type& nodes_x, Vec_Type& nodes_y, Vec_Type& nodes_z, Vec_Type& weights_x, Vec_Type& weights_y, Vec_Type& weights_z, const limit_type& limits){

	//-----Do the Gauss-Chebychev Integration-----
	double sum = 0.0;

	//-----Transform the Upper / Lower Limits to [-1, 1]-----
	auto f = [&function](const double& xpara, const double& ypara, const double& zpara, const double& a, const double& b, const double& c, const double& d, const double& e, const double& f) {

		double x, y, z= 0;
		x = (b - a) / 2.0 * xpara + (b + a) / 2.0;
		y = (d - c) / 2.0 * ypara + (d + c) / 2.0;
		z = (f - e) / 2.0 * zpara + (f + e) / 2.0;

		return ((b - a) / 2.0 * (d - c) / 2.0 * (f - e) / 2.0 * function(x, y, z));
	};

	//-----Perform the Gauss-Chebychev-Quadrature Summation-----
	#pragma omp parallel for reduction(+:sum)
	for(int k = 0; k < (int)nodes_z.size(); k++){
		for(int j = 0; j < (int)nodes_y.size(); j++){
			for(int i = 0; i < (int)nodes_x.size(); i++){
				sum += weights_x[i] * weights_y[j] * weights_z[k] * f(nodes_x[i], nodes_y[j], nodes_z[k], limits[0], limits[1], limits[2], limits[3], limits[4], limits[5]) * sqrt(1 - nodes_x[i] * nodes_x[i]) * sqrt(1 - nodes_y[j] * nodes_y[j]) * sqrt(1 - nodes_z[k] * nodes_z[k]);
			}
		}
	}

	return sum;
}

int main() {

	//-----Number of Nodes & Pi-----
	double pi = 3.1415926535897932384626434;
	double n = 500;

	//-----Initialize the Integration Grid-----
    vector<double> nodes(n);
    vector<double> weights(n);

    for(int i = 1; i <= n; i++){
        nodes[i-1] = cos((2.0 * i - 1.0)/(2.0*n) * pi);
        weights[i-1] = pi / n;
    }

    //-----Test Lambda Functions-----
	auto test_func = [](double x) {return exp(- x * x);};
	auto test_func_2D = [](double x, double y) {return exp(- x * y);};
	auto test_func_3D = [](double x, double y, double z) {return exp(- x * y * z);};

	vector<double> limits3D = {0,10,0,5,0,10};

	//-----Integrate the Functions-----
	double sum   = Gauss_Chebychev_Quadrature(test_func, nodes, weights, 0, 10);
	double sum2D = Gauss_Chebychev_Quadrature_2D(test_func_2D, nodes, nodes, weights, weights, 0, 10, 0 ,5);
	double sum3D = Gauss_Chebychev_Quadrature_3D(test_func_3D, nodes, nodes, nodes, weights, weights, weights, limits3D);

	cout << "Integral Result with Quadrature: " << sum << "  " << sum2D << "  " << sum3D << endl;

	return 0;
}
