///////////////////////////////////////////////
//  ublas.hpp
//	Created by Dominik Lindorfer 11.10.2016
//	Implements a simple Boost::Vector Template
///////////////////////////////////////////////

#ifndef FEMSchroed_ublas_hpp
#define FEMSchroed_ublas_hpp


#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/assign.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/assignment.hpp>


template <typename T>
using Vektor = boost::numeric::ublas::vector<T>;

template <typename T>
using Matrix = boost::numeric::ublas::matrix<T>;

namespace ublas = boost::numeric::ublas;

//*****Push_Back Routine for UBLAS Vektor*****

template <typename T, typename A>
void push_back(Vektor<T>& vec, A item) {

	vec.resize((int)vec.size() + 1);
	vec(vec.size() - 1) = item;
}

#endif
