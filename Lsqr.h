#pragma once
#include<iostream>
#include<math.h>
#include "Array.h"



class Lsqr {
    public:
	double c, s, r;
	double conlim = 100000000, damp = 0, btol = 0.00000001, atol = 0.00000001, eps = 0.00000000001;
	//Array x0;
	int iter_lim=0;

	void _sym_ortho(double a, double b);
	
	Array& solve(Applyble& A, Array& b, Array& x);
};

