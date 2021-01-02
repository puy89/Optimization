#include "Func.h"
#include "Misc/Lsqr.h"
#include <iostream>

Array& Func::operator()(Array& x) {
	Array& res = *(new Array(shape, dims));
	return (*this)(x, res);
}

Array& Func::newton_raphson(Array &tar, int iters){
	
	Lsqr lsqr;
	Array &f_x = (*this)(tar);
	Array &df_x = diff(tar);
	Array sol(tar.shape, tar.dims);
	if (tar.size > 1){
		tar -= lsqr.solve(df_x, f_x, sol);
		for (int i = 1; i < iters; i++)
			tar -= lsqr.solve(diff(tar, df_x), (*this)(tar, f_x), sol);
	}
	else{
		tar -= f_x /= df_x;
		for (int i = 1; i < iters; i++)
			tar -= (*this)(tar, f_x)/=diff(tar, df_x);
	}

	delete &f_x;
	delete &df_x;
	return tar;
}

Array& Func::newton(Array& tar, int iters) {
	return deriv_func().newton_raphson(tar, iters);
}

Array& Func::gd(Array& tar, int iters, double step) {
	Array& df_x = diff(tar);
	Array sol(tar.shape, tar.dims);
	for (int i=0; i < iters - 1; i++)
		tar -= (diff(tar, df_x) *= step);
	delete& df_x;
	return tar;
}


Plus::Plus(Func& func1, Func& func2) : func1(func1), func2(func2){
}


Array& Plus::operator()(Array &x, Array &tar){
	func1(x, tar);
	tar += func2(x);
	return tar;
}

Func& Plus::get_deriv_func() {
	return *(new Plus(func1.deriv_func(), func2.deriv_func()));
}


Func& operator+(Func &func1, Func &func2) {
	return *(new Plus(func1, func2));
}

Minus::Minus(Func& func1, Func& func2) : func1(func1), func2(func2) {
}


Array& Minus::operator()(Array& x, Array& tar) {
	func1(x, tar);
	tar -= func2(x);
	return tar;
}

Func& Minus::get_deriv_func() {
	return *(new Minus(func1.deriv_func(), func2.deriv_func()));
}



Func& operator-(Func& func1, Func& func2) {
	return *(new Minus(func1, func2));
}


Prod::Prod(Func& func1, Func& func2) : func1(func1), func2(func2) {
}


Array& Prod::operator()(Array& x, Array& tar) {
	func1(x, tar);
	tar *= func2(x);
	return tar;
}



Func& Prod::get_deriv_func() {

	return *(new Plus(*(new Prod(func1.deriv_func(), func2)),
		              *(new Prod(func1, func2.deriv_func()))));
}


Func& operator*(Func& func1, Func& func2) {
	
	return *(new Prod(func1, func2));
}

Const::Const(Array &val):val(val){
	is_const = true;
}

Const::Const(double val) : val(*(new Array(val))){
	is_const = true;
}


Array& Const::operator()(Array& x, Array& tar) {
	tar = val;
	return tar;

}
Func& Const::get_deriv_func() {
	return *(new Const(*(new Array(ishape, idims))));
}

Var::Var(int size){
	ishape = shape = new int[1]{1};
	idims = dims = 1;
	this->size = size;
}

Var::~Var(){
	delete shape;
}

Array& Var::operator()(Array& x, Array& tar){
	tar = x;
	return tar;
}

Func& Var::get_deriv_func(){
	return *(new Const(eye(size)));
}


Poly::Poly(Func &func, double* coeffs, int deg):func(func){
	this->coeffs = new double[deg+1];
	this->deg = deg;
	for (int i = 0; i <= deg; i++)
		this->coeffs[i] = coeffs[i];
}

Poly::Poly(Func &func, double *coeffs, int deg, bool copy) :func(func) {
	this->deg = deg;
	if (copy) {
		this->coeffs = new double[deg + 1];
		for (int i = 0; i <= deg; i++)
			this->coeffs[i] = coeffs[i];
	}
	else
		this->coeffs = coeffs;
}


Poly::~Poly() {
	delete[] coeffs;
}

Array& Poly::operator()(Array& x, Array& tar){
	tar = coeffs[0];
	Array &funcx = func(x);
	Array power(funcx);
	//std::cout << "before"<< " " << x.toString() << '\n';
	for (int i = 1; i <= deg; i++) {
		Array mul(power);
		tar += (mul*=coeffs[i]);
		//std::cout << "iter " << i << " " << tar.toString()<<'\n';
		power *= funcx;
	}
	delete &funcx;
	return tar;
}

Func& Poly::get_deriv_func() {
	double *dcoeffs = new double[deg];
	for (int i = 0; i < deg; i++)
		dcoeffs[i] = (i+1)*coeffs[i+1];
	return *(new Poly(func, dcoeffs, deg-1)) * func.deriv_func();
}

Norm2::Norm2(Func& func):func(func){
}

Norm2::~Norm2(){
}

Array& Norm2::operator()(Array& x, Array& tar) {
	Array& funcx = func(x);
	tar.arr[0] = funcx^funcx;
	tar /= 2;
	delete &funcx;
	return tar;
}

Func& Norm2::get_deriv_func(){
	return func.deriv_func(); 
}

DiffNorm2::DiffNorm2(Func &func1, Func &func2) :func1(func1), func2(func2){
}

DiffNorm2::~DiffNorm2() {
}

Array& DiffNorm2::operator()(Array& x, Array& tar) {
	Array& func1x = func1(x);
	Array& func2x = func2(x);
	func1x -= func2x;
	tar.arr[0] = func1x ^ func2x;
	tar /= 2;
	delete& func1x;
	delete& func2x;
	return tar;
}

Func& DiffNorm2::get_deriv_func() {
	return func1 - func2;
}

Array& DiffNorm2::gauss_newton(Array& tar, int iters) {
	Lsqr lsqr;
	Array& f_x = func1(tar);
	Array& f2_x = func2(tar);
	f_x -= f2_x;
	Array& df_x = diff(tar);
	Array sol(tar.shape, tar.dims);
	tar -= lsqr.solve(diff(tar, df_x), f_x, sol);
	std::cout << "tar = " << tar.toString() << "\n";
	for (int i = 1; i < iters; i++)
		tar += lsqr.solve(diff(tar, df_x), (*this)(tar, func1(tar, f_x)-=f2_x), sol);
	delete& f_x;
	delete& df_x;
	delete& f2_x;
	return tar;
}

Dot::Dot(Func &func1, Func &func2):func1(func1), func2(func2){
}

Dot::~Dot() {
}

Array& Dot::operator()(Array& x, Array& tar) {
	Array& func1x = func1(x);
	Array& func2x = func2(x);
	tar.arr[0] = func2x ^ func1x;
	delete& func1x;
	delete& func2x;
	return tar;
}

Func& Dot::get_deriv_func(){
	return func1^func2.deriv_func() + func1.deriv_func()^func2;
}

Func& operator^(Func& func1, Func& func2) {
	return *(new Dot(func1, func2));
}