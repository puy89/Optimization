#pragma once
#include "Misc/Array.h"
class Plus;
class Minus;
class UMinus;
class Prod;
class Div;
class Const;

class Func
{
protected:
	int *ishape, idims, *shape, dims;
	Func* dfunc = nullptr;
	bool is_const = 0;
public:
	virtual Array& operator()(Array& x);
	virtual Array& operator()(Array& x, Array& tar) = 0;
	virtual Func& get_deriv_func() = 0;
	Func& deriv_func() {
		dfunc = (dfunc != nullptr) ? dfunc : &get_deriv_func();
		return *dfunc;
	}

	virtual Array& diff(Array& x) {
		return deriv_func()(x);
	}

	virtual Array& diff(Array& x, Array& tar) {
		return deriv_func()(x, tar);
	}

	Array& newton_raphson(Array &tar, int iters);
	Array& newton(Array& tar, int iters);
	Array& gd(Array& tar, int iters, double step);

	friend Func& operator+(Func &func1, Func &func2);
	friend Func& operator-(Func &func1, Func& func2);
	friend Func& operator*(Func &func1, Func& func2);
	//friend Func operator^(Func& func);
	//friend Func operator/(Func& func);
};

Func& operator+(Func& func1, Func& func2);
Func& operator-(Func& func1, Func& func2);
Func& operator*(Func& func1, Func& func2);
Func& operator^(Func& func1, Func& func2);


class Plus: public Func{
	Func& func1, & func2;
public:
	//virtual Array& operator()(Array& x);
	Plus(Func& func1, Func& func2);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class Minus : public Func {
	Func& func1, & func2;
public:
	//virtual Array& operator()(Array& x);
	Minus(Func& func1, Func& func2);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class UMinus : public Func {
	virtual Array& operator()(Array& x);
	virtual Func& get_deriv_func() = 0;
};

class Prod: public Func {
	Func& func1, & func2;
public:
	//virtual Array& operator()(Array& x);
	Prod(Func& func1, Func& func2);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class Div : public Func {
public:
	Func& func1, &func2;
	//virtual Array& operator()(Array& x);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func() = 0;
};

class Const:public Func {
	Array &val;
public:
	Const(Array& val);
	Const(Array& val, int* shape, int dim, int ishape, int idim);
	Const(double val);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class Var :public Func{
	int size;
public:
	Var(int size);
	~Var();
	//virtual Array& operator()(Array& x);
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class Poly :public Func {
	Func& func;
	double *coeffs;
	int deg;
public:
	//virtual Array& operator()(Array& x);
	Poly(Func& func, double* coeffs, int deg);
	Poly(Func& func, double* coeffs, int deg, bool copy);
	~Poly();
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};


class Norm2 :public Func {
	Func& func;
public:
	//Array& operator()(Array& x);
	Norm2(Func& func);
	~Norm2();
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};

class DiffNorm2 :public Func {
	Func &func1, &func2;
public:
	//Array& operator()(Array &x);
	DiffNorm2(Func &func1, Func &func2);
	~DiffNorm2();
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
	Array& gauss_newton(Array& tar, int iters);
};

class Dot :public Func {
	Func &func1, &func2;
public:
	//virtual Array& operator()(Array& x);
	Dot(Func &func1, Func &func2);
	~Dot();
	virtual Array& operator()(Array& x, Array& tar);
	virtual Func& get_deriv_func();
};