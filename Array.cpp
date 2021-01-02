#include "Array.h"
#include "Lsqr.h"
#include <iostream>
#include <sstream>
#include <fstream>


double Applyble::test(Array vec) {
	/*
	This function tests if applyT actually perform transpose of apply
	Should be close to zero
	*/
	Array &Ax = zeros(get_shape()[0]);
	Array &ATAx = zeros(get_shape()[1]);
	apply(Ax, vec);
	applyT(ATAx, Ax);
	delete &Ax;
	delete &ATAx;
	return (Ax^Ax) - (vec^ATAx);
}

int prod_all(int shape[], int dims) {
	int rv = 1;
	for (int i = 0; i < dims; i++)
		rv *= shape[i];
	return rv;
}


void Array::init_shape(int shape[], int dims, int size) {
	this->size = size;
	this->dims = dims;
	this->shape = new int[dims];
	for (int i = 0; i < dims; i++)
		this->shape[i] = shape[i];
	int stride = 1;
	strides = new int[dims];
	for (int i = dims - 1; i >= 0; i--) {
		strides[i] = stride;
		stride *= shape[i];
	}
	if (stride != size)
		std::cerr << "error:"<<"size = "<<size<<" stride="<<stride<<" strides = "<<" "<<strides[0]<<" "<<strides[1]<<"\n";
}

Array::~Array() {
	if(!sub_arr)
	    delete arr;
	delete shape;
	//std::cout << "delete";
}

Array::Array(const char* filename) {
	std::fstream file;
	file.open(filename);
	std::vector<double> vals;
	char s[21];
	if (!file)
		return;
	double val;
	while (file >> val)
		vals.push_back(val);
	//std::cout << (file>>val);
	file.close();
	size = vals.size();
	arr = new double[size];
	for (int i = 0; i < size; i++)
		arr[i] = vals[i];
	shape = new int[1]{ size };
	init_shape(shape, 1, size);
	
}


Array::Array(double arr[], int shape[], int dims) {
	init_shape(shape, dims, prod_all(shape, dims));
	this->arr = new double[size];
	for (int i = 0; i < size; i++)
		this->arr[i] = arr[i];
}

Array::Array(double arr[], int size) {
	
	init_shape(new int[1]{size}, 1, size);
	size = prod_all(shape, dims);
	this->arr = new double[size];
	for (int i = 0; i < size; i++)
		this->arr[i] = arr[i];
}

Array::Array(Array& other){
	init_shape(other.shape, other.dims, other.size);
	this->arr = new double[size];
	for (int i = 0; i < size; i++)
		this->arr[i] = other.arr[i];

}

Array::Array(double val, int shape[], int dims) {
	size = prod_all(shape, dims);
	init_shape(shape, dims, size);
	this->arr = new double[size];
	for (int i = 0; i < size; i++)
		this->arr[i] = val;

}

Array::Array(int shape[], int dims) {
	size = prod_all(shape, dims);
	init_shape(shape, dims, size);
	this->arr = new double[size];
}

Array::Array(double val) {
	arr = new double[1]{val};
	shape = new int[1]{1};
	size = 1;
	dims = 1;
}



Array& zeros(int* shape, int dims) {
	return *(new Array(0., shape, dims));
}

Array& zeros(int size) {
	Array& res = *(new Array());
	res.arr = new double[size];
	res.init_shape(new int[1]{ size }, 1, size);
	for (int i = 0; i < size; i++)
		res.arr[i] = 0;
	return res;
}


Array& ones(int* shape, int dims) {
	return *(new Array(1., shape, dims));
}

Array& ones(int size) {
	Array& res = *(new Array());
	res.arr = new double[size];
	res.init_shape(new int[1]{ size }, 1, size);
	for (int i = 0; i < size; i++)
		res.arr[i] = 1;
	return res;
}

Array::Array(double **arr, int shape[]) {
	init_shape(shape, 2, shape[0] * shape[1]);
	this->arr = new double[size];
	for (int i = 0, k = 0; i < shape[0]; i++)
		for (int j = 0; j < shape[1]; j++, k++)
			this->arr[k] = arr[i][j];
}


Array& eye(int n) {
	Array &res = *(new Array());
	res.size = n * n;
	res.arr = new double[res.size];
	res.init_shape(new int[2]{n, n}, 2, res.size);
	for(int i = 0; i < res.size; i++)
		res.arr[i] = 0;
	for (int i = 0; i < res.size; i+=n+1)
		res.arr[i] = 1;
	return res;
}


Array::Array(int size) {
	
	arr = new double[size];
	init_shape(new int[1]{ size }, 1, size);
	for (int i = 0; i < size; i++)
		arr[i] = 0;
}




Array& range(int size) {
	Array &res = *(new Array());
	res.arr = new double[size];
	res.init_shape(new int[1]{ size }, 1, size);
	//std::cout <<"range size = " <<size;
	for (int i = 0; i < size; i++)
		res.arr[i] = i;
	return res;
}


std::string Array::toString() {
	int enter = (dims < 2) ? shape[0] : shape[1];
	std::ostringstream s;
	//std::cout <<"shape"<<shape[0]<<" "<<shape[1]<<"\n";
	for (int i = 0; i<size; i++) {
		
		s<<arr[i]<<" ";
		if ((i + 1) % enter == 0)
			s<<"\n";
	}
	return s.str();
}


Array& Array::reshape(int shape[], int dims){
	Array &res = *(new Array());
	res.arr = arr;
	res.init_shape(shape, dims, size);
	res.sub_arr = true;
	return res;
}

Array& Array::ravel() {
	Array& res = *(new Array());
	res.arr = arr;
	res.init_shape(new int[]{size}, 1, size);
	res.sub_arr = true;
	return res;
}

Array& Array::operator=(Array& other) {
	for (int i = 0; i < size; i++)
		this->arr[i] = other.arr[i];
	return (*this);
}

Array& Array::operator=(double scalar) {
	for (int i = 0; i < size; i++)
		this->arr[i] = scalar;
	return (*this);
}


Array& Array::operator+=(Array &other) {
	for (int i = 0; i < size; i++)
		this->arr[i] += other.arr[i];
	return (*this);
}

Array& Array::operator-=(Array& other) {
	for (int i = 0; i < size; i++)
		this->arr[i] -= other.arr[i];
	return (*this);
}

Array& Array::iminus_row(Array& other) {
	for(int i = 0, k=0; i < shape[0]; i++)
		for (int j = 0; j < shape[1]; j++, k++)
		   this->arr[k] -= other.arr[j];
	return (*this);
}


Array& Array::operator*=(Array& other) {
	for (int i = 0; i < size; i++) 
		this->arr[i] *= other.arr[i];
	return (*this);
}

Array& Array::operator/=(Array& other) {
	for (int i = 0; i < size; i++)
		this->arr[i] /= other.arr[i];
	return (*this);
}

Array& Array::operator+=(double scalar) {
	for (int i = 0; i < size; i++)
		this->arr[i] += scalar;
	return (*this);
}


Array& Array::operator-=(double scalar) {
	for (int i = 0; i < size; i++)
		this->arr[i] -= scalar;
	return (*this);
}

Array& Array::operator*=(double scalar) {
	for (int i = 0; i < size; i++)
		this->arr[i] *= scalar;
	return (*this);
}

Array& Array::operator/=(double scalar){
	for (int i = 0; i < size; i++)
		this->arr[i] /= scalar;
	return (*this);
}

double Array::operator^(Array& other) {
	double res = 0;
	for (int i = 0; i < size; i++)
		res += this->arr[i]*other.arr[i];
	return res;
}

Array& Array::iplus_weighted(Array &other, double w) {
	for (int i = 0; i < size; i++)
		arr[i] += other.arr[i] * w;
	return (*this);
}


Array& Array::iuminus() {
	for (int i = 0; i < size; i++)
		arr[i] = -arr[i];
	return (*this);
}

Array& Array::iexp() {
	for (int i = 0; i < size; i++)
		arr[i] = exp(arr[i]);
	return (*this);
}

Array& Array::ilog() {
	for (int i = 0; i < size; i++)
		arr[i] = log(arr[i]);
	return (*this);
}


double Array::sum4() {
	double res = 0;
	for (int i = 0; i < size;i+=4)
		res += arr[i]+arr[i+1]+arr[i+2]+arr[i+3];

	return res;
}

double Array::sump() {
	double res = 0;
	for (int i = 0; i < size;)
		res += arr[i++] + arr[i++] + arr[i++] + arr[i++];
	return res;
}

double Array::sum() {
	double res = 0;
	for (int i = 0; i < size;i++)
		res += arr[i];
	return res;
}

double Array::norm2() {
	double res = 0;
	for (int i = 0; i < size; i++)
		res += arr[i]*arr[i];
	return res;
}

double Array::norm() {
	return sqrt(norm2());
}

int Array::argmax() {
	double max = arr[0];
	int arg = 0;
	for (int k = 1; k < size; k++)
		if (arr[k] > max) {
			max = arr[k];
			arg = k;
		}
	return arg;
}

Array& Array::operator&(Array& other) {
	Array &res = *(new Array());
	res.arr = new double[this->size+other.size];
	res.shape = new int[1]{this->size + other.size};
	for (int k = 0; k < this->size; k++)
		res.arr[k] = this->arr[k];
	for (int i = 0, k = this->size; i < other.size; k++, i++)
		res.arr[k] = other.arr[i];
	return res;
}

Array& Array::operator[](int i){
	double *row_arr = arr + i*strides[0];
	Array &res = *(new Array());
	res.arr = row_arr;
	res.sub_arr = true;
	res.size = strides[0];
	res.shape = new int[1]{strides[0]};
	res.init_shape(res.shape, 1, res.size);
	return res;
}

Array& Array::crop(int a1, int b1, int a2, int b2) {
	int *new_shape = new int[] {b1 - a1, b2 - a2};
	Array &res = *(new Array(new_shape, 2));
	for (int i = 0, raveli = a1 * strides[0] + a2 * strides[1], k = 0; i < new_shape[0]; i++, raveli += this->strides[0])
		for (int j = 0, l = raveli; j < new_shape[1]; j++, k++, l++)
			res.arr[k] = arr[l];
	return res;
}

Array& Array::crop_iplus(int a1, int b1, int a2, int b2, Array &other) {
	int m=b1 - a1, n=b2 - a2;
	for (int i = 0, raveli = a1 * strides[0] + a2 * strides[1], k = 0; i < m; i++, raveli += strides[0])
		for (int j = 0, l = raveli; j < n; j++, k++, l++) {
			arr[l] += other.arr[k];
			//std::cout << "i = " << i << "j = " << j << "l = " << l << "k = " << k<<'\n';
		}
	return *this;
}

Array& Array::crop_set(int a1, int b1, int a2, int b2, double scalar) {
	int m = b1 - a1, n = b2 - a2;
	for (int i = 0, raveli = a1 * strides[0] + a2 * strides[1], k = 0; i < m; i++, raveli += strides[0])
		for (int j = 0, l = raveli; j < n; j++, k++, l++) {
			arr[l] = scalar;
			//std::cout << "i = " << i << "j = " << j << "l = " << l << "k = " << k<<'\n';
		}
	return *this;
}

Array& Array::sub(int a, int b) {
	double* row_arr = arr + a;
	Array& res = *(new Array());
	res.arr = row_arr;
	res.sub_arr = true;
	res.size = b-a;
	res.shape = new int[1]{ res.size};
	res.init_shape(res.shape, 1, res.size);
	return res;
}

Array& Array::set(int a, int b, Array &other) {
	for (int i = a, k=0; i < b; i++, k++)
		this->arr[i] = other.arr[k];
	return *this;
}

Array& Array::set_iplus(int a, int b, Array& other) {
	for (int i = a, k = 0; i < b; i++, k++)
		this->arr[i] += other.arr[k];
	return *this;
}

Array& Array::apply(Array &tar, Array &vec) {
	for (int i = 0, k = 0; i < shape[0]; i++)
		for (int p = 0; p < shape[1]; p++, k++) {
			tar.arr[i] += this->arr[k] * vec.arr[p];
			//std::cout << " i = " << i << " p = " << p << " k = " << k<<"\n";
		}
	return tar;
}


Array& Array::applyT(Array &tar, Array &vec){
	Array &trans = t(); //fix
	return trans.apply(tar, vec);	
}

Array& Array::t(){
	Array &res = zeros(new int[2]{ this->shape[1], this->shape[0] }, 2);
	for (int i = 0, raveledi = 0, k = 0; i < res.shape[0]; i++, raveledi += this->strides[1]) {
		for (int j = 0, raveledj = 0; j < res.shape[1]; j++, k++, raveledj += this->strides[0]) {
			res.arr[k] = this->arr[raveledi + raveledj];
		}
	}
	return res;
}

Array& Array::inv() {
	Array& res = *(new Array(shape, dims));
	Lsqr lsqr;
	Array vec(shape[1]);
	Array &e = eye(shape[0]);
	for (int j = 0; j < shape[1]; j++){
		Array &r = e[j];
		lsqr.solve(*this, r, vec);
		for (int i = 0, raveledi = j; i < shape[0]; i++, raveledi += strides[0])
			res.arr[raveledi] = vec.arr[i];
		delete &r;
	}
    delete &e;
	return res;
}

void Array::save(const char* filename){
	std::ofstream s(filename);
	for (int i = 0; i < size; i++)
		s << arr[i] << ' ';
	s.close();
}

