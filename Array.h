#pragma once
#include<string>
#include<vector>

class Array;
class Applyble {
public:
	virtual int* get_shape() = 0;
	virtual Array& apply(Array& tar, Array& arr) = 0;

	virtual Array& applyT(Array& tar, Array& arr) = 0;

	double test(Array vec);
};

class Array: public Applyble{

private:
	Array() {
	}
    public:
		double *arr=nullptr;
		int *shape=nullptr;
		int dims=0;
		int *strides=nullptr;
		//char seps = "";
		int size=0;
		bool sub_arr = false;

		void init_shape(int shape[], int dims, int size);


		~Array();

		Array(Array& other);

		Array(const char *filename);

		Array(double arr[], int size);

		Array(int size);

		Array(double arr[], int shape[], int dims);

		Array(double** arr, int shape[]);

		Array(double val, int shape[], int dims);
		
		Array(int shape[], int dims);

		Array(double val);

		std::string toString();

		int* get_shape() {
			return shape;
		}



		Array(double val, int shape[]);

	    Array& reshape(int shape[], int dims);

		Array& ravel();

		Array& operator=(Array& other);

		Array& operator=(double scalar);

		Array& operator+=(Array& other);
	
		Array& operator-=(Array& other);


		Array& operator*=(Array& other);
		
		double operator^(Array& other);

		Array& operator/=(Array& other);

		Array& operator+=(double scalar);

		Array& operator-=(double scalar);

		Array& iminus_row(Array &vector);

		Array& operator*=(double scalar);

		Array& operator/=(double scalar);

		Array& operator&(Array &other);
		
		Array& operator[](int i);

		Array& crop(int a1, int b1, int a2, int b2);

		Array& crop_iplus(int a1, int b1, int a2, int b2, Array &other);

		Array& crop_set(int a1, int b1, int a2, int b2, double scalar);

		Array& sub(int a, int b);

		Array& set(int a, int b, Array &other);

		Array& set_iplus(int a, int b, Array& other);


		Array& iuminus();
		
		Array& iplus_weighted(Array &other, double w);

		Array& iexp();

		Array& ilog();

		double sum4();

		double sump();

		double sum();

		double norm2();

		double norm();

		int argmax();

		Array& inv();


		Array& t();

		Array& apply(Array& tar, Array& vec);

		Array& applyT(Array& tar, Array& vec);

		void save(const char* filename);

		friend Array& zeros(int size);

		friend Array& zeros(int* shape, int dims);

		friend Array& ones(int size);

		friend Array& ones(int* shape, int dims);

		friend Array& range(int size);





		friend Array& eye(int n);


};


Array& zeros(int size);

Array& zeros(int* size, int dims);

Array& ones(int size);

Array& ones(int* size, int dims);

Array& range(int size);

Array& eye(int n);



