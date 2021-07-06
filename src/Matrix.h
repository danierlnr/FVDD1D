#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <string>
#include <sstream>
#include "DVector.h"
#include <iomanip>
#include <fstream>

class Matrix
{
public:
	Matrix();
	Matrix(size_t,size_t);
	virtual ~Matrix();
	//virtual void   set(double,size_t,size_t){}
	virtual void   set(size_t,size_t,double){}
	virtual double get(size_t,size_t)const{return 0;}
	size_t nrows() const {return m_n;};
	size_t ncols() const {return m_m;};
	//virtual std::string to_string()const{return "";}
	std::string to_string()const;
    //virtual double get(size_t,size_t)const;
	
protected:
	size_t m_n;
	size_t m_m;
};


class DenseMatrix : public Matrix
{
public:
	DenseMatrix();
	DenseMatrix(size_t,size_t);
	DenseMatrix(const double&,size_t,size_t);
	DenseMatrix(const DenseMatrix&); 
	virtual ~DenseMatrix();
	size_t nrows() const {return m_n;};
	size_t ncols() const {return m_m;};
	double* operator[](const size_t i)     {return m_v[i];}
	double* operator[](const size_t i)const{return m_v[i];}
	//void   set(double,size_t,size_t);
	void   set(size_t,size_t,double);
	double get(size_t,size_t)const;
	DenseMatrix& operator=(const double &a);
	void resize(size_t,size_t);
	void resize(const double&,size_t,size_t);
	//std::string to_string()const;
	
private:
	double **m_v;
};

class BandMatrix : public Matrix
{
friend class LSolver;
friend DVector operator*(const BandMatrix&,const DVector&);
public:
	BandMatrix();
	BandMatrix(size_t);
	BandMatrix(const BandMatrix&);
	//BandMatrix(size_t,size_t);
	virtual ~BandMatrix();
	//void   set(double,size_t,size_t);
	void   set(size_t,size_t,double);
	double get(size_t,size_t)const;
	//std::string to_string()const;
	
	void compress(); //(private?)
	size_t m1()const{return m_m1;}
	size_t m2()const{return m_m2;}
    
    // Traspuesta por vector
    int gaxpyt(const DVector&, DVector&);    

private:
	size_t m_m1;
	size_t m_m2;
	bool m_is_triplet;
	DenseMatrix m_dm;
	size_t* m_i;
	size_t* m_j;
	double* m_x;
	size_t  m_nx;
	size_t  m_nn;
	
	void bmrealloc();
	//void compress();
	void get_band();
	
	
};

class SPRMatrix : public Matrix
{
	friend class Solver;
public:
	SPRMatrix(void);
	SPRMatrix(int);
	SPRMatrix(const SPRMatrix&);
	//SPRMatrix(int , int , int , int , int );
	virtual ~SPRMatrix(void);

	int load(char*);
	//void print();
	void printm();
	//std::string toString()const;
	//double get(int, int)const;
	void show();
	
	SPRMatrix& operator=(const SPRMatrix&);
	void print();
	void set(int, int, double);

	void compress();

	int gaxpy(const DVector&, DVector&);
	int gaxpyt(const DVector&, DVector&);
	//int load(char* filename);
	int dupl();
	void sort1(size_t n, size_t* arr, double* brr);
	void get(size_t**,double**);
	//double get(int,int);
	double get(size_t,size_t)const;



protected:
	int     m_nzmax ;	   // numero maximo de elementos de matriz
    int    *m_p ;	       // punteros a las columnas (compressed, n+1 elementos)
					       // o indice de columna (triplet, nzmax elementos)
    int    *m_j ;	       // indices de columna, nzmax elementos
    double *m_x ;	       // valores de los elementos de matriz, nzmax elementos
    int     m_nz ;	       // # de elementos si esta en forma triplet, o -1 si esta compressed
					       // (en este �ltimo caso el # de elementos lo da p[n])
    int    *m_c ;
	size_t* m_ija;
	double* m_sa;


	inline bool TRIPLET()const{return m_nz >=  0;}
	inline bool COMPRES()const{return m_nz == -1;}

	double sp_cumsum(int*,int*,int);
	int sp_compress();

	int MIN(int a, int b){return a<b?a:b;};
	int MAX(int a, int b){return a<b?b:a;};

};

class SPMatrix
{
	friend class Solver;
public:
	SPMatrix(void);
	SPMatrix(int);
	SPMatrix(const SPMatrix&);
	SPMatrix(int , int , int , int , int );
	virtual ~SPMatrix(void);

	int load(char*);
	//void print();
	void printm();
	std::string toString()const;
	double get(int, int)const;
	void show();
	void save();
	
	SPMatrix& operator=(const SPMatrix&);
	void print();
	void set(int, int, double);
	void addDiag(size_t, double);

	void compress();

	int gaxpy(const DVector&, DVector&);
	int gaxpyt(const DVector&, DVector&);
	//int load(char* filename);
	int dupl();


protected:
	int nzmax ;	    // numero maximo de elementos de matriz
    int     m ;	    // numero de filas 
    int     n ;	    // numero de columnas
    int    *p ;	    // punteros a las columnas (compressed, n+1 elementos) 
					// o indice de columna (triplet, nzmax elementos)
    int    *i ;	    // indices de fila, nzmax elementos
    double *x ;	    // valores de los elementos de matriz, nzmax elementos
    int    nz ;	    // # de elementos si esta en forma triplet, o -1 si esta compressed
					// (en este �ltimo caso el # de elementos lo da p[n])

	inline bool TRIPLET()const{return nz >=  0;}
	inline bool COMPRES()const{return nz == -1;}

	double sp_cumsum(int*,int*,int);
	int sp_compress();


};


class ivector
{
public:
	ivector(int n): m_n(n)
	{
		m_p = new(std::nothrow) int[m_n];
	}

	ivector(int n, int val = 0): m_n(n)
	{
		m_p = new(std::nothrow) int[m_n];

		if(m_p != 0)
		{
			for(int i = 0; i < m_n; i++)
				m_p[i] = val;
		}

	}

	virtual ~ivector()
	{
		if(m_p != 0)
		{
			delete[] m_p;
		}
	}

	int& operator[](int i){return m_p[i];};
	int  operator[](int i)const{return m_p[i];};

	int resize(int nnew)
	{
		int* pnew = new(std::nothrow) int[nnew];

		if(pnew == 0)
		{
			return 0; // Error
		}

		for(int i = 0; i < MIN(m_n, nnew); i++)
			pnew[i] = m_p[i];

		delete[] m_p;

		m_p = pnew;
		m_n = nnew;

		return 1;
	}

	void assign(ivector& p)
	{
		if(m_p != 0)
			delete m_p;

		// Toma los datos de p
		m_n = p.m_n;
		m_p = p.m_p;

		// El otro vector, p, queda vacío
		p.m_n = 0;
		p.m_p = 0;
	}

protected:
	int* m_p;
	int  m_n;
	int MIN(int a, int b){return a<b?a:b;};

};

#endif /* MATRIX_H */ 
