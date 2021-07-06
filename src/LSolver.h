#ifndef LSOLVER_H
#define LSOLVER_H

#include <cmath>
#include <vector>
#include "nr.h"
#include "Matrix.h"
using namespace std;

// SPRMatrix
#define NR_END 1
#define FREE_ARG char*


class LSolver
{
public:
	LSolver();
	virtual ~LSolver();
	
	//void lu(Mat_DP &a, Vec_DP &b);
	
	//DVector solve(const BandMatrix&, const DVector&);
	//DVector solve(SPRMatrix&, const DVector&);	
	//DVector solve(SPMatrix&, const DVector&);	

	void solve(const BandMatrix&, const DVector&, DVector&);
	void solve(SPRMatrix&, const DVector&, DVector&);	
	void solve(SPMatrix&, const DVector&, DVector&);	
	
private:
	//void lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b);
	//void ludcmp(Mat_IO_DP&, Vec_O_INT&, DP&);
	void nrerror(const string error_text);
	
	void LU(const BandMatrix&,DenseMatrix&,DenseMatrix&,std::vector<size_t>&,double&);
	DVector back_subs(const DenseMatrix& al,const DenseMatrix& au,const DVector& b,std::vector<size_t>& indx,double& d);

    // SPRMatrix    
    size_t* m_ija;
	double* m_sa;
	void free_dvector(double *v, long nl, long nh);
	void free_stvector(size_t *v, long nl, long nh);
	double *nrdvector(long nl, long nh);
	size_t *stvector(long nl, long nh);
	void linbcg(unsigned long n, double b[], double x[], int itol, double tol,int itmax, int *iter, double *err);
	double snrm(unsigned long n, double sx[], int itol);
	void asolve(unsigned long n, double b[], double x[], int itrnsp);
	void atimes(unsigned long n, double x[], double r[], int itrnsp);
	void dsprstx(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);
	void dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n);

};

#endif /* LSOLVER_H */ 
