#include "LSolver.h"


LSolver::LSolver()
{
	
}


LSolver::~LSolver()
{
	
}

void LSolver::nrerror(const string error_text)
// Numerical Recipes standard error handler
{
	cerr << "Numerical Recipes run-time error..." << endl;
	cerr << error_text << endl;
	cerr << "...now exiting to system..." << endl;
	//exit(1);
}

void LSolver::solve(const BandMatrix& a, const DVector& b, DVector& x)
//DVector LSolver::solve(const BandMatrix& a, const DVector& b)
{
    //DVector x(0.0,b.size());
    
    DenseMatrix al(a.nrows(),a.m1()), au(a.m_dm);
    
    std::vector<size_t> indx(a.nrows());
    
    double d;
    //std::cout << "LU" << std::endl;
    LU(a,al,au,indx,d);
    
    //std::cout << "BkS" << std::endl;
    x = back_subs(al,au,b,indx,d);

    //return x;
}

void LSolver::solve(SPMatrix& a, const DVector& b, DVector& x)
//DVector LSolver::solve(SPMatrix& a, const DVector& b)
{
    //DVector x(b);

    //return x;
}


void LSolver::LU(const BandMatrix& a,DenseMatrix& al,DenseMatrix& au,std::vector<size_t>& indx,double& d)
{
	size_t i,j,k,l;
	const double tiny = 1.0e-40;

	
	size_t mm;
	double dum;
	
	size_t n = a.nrows();

	mm = a.m1() + a.m2() + 1;
	l  = a.m1();
	
	for (i = 0; i < a.m1(); i++)
	{
		for(j = a.m1() - i; j < mm;j++) 
			au[i][j-l] = au[i][j];
		l--;
		for(j=mm - l - 1; j < mm;j++) 
			au[i][j] = 0.0;
	}

	//*d=1.0;
	d = 1.0;
	
	l=a.m1();
	
	for (k = 0;k < n;k++) 
	{
		dum = au[k][0];
		i = k;
		
		if (l < n) 
			l++;
		
		for (j = k + 1; j < l; j++) 
		{
			if (fabs(au[j][0]) > fabs(dum))
			{
				dum = au[j][0];
				i = j;
			}
		}
		
		indx[k] = i + 1;
		
		if (dum == 0.0) 
			au[k][0]=tiny;
			
		// Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in
        // some applications)

		if (i != k) 
		{
			d = -d;

			for (j = 0; j < mm;j++) 
				swap(au[k][j],au[i][j]);
		}
		
		// Do the elimination
		for (i = k + 1; i < l; i++)
		{
			dum = au[i][0]/au[k][0];
			al[k][i-k-1] = dum;

			for (j = 1; j < mm;j++) 
				au[i][j-1]=au[i][j]-dum*au[k][j];
			
			au[i][mm-1]=0.0;
		}
	}

}

DVector LSolver::back_subs(const DenseMatrix& al,const DenseMatrix& au,const DVector& b,std::vector<size_t>& indx,double& d)
{
	DVector x(b);

    size_t j,k,l,mm,n;
    double dum;
    //mm=m1+m2+1;
    //l=m1;
	n  = au.nrows();
	mm = au.ncols();
	l  = al.ncols();
    
    for (k=0;k<n;k++) 
        x[k] = b[k];
    //Forward substitution, unscrambling the permuted rows as we go.
    //std::cout << "Forward substitution" << std::endl;
    for (k=0;k<n;k++) 
    {
        j=indx[k]-1;
        if (j!=k) 
            SWAP(x[k],x[j]);
        if(l<n) 
            l++;
        for (j=k+1;j<l;j++) 
            x[j] -= al[k][j-k-1]*x[k];
    }

    // Backsubstitution.
    //std::cout << "Backsubstitution" << std::endl;
    
    l = 1;
    //std::cout << "n: " << n << std::endl;
    for(int i = n - 1; i >= 0; i--) 
    {
        //std::cout << "i: " << i << std::endl;
        dum=x[i];
        for (k = 1; k < l; k++) 
            dum -= au[i][k]*x[k+i];
            
        x[i]=dum/au[i][0];
        if (l<mm) 
            l++;
    }
    //std::cout << "returning" << std::endl;
	return x;

}


/*
void LSolver::lu(Mat_DP &a, Vec_DP &b)
{
	int n = b.size();
    Vec_INT indx(n);
    DP p;
    ludcmp(a,indx,p);
    lubksb(a,indx,b);


}

void LSolver::lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
	int i,ii=0,ip,j;
	DP sum;

	int n=a.nrows();
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void LSolver::ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d)
{
	const DP TINY=1.0e-20;
	int i,imax,j,k;
	DP big,dum,sum,temp;

	int n=a.nrows();
	Vec_DP vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}
*/

void LSolver::solve(SPRMatrix& sm, const DVector&  rhs, DVector& x)
//DVector LSolver::solve(SPRMatrix& sm, const DVector&  rhs)
{
	size_t* ija;
	double* sa;

	sm.get(&ija,&sa);

	size_t n  = sm.nrows();

	size_t nz = ija[n];

	m_ija = stvector(1,nz);
	m_sa  = nrdvector(1,nz);

	for(size_t i = 0; i < nz; i++)
	{
		m_ija[i+1] = ija[i] + 1;
		m_sa[i+1]  = sa[i];

		//std::cout << "ija[" << i+1 <<"] = " << m_ija[i+1] << std::endl;
		//std::cout << "sa[" << i+1 <<"] = "  << m_sa[i+1]  << std::endl;

	}
	/*
	for(size_t i = 0; i < nz; i++)
	{
		std::cout << "sa[" << i+1 <<"] = "  << m_sa[i+1]  << std::endl;

	}
	for(size_t i = 0; i < nz; i++)
	{
		std::cout << "ija[" << i+1 <<"] = " << m_ija[i+1] << std::endl;

	}
    */
	double* x0 = nrdvector(1,n);
	double* b  = nrdvector(1,n);

	for(size_t i = 0; i < n; i++)
	{
		b[i+1]  = rhs[i];
		x0[i+1] = x[i];
	}

	double err;
	int itol = 1, itmax = 4000;
	int iter;
	double tol = 1e-6;

	linbcg(n,b,x0,itol,tol,itmax,&iter, &err);
    
    //DVector ans(rhs);

	for(size_t i = 0; i < n; i++)
	{
		//rhs[i] = x[i+1];
		x[i] = x0[i+1];
	}

	free_dvector(x0,1,n);
	free_dvector(b,1,n);

	if(ija != 0)
		delete[] ija;
	if(sa != 0)
		delete[] sa;
        
    //return ans;    

}

// Funciones para resolver con SPRMatrix
#define EPS 1.0e-14

void LSolver::linbcg(unsigned long n, double b[], double x[], int itol, double tol,int itmax, int *iter, double *err)
{
	//void asolve(unsigned long n, double b[], double x[], int itrnsp);
	//void atimes(unsigned long n, double x[], double r[], int itrnsp);
	//double snrm(unsigned long n, double sx[], int itol);
	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p  = nrdvector(1,n);
	pp = nrdvector(1,n);
	r  = nrdvector(1,n);
	rr = nrdvector(1,n);
	z  = nrdvector(1,n);
	zz = nrdvector(1,n);

	*iter=0;
	atimes(n,x,r,0);
	for (j=1;j<=n;j++)
	{
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	znrm=1.0;
	if (itol == 1)
		bnrm=snrm(n,b,itol);
	else if (itol == 2)
	{
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
	}
	else if (itol == 3 || itol == 4)
	{
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0);
		znrm=snrm(n,z,itol);
	}
	else nrerror("illegal itol in linbcg");
	asolve(n,r,z,0);
	while (*iter <= itmax)
	{
		++(*iter);
		zm1nrm=znrm;
		asolve(n,rr,zz,1);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(n,p,z,0);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(n,pp,zz,1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(n,r,z,0);
		if (itol == 1 || itol == 2) {
			znrm=1.0;
			*err=snrm(n,r,itol)/bnrm;
		} else if (itol == 3 || itol == 4) {
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		printf("iter=%4d err=%12.6f\n",*iter,*err);
	if (*err <= tol) break;
	}

	free_dvector(p,1,n);
	free_dvector(pp,1,n);
	free_dvector(r,1,n);
	free_dvector(rr,1,n);
	free_dvector(z,1,n);
	free_dvector(zz,1,n);
}
/*
void Solver::nrerror(char* error_text)
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
*/

void LSolver::dsprsax(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,k;

	if (ija[1] != n+2)
		return;
		//nrerror("dsprsax: mismatched nrdvector and matrix");
	for (i=1;i<=n;i++) {
		b[i]=sa[i]*x[i];
		for (k=ija[i];k<=ija[i+1]-1;k++) b[i] += sa[k]*x[ija[k]];
	}
}

void LSolver::dsprstx(double sa[], unsigned long ija[], double x[], double b[], unsigned long n)
{
	void nrerror(char error_text[]);
	unsigned long i,j,k;
	if (ija[1] != n+2)
		return;
		//nrerror("mismatched vector and matrix in dsprstx");
	for (i=1;i<=n;i++) b[i]=sa[i]*x[i];
	for (i=1;i<=n;i++) {
		for (k=ija[i];k<=ija[i+1]-1;k++) {
			j=ija[k];
			b[j] += sa[k]*x[i];
		}
	}
}

//extern unsigned long ija[];
//extern double sa[];

void LSolver::asolve(unsigned long n, double b[], double x[], int itrnsp)
{
	unsigned long i;

	for(i=1;i<=n;i++)
		x[i] = (m_sa[i] != 0.0 ? b[i]/m_sa[i] : b[i]);
}

void LSolver::atimes(unsigned long n, double x[], double r[], int itrnsp)
{
	//void dsprsax(double sa[], unsigned long ija[], double x[], double b[],unsigned long n);
	//void dsprstx(double sa[], unsigned long ija[], double x[], double b[],unsigned long n);

	if (itrnsp)
		dsprstx(m_sa,m_ija,x,r,n);
	else
		dsprsax(m_sa,m_ija,x,r,n);
}

double LSolver::snrm(unsigned long n, double sx[], int itol)
{
	unsigned long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}

double* LSolver::nrdvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in nrdvector()");
	return v-nl+NR_END;
}

void LSolver::free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with nrdvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

size_t* LSolver::stvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	size_t *v;

	v=(size_t *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(size_t)));
	if (!v) nrerror("allocation failure in stvector()");
	return v-nl+NR_END;
}

void LSolver::free_stvector(size_t *v, long nl, long nh)
/* free a double vector allocated with nrdvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


