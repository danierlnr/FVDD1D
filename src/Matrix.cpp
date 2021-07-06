#include "Matrix.h"


Matrix::Matrix():
m_n(0),m_m(0)
{
	
}

Matrix::Matrix(size_t n,size_t m):
m_n(n),m_m(m)
{
	
}


Matrix::~Matrix()
{
	
}

std::string Matrix::to_string()const
{
    std::ostringstream oss;
    
    oss << std::setprecision(2);

	for (size_t i = 0; i < m_n; i++)
    {
		for (size_t j = 0; j < m_m; j++)
        {
            //std::cout << setw(5) << i << setw(13) << x[i] << setw(13) << f[i] << endl;
            oss << std::setw(7) << get(i,j);
			//m_v[i][j] = a;
        }
        
        oss << std::endl;
    }
    
    return oss.str();
}


//---------------( DenseMatrix)-------------------

DenseMatrix::DenseMatrix():
Matrix(),m_v(0)
{
	
}

DenseMatrix::DenseMatrix(size_t n,size_t m):
Matrix(n,m),m_v(new double*[n])
{
	m_v[0] = new double[m_m*m_n];
	for (size_t i=1; i< n; i++)
		m_v[i] = m_v[i-1] + m;
	
}

DenseMatrix::DenseMatrix(const double& a, size_t n,size_t m):
Matrix(n,m),m_v(new double*[n])
{
	size_t i,j;
	m_v[0] = new double[m_m*m_n];
	for (i=1; i< m_n; i++)
		m_v[i] = m_v[i-1] + m_m;
	for (i=0; i< n; i++)
		for (j=0; j<m_m; j++)
			m_v[i][j] = a;
}

DenseMatrix::DenseMatrix(const DenseMatrix& rhs): 
Matrix(rhs.m_n,rhs.m_m), m_v(new double*[m_n])
{
	size_t i,j;
	m_v[0] = new double[m_m*m_n];
	for (i=1; i< m_n; i++)
		m_v[i] = m_v[i-1] + m_m;
	for (i=0; i< m_n; i++)
		for (j=0; j<m_m; j++)
			m_v[i][j] = rhs[i][j];
}


DenseMatrix::~DenseMatrix()
{
	if (m_v != 0) {
		delete[] (m_v[0]);
		delete[] (m_v);
	}
	
}

inline void DenseMatrix::set(size_t i,size_t j,double a)
{
	m_v[i][j] = a;
}

inline double DenseMatrix::get(size_t i,size_t j)const
{
	return m_v[i][j];
}

DenseMatrix& DenseMatrix::operator=(const double &a)
{
	for(size_t i=0; i< m_n; i++)
		for(size_t j=0; j<m_m; j++)
			m_v[i][j] = a;
	return *this;
}

void DenseMatrix::resize(size_t n,size_t m)
{
	m_n = n;
	m_m = m;

	if (m_v != 0) {
		delete[] (m_v[0]);
		delete[] (m_v);
	}

	m_v[0] = new double[m_m*m_n];
	for (size_t i=1; i< n; i++)
		m_v[i] = m_v[i-1] + m;
	
}

void DenseMatrix::resize(const double& a, size_t n,size_t m)
{
	m_n = n;
	m_m = m;

	if (m_v != 0) {
		delete[] (m_v[0]);
		delete[] (m_v);
	}
	
	size_t i,j;
	
	m_v = new double* [m_n];
	
	m_v[0] = new double[m_m*m_n];
	for (i=1; i< n; i++)
		m_v[i] = m_v[i-1] + m;
	
	for (i=0; i< n; i++)
		for (j=0; j<m_m; j++)
			m_v[i][j] = a;
	
}
/*
std::string DenseMatrix::to_string()const
{
    std::ostringstream oss;
    
   	for(size_t i = 0; i < m_n; i++)
	{
		for(size_t j = 0; j < m_m; j++)
			oss << m_v[i][j] << "  ";
		oss << std::endl;
	}

    return oss.str();
}
*/
//---------------( BandMatrix)-------------------
BandMatrix::BandMatrix():
Matrix(),m_m1(0),m_m2(0),m_is_triplet(true),m_i(0),m_j(0),m_x(0),m_nx(0),m_nn(0)
{
	
}
/*
BandMatrix::BandMatrix(size_t n,size_t m):
m_n(n),m_m(m)
{
	
}
*/
BandMatrix::BandMatrix(size_t n):
Matrix(n,n),m_m1(0),m_m2(0),m_is_triplet(true),m_nx(0),m_nn(n)
{
	m_i = new size_t[m_nn];
	m_j = new size_t[m_nn];
	m_x = new double[m_nn];
	
}

BandMatrix::BandMatrix(const BandMatrix& rhs)
{
    m_is_triplet = rhs.m_is_triplet;
	m_nx = rhs.m_nx;
	m_nn = rhs.m_nx; // solo los necesarios
	m_m1 = rhs.m_m1;
	m_m2 = rhs.m_m2;
    
    if(m_is_triplet)
    {
    	m_i = new size_t[m_nn];
    	m_j = new size_t[m_nn];
    	m_x = new double[m_nn];
    	
    	for(size_t i = 0; i < m_nx; i++)
    	{
    	    m_i[i] = rhs.m_i[i];
    	    m_j[i] = rhs.m_j[i];
    	    m_x[i] = rhs.m_x[i];
    	}
    
    }
    
    else
    {
        m_dm = rhs.m_dm;
        
	    m_i = 0;
	    m_j = 0;
	    m_x = 0;
    }

}


BandMatrix::~BandMatrix()
{
	if(m_i != 0)
		delete[] m_i;
	if(m_j != 0)
		delete[] m_j;
	if(m_x != 0)
		delete[] m_x;
	
}

void BandMatrix::set(size_t i,size_t j,double x)
{
    if(m_is_triplet)
    {
	    m_nx++;
	    //std::cout << "     nx = " << m_nx << std::endl;
	
	    if(m_nx > m_nn)
		    bmrealloc();
		
	    m_i[m_nx-1] = i;	
	    m_j[m_nx-1] = j;	
	    m_x[m_nx-1] = x;	
    }
    
    else
    {
		m_dm[i][ m_m1 + ( j - i )] = x;    
    }
    
}

int BandMatrix::gaxpyt(const DVector& x, DVector& y)    
{
    //Implementacion no eficiente
    double sum=0;
    for (int i=0;i<m_n;i++) 
	{
		sum=0.0;
		for (int j=0;j<m_m;j++) 
		    //sum += fjac[j][i]*fvec[j];
		    sum += get(j,i)*x[j];
		y[i]=sum;
	}

    return 1;
}

double BandMatrix::get(size_t i,size_t j)const
{
	if(m_is_triplet)
		return 0;
		
	if(i <= j + m_m1 && j <= i + m_m2)
		return m_dm[i][ m_m1 + ( j - i )] ;
	
	else return 0;
	
}

void BandMatrix::get_band()
{
	m_m1 = 0;
	m_m2 = 0;
	
	for(size_t i = 0; i < m_nx; i++)
	{
		if(m_i[i] > m_j[i])
		{
		    if(m_i[i] - m_j[i] > m_m1)
		        m_m1 = m_i[i] - m_j[i];
		}
		else if(m_j[i] > m_i[i])
		{
		    if(m_j[i] - m_i[i] > m_m2)
		        m_m2 = m_j[i] - m_i[i];
		}

	}

}


void BandMatrix::compress()
{
    if(!m_is_triplet) return;
        
	get_band();
	
	m_dm.resize(0.0,m_n, m_m1 + m_m2 + 1);
	
	//std::cout << "m1 = " << m_m1 << std::endl;
	//std::cout << "m2 = " << m_m2 << std::endl;
    
    //m_nx
	
	for(size_t i = 0; i < m_nx; i++)
	{
	    //std::cout << "   compressing : " << i << "(" << m_nx << ")" << std::endl;
		m_dm[m_i[i]][ m_m1 + ( m_j[i] - m_i[i] )] += m_x[i];
	}

	if(m_i != 0)
	{
		delete[] m_i;
		m_i = 0;
	}
	if(m_j != 0)
	{
		delete[] m_j;
		m_j = 0;
	}
	if(m_x != 0)
	{
		delete[] m_x;
		m_x = 0;
	}
	
	m_is_triplet = false;

}

void BandMatrix::bmrealloc()
{
	m_nn *= 2;
	
	size_t* iaux;
	size_t* jaux;
	double* xaux;
	
	iaux = new size_t[m_nn];
	for(size_t i = 0; i < m_nx; i++)
	    iaux[i] = m_i[i];
	delete[] m_i;
	m_i = iaux;

	jaux = new size_t[m_nn];
	for(size_t i = 0; i < m_nx; i++)
	    jaux[i] = m_j[i];
	delete[] m_j;
	m_j = jaux;

	xaux = new double[m_nn];
	for(size_t i = 0; i < m_nx; i++)
	    xaux[i] = m_x[i];
	delete[] m_x;
	m_x = xaux;


}
/*
std::string BandMatrix::to_string()const
{
    std::ostringstream oss;
    
   	for(size_t i = 0; i < m_n; i++)
	{
		for(size_t j = 0; j < m_m; j++)
			oss << get(i,j) << "  ";
		oss << std::endl;
	}

    return oss.str();
}
*/
DVector operator*(const BandMatrix& bm,const DVector& v)
{
    DVector u(0.0,v.size());

    size_t i,j,n = bm.nrows();
    
    for (i = 0; i < n; i++) 
    {
        /*
        for (j = 0; j < n; j++) 
        {
            u[i] += bm.get(i,j)*v[j];
        }
        */
	    std::cout << "i = " << i << std::endl;

        if(i < bm.m_m1)
        {
            for (j = 0; j < i + bm.m_m2 + 1; j++) 
                //u[i] += bm.m_dm[i][j]*v[j];
                u[i] += bm.get(i,j)*v[j];
        }

        else if(i + bm.m_m2 > n)
        {
            for (j = i - bm.m_m1; j < n; j++) 
                u[i] += bm.get(i,j)*v[j];
        }

        //if(i - bm.m_m1 > 0 && i + bm.m_m2 + 1 < n)
        else
        {
            for (j = i - bm.m_m1; j < i + bm.m_m2 + 1; j++) 
                u[i] += bm.get(i,j)*v[j];
        }
        
    }

    return u;
}



//---------------( SparseRowMatrix)-------------------
SPRMatrix::SPRMatrix(void):Matrix(),m_nzmax(1),m_p(0),m_j(0),m_x(0),m_nz(0),m_c(0),m_ija(0),m_sa(0)
{

	m_p = new(std::nothrow) int[m_nzmax];
	m_j = new(std::nothrow) int[m_nzmax];
	m_x = new(std::nothrow) double[m_nzmax];
	m_c = new(std::nothrow) int[m_nzmax];

	m_p[m_nz] = 0;
	m_j[m_nz] = 0;
	m_x[m_nz] = 0.0;
	m_c[m_nz] = 0;
	
}


SPRMatrix::SPRMatrix(int n):
Matrix(n,n),m_nzmax(n),m_p(0),m_j(0),m_x(0),m_nz(0),m_c(0),m_ija(0),m_sa(0)
{
	//this->m = m;
	//this->n = m;
	//this->nzmax = m;
	
	//nz = 0;
		
	m_p = new(std::nothrow) int[m_nzmax];
	m_j = new(std::nothrow) int[m_nzmax];
    m_x = new(std::nothrow) double[m_nzmax];
	m_c = new(std::nothrow) int[n];
	for(int k = 0; k < n; k++)
	{
		m_c[k] = 0;
	}
}

SPRMatrix::SPRMatrix(const SPRMatrix& csm):
Matrix(csm.m_n,csm.m_m),m_nzmax(csm.m_nzmax),m_p(0),m_j(0),m_x(0),m_nz(csm.m_nz),m_c(0),m_ija(0),m_sa(0)
{
	//m     = csm.m;
	//n     = csm.n;
	//nz    = csm.nz;
	//m_nzmax = csm.m_nzmax;

	if(csm.TRIPLET())
	{
		m_p = new(std::nothrow) int[m_nzmax];
		m_j = new(std::nothrow) int[m_nzmax];
		m_x = new(std::nothrow) double[m_nzmax];
		m_c = new(std::nothrow) int[m_n];

		for(int k = 0; k < m_nz; k++)
		{
			m_j[k] = csm.m_j[k];
			m_p[k] = csm.m_p[k];
			m_x[k] = csm.m_x[k];
		}

		for(size_t k = 0; k < m_n; k++)
		{
			m_c[k] = csm.m_c[k];
		}
	}

	else
	{
		m_p = new(std::nothrow) int[m_m + 1];
		m_j = new(std::nothrow) int[m_nzmax];
		m_x = new(std::nothrow) double[m_nzmax];
		m_c = new(std::nothrow) int[m_n];

		for(size_t k = 0; k < m_m + 1; k++)
			m_p[k] = csm.m_p[k];

		for(int k = 0; k < m_p[m_m]; k++)
		{
			m_j[k] = csm.m_j[k];
			m_x[k] = csm.m_x[k];
			m_c[k] = csm.m_c[k];
		}
	}

}
/*
SPRMatrix::SPRMatrix(int m, int n, int m_nzmax, int values, int triplet)
{

}
*/

SPRMatrix::~SPRMatrix(void)
{
	if(m_p != 0)
		delete[] m_p;
	if(m_j != 0)
		delete[] m_j;
	if(m_x != 0)
		delete[] m_x;
	if(m_c != 0)
		delete[] m_c;
	if(m_ija != 0)
		delete[] m_ija;
	if(m_sa != 0)
		delete[] m_sa;
		
}

int SPRMatrix::gaxpy(const DVector& x, DVector& y)
{
	unsigned long i,k;
	size_t n = x.size();

    double yi=0;
    for (i=0;i<n;i++) 
    {
        yi=y.get(i);
        yi+=m_sa[i]*x[i];
        for (k = m_ija[i]; k <= m_ija[i+1]-1; k++)      //  Loop over off-diagonal terms.
             yi += m_sa[k]*x[m_ija[k]];
        y.set(i,yi);     
    }
    
    return (1) ;
}



int SPRMatrix::gaxpyt(const DVector& x, DVector& y)
{
	unsigned long i,j,k;
	size_t n = x.size();

	//if (ija[1] != n+2) nrerror("mismatched vector and matrix in sprstx");
	for (i=1;i<n;i++)
		y[i]=m_sa[i]*x[i];
	for (i=1;i<n;i++)
	{
		for (k=m_ija[i];k<=m_ija[i+1]-1;k++)
		{
			j=m_ija[k];
			y[j] += m_sa[k]*x[i];
		}
	}

    return (1) ;
}


void SPRMatrix::get(size_t** ija,double** sa)
{
	compress();

	size_t n = m_ija[m_n];

	(*ija) = new size_t[n];
	(*sa)  = new double[n];

	for(size_t i = 0; i < n; i++)
	{
		(*ija)[i] = m_ija[i];
		(*sa)[i]  = m_sa[i];
	}

}

void SPRMatrix::show()
{
		for(int i = 0; i < m_nz; i++)
		{
			std::cout << "i = " << m_p[i] << " , j = " << m_j[i] << " aij = " << m_x[i]<< std::endl;
		}
}

double SPRMatrix::get(size_t k,size_t l)const
{
	if(TRIPLET())
	{
		int q = 0;

		while(q < m_nz)
		{
			if(m_p[q] == (int)k && m_j[q] == (int)l)
				return m_x[q];

			q++;
		}
	}

	else
	{
        if(k==l)
            return m_sa[k];
        else
        {
            int ic = m_ija[k];
            int fc = m_ija[k+1];
            
            for(int j = ic; j < fc; j++)
            {
                if(m_ija[j]==l)
                {
                    return m_sa[j];
                    continue;
                }
            }
                
            return 0;
        }    
        /*
		int q = m_p[l];

		while(q < m_p[l + 1])
		{
			if(m_j[q] == k)
				return m_x[q];

			q++;
		}
        */ 
	}

	return 0;
}


void SPRMatrix::set(int k, int l, double akl)
{
	if(m_nz < m_nzmax)
	{
	    bool ais=false;
	    int jj = -1; // con este alcanza, no hace falta ais
	    for(int ii=0; ii < m_nz; ii++) // salir de lazo si true
	    {
	        if(m_p[ii] == k && m_j[ii] == l)
	        {
	            ais = true;
	            jj = ii; 
				//std::cout << "i = " << m_p[jj] << " , j = " << m_j[jj] << " ya estaba" << std::endl;
	        }
		}
		
		if(!ais)
		{
		    m_p[m_nz] = k;
		    m_j[m_nz] = l;
		    m_x[m_nz] = akl;
		    if(k != l)
		        m_c[k]++ ;
		    m_nz++;
		}    
		else
		{
		    //std::cout << "i = " << m_p[jj] << " , j = " << m_j[jj] << " repetido" << std::endl;
		    m_x[jj] += akl;
		}    

	}

	else
	{
		int* p2 = new(std::nothrow) int[2*m_nzmax];
		for(int j = 0; j < m_nz; j++)
			p2[j] = m_p[j];
		delete[] m_p;
		m_p = p2;

		int* i2 = new(std::nothrow) int[2*m_nzmax];
		for(int j = 0; j < m_nz; j++)
			i2[j] = m_j[j];
		delete[] m_j;
		m_j = i2;

		double* x2 = new(std::nothrow) double[2*m_nzmax];
		for(int j = 0; j < m_nz; j++)
			x2[j] = m_x[j];
		delete[] m_x;
		m_x = x2;

		m_nzmax = 2*m_nzmax;
        /*
		m_p[m_nz] = k;
		m_j[m_nz] = l;
		m_x[m_nz] = akl;
		if(k != l)
		    m_c[k]++ ;

		m_nz++;
		*/
	    bool ais=false;
	    int jj = -1; // con este alcanza, no hace falta ais
	    for(int ii=0; ii < m_nz; ii++) // salir de lazo si true
	    {
	        if(m_p[ii] == k && m_j[ii] == l)
	        {
	            ais = true;
	            jj = ii; 
	        }
		}
		
		if(!ais)
		{
		    m_p[m_nz] = k;
		    m_j[m_nz] = l;
		    m_x[m_nz] = akl;
		    if(k != l)
		        m_c[k]++ ;
		    m_nz++;
		}    
		else
		{
		    //std::cout << "i = " << m_p[jj] << " , j = " << m_j[jj] << " repetido" <<  std::endl;
		    m_x[jj] += akl;
		}    
	}


}
/*
void SPRMatrix::print()
{
	for(int i = 0; i < m_n; i++)
	{
	    nondiag += m_c[k];
	}

}
*/
void SPRMatrix::compress()
{
	if(COMPRES())
	    return;
	    
	int nondiag = 0;
    //std::cout<< "Compress.." << std::endl;
	for(size_t k = 0; k < m_n; k++)
	{
	    //std::cout<< "c[" << k << "] = " << m_c[k] << std::endl;
	    nondiag += m_c[k];
	}
    /*
	for(int k = 0; k < m_n; k++)
	{
	    std::cout<< "c[" << k << "] = " << m_c[k] << std::endl;
	}
    */
	int nm = m_n + nondiag + 1;
	//int* ija = new int[nm];
	//double* sa = new double[nm];
	m_ija = new size_t[nm];
	m_sa  = new double[nm];
	//for(int k = 0; k < m_n; k++)
	for(int k = 0; k < nm; k++)
	{
		m_ija[k] = 1;
		m_sa[k]  = 0.0;
	}

	m_ija[0] = m_n+1;
	for(size_t k = 1; k < m_n; k++)
		m_ija[k] = m_ija[k-1] + m_c[k-1];

	m_ija[m_n] = nm;
	
	//std::ostringstream oss;
	
	//oss << m_nz << std::endl;

	for(int k = 0; k < m_nz; k++)
	{
		//oss << "p[ " << k << " ] =  " << m_p[k] << "   ;   j[ " << k << " ] =  " << m_j[k] << std::endl;
		if(m_p[k] == m_j[k])
		{
			m_sa[m_p[k]] += m_x[k];
		}
        // m_ija[m_p[k]] apunta al inicio de la fila m_p[k] en m_sa
		else
		{
			int j0 = m_ija[m_p[k]];

			m_ija[j0 + m_c[m_p[k]] - 1] = m_j[k];
			m_sa[j0 + m_c[m_p[k]] - 1] = +m_x[k];
			m_c[m_p[k]] -= 1;
		}

	}
    
   	m_nz = -1;

	/*
	for(int k = 0; k < m_n; k++)
	{
        oss << "m_ija[ " << k << " ] =  " << m_ija[k] << "   ;   m_sa[ " << k << " ] =  " << m_sa[k] << std::endl;
        std::cout<< "k = " << k << std::endl;
		sort1(m_ija[k + 1] - m_ija[k],&m_ija[m_ija[k]],&m_sa[m_ija[k]] );
	}
	
 	std::ofstream fout("compress.txt"); 
    fout << oss.str();
    fout.close();
    */
    
    /*
    std::ostringstream oss;
	for(int k = 0; k < nm; k++)
	{
	    //std::cout<< "+m_ija[" << k << "] = " << m_ija[k] << std::endl;
        oss << "ija[ " << k << " ] =  " << m_ija[k] << "   ;   sa[ " << k << " ] =  " << m_sa[k] << std::endl;
	}
	std::ofstream fout("da.txt"); 
    fout << oss.str();
    fout.close();
    */

	/*
	for(int k = 0; k < nm; k++)
	{
	    std::cout<< "+m_ija[" << k << "] = " << m_ija[k] << std::endl;
	}
	for(int k = 0; k < nm; k++)
	{
	    std::cout<< "+m_sa[" << k << "] = " << m_sa[k] << std::endl;
	}

	for(int k = 0; k < m_n; k++)
	{
	    std::cout<< "c[" << k << "] = " << m_c[k] << std::endl;
	}
	*/
	
	//return 1;
}

void SPRMatrix::sort1(size_t n, size_t arr[], double brr[])
//Sorts an array arr[1..n] into ascending numerical order, by straight insertion, while making
//the corresponding rearrangement of the array brr[1..n] .
{
	int i;
    size_t a;
    double b;
    for(size_t j = 1; j < n; j++)
    {
        a = arr[j];
        b = brr[j];
        i = j-1;
        while (i > -1 && arr[i] > a)
        {
            arr[i+1]=arr[i];
            brr[i+1]=brr[i];
            i--;
        }
        arr[i+1]=a;
        brr[i+1]=b;
    }
}


////////////////////////////////////////
// SPMatrix
////////////////////////////////////////
SPMatrix::SPMatrix(void)
{

}


SPMatrix::SPMatrix(int m)
{
	this->m = m;
	this->n = m;
	this->nzmax = 2*m;
	
	nz = 0;
		
	p = new(std::nothrow) int[nzmax];
	i = new(std::nothrow) int[nzmax];
	x = new(std::nothrow) double[nzmax];
    /*
    for(int j = 0; j < nz; j++)
    {
        i[j]=j;
        p[j]=j;
        x[j]=0;
    }
    */
	
}

void SPMatrix::compress()
{
	sp_compress();
}

void SPMatrix::show()
{
		for(int j = 0; j < nz; j++)
		{
			std::cout << "i = " << i[j] << " , j = " << p[j] << " aij = " << x[j]<< std::endl;
		}
}

void SPMatrix::save()
{
    std::ostringstream oss;
    
    for(int j = 0; j < nz; j++)
	{
		oss << "i = " << i[j] << " , j = " << p[j] << " aij = " << x[j]<< std::endl;
	}
    
    std::ofstream fout("mm.txt");
    
    fout << oss.str();
    
    fout.close();

}



void SPMatrix::addDiag(size_t ij, double aij)
{
    x[ij] += aij;
}


void SPMatrix::set(int k, int l, double akl)
{
	// Este for es provisorio, no dbería esta rutina trabajar asi
	/*
	for(int k = 0; k < nz; k++)
	{
		if((i[k] == k)||(p[k] == l))
		{
			//i[k] = k;
			//p[k] = l;
			x[k] += akl;

			return;
		}

	}
	*/
	if(nz < nzmax)
	{
		i[nz] = k;
		p[nz] = l;
		x[nz] = akl;
		
		nz++;
	}
	
	else
	{
		int* p2 = new(std::nothrow) int[2*nzmax];
		for(int j = 0; j < nz; j++)
			p2[j] = p[j];
		delete[] p;
		p = p2;	
		
		int* i2 = new(std::nothrow) int[2*nzmax];
		for(int j = 0; j < nz; j++)
			i2[j] = i[j];
		delete[] i;
		i = i2;	
		
		double* x2 = new(std::nothrow) double[2*nzmax];
		for(int j = 0; j < nz; j++)
			x2[j] = x[j];
		delete[] x;
		x = x2;	
		
		nzmax = 2*nzmax;
		
		i[nz] = k;
		p[nz] = l;
		x[nz] = akl;
		
		nz++;
	}
	
	
}

SPMatrix::SPMatrix(const SPMatrix& csm)
{
	m     = csm.m;
	n     = csm.n;
	nz    = csm.nz;
	nzmax = csm.nzmax;

	if(csm.TRIPLET())
	{
		p = new(std::nothrow) int[nzmax];
		i = new(std::nothrow) int[nzmax];
		x = new(std::nothrow) double[nzmax];

		for(int k = 0; k < nz; k++)
		{
			i[k] = csm.i[k];
			p[k] = csm.p[k];
			x[k] = csm.x[k];
		}
	}
	
	else
	{
		p = new(std::nothrow) int[n+1];
		i = new(std::nothrow) int[nzmax];
		x = new(std::nothrow) double[nzmax];

		for(int k = 0; k < n + 1; k++)
			p[k] = csm.p[k];

		for(int k = 0; k < p[n]; k++)
		{
			i[k] = csm.i[k];
			x[k] = csm.x[k];
		}
	}

}

SPMatrix::SPMatrix(int m, int n, int nzmax, int values, int triplet)
{

}

SPMatrix::~SPMatrix(void)
{
	delete[] p;
	delete[] i;
	delete[] x;
}

// (Revisar) Vector xx or traspuesta de matriz(?)
int SPMatrix::gaxpy(const DVector& xx, DVector& y)
{
    int pp, j;

    if(!COMPRES())
    	compress();

    for (j = 0 ; j < m ; j++)
    {
		for (pp = p[j] ; pp < p[j+1] ; pp++)
		{
			y [i[pp]] += x[pp] * xx[j] ;
		}
    }

    return (1) ;
}

int SPMatrix::gaxpyt(const DVector& xx, DVector& y)
{
    int pp, j;

    if(!COMPRES())
    	compress();

    for (j = 0 ; j < m ; j++)
    {
		for (pp = p[j] ; pp < p[j+1] ; pp++)
		{
			y [j] += x[pp] * xx[i[pp]] ;
		}
    }

    return (1) ;
}

int SPMatrix::dupl()
{
	if(TRIPLET())
		return 0;

	//int* w = sp_icalloc(m,-1);
	ivector w(n,-1);
	//if(!w)
	//	return 0;
	/*
	for(int l = m_n; l < m_m; l++)
	{
		w[l] = -1;
	}
	*/
	int nz2 = 0, q;

	for(int j = 0; j < m; j++)
	{
		q = nz2;
		for(int k = p[j]; k < p[j+1]; k++)
		{
			int l = i[k];

			if(w[l] >= q)
			{
				x[w[l]] += x[k];
			}

			else
			{
				w[l] = nz2;
				i[nz2] = l;
				x[nz2] = x[k];
				nz2++;
			}

			p[j] = q;
		}
	}

	p[m] = nz2;


	return 1;
}
SPMatrix& SPMatrix::operator=(const SPMatrix& csm)
{
	if(&csm == this)
		return *this;

	delete[] p;
	delete[] i;
	delete[] x;
	
	m     = csm.m;
	n     = csm.n;
	nz    = csm.nz;
	nzmax = csm.nzmax;

	if(csm.TRIPLET())
	{
		p = new(std::nothrow) int[nzmax];
		i = new(std::nothrow) int[nzmax];
		x = new(std::nothrow) double[nzmax];

		for(int k = 0; k < nz; k++)
		{
			i[k] = csm.i[k];
			p[k] = csm.p[k];
			x[k] = csm.x[k];
		}
	}
	
	else
	{
		p = new(std::nothrow) int[n+1];
		i = new(std::nothrow) int[nzmax];
		x = new(std::nothrow) double[nzmax];

		for(int k = 0; k < n + 1; k++)
			p[k] = csm.p[k];

		for(int k = 0; k < p[n]; k++)
		{
			i[k] = csm.i[k];
			x[k] = csm.x[k];
		}
	}

	return *this;
}
	

std::string SPMatrix::toString()const
{
	std::ostringstream ost;

	for(int k = 0; k < m; k++)
	{
		for(int l = 0; l < n; l++)
		{
			ost.width(7);
			ost.precision(3);
            
            if(get(k,l) != 0.0 || k==l)
            {
                ost << "M(" << k <<  "," << l << ") = " << get(k,l) << "  ;  ";// << std::endl;
            }
			
			//ost << get(k,l);
		}
		
		ost << std::endl;
	}

	return ost.str();
}

void SPMatrix::print()
{
	if(TRIPLET())
	{
		for(int k = 0; k < nz; k++)
			std::cout << i[k] << "  " << p[k] << "  " <<  x[k] << std::endl;
	}
	else
	{
		for(int k = 0; k < n + 1; k++)
			std::cout << i[k] << "  " << x[k] << "  " <<  p[k] << std::endl;
		for(int k = n + 1; k < p[n]; k++)
			std::cout << i[k] << "  " << x[k] << std::endl;
	}
}


double SPMatrix::get(int k,int l)const
{
	if(TRIPLET())
	{
		int q = 0;

		while(q < nz)
		{
			if(i[q] == k && p[q] == l)
				return x[q];

			q++;
		}
	}

	else
	{
		int q = p[l];

		while(q < p[l + 1])
		{
			if(i[q] == k)
				return x[q];

			q++;
		}
	}

	return 0;
}

int SPMatrix::sp_compress()
{
	if(COMPRES())
		return 1;

	// Las dimensiones m y n de la matriz en forma compress
	// son las misma, obviamente. Se eliminarán elementos
	// no usados, de modo que nzmax = nz. Al final se pone
	// nz = -1 (expresando asi la forma compress)

	// Contendrá la cantidad de elementos de cada columna
	//int* w = sp_icalloc(n);
	int* w = new int[n];
	if(!w)
		return 0;
	for(int k = 0; k < n; k++)
		w[k] = 0;

	// Cuenta la cantidad de elementos que tiene cada columna
	for(int k = 0; k < nz; k++)
		w[p[k]]++;
	//for(int k = 0; k < n; k++)
	//	std::cout << w[k] << std::endl;

	// El nuevo vector p
	//int* p2 = sp_icalloc(n+1);
	int* p2 = new int[n+1];
	if(!p2)
	{
		delete[] w;
		return 0;
	}

	// Ahora p2 contiene los apuntadores a columna
	sp_cumsum(p2,w,n);
	//printv(p2,n+1);

	// Los nuevos vectores x e i (solo los elementos no
	// nulos: nzamax = nz en la matriz en forma compress
	//int* i2 = sp_icalloc(nz);
	int* i2 = new int[nz];
	if(!i2)
	{
		delete[] w;
		delete[] p2;
		return 0;
	}

	double* x2;
	if(x)
	{
		//x2 = sp_dcalloc(nz);
		x2 = new double[nz];
		if(!x2)
		{
			delete[] w;
			delete[] p2;
			delete[] i2;
			return 0;
		}
	}

	int q;
	for(int k = 0; k < nz; k++)
	{
		// p[k] tiene un valor de columna ->
		// -> w[p[k]] tiene la posición de inicio de esa columna,
		// si es que esa columna aun no fue usada. Si lo fue una vez, por ej.,
		// contiene el lugar siguiente al del inicio, debido al posincremento
		// (para eso esta el posincremento). Así, q tiene la posición de i2 y
		// x2 en que debe almacenarse el valor.
		q = w[p[k]]++;
		i2[q] = i[k];
		if(x)
			x2[q] = x[k];
	}

	delete[] w;

	// Actualizo los vectores, ahora en forma compress
	delete[] p;
	delete[] i;
	delete[] x;

	p = p2;
	i = i2;
	x = x2;

	nzmax = nz;
	nz    = -1;

	return 1;
}
// Suma acumulativa: p[i] = c[0] + .... + c[i-1]
// Se pierde c: es sobreescrito por p.
double SPMatrix::sp_cumsum(int* p, int* c, int n)
{
	int nz = 0;
	double nz2 = 0;

	if(!p || !c)
		return -1;

	for(int i = 0; i < n; i++)
	{
		p[i] = nz;   // Suma acumulada en el ciclo anterior (i-1)
		nz  += c[i]; //Acumula el nuevo elemento
		nz2 += c[i]; //Resguardo por si overflow
		c[i] = p[i]; //Para no volver a sumar lo ya sumado;
	}

	p[n] = nz;

	return nz2;

}



