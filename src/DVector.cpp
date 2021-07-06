#include "DVector.h"


DVector::DVector():
m_n(0),m_v(0)
{
	
}

DVector::DVector(size_t n):
m_n(n)
{
	m_v = new double[m_n];
}

DVector::DVector(const double& a,size_t n):
m_n(n)
{
	m_v = new double[m_n];

	for (size_t i = 0; i < m_n; i++)
		m_v[i] = a;

}

DVector::DVector(const DVector& rhs):
m_n(rhs.m_n)
{
	m_v = new double[m_n];

	//std::cout << "DVector::DVector(const DVector& rhs)" << std::endl;
	
	for (size_t i = 0; i < m_n; i++)
		m_v[i] = rhs[i];
	
}


DVector::~DVector()
{
	if(m_v!=0)
	    delete[] m_v;
}

std::string DVector::to_string()const
{
    std::ostringstream oss;
    
   	for(size_t i = 0; i < m_n; i++)
	{
		oss << m_v[i] << std::endl;
	}

    return oss.str();
}

double DVector::norm()const
{
    double ans = 0.0;
    
   	for(size_t i = 0; i < m_n; i++)
	{
		ans += m_v[i]*m_v[i];
	}
	
	return sqrt(ans);
}

DVector& DVector::operator = (const DVector& rhs)
{
	if (this != &rhs)
	{
		if (m_n != rhs.m_n) 
		{
			if (m_v != 0) 
			    delete [] (m_v);
			m_n = rhs.m_n;
			m_v = new double[m_n];
		}
		for (size_t i = 0; i < m_n; i++)
			m_v[i] = rhs[i];
	}
	
	return *this;

}

DVector& DVector::operator = (const double& a)
{
	for (size_t i = 0; i < m_n; i++)
		m_v[i] = a;
	
	return *this;

}


void DVector::operator*=(const double a)
{
	for (size_t i = 0; i < m_n;i++)
		m_v[i] *= a;
}

void DVector::operator/=(const double a)
{
	for (size_t i = 0; i < m_n;i++)
		m_v[i] /= a;
}



DVector operator+(const DVector& u,const DVector& v)
{
	DVector w(u);

	size_t n = w.size();
		
	for(size_t i = 0; i < n; i++)
		w.m_v[i] += v.m_v[i];

	return w;
}

DVector operator-(const DVector& u,const DVector& v)
{
	DVector w(u);

	size_t n = w.size();
		
	for(size_t i = 0; i < n; i++)
		w.m_v[i] -= v.m_v[i];

	return w;
}

double operator*(const DVector& u,const DVector& v)
{
	double ans = 0.0;

	size_t n = u.size();
		
	for(size_t i = 0; i < n; i++)
		ans += u.m_v[i]*v.m_v[i];

	return ans;
}

DVector operator*(const DVector& u,const double& a)
{
	DVector w(u);

	size_t n = w.size();
		
	for(size_t i = 0; i < n; i++)
		w.m_v[i] *= a;

	return w;
}

