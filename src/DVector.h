#ifndef DVECTOR_H
#define DVECTOR_H
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

class DVector
{
friend DVector operator+(const DVector&,const DVector&);
friend DVector operator-(const DVector&,const DVector&);
friend double  operator*(const DVector&,const DVector&);
friend DVector operator*(const DVector&,const double&);

public:
	DVector();
	DVector(size_t);
	DVector(const double&, size_t);
	DVector(const DVector&);
	virtual ~DVector();
	void   set(size_t i,double v) {m_v[i]=v;}
	double get(size_t i)const{return m_v[i];}
	size_t size() const {return m_n;};
	std::string to_string()const;//{return "";}
	double& operator[](const size_t i){return m_v[i];}
	double& operator[](const size_t i)const{return m_v[i];}
	DVector& operator = (const DVector&);
	DVector& operator = (const double&);
	double norm()const;
	void operator *= (double);
	void operator /= (double);
protected:
	size_t  m_n;
    double* m_v;
};

#endif /* DVECTOR_H */ 
