#ifndef NLSOLVER_H
#define NLSOLVER_H

#include <cmath>
#include <limits>
#include "nr.h"
#include "Matrix.h"
#include "LSolver.h"
//#include "NLSolverImp.h"

//using namespace std;

template <class T>
class NLSolver
{
public:
	NLSolver(T* vecfunc):m_vecfunc(vecfunc){}
	virtual ~NLSolver(){}
	//void newt(Vec_IO_DP &x, bool &check);
	void newt(DVector &x, bool &check)
    {
        const int MAXITS=8;
        const double TOLF=1.0e-8,TOLMIN=1.0e-12,STPMX=100.0;
        const double TOLX=numeric_limits<DP>::epsilon();
        int i,its;
        DP den,f,fold,stpmax,sum,temp,test;

        int n=x.size();
        Vec_INT indx(n);
        //Vec_DP g(n),p(n),xold(n);
        DVector g(n),p(n),xold(n);

        //fvec_p=new Vec_DP(n);
        fvec_p=new DVector(n);

        //Vec_DP &fvec=*fvec_p;
        DVector &fvec=*fvec_p;



        f=fmin(x);
        test=0.0;
        for (i=0;i<n;i++)
            if (fabs(fvec[i]) > test) 
                test=fabs(fvec[i]);

        if (test < 0.01*TOLF) 
        {
            check=false;
            delete fvec_p;
            return;
        }
        
        //for(int ii=0;ii<n;ii++)
        //    std::cout << fvec[ii] << std::endl;
        sum=0.0;
        for (i=0;i<n;i++) 
            sum += SQR(x[i]);

        stpmax=STPMX*MAX(sqrt(sum),DP(n));

        for (its=0;its<MAXITS;its++) 
        {
            std::cout << "newt it : " << its << std::endl;
            //for(int ii=0;ii<n;ii++)
            //    std::cout << fvec[ii] << std::endl;

            BandMatrix fjac(n);
            //SPRMatrix fjac(n);
            //fdjac(x,fvec,fjac,vecfunc);
            
            //fdjac(x,fvec,fjac);
            m_vecfunc->J(x,fjac);

            //std::cout << "Jacobian" << std::endl;
            //std::cout << fjac.to_string() << std::endl;
            // Simple matriz*vector
            /*
            for (i=0;i<n;i++) 
            {
                sum=0.0;
                for (j=0;j<n;j++) 
                    //sum += fjac[j][i]*fvec[j];
                    sum += fjac.get(j,i)*fvec[j];
                g[i]=sum;
            }
            */ 
            fjac.compress();
            fjac.gaxpyt(fvec,g);
            
            //for (i=0;i<n;i++) 
            //    xold[i]=x[i];
            xold=x;    
                
            fold=f;
            
            for (i=0;i<n;i++) 
                p[i] = -fvec[i];
            //p=-1*fvec;

            LSolver ls;
            DVector p2(p.size());
            for (size_t ii=0;ii<p.size();ii++)
                p2[ii]=p[ii]; 
            
            //DVector xx = ls.solve(fjac,p2);
            ls.solve(fjac,p2,x);
            
            //std::cout << fjac.to_string() << std::endl;
            
            for (size_t ii=0;ii<p.size();ii++)
            {
                p[ii]=x[ii]; 
            }
                
            //std::cout << "fjac:"          << std::endl;    
            //std::cout << fjac.to_string() << std::endl;    
            std::cout << "p2:"            << std::endl;    
            std::cout << p2.to_string()   << std::endl;    
            std::cout << "x:"             << std::endl;    
            std::cout << x.to_string()    << std::endl;   
             
            //lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin);
            
            lnsrch(xold,fold,g,p,x,f,stpmax,check);
            /*
            std::cout << ">x:"             << std::endl;    		
            for (int ii=0;ii<p.size();ii++)
            {
                std::cout << x[ii]   << std::endl;    
            }
            */

            test=0.0;
            //std::cout << "**1:"             << std::endl;    

            for (i=0;i<n;i++)
                if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
            if (test < TOLF) {
            //std::cout << "**1a:"             << std::endl;    
                check=false;
                delete fvec_p;
                return;
            }
            if (check) {
                test=0.0;
                den=MAX(f,0.5*n);
                for (i=0;i<n;i++) {
                    temp=fabs(g[i])*MAX(fabs(x[i]),1.0)/den;
                    if (temp > test) test=temp;
                }
                //std::cout << "**1b:"             << std::endl;    
                check=(test < TOLMIN);
                delete fvec_p;
                return;
            }
            //std::cout << "**2:"             << std::endl;    
            test=0.0;
            for (i=0;i<n;i++) 
            {
                std::cout << "i:" << i << std::endl;    
                temp=(fabs(x[i]-xold[i]))/MAX(fabs(x[i]),1.0);
                if (temp > test) 
                    test=temp;
            }
            if (test < TOLX) {
                delete fvec_p;
                return;
            }
        }
        //std::cout << "**3:"             << std::endl;    	
        nrerror("MAXITS exceeded in newt");
    }
    
	
private:
    T* m_vecfunc;
	//Vec_DP *fvec_p;
	DVector *fvec_p;
	//void (*nrfuncv)(Vec_I_DP &v, Vec_O_DP &f);
	void (*nrfuncv)(const DVector &v, DVector &f);

    void fdjac(const DVector& x, DVector& fvec, BandMatrix& df)
    {
        //std::cout << "  fdjac...."             << std::endl;    	
        const DP EPS=1.0e-8;
        int i,j;
        DP h,temp;

        int n=x.size();
        //Vec_DP f(n);
        DVector f(n);

        for (j=0;j<n;j++) 
        {
            //std::cout << "      j = " << j << "(" << n << ")" << std::endl;    	
            temp=x[j];
            h=EPS*fabs(temp);
            if (h == 0.0) 
                h=EPS;
            x[j]=temp+h;
            h=x[j]-temp;
            //vecfunc(x,f);
            m_vecfunc->F(x,f);
            x[j]=temp;
            for (i=0;i<n;i++){
                //df.set((f[i]-fvec[i])/h,i,j);
                df.set(i,j,(f[i]-fvec[i])/h);
            }
        }
        
        std::cout << " Compress..." << std::endl;    	
        df.compress();
    }

	void fdjac(const DVector &x, DVector &fvec, SPRMatrix &df)
    {
        const DP EPS=1.0e-8;
        int i,j;
        DP h,temp;

        int n=x.size();
        //Vec_DP f(n);
        DVector f(n);

        for (j=0;j<n;j++) 
        {
            std::cout << "      j = " << j  << std::endl;    	
            temp=x[j];
            h=EPS*fabs(temp);
            if (h == 0.0) 
                h=EPS;
            x[j]=temp+h;
            h=x[j]-temp;
            //vecfunc(x,f);
            m_vecfunc->F(x,f);
            x[j]=temp;
            for (i=0;i<n;i++){
                //df.set((f[i]-fvec[i])/h,i,j);
                df.set(i,j,(f[i]-fvec[i])/h);
                //std::cout << "Jacobian" << std::endl;
                //std::cout << df.to_string() << std::endl;
                //std::cout << (f[i]-fvec[i])/h << std::endl;
            }
        }
        
        df.compress();
    }

	//DP fmin(NRVec<double> &x);
	DP fmin(const DVector &x)
    {
        int i;
        DP sum;

        //Vec_DP &fvec=*fvec_p;
        DVector &fvec=*fvec_p;
        m_vecfunc->F(x,fvec);
        int n=x.size();
        for (sum=0.0,i=0;i<n;i++) sum += SQR(fvec[i]);
        return 0.5*sum;
    }

    void lnsrch(const DVector &xold, const DP fold, DVector &g, DVector &p,DVector &x, DP &f, const DP stpmax, bool &check)
    {
        const double ALF=1.0e-4, TOLX=numeric_limits<DP>::epsilon();
        int i;
        double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
        double rhs1,rhs2,slope,sum,temp,test,tmplam;

        int n=xold.size();
        check=false;
        sum=0.0;
        for (i=0;i<n;i++) sum += p[i]*p[i];
        sum=sqrt(sum);
        if (sum > stpmax)
            for (i=0;i<n;i++) p[i] *= stpmax/sum;
        slope=0.0;
        for (i=0;i<n;i++)
            slope += g[i]*p[i];
        if (slope >= 0.0) nrerror("Roundoff problem in l*nsrch.");
        test=0.0;
        for (i=0;i<n;i++) {
            temp=fabs(p[i])/MAX(fabs(xold[i]),1.0);
            if (temp > test) test=temp;
        }
        alamin=TOLX/test;
        alam=1.0;
        for (;;) {
            for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
            //f=func(x);
            f=fmin(x);
            if (alam < alamin) {
                for (i=0;i<n;i++) x[i]=xold[i];
                check=true;
                return;
            } else if (f <= fold+ALF*alam*slope) return;
            else {
                if (alam == 1.0)
                    tmplam = -slope/(2.0*(f-fold-slope));
                else {
                    rhs1=f-fold-alam*slope;
                    rhs2=f2-fold-alam2*slope;
                    a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                    b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                    if (a == 0.0) tmplam = -slope/(2.0*b);
                    else {
                        disc=b*b-3.0*a*slope;
                        if (disc < 0.0) tmplam=0.5*alam;
                        else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
                        else tmplam=-slope/(b+sqrt(disc));
                    }
                    if (tmplam>0.5*alam)
                        tmplam=0.5*alam;
                }
            }
            alam2=alam;
            f2 = f;
            alam=MAX(tmplam,0.1*alam);
        }
    }
    
	void nrerror(const string error_text)
    {
        cerr << "Numerical Recipes run-time error..." << endl;
        cerr << error_text << endl;
        cerr << "...now exiting to system..." << endl;
    }
    

};

#endif /* NLSOLVER_H */ 
