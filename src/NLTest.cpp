
#include "NLTest.h"

NLTest::NLTest():
pi(3.14159)
{
    //NLSolver nls;

    DVector x0(3);
    
    x0.set(0, 0.1);
    x0.set(1, 0.1);
    x0.set(2,-0.2);

    bool check;
    
    //nls.newt4(x0,check,this);
    //nls.mnewt4(x0,check,this);
    NLSolver<NLTest> nls(this);
    //nls.newt(x,check,funcv);
    nls.newt(x0,check);
    

    DVector f(3);
    
    std::cout << "==========" << std::endl;
    F(x0,f);
    //x0.show();
    //f.show();
    std::cout << "x:"          << std::endl;    
    std::cout << x0.to_string() << std::endl;    
    std::cout << "f:"          << std::endl;    
    std::cout << f.to_string() << std::endl;    

}


NLTest::~NLTest()
{

}


void NLTest::operator()(const DVector& x, DVector& F)
{
    double x0 = x.get(0);
    double x1 = x.get(1);
    double x2 = x.get(2);
    
    double F0 = 3*x0 - cos(x1*x2) - 0.5;
    double F1 = x0*x0 - 81*(x1+0.1)*(x1+0.1) + sin(x2) + 1.06;
    double F2 = exp(-x0*x1) + 20*x2 + (10*pi-3.0)/2.0;
    
    F.set(0,F0);
    F.set(1,F1);
    F.set(2,F2);
    
}

void NLTest::F(const DVector& x, DVector& F)
{
    double x0 = x.get(0);
    double x1 = x.get(1);
    double x2 = x.get(2);
    
    double F0 = 3*x0 - cos(x1*x2) - 0.5;
    double F1 = x0*x0 - 81*(x1+0.1)*(x1+0.1) + sin(x2) + 1.06;
    double F2 = exp(-x0*x1) + 20*x2 + (10*pi-3.0)/3.0;
    
    F.set(0,F0);
    F.set(1,F1);
    F.set(2,F2);
    
}

int NLTest::J(DVector& x, SPRMatrix& J)
{
    double x0 = x.get(0);
    double x1 = x.get(1);
    double x2 = x.get(2);
    
    //double F0 = 3*x0 - cos(x1*x2) - 0.5;
    double J00 = 3;
    double J01 = x2*sin(x1*x2);
    double J02 = x1*sin(x1*x2);


    //double F1 = x0*x0 - 81*(x1+0.1)*(x1+0.1) + sin(x2) + 1.06;
    double J10 =  2*x0;
    double J11 = -162*(x1+0.1);
    double J12 =  cos(x2);
    
    
    //double F2 = exp(-x0*x1) + 20*x2 + (10*pi-3.0)/2.0;
    double J20 =  -x1*exp(-x0*x1);
    double J21 =  -x0*exp(-x0*x1);
    double J22 =   20;
    
    J.set(0,0,J00);
    J.set(0,1,J01);
    J.set(0,2,J02);

    J.set(1,0,J10);
    J.set(1,1,J11);
    J.set(1,2,J12);

    J.set(2,0,J20);
    J.set(2,1,J21);
    J.set(2,2,J22);

    return 1;
}

int NLTest::J(DVector& x, BandMatrix& J)
{
    double x0 = x.get(0);
    double x1 = x.get(1);
    double x2 = x.get(2);
    
    //double F0 = 3*x0 - cos(x1*x2) - 0.5;
    double J00 = 3;
    double J01 = x2*sin(x1*x2);
    double J02 = x1*sin(x1*x2);


    //double F1 = x0*x0 - 81*(x1+0.1)*(x1+0.1) + sin(x2) + 1.06;
    double J10 =  2*x0;
    double J11 = -162*(x1+0.1);
    double J12 =  cos(x2);
    
    
    //double F2 = exp(-x0*x1) + 20*x2 + (10*pi-3.0)/2.0;
    double J20 =  -x1*exp(-x0*x1);
    double J21 =  -x0*exp(-x0*x1);
    double J22 =   20;
    
    J.set(0,0,J00);
    J.set(0,1,J01);
    J.set(0,2,J02);

    J.set(1,0,J10);
    J.set(1,1,J11);
    J.set(1,2,J12);

    J.set(2,0,J20);
    J.set(2,1,J21);
    J.set(2,2,J22);

    return 1;
}

