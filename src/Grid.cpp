/*
 * Grid.cpp
 *
 *  Created on: Sep 3, 2016
 *      Author: daniel
 */


#include "Grid.h"

Grid::Grid()
{

}

Grid::Grid(double xmin,double xmax,size_t n)
{
	double dx = (xmax-xmin)/n;

	for(size_t i = 0; i < n; i++)
	{
		double x1 = xmin + i*dx;
		double x2 = x1 + dx;

		ControlVolume* new_cv;
		if(i == 0)
		{
			new_cv = new ControlVolume(x1,x2,1);
			new_cv->set_bc(1.00);
		}

		else if(i == n-1)
		{
			new_cv = new ControlVolume(x1,x2,2);
			new_cv->set_bc(0.00);
		}
		else
		{
			new_cv = new ControlVolume(x1,x2);
		}

        new_cv->setCVNum(i);

		m_cont_vols.push_back(new_cv);

	}

}

Grid::~Grid()
{

}

void Grid::set_vals(const DVector& v)
{
	for(size_t i = 0; i < m_cont_vols.size(); i++)
	{
		m_cont_vols[i]->set_val(v[i]);
		//std::cout << v[i] << std::endl;
	}

}

void Grid::results(const std::string& filename)
{
	std::ofstream fout(filename.c_str());
	
	//std::cout << "n=" << m_cont_vols.size() << std::endl;

	fout << "import matplotlib.pyplot as plt" << std::endl;
	fout << "x=[]" << std::endl;
	fout << "y=[]" << std::endl;

	//int N = x.size();
	fout << "x.append(" << m_cont_vols[0]->xmin()   << ")" << std::endl;
	fout << "y.append(" << m_cont_vols[0]->get_bc() << ")" << std::endl;

	for(int i = 0; i < m_cont_vols.size(); i++)
	{
		//std::cout << "x = " << m_cont_vols[i]->xnod() << std::endl;
		fout << "x.append(" << m_cont_vols[i]->xc() << ")" << std::endl;
		fout << "y.append(" << m_cont_vols[i]->get_val() << ")" << std::endl;
	}

	fout << "x.append(" << m_cont_vols[m_cont_vols.size()-1]->xmax()   << ")" << std::endl;
	fout << "y.append(" << m_cont_vols[m_cont_vols.size()-1]->get_bc() << ")" << std::endl;

	fout << "plt.xlabel('x')" << std::endl;
	fout << "plt.ylabel('y')" << std::endl;
	fout << "plt.plot(x,y)" << std::endl;
	fout << "plt.show()" << std::endl;

	fout.close();

	fout << "System('python out.py')" << std::endl;

}

//-----------------------------------------------------------

ControlVolume::ControlVolume(double xmin,double xmax,int has_bc)
{
	m_xmin    = xmin;
	m_xmax    = xmax;
	m_xc      = (xmin + xmax)/2.0;
	m_vc      = 0;
	m_has_bc  = has_bc;
	m_bc_xmin = 0;
	m_bc_xmax = 0;

}

ControlVolume::ControlVolume()
{
	m_xmin    = 0;
	m_xmax    = 0;
	m_xc       = 0;
	m_vc      = 0;
	m_has_bc  = 0;
	m_bc_xmin = 0;
	m_bc_xmax = 0;

}

ControlVolume::~ControlVolume()
{

}

void ControlVolume::set_bc(double bc)
{
    if(m_has_bc == 1)
    	m_bc_xmin = bc;
    else if(m_has_bc == 2)
    	m_bc_xmax = bc;
    else // Error
    	;
}

//-------------------------------------------------------

