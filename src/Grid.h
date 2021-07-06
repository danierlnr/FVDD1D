/*
 * Grid.h
 *
 *  Created on: Sep 3, 2016
 *      Author: daniel
 */

#ifndef SRC_GRID_H_
#define SRC_GRID_H_

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "DVector.h"

class ControlVolume
{
public:
	ControlVolume();
	ControlVolume(double, double ,int=0);
	virtual ~ControlVolume();
	void set_bc(double);
	double xmin()const{return m_xmin;}
	double xmax()const{return m_xmax;}
	double xc(){return m_xc;}
    void setCVNum(int n){m_cvn=n;}
    int  getCVNum(){return m_cvn;}
	void set_val(const double val){m_vc = val;}
	double get_val(){return m_vc;}
	double get_bc()
	{
		if(m_has_bc == 1)
			return m_bc_xmin;
		else if(m_has_bc == 2)
			return m_bc_xmax;
	}
protected:
	double m_xmin;
	double m_xmax;
	double m_xc;
	double m_vc;
	int    m_has_bc;
	double m_bc_xmin;
	double m_bc_xmax;
    int m_cvn;
};

class Node
{
public:
	Node();
	Node(double x){m_x=x;}
	virtual ~Node();
protected:
	double m_x;
};


class Grid
{
public:
	Grid();
	Grid(double,double,size_t);
	virtual ~Grid();
	double xmin(size_t i)const{return m_cont_vols[i]->xmin();}
	double xmax(size_t i)const{return m_cont_vols[i]->xmax();}
	double xc(size_t i){return m_cont_vols[i]->xc();}
	size_t size(){return m_cont_vols.size();}
	void set_vals(const DVector&);
	void results(const std::string& nomArch="out.py");


protected:
	std::vector<ControlVolume*> m_cont_vols;
};

#endif /* SRC_GRID_H_ */
