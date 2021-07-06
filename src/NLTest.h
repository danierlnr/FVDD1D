/*
 * NodoVF.h
 *
 *  Created on: Nov 20, 2013
 *      Author: daniel
 */

#ifndef NLTEST_H_
#define NLTEST_H_

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

#include "NLSolver.h"

class NLTest;
template class NLSolver<NLTest>;

class NLTest
{
public:
	NLTest();
	virtual ~NLTest();

    void F(const DVector& x, DVector& F);
    int J(DVector& x, SPRMatrix& J);
    int J(DVector& x, BandMatrix&  J);
    void operator()(const DVector& x, DVector& f);

protected:
    const double pi;

};

#endif /* NLTEST_H_ */
