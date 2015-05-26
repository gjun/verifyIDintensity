/*
 *   Author: Goo Jun (goo.jun@uth.tmc.edu)
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __COMMON_H__
#define __COMMON_H__

#include <math.h>

const double pi=3.14159265;
double min_af;

double det(double *M)
{
	return (M[0]*M[3]-M[1]*M[2]);
}

double normpdf2(float nA, float nB, double m0, double m1, double* cov)
{
	double val,  D, prc[4];
	double y[2];
	y[0] = (double)nA;
	y[1] = (double)nB;

	D = det(cov);
	prc[0] = cov[3]/D;
	prc[3] = cov[0]/D;
	prc[1] = prc[2] = -1.0*cov[1]/D;

	val= -0.5*((y[0]-m0)*((y[0]-m0)*prc[0] + (y[1]-m1)*prc[1]) + (y[1]-m1)*((y[0]-m0)*prc[2]+(y[1]-m1)*prc[3]));
	val= 0.5*exp(val)/(pi*sqrt(D));

	return val;
}

double lognormpdf2(float nA, float nB, double m0, double m1, double* cov)
{
	double val,  D, prc[4];
	double y[2];
	y[0] = (double)nA;
	y[1] = (double)nB;


	D = det(cov);
	prc[0] = cov[3]/D;
	prc[3] = cov[0]/D;
	prc[1] = prc[2] = -1.0*cov[1]/D;

	val = -0.5*((y[0]-m0)*((y[0]-m0)*prc[0] + (y[1]-m1)*prc[1]) + (y[1]-m1)*((y[0]-m0)*prc[2]+(y[1]-m1)*prc[3]));
	val = val-log(2*pi)-0.5*log(D);

	return val;
}


#endif
