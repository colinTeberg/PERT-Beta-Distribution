/*
	Author: Colin Teberg
	Contact: cteberg@nd.edu, cdt7654@truman.edu, cteberg@gmail.com
	Created: October 2nd, 2012 @ Truman State University
	Purpose: Generate a beta distribution for PERT. Created to be used with
	the Michael Lab Lymphatic Filarisis simulation.
	
	Copyright (C) 2012 Colin Teberg

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
	
*/
#ifndef PERT_H_INCLUDED
#define PERT_H_INCLUDED

#define N 1000000.0 //accuracy of the integral

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

typedef struct {
	unsigned int alpha; //parameter
	unsigned int beta; //parameter
	
	double pessi; //the worst case scenario
	double opti; //the best case scenario
	double interval; //step between a and b
	
	double *pdf; // E(x)
	double *cdf;
}pert;

/*
	Allocates memory for each of the arrays in pert based on size
*/
void initPert(unsigned int alpha, unsigned int beta, double interval,
							double pessi, double opti, pert *p);

/*
	Frees allocated memory in pert
*/
void freePert(pert *p);

/*
	Calculate the pdf of the distribution 
*/
void fillDist(pert *p);


/*
	Calculate the integral of the beta distribution
*/
double integrate(pert *p, double a, double b);

/*
	Calculate the factorial of some number minus 1, the gamma function
*/
double gamma(double num);


/*
	Prints the values of the pdf to the screen
*/
void dumpPdf(pert *p);

/*
	Prints the values of the cdf to the screen
*/
void dumpCdf(pert *p);

#endif