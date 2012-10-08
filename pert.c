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

#include "pert.h"

/*
	Allocates memory for each of the arrays in pert based on size
*/

void initPert(unsigned int alpha, unsigned int beta, double interval,
							double pessi, double opti, pert *p)
{
	p->alpha = alpha;
	p->beta = beta;
	p->interval = interval;
	p->pessi = pessi;
	p->opti = opti;
	
	//memory allocation
	p->pdf = (double*)malloc((opti - pessi)/interval * sizeof(double));
	p->cdf = (double*)malloc((opti - pessi)/interval * sizeof(double));
	
	//now fill the pert distribution
	fillDist(p);
}

/*
	Frees allocated memory in pert
*/
void freePert(pert *p)
{
	free(p->pdf);
	free(p->cdf);
}

/*
	This is where the meat and potatoes will be, TBD	
*/
void fillDist(pert *p)
{
	double prev_cdf = 0; //keep track of the previous pdf
	for(unsigned int i = 0; i < (p->opti- p->pessi)/p->interval; i++) {
		double a = p->pessi + (p->interval * i);
		double b = p->pessi + (p->interval * (i+1));
		p->pdf[i] = integrate(p, a, b);
		p->cdf[i] = p->pdf[i] + prev_cdf;
		prev_cdf += p->pdf[i];
	}
}

/*
	Calculate the integral of the beta distribution
*/
double integrate(pert *p, double a, double b)
{
	double sum = 0;

	for(double x = a; x < b; x += (b - a)/N) {
		double y = (1.0/(p->opti - p->pessi));
		y *= gamma(p->alpha + p->beta)/(gamma(p->alpha) * gamma(p->beta));
		y *= pow( (x - p->pessi)/(p->opti - p->pessi), p->alpha - 1);
		y *= pow( (p->opti - x)/(p->opti - p->pessi), p->beta - 1);
		sum += y * ((b - a) / N);
	}

	return sum;
}

/*
	Calculate the factorial of some number minus 1, the gamma function
*/
double gamma(double number) {
	double sum = 1;
	
	for(unsigned int i = 1; i < number; i++) {
		sum *= i;
	}
	
	return sum;
}

/*
	Prints the values of the pdf to the screen
*/
void dumpPdf(pert *p)
{
	double sum = 0;
	
	for(unsigned int i = 0; i < (p->opti- p->pessi)/p->interval; i++) {
		double a = p->pessi + (p->interval * i);
		double b = p->pessi + (p->interval * (i+1));
		fprintf(stderr, "%lf to %lf is %lf\n", a, b, p->pdf[i]);
		sum += p->pdf[i];
	}
	
	fprintf(stderr, "Which sums to: %lf\n", sum);
}


/*
	Prints the values of the cdf to the screen
*/
void dumpCdf(pert *p)
{
	double sum = 0;
	
	for(unsigned int i = 0; i < (p->opti- p->pessi)/p->interval; i++) {
		double a = p->pessi + (p->interval * i);
		double b = p->pessi + (p->interval * (i+1));
		fprintf(stderr, "%lf to %lf is %lf\n", a, b, p->cdf[i]);
		sum += p->pdf[i];
	}
	
	fprintf(stderr, "Which sums to: %lf\n", sum);
}