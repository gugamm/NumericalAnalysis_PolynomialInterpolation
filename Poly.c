/**************************************
* NOME : GUSTAVO MARQUES MARTINS
* MATRICULA : 1310630
* LAB 6 : INTERPOLACAO DE POLINOMIOS
* ULTIMA MODIFICAÇÃO : 11/10/2015
**************************************/

#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Poly.h"


/***************************
* Internal type definition *
****************************/

struct sample {
	int n;     /* Numero de amostras */
	double *x; /* Valores x das amostras */
	double *y; /* Valores y das amostras */
};

/*********************************
* INTERNAL FUNCTIONS DEFINITIONS *
**********************************/

static double *calculateLagrangeNumerators(Sample *s, double x);

/**********************
* Calcula f[x0...xn]
* a precisa ser 0
* b >= a
***********************/
static double calculateNewtonB(Sample *s, int a, int b);

/**************************************
* INTERFACE FUNCTIONS IMPLEMENTATIONS *
***************************************/

/****************************************************************
* Determinar as n amostras para a aproximação de uma função f
* dentro do intervalo [a,b]
* Exemplo : para n = 1(1 amostra)
* (a+b)/2 + (b-a)/2 * cos( PI/2 )
*****************************************************************/
Sample *chebyshev(int n, double a, double b, double(*f) (double x)) {
	Sample *s;
	int i;

	//CreateSample
	s = plCreateSample(n);
	if (s == NULL)
		return NULL;

	for (i = 0; i < n; i++) {
		s->x[i] = (a + b) / 2.0 + (b - a) / 2.0 * cos((2 * i + 1) * M_PI / (2 * (n-1) + 2));
		s->y[i] = f(s->x[i]);
	}

	return s;
}

/****************************************************************************
* A funcao recebe um conjunto de amostras e deve retornar os denominadores
* da funcao Lk, ou seja (xi - xj) onde xj != xi
*****************************************************************************/
double *lagrangeCompute(Sample *s) {
	double *den;
	int n, i, j;
	double *x;
	double *y;
	double value = 1.0;

	//inicia variaveis
	n = plGetN(s);
	x = plGetX(s);
	y = plGetY(s);

	//Cria o vetor de denominador
	den = (double *)calloc(n,sizeof(double));
	if (den == NULL)
		return NULL;

	for (i = 0; i < n; i++) { //Li(X)
		for (j = 0; j < n; j++) { //II (xi - xj) -> i != j
			if (j == i)
				continue;
			value *= (x[i] - x[j]);
		}
		den[i] = value;
		value = 1.0;
	}
	return den;
}

/***********************************************************************************
* Essa funcao deve avaliar o polinomio interpolante de Lagrange em um ponto x dado
************************************************************************************/
double lagrangeEval(Sample *s, double *den, double x) {
	double value = 0.0;
	int i;
	int n = plGetN(s);
	double *y = plGetY(s);
	double *numerators = calculateLagrangeNumerators(s, x);

	if (numerators == NULL) {
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < n; i++) {
		value += y[i] * (numerators[i] / den[i]);
	}
	
	free(numerators);
	return value;
}

/****************************************************************
* Essa funcao calcula os fatores "b"s do polinomio de newton
* Para isso, a funcao utiliza-ra os valores das amostras para calcular
* Bk = f[x0...xk]
* A funcao retorna um vetor alocado dinamicamente
****************************************************************/
double *newtonCompute(Sample *s) {
	double *bs;
	int n, i, j;

	n = plGetN(s);
	bs = (double *)calloc(n, sizeof(double));

	if (bs == NULL)
		return NULL;

	for (i = 0; i < n; i++) {
		bs[i] = calculateNewtonB(s, 0, i);
	}
	return bs;
}

/*************************************************************
* Essa funcao deve avaliar o polinomio interpolante de Newton
* em um ponto x dado
*************************************************************/
double newtonEval(Sample *s, double *coef, double x) {
	double value = 0.0;
	int n, i, coefCount, j;
	double *x_vet;
	double coefValue = 1.0;

	n = plGetN(s);
	x_vet = plGetX(s);
	coefCount = 0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < coefCount; j++) {
			coefValue *= (x - x_vet[j]); //Calculates (x - x0)...(x-xn-2)
		}
		value += coef[i] * coefValue;
		coefValue = 1.0;
		coefCount++;
	}
	return value;
}

double *plGetX(Sample *s) {
	return s->x;
}

double *plGetY(Sample *s) {
	return s->y;
}

int plGetN(Sample *s) {
	return s->n;
}

/*************************************
* INTERNAL FUNCTIONS IMPLEMENTATIONS *
**************************************/

Sample *plCreateSample(int n) {
	Sample *newSample;

	newSample = (Sample *)malloc(sizeof(struct sample) * 1);
	if (newSample == NULL)
		return newSample;

	newSample->x = (double *)calloc(n, sizeof(double));
	if (newSample->x == NULL) {
		free(newSample);
		return NULL;
	}

	newSample->y = (double *)calloc(n, sizeof(double));
	if (newSample->y == NULL) {
		free(newSample->x);
		free(newSample);
		return NULL;
	}

	newSample->n = n;
	return newSample;
}

void plFreeSample(Sample *s) {
	if (s != NULL) {
		free(s->x);
		free(s->y);
		free(s);
	}
}

static double *calculateLagrangeNumerators(Sample *s, double x) {
	double *numerators;
	int i, j;
	int n = plGetN(s);
	double *x_vet = plGetX(s);
	double value = 1.0;

	numerators = (double *)calloc(n, sizeof(double));
	if (numerators == NULL) {
		return NULL;
	}

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				continue;
			value *= (x - x_vet[j]);
		}
		numerators[i] = value;
		value = 1.0;
	}

	return numerators;
}

static double calculateNewtonB(Sample *s, int a, int b) {
	double *y = plGetY(s);
	double *x = plGetX(s);
	if (a == b)
		return y[a];

	return (calculateNewtonB(s, a + 1, b) - calculateNewtonB(s, a, b - 1)) / (x[b] - x[a]);
}