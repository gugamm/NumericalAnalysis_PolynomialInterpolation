/**************************************
* NOME : GUSTAVO MARQUES MARTINS
* MATRICULA : 1310630
* LAB 6 : INTERPOLACAO DE POLINOMIOS
* ULTIMA MODIFICAÇÃO : 11/10/2015
**************************************/

#include <conio.h>
#include <stdio.h>
#define _USE_MATH_DEFINES //Para usar o valor de PI --> M_PI
#include <math.h>
#include <windows.h>
#include "Poly.h"

static void printRed(const char *text);
static void questaoABC();

int main(void) {
	questaoABC();
	printf("Press any key to exit...");
	_getch();
	return 0;
}

static void printRed(const char *text) {
	HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_SCREEN_BUFFER_INFO consoleInfo;
	WORD saved_attributes;

	/* Save current attributes */
	GetConsoleScreenBufferInfo(hConsole, &consoleInfo);
	saved_attributes = consoleInfo.wAttributes;

	SetConsoleTextAttribute(hConsole, FOREGROUND_RED);
	printf("%s",text);

	/* Restore original attributes */
	SetConsoleTextAttribute(hConsole, saved_attributes);
}

static void questaoABC() {
	Sample *s;
	int n, i;
	double *x;
	double *y;
	double *den;
	double *coefs;
	double points[] = { 0.0,0.5,0.2,0.3,0.7,0.9,0.987,1.0,1.12,1.3 };
	int pointCount = 10;
	double lagrange, newton, dif;
	/*QUESTAO A*/
	s = chebyshev(6, 0, M_PI / 2, cos);
	n = plGetN(s);
	x = plGetX(s);
	y = plGetY(s);
	printRed("====CHEBYSHEV VALUES====\n");
	for (i = 0; i < n; i++) {
		printf("X[%d] : %f\n"
			"Y[%d] : %f\n"
			"\n", i, x[i], i, y[i]);
	}
	/*QUESTAO B*/
	printRed("====CALCULANDO POLINOMIO INTERPOLANTE LAGRANGE====\n");
	den = lagrangeCompute(s);
	for (i = 0; i < n; i++) {
		printf("den[%d] : %f\n", i, den[i]);
	}
	printf("\n");
	/*QUESTAO C*/
	printRed("====CALCULANDO POLINOMIO INTERPOLANTE NEWTON====\n");
	coefs = newtonCompute(s);
	for (i = 0; i < n; i++) {
		printf("coefs[%d] : %f\n", i, coefs[i]);
	}
	/*QUESTAO D*/
	printRed("\n====COMPARANDO VALORES DOS DOIS POLINOMIOS====\n");
	for (i = 0; i < pointCount; i++) {
		lagrange = lagrangeEval(s, den, points[i]);
		newton = newtonEval(s, coefs, points[i]);
		dif = lagrange - newton;
		printf("Ponto    : %f\n"
			"Lagrange : %.20f\n"
			"Newton   : %.20f\n" , points[i], lagrange, newton);
		if (dif) {
			printf("Sao diferentes\n");
			printf("Diferenca(L - N) : %.20f\n", dif);
		}
		else {
			printf("Sao iguais\n");
		}
		printf("\n");
	}
	printf("\n");
	free(den);
	free(coefs);
	plFreeSample(s);
}
