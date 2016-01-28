/**************************************
* NOME : GUSTAVO MARQUES MARTINS
* LAB 6 : INTERPOLACAO DE POLINOMIOS
* ULTIMA MODIFICAÇÃO : 11/10/2015
**************************************/

typedef struct sample Sample;

/******************************************************************************
*  Create a sample. Alocate memory for arrays for xs and ys. In case of fail,
*  A null pointer is returned
******************************************************************************/
Sample *plCreateSample(int n);

double *plGetX(Sample *s);

double *plGetY(Sample *s);

int plGetN(Sample *s);

/*Free a sample*/
void plFreeSample(Sample *s);

/****************************************************************
* Determinar as n amostras para a aproximação de uma função f
* dentro do intervalo [a,b]
*****************************************************************/
Sample *chebyshev(int n, double a, double b, double(*f) (double x));

/****************************************************************************
* A funcao recebe um conjunto de amostras e deve retornar os denominadores
* da funcao Lk, ou seja (xi - xj) onde xj != xi
*****************************************************************************/
double *lagrangeCompute(Sample *s);


/***************************************************************
* Essa funcao deve avaliar o polinomio interpolante de Lagrange 
* em um ponto x dado
****************************************************************/
double lagrangeEval(Sample *s, double *den, double x);

/****************************************************************
* Essa funcao calcula os fatores "b"s do polinomio de newton
* Para isso, a funcao utiliza-ra os valores das amostras para calcular
* Bk = f[x0...xk]
* A funcao retorna um vetor alocado dinamicamente
****************************************************************/
double *newtonCompute(Sample *s);

/*************************************************************
* Essa funcao deve avaliar o polinomio interpolante de Newton
* em um ponto x dado
*************************************************************/
double newtonEval(Sample *s, double *coef, double x);
