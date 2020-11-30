#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SIR_Model_funcoes.h"

#define DEBUG 0

/** 
 * Setando parâmetros básicos caso não encontre o arquivo na pasta
*/
double b = 520 / 365;
double k = 1 / 7;
double S0 = 60;
double I0 = 6;
double h = 70;

/** 
 * Setando variáveis e taxas de mudança
*/
double t, S, I, R, Pop[3];
double dPop[3];

/**
 * Inicialização das equações e Método de Runge-Kutta
*/
void Diff(double Pop[3])
{
	double tmpS, tmpI, tmpR;

	tmpS = Pop[0];
	tmpI = Pop[1];
	tmpR = Pop[2];

	dPop[0] = -b * tmpS * tmpI;
	dPop[1] = b * tmpS * tmpI - k * tmpI;
	dPop[2] = k * tmpI;

	return;
}

/**
 * Função Principal
*/
int main()
{
	FILE *arq;
	arq = fopen("parametros.csv", "r");

	if (arq == NULL)
	{
		printf("Problemas na criação do arquivo\n");
		return 0;
	}
	else
	{
		fscanf(arq, "%lf,%lf,%lf,%lf,%lf", &b, &k, &S0, &I0, &h);
	}

	double step, Every;

	checkpoint();

	S = S0;
	I = I0;
	R = 1 - S - I;

	/**
	 * Escala de tempo adequada
	*/
	step = 0.01 / ((b + k) * S0);

	Every = pow(10, floor(log10((1.0 / ((b + k) * S0)))));
	while (h / Every > 10000)
	{
		Every *= 10.0;
	}

	if (DEBUG)
		printf("Usando um tempo provisório %lf and apresentando toda data %lf\n\n", step, Every);

	if ((arq = fopen("SIR_Model.csv", "w")) == NULL)
	{
		printf("O arquivo não pode ser aberto para que seja escrito\n");
		exit(1);
	}

	/** 
	 * Iteração principal
	*/
	t = 0;

	saida(arq);

	do
	{
		Runge_Kutta(step);
		t += step;

		if (floor(t / Every) > floor((t - step) / Every))
		{
			saida(arq);
		}
		else
		{
			if (DEBUG < 0)
			{
				saida(stdout);
			}
		}
	} while (t < h);

	saida(arq);
	fclose(arq);
}

/**
 * Funções
*/
void leitura(FILE *arq)
{
	char str[200];
	fscanf(arq, "%s", str);
	fscanf(arq, "%lf", &b);

	fscanf(arq, "%s", str);
	fscanf(arq, "%lf", &k);

	fscanf(arq, "%s", str);
	fscanf(arq, "%lf", &S0);

	fscanf(arq, "%s", str);
	fscanf(arq, "%lf", &I0);

	fscanf(arq, "%s", str);
	fscanf(arq, "%lf", &h);

	fclose(arq);
}

void checkpoint()
{
	if (S0 <= 0)
	{
		printf("O número de suscetíveis (%lf) é menor ou igual a zero\n", S0);
		exit(1);
	}

	if (I0 <= 0)
	{
		printf("O número inicial de infectados (%lf) é menor ou igual a zero\n", I0);
		exit(1);
	}

	if (b <= 0)
	{
		printf("A taxa de transmissão (%lf) é menor ou igual a zero\n", b);
		exit(1);
	}

	if (k <= 0)
	{
		printf("A taxa de recuperação (%lf) é menor ou igual a zero\n", k);
		exit(1);
	}

	if (h <= 0)
	{
		printf("O intervalo (%lf) é menor ou igual a zero\n", h);
		exit(1);
	}

	if (S0 + I0 > 1)
	{
		printf("O nível inicial de suspeitos e infectados (%lf + %lf = %lf) é maior do que 1\n", S0, I0, S0 + I0);
	}

	if (b < k)
	{
		printf("Taxa de reprodução básica (R_0 = %lf) é menor do que 1\n", b / k);
	}
}

void saida(FILE *arq)
{
	if (arq != stdout)
	{
		fprintf(arq, "%lf,%lf,%lf,%lf\n", S, I, R, t);
	}

	if (DEBUG)
	{
		printf("%lf,%lf,%lf,%lf\n", S, I, R, t);
	}
}

void Runge_Kutta(double step)
{
	int i;
	double dPop1[3], dPop2[3], dPop3[3], dPop4[3];
	double tmpPop[3], inicialPop[3];

	inicialPop[0] = S;
	inicialPop[1] = I;
	inicialPop[2] = R;

	Diff(inicialPop);

	for (i = 0; i < 3; i++)
	{
		dPop1[i] = dPop[i];
		tmpPop[i] = inicialPop[i] + step * dPop1[i] / 2;
	}
	Diff(tmpPop);
	for (i = 0; i < 3; i++)
	{
		dPop2[i] = dPop[i];
		tmpPop[i] = inicialPop[i] + step * dPop2[i] / 2;
	}
	Diff(tmpPop);
	for (i = 0; i < 3; i++)
	{
		dPop3[i] = dPop[i];
		tmpPop[i] = inicialPop[i] + step * dPop3[i];
	}
	Diff(tmpPop);
	for (i = 0; i < 3; i++)
	{
		dPop4[i] = dPop[i];
		tmpPop[i] = inicialPop[i] + (dPop1[i] / 6 + dPop2[i] / 3 + dPop3[i] / 3 + dPop4[i] / 6) * step;
	}

	S = tmpPop[0];
	I = tmpPop[1];
	R = tmpPop[2];

	return;
}
