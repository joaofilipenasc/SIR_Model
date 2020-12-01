#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "SIR_Model_funcoes.h"

#define DEBUG 0

/**
O Modelo SIR é um modelo compartimental que simplifica a modelagem matemática de doenças infecciosas. 
A população é atribuída a compartimentos com etiquetas - por exemplo, S, I ou R - e elas podem progredir entre compartimentos. 
Um dos modelos mais simples e populares de SIR é o Kermack-McKendrick model. 

fonte: https://mathworld.wolfram.com/SIRModel.html

LEGENDA

S:número de indivíduos suscetíveis (que ainda não estão contaminados);
I:número de indivíduos infectados (capazes de infectar indivíduos S);
R:número de indivíduos removidos (que se recuperaram, tornaram-se imunes ou faleceram).
h:pequeno intervalo de tempo (em horas);
b:facilidade de contágio de um indivíduo;
k:probabilidade que um indivíduo se recupere;
t:instantes de tempo nos quais o modelo é simulado (em horas);
Nb:número de pessoas suscetíveis que se infectaram em um intervalo de tempo Tb;
Tb:intervalo de tempo de indivíduos suscetíveis que se infectaram;
Sb:número de pessoas suscetíveis no início da observação;
Ib:número de pessoas infectadas no início da observação;
mk:quantos indivíduos se recuperaram; 
nk:total de indivíduos;
Tk: intervalo de tempo de indivíduos recuperados;
S0:número de indivíduos suscetíveis no início da análise;
I0:número de indivíduos infectados no início da análise;
R0:número de indivíduos recuperados no início da análise;
*/

/** 
 * Parâmetros básicos
*/
double b;
double k;
double h;
double Nb, Tb, Sb, Ib;
double mk, nk, Tk;
double S0;
double I0;
double R0;

/** 
 * Variáveis
*/
double t, S, I, R, Populacao[3];
double dPopulacao[3];

/**
 * Função com as equações de S, I e R
*/
void Calc(double Populacao[3])
{
	double tempS, tempI, tempR;

	tempS = Populacao[0];
	tempI = Populacao[1];
	tempR = Populacao[2];

	dPopulacao[0] = -b * tempS * tempI;
	dPopulacao[1] = b * tempS * tempI - k * tempI;
	dPopulacao[2] = k * tempI;

	return;
}

/**
 * Função Principal
*/
int main()
{
	/**
	* Abertura do arquivo SIR_Model_parametros.csv e leitura dos dados
	*/
	FILE *arquivo;
	arquivo = fopen("SIR_Model_parametros.csv", "r");

	if (arquivo == NULL)
	{
		printf("Problemas na criação do arquivo!\n");
		return 0;
	}
	else
	{
		fscanf(arquivo, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &S0, &I0, &R0, &h, &Nb, &Tb, &Sb, &Ib, &mk, &nk, &Tk, &t);
	}

	double step, Every;
  	
  	/**
  	* Equações do b e k, importantes para o cálculo das equações principais (linhas 66, 67 e 68)
  	*/
  	b = (Nb)/(Tb*Sb*Ib);
  	k = (mk)/(nk*Tk);

  	/**
 	* Verificar os dados puxados do arquivo SIR_Model_parametros.csv para o programa
	*/
	checkpoint();

	S = S0;
	I = I0;
	R = 0;

	/**
	 * Adequação da escala de tempo caso os valores de h e t fornecidos no arquivo SIR_Model_parametros não estejam em horas, para HORAS
	*/
	step = 0.01 / ((b + k) * S0);

	Every = pow(10, floor(log10((1.0 / ((b + k) * S0)))));
	while (h / Every > 10000)
	{
		Every *= 10.0;
	}

	if (DEBUG)
	{
		printf("Usando um tempo provisório %lf e saída de dados %lf", step, Every);
	}

	/**
	* Arquivo de saída
	*/

	if ((arquivo = fopen("SIR_Model_final.csv", "w")) == NULL)
	{
		printf("O arquivo não pode ser aberto para que seja escrito\n");
		exit(1);
	}

	/** 
	 * Iteração principal
	*/
	t = 0;
   
	saida(arquivo);

	do
	{
		Runge_Kutta(step);
		t += step;

		if (floor(t / Every) > floor((t - step) / Every))
		{
			saida(arquivo);
		}
		else
		{
			if (DEBUG < 0)
			{
				saida(stdout);
			}
		}
	} while (t < h);

	saida(arquivo);
	fclose(arquivo);
}

/**
 * Funções
*/
void leitura(FILE *arquivo)
{
	char str[200];
	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &S0);

	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &I0);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &R0);

	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &h);

	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &Nb);

	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &Tb);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &Sb);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &Ib);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &mk);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &nk);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &Tk);

  	fscanf(arquivo, "%s", str);
	fscanf(arquivo, "%lf", &t);

	fclose(arquivo);
}

void checkpoint()
{
	if (S0 <= 0)
	{
		printf("O número de suscetíveis (%.lf) é menor ou igual a zero\n", S0);
		exit(1);
	}

	if (I0 <= 0)
	{
		printf("O número inicial de infectados (%.lf) é menor ou igual a zero\n", I0);
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
		printf("O nível inicial de suspeitos e infectados (%.lf + %.lf = %.lf) é maior do que 1\n", S0, I0, S0 + I0);
	}

	if (b < k)
	{
		printf("Taxa de reprodução básica (R_0 = %lf) é menor do que 1\n", b / k);
	}
}

void saida(FILE *arquivo)
{
	if (arquivo != stdout)
	{
		fprintf(arquivo, "%.lf,%.lf,%.lf,%.lf\n", S, I, R, t);
	}

	if (DEBUG)
	{
		printf("%.lf,%.lf,%.lf,%.lf\n", S, I, R, t);
	}
}

/**
* Utilização do Método de Runge–Kutta como método de cálculo iterativo entre os arrays de população
*/

void Runge_Kutta(double step)
{
	int i;
	double dPopulacao1[3], dPopulacao2[3], dPopulacao3[3], dPopulacao4[3];
	double tempPopulacao[3], inicialPopulacao[3];

	inicialPopulacao[0] = S;
	inicialPopulacao[1] = I;
	inicialPopulacao[2] = R;

	Calc(inicialPopulacao);

	for (i = 0; i < 3; i++)
	{
		dPopulacao1[i] = dPopulacao[i];
		tempPopulacao[i] = inicialPopulacao[i] + step * dPopulacao1[i] / 2;
	}
	Calc(tempPopulacao);
	for (i = 0; i < 3; i++)
	{
		dPopulacao2[i] = dPopulacao[i];
		tempPopulacao[i] = inicialPopulacao[i] + step * dPopulacao2[i] / 2;
	}
	Calc(tempPopulacao);
	for (i = 0; i < 3; i++)
	{
		dPopulacao3[i] = dPopulacao[i];
		tempPopulacao[i] = inicialPopulacao[i] + step * dPopulacao3[i];
	}
	Calc(tempPopulacao);
	for (i = 0; i < 3; i++)
	{
		dPopulacao4[i] = dPopulacao[i];
		tempPopulacao[i] = inicialPopulacao[i] + (dPopulacao1[i] / 6 + dPopulacao2[i] / 3 + dPopulacao3[i] / 3 + dPopulacao4[i] / 6) * step;
	}

	S = tempPopulacao[0];
	I = tempPopulacao[1];
	R = tempPopulacao[2];

	return;
}
