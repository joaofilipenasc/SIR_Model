Modelo SIR para simulação epidemiológica na linguagem C

Aluno:
João Filipe do Nascimento e Silva

Para gerar o arquivo SIR_Model_final.csv, deve-se preencher o arquivo "SIR_Model_parametros.csv" com os seguintes dados, respeitando os espaçamentos e as vírgulas:
S,I,R,h,Nb,Tb,Sb,Ib,mk,nk,Tk,t

Observação: um arquivo .csv com esses dados preenchidos, chamado de "SIR_Model_parametros.csv", já se encontra na pasta do projeto.

Comando para compilar e rodar o programa no terminal do Linux:
gcc -o main  main.c  -lm
./main SIR_Model_parametros.csv

Observação:
O programa foi feito no site Repl.it, portanto os arquivos com os códigos também podem ser simulados nessa plataforma.
Link: https://repl.it/join/raewqjam-joaofilipenasc

O que pode ser feito no Checkpoint 2?
- Separar as funções em outro arquivo, modularizando ainda mais o projeto.
- Mudar os parâmetros dos dados fornecidos no arquivo SIR_Model_parametros.csv, gerando assim gráficos a partir de diferentes cenários.

Foi plotado um gráfico no calc libreoffice apenas teste e ele se chama "SIR_Model_grafico.jpg".

