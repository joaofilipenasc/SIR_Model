Modelo SIR para simulação epidemiológica na linguagem C

Aluno:
João Filipe do Nascimento e Silva

Para gerar o "arquivo SIR_Model_final.csv", deve-se preencher o arquivo "SIR_Model_parametros.csv" com os seguintes dados, respeitando os espaçamentos e as vírgulas:
S,I,R,h,Nb,Tb,Sb,Ib,mk,nk,Tk,t

Comando para compilar e rodar o programa no terminal do Linux:
gcc -o main  main.c  -lm
./main SIR_Model_parametros.csv
