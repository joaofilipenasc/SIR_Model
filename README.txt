Modelo SIR para simulação epidemiológica na linguagem C

Aluno:
João Filipe do Nascimento e Silva

Para gerar o arquivo SIR_Model.csv, deve-se preencher o arquivo "parametros.csv" com os seguintes dados, respeitando os espaçamentos e as vírgulas:
b,k,S,I,t

LEGENDA
b: facilidade de contágio de um indivíduo;
k: probabilidade que um indivíduo se recupere; 
S: número de indivíduos suscetíveis (que ainda não estão contaminados);
I: número de indivíduos infectados (capazes de infectar indivíduos S);
t: instantes de tempo nos quais o modelo é simulado (em horas).

Comando de compilar e rodar o programa no Linux:
make all