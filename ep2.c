#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>

/*********************** VARIAVEIS GLOBAIS *********************************/
int q;
mpf_t x;
mpf_t precisao;
char criterio;
char opc;
mpf_t cosseno;

mpf_t somaRodada; // idem
int stop;
unsigned long int ind_sum;
int contador_rodadas;

pthread_barrier_t bar;
sem_t mutex;

/*********************** PROTOTIPOS DAS FUNCOES *********************************/
void parserEntrada(int, char**);
void potencia(mpf_t, mpf_t, int);
void fatorial(mpf_t, unsigned long int);
void calculaSequencial(void);
void* threadF(void*);
void* threadM(void*);
void calculaCoeficiente(mpf_t, unsigned long int);

/********************************************************************************/

int main(int argc, char** argv){
	int i;
	mpf_set_default_prec(1000000);
	mpf_init(cosseno);

	parserEntrada(argc, argv);

	if(opc == 's'){
    	calculaSequencial();
	}

	else{
		pthread_t parcela[q];

		if(sem_init(&mutex, 0, 1)){
			fprintf(stderr, "erroR creating semaphore\n");
			exit(0);
		}
		
		if(pthread_barrier_init(&bar ,NULL, q)){
			fprintf(stderr, "erroR initializing barrier\n");
			exit(0);
		}

		ind_sum = 0;
		stop = 0;
		mpf_init(somaRodada);
		
		if(criterio == 'f'){
			for(i = 0; i < q; i++){
			 	if(pthread_create(&parcela[i], NULL, &threadF, (void*) i)){
					fprintf(stderr, "erroR creating threadsF\n");
					exit(0);
			 	}
		 	}
		}

		else{
			for(i = 0; i < q; i++){
			 	if(pthread_create(&parcela[i], NULL, &threadM, (void*) i)){
					fprintf(stderr, "erroR creating threadsM\n");
					exit(0);
			 	}
		 	}	
		}
		
		for(i = 0; i < q; i++){
			if(pthread_join(parcela[i], NULL)){
				fprintf(stderr, "erroR joining threads\n");
				exit(0);	
			}
		}

		printf("\n---------------------------------------------------------\n");
		gmp_printf("Cosseno de %Ff eh igual a: \n%.1000Ff\n", x, cosseno);
		printf("\n# de Rodadas: %ld", ind_sum/q);
		printf("\n---------------------------------------------------------");
	}

	mpf_clear(x);
	mpf_clear(cosseno);
	mpf_clear(precisao);
	pthread_exit(NULL);

	return 0;
}

void* threadF(void* idThread){
	unsigned long int i;
	mpf_t termo;
	mpf_t pot_x;
	mpf_t erro;
	int id = (int) idThread;

	mpf_init(termo);
	mpf_init(pot_x);
	mpf_init_set(erro, precisao);

	while(!stop){
		sem_wait(&mutex);
		i = ind_sum;
		ind_sum++;
		sem_post(&mutex);

		/********** calcula termo  ************/
		potencia(pot_x, x, 2 * i);
		calculaCoeficiente(termo, i);
		mpf_mul(termo, termo, pot_x);
		/**************************************/

		sem_wait(&mutex);
		mpf_add(somaRodada, somaRodada, termo);
		mpf_add(cosseno, cosseno, termo);
		sem_post(&mutex);
		
		if(opc == 'd'){
			sem_wait(&mutex);
			printf("\nNº Rodada: %ld;		Thread_ID: %d", i/q, id);
			sem_post(&mutex);
		}

		/* barreira 1 */
		pthread_barrier_wait(&bar);
		//aqui todos somaram em somaRodada
		
		if(id == 0){
			mpf_abs(somaRodada, somaRodada);

			if(mpf_cmp(somaRodada, precisao) < 0){
				stop = 1;
			}

			mpf_set_ui(somaRodada, 0);

			if(opc == 'd'){
				gmp_printf("\nCosseno de %Ff:\n %.60Ff\n", x, cosseno);
			}
		}

		/* barreira 2 */
		pthread_barrier_wait(&bar);
	}

	return NULL;
}

void* threadM(void* idThread){
	unsigned long int i;
	mpf_t termo;
	mpf_t pot_x;
	int id = (int) idThread;

	mpf_init(termo);
	mpf_init(pot_x);
	
	while(!stop){
		sem_wait(&mutex);
		i = ind_sum;
		ind_sum++;
		sem_post(&mutex);

		/********** calcula termo  ************/
		potencia(pot_x, x, 2 * i);
		calculaCoeficiente(termo, i);
		mpf_mul(termo, termo, pot_x);
		/**************************************/
		
		sem_wait(&mutex);
		mpf_add(cosseno, cosseno, termo);
		sem_post(&mutex);

		mpf_abs(termo, termo);
		
		if(mpf_cmp(termo, precisao) < 0){
			stop = 1; // thread diz: "hora de terminar, termo pequeno"
		}
		
		if(opc == 'd'){
			sem_wait(&mutex);
			printf("\nNº Rodada: %ld;		Thread_ID: %d", i/q, id);
			sem_post(&mutex);
		}

		// TODAS threads chegam na barreira com o termo já calculado
		pthread_barrier_wait(&bar); 

		if(opc == 'd' && id == 0){
			sem_wait(&mutex);
			gmp_printf("\nCosseno de %Ff:\n %.60Ff\n", x, cosseno);
			sem_post(&mutex);
		}
	}

	return NULL;
}

void calculaCoeficiente(mpf_t cf, unsigned long int i){
	mpf_t fact;
	mpf_init(fact);

	fatorial(fact, 2 * i);
	mpf_ui_div(cf, 1, fact);

	if(i % 2 == 1){
		mpf_neg(cf, cf);
	}
}

void calculaSequencial(void){
	unsigned long int i;
	mpf_t term;
	mpf_t pot_x;
	mpf_t termAnt;
	mpf_t erro;	

	mpf_init_set_ui(term, 1);
	mpf_init_set_ui(pot_x, 1);
	mpf_init_set_ui(erro, 1);
	mpf_init_set_ui(termAnt, 2);

	for(i = 0; ; i++){
		potencia(pot_x, x, 2 * i);    						// pot_x = potencia(x, 2*i)
		calculaCoeficiente(term, i); 						// calcula coeficiente da parcela i
		mpf_mul(term, term, pot_x);							// termo final

		mpf_add(cosseno, cosseno, term);

		if(criterio == 'f'){
			mpf_sub(erro, term, termAnt);
			mpf_set(termAnt, term);
		}

		else if(criterio == 'm'){
			mpf_set(erro, term);	
		}

		mpf_abs(erro, erro);

		if(mpf_cmp(erro, precisao) < 0){
			break;	
		}	
		
		gmp_printf ("\nRodada: %ld Cosseno:\n%.50Ff\n", i, cosseno);
	}

	printf("\n---------------------------------------------------------\n");
	gmp_printf("Cosseno de %Ff eh igual a: \n%.100000Ff\n", x, cosseno);
	printf("\n# de termos calculados: %ld", i+1);
	printf("\n---------------------------------------------------------");

	mpf_set_si(cosseno, 0);
	
	mpf_clear(term);
	mpf_clear(pot_x);
	mpf_clear(termAnt);
	mpf_clear(erro);
}

void parserEntrada(int argc, char** argv){
	if(argc < 5 || argc > 6){
		printf("Formato esperado: ./ep2 <arg1> <arg2> <arg3> <arg4> [arg5]\n");
		printf("arg1 := #Threads: unsigned_int|0\n");
		printf("arg2 := criterio_de_parada: (f|m)\n");
		printf("arg3 := precisao: unsigned_int\n");
		printf("arg4 := valor_de_x: float\n");
		printf("arg5(opc) := opcao_debugger/sequencial: (d|s)\n\n");
		exit(-1);
	}

	else{
		q = atoi(argv[1]);
		int expoente;

		if(q == 0){
			q = sysconf(_SC_NPROCESSORS_ONLN);		// descobre a qtde de nucleos do computador
		}

		criterio = argv[2][0];
		
		if(criterio != 'f' && criterio != 'm'){
			printf("criterio f ou m!\n");
			exit(0);
		}

		expoente = atoi(argv[3]);

		if(expoente >= 0){
			mpf_init_set_d(precisao, 0.1);
		}

		else{
			expoente = -expoente; // expoente tem que ser positivo
			mpf_init_set_d(precisao, 10);	
		}

		potencia(precisao, precisao, expoente);
		
		mpf_init_set_str(x, argv[4], 10);

		opc = '0';

		if(argc == 6){
			opc = argv[5][0];
		}

	}
}

/*############################# FUNÇÕES MATEMATICAS #############################*/

void fatorial(mpf_t fact, unsigned long int n){
	if(n == 0){
		mpf_set_ui(fact, 1);
	}
	else{
		fatorial(fact, n-1);
		mpf_mul_ui(fact, fact, n);
	}
}

void potencia(mpf_t result, mpf_t base, int expo){
	int i;
	mpf_t tmp;

	mpf_init_set_d(tmp, 1);

	for(i = 0; i < expo; i++){
		mpf_mul(tmp, tmp, base);
	}

	mpf_set(result, tmp);
	mpf_clear(tmp);
}

/*###############################################################################*/