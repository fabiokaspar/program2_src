#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>
#include <pthread.h>
#include <semaphore.h>
#include <time.h>

/*********************** VARIAVEIS GLOBAIS *********************************/
int q;
char criterio;
char opc;
time_t inicio;

mpf_t boxTermAnt;
int stop;
unsigned long int ind_sum;

mpf_t precisao;
mpf_t x;
mpf_t cosseno;

pthread_barrier_t bar;
sem_t mutex;

/*********************** PROTOTIPOS DAS FUNCOES *********************************/
void parserEntrada(int, char**);
void potencia(mpf_t, mpf_t, int);
void fatorial(mpf_t, unsigned long int);
void calculaSequencial(void);
void* threadParcela(void*);
void calculaCoeficiente(mpf_t, unsigned long int);

/********************************************************************************/

int main(int argc, char** argv){
	int i;
	mpf_init(precisao);
	mpf_init(x);
	mpf_init(cosseno);
	mpf_set_default_prec(1000000);
	
	parserEntrada(argc, argv);

	if(opc == 's'){
		inicio = time(NULL);
    	calculaSequencial();
	}

	else{
		
		pthread_t parcela[q];

		if(sem_init(&mutex, 0, 1)){
			fprintf(stderr, "ERROR creating semaphore\n");
			exit(0);
		}
		
		if(pthread_barrier_init(&bar ,NULL, q)){
			fprintf(stderr, "ERROR initializing barrier\n");
			exit(0);
		}

		ind_sum = 0;
		stop = 0;
		mpf_init_set_si(boxTermAnt, -2);

		inicio = time(NULL);
		
		for(i = 0; i < q; i++){
		 	if(pthread_create(&parcela[i], NULL, &threadParcela, (void*) i)){
				fprintf(stderr, "ERROR creating threads\n");
				exit(0);
		 	}
		}

		for(i = 0; i < q; i++){
			if(pthread_join(parcela[i], NULL)){
				fprintf(stderr, "ERROR joining threads\n");
				exit(0);	
			}
		}

		printf("\n---------------------------------------------------------\n");
		gmp_printf("Cosseno de %Ff eh igual a: \n%.50Ff\n", x, cosseno);
		printf("\n# de termos calculados: %ld", ind_sum);
		printf("\n---------------------------------------------------------\n\n");
	}

	printf("Tempo estimado: %.2f seg\n", (float)(time(NULL)-inicio));

	mpf_clear(x);
	mpf_clear(cosseno);
	mpf_clear(precisao);
	pthread_exit(NULL);

	return 0;
}

void* threadParcela(void* id){
	unsigned long int i;
	mpf_t parcela;
	mpf_t pot_x;
	mpf_t termDif;

	mpf_init(parcela);
	mpf_init(pot_x);
	mpf_init(termDif);

	while(!stop){
		sem_wait(&mutex);
		i = ind_sum;
		ind_sum++;
		sem_post(&mutex);

		potencia(pot_x, x, 2 * i);
		calculaCoeficiente(parcela, i);
		mpf_mul(parcela, parcela, pot_x);

		if(criterio == 'f'){
			sem_wait(&mutex);
			mpf_sub(termDif, parcela, boxTermAnt);
			mpf_set(boxTermAnt, parcela);
			sem_post(&mutex);
		}

		mpf_abs(termDif, termDif);

		if(mpf_cmp(termDif, precisao) < 0){
			stop = 1;
		}
		
		sem_wait(&mutex);
		mpf_add(cosseno, cosseno, parcela);
		//gmp_printf("\nCosseno de %Ff eh igual a: \n%.30Ff\n", x, cosseno);
		sem_post(&mutex);

		pthread_barrier_wait(&bar);

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
	mpf_t termDif;	

	mpf_init_set_ui(term, 1);
	mpf_init_set_ui(pot_x, 1);
	mpf_init_set_ui(termDif, 1);
	mpf_init_set_ui(termAnt, 2);

	for(i = 0; ; i++){
		potencia(pot_x, x, 2 * i);    						// pot_x = potencia(x, 2*i)
		calculaCoeficiente(term, i); 						// calcula coeficiente da parcela i
		mpf_mul(term, term, pot_x);							// termo final

		mpf_add(cosseno, cosseno, term);

		if(criterio == 'f'){
			mpf_sub(termDif, term, termAnt);
			mpf_set(termAnt, term);
		}

		else if(criterio == 'm'){
			mpf_set(termDif, term);	
		}

		mpf_abs(termDif, termDif);

		if(mpf_cmp(termDif, precisao) < 0){
			//printf("RODADA TERMINO: %ld\n", i+1);
			//gmp_printf("\npenultimo termo =  %.300Ff\n", termAnt);
			//gmp_printf("\nultimo    termo =  %.300Ff\n", term);
			break;	
		}	
		
		gmp_printf ("\nCosseno:\n%.50Ff\n", cosseno);
	}

	printf("\n---------------------------------------------------------\n");
	gmp_printf("Cosseno de %Ff eh igual a: \n%.50Ff\n", x, cosseno);
	printf("\n# de termos calculados: %ld", i+1);
	printf("\n---------------------------------------------------------\n\n");

	mpf_set_si(cosseno, 0);
	
	mpf_clear(term);
	mpf_clear(pot_x);
	mpf_clear(termAnt);
	mpf_clear(termDif);
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
		if(q == 0){
			q = sysconf(_SC_NPROCESSORS_ONLN);		// descobre a qtde de nucleos do computador
		}

		criterio = argv[2][0];

		mpf_set_d(precisao, 0.1);
		potencia(precisao, precisao, atoi(argv[3]));
		//precisao = atoi(argv[3]);
		mpf_init_set_str(x, argv[4], 10);

		opc = 0;

		if(argc == 6){
			opc = argv[5][0];
		}

		//gmp_printf("q = %d    criterio = %c    precisao = %.17Ff    x = %Ff    opc = %c \n", q, criterio, precisao, x, opc);
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