#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <gmp.h>

/*********************** VARIAVEIS GLOBAIS *********************************/
int q;
char criterio;
char opc;

mpf_t precisao;
mpf_t x;
mpf_t cosseno;


/*********************** PROTOTIPOS DAS FUNCOES *********************************/
void parserEntrada(int, char**);
void potencia(mpf_t, mpf_t, int);
void fatorial(mpz_t, unsigned long int);
void calculaSequencial(void);
/********************************************************************************/

int main(int argc, char** argv){
	mpf_init(precisao);
	mpf_init(x);
	mpf_init(cosseno);
	mpf_set_default_prec(1000000);
	
	parserEntrada(argc, argv);

	if(opc == 's'){
    	calculaSequencial();
	}

	mpf_clear(x);
	mpf_clear(cosseno);
	mpf_clear(precisao);
	return 0;
}


/*********************** FUNCOES QUE FALTAM IMPLEMENTAR *********************************/

void calculaSequencial(void){
	unsigned long int i;
	mpf_t term;
	mpf_t pot_x;
	mpf_t termAnt;
	mpf_t termDif;	
	mpz_t fact;

	mpf_init_set_ui(term, 1);
	mpf_init_set_ui(pot_x, 1);
	mpf_init_set_ui(termAnt, 2);
	mpf_init_set_ui(termDif, 1);
	mpz_init_set_ui(fact, 1);	
	
	for(i = 0; ; i++){
		potencia(pot_x, x, 2*i);    						// pot_x = potencia(x, 2*i)

		if(i % 2 == 1){										// pot_x = potencia(-1, i) * potencia(x, 2*i) 
			mpf_neg(pot_x, pot_x);
		}
											
		// gmp_printf ("%lu) Mull:  %2.20Ff\n", i, aux);

		fatorial(fact, 2*i);	          					// fact = fatorial(2*i)
		mpf_set_z(term, fact);								// term = (mpf_t)fatorial(2*i) 
		mpf_div(term, pot_x, term);							// term = potencia(-1, i)*potencia(x, 2*i)/fatorial(2*i)
		// gmp_printf ("%lu) Div aux :  %2.20Ff\n", i, aux);
		// gmp_printf ("%lu) Div aux1:  %2.20Ff\n", i, aux1);
		
		mpf_add(cosseno, cosseno, term);

		if(criterio == 'f'){
			mpf_sub(termDif, term, termAnt);
			//mpf_set(termAnt, term);
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
		
		mpf_set(termAnt, term);	
		
		gmp_printf ("\nCosseno:\n%.50Ff\n", cosseno);
	}

	printf("\n---------------------------------------------------------\n");
	gmp_printf("Cosseno de %Ff eh igual a: \n%.50Ff\n", x, cosseno);
	printf("\n# de termos calculados: %ld", i+1);
	printf("\n---------------------------------------------------------\n\n");

	mpf_set_si(cosseno, 0);
	
	mpf_clear(term);
	mpf_clear(pot_x);
	mpz_clear(fact);
	mpf_clear(termAnt);
	mpf_clear(termDif);
}

/*********************** FUNCOES PRONTAS *********************************/
void parserEntrada(int argc, char** argv){
	if(argc < 5 || argc > 6){
		printf("Erro na quantidade dos parametros.\n");
		printf("Formato: ./ep2  <q (int)>  <f|m>  <exp (int)>  <x>  [d|s]\n");
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
		
		if(argc == 6){
			opc = argv[5][0];
		}

		//gmp_printf("q = %d    criterio = %c    precisao = %.17Ff    x = %Ff    opc = %c \n", q, criterio, precisao, x, opc);
	}
}


/*############################# FUNÇÕES MATEMATICAS #############################*/

void fatorial(mpz_t fact, unsigned long int n){
	if(n == 0){
		mpz_set_si(fact, 1);
	}
	else{
		fatorial(fact, n-1);
		mpz_mul_ui(fact, fact, n);
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

/*
double potencia(double base, int exp){
	int i;
	double pot = 1;

	for(i = 0; i < exp; i++){
		pot *= base;
	}

	return pot;
}

void fatorial (mpz_t result, unsigned long n){
  unsigned long i;
  mpz_init_set_str (result, "1", 0);

  for (i = 1; i <= n; i++)
    mpz_mul_ui (result, result, i);
}
*/
/*###############################################################################*/