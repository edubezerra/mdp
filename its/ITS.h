#include <stdio.h>

#define TABU_COEF                    4
#define TABU_TIME1                  20
#define TABU_TIME2                  20
#define ITERATIONS_FIXED_BOUND   10000
#define PER_COEF                   0.1
#define MIN_PER_COUNT               10
#define LIST_SIZE                    5 

#define NEG_LARGE_LONG        -1000000
#define POS_LARGE_INT         30000
#define	CLK_TCK		CLOCKS_PER_SEC
#define ALS(X,Y,Z) if ((X=(Y *)calloc(Z,sizeof(Y)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALI(X,Z) if ((X=(int *)calloc(Z,sizeof(int)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALF(X,Z) if ((X=(double *)calloc(Z,sizeof(double)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALM(X,Z) if ((X=(int **)calloc(Z,sizeof(int *)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}
#define ALMF(X,Z) if ((X=(double **)calloc(Z,sizeof(double *)))==NULL) \
       {fprintf(out,"  failure in memory allocation\n");exit(0);}

#define weight(X,Y) *(*(pweight+X)+Y)
#define tabu(X,Y) *(*(ptabu+X)+Y)
#define sol(Y) *(pres->sol+Y)

typedef struct
     {int sol;            /*  */
      int best_sol;       /*  */
      int val;            /*  */
      int t;              /*  */
      int cand1;          /*  */
      int cand2;          /*  */
      int performance;    /*  */
      double cl;           /*  */
      double d;            /*  */
      double s1;           /*  */
      double sf;           /*  */
     }Solution;

typedef struct
     {int *sol;           /* solution obtained                                */
      double value;       /* its value                                        */
      long time_to_opt;   /* time to solution, secs                           */
      long total_time;    /* total time, secs                                 */
      long characts[10];  /* some characteristics:                            */
                          /*    characts[0] - graph order                     */
                          /*    characts[1] - time limit                      */
                          /*    characts[2] - number of starts executed       */
                          /*    characts[3] - number of solution improvements */
                          /*    characts[4] - start no. for the last          */
                          /*                  improvement                     */
                          /*    characts[5] - lower bound on subgraph's size  */
                          /*    characts[6] - upper bound on subgraph's size  */
                          /*    characts[7] - subgraph's size                 */
     }Results;
