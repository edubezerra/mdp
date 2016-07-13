#include <memory.h>
//#include <process.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define	CLK_TCK		CLOCKS_PER_SEC

typedef struct {
	int *sol; /* solution obtained                                */
	double value; /* its value                                        */
	long time_to_opt; /* time to solution, secs                           */
	long total_time; /* total time, secs                                 */
	long characts[10]; /* some characteristics:                            */
	/*    characts[0] - graph order                     */
	/*    characts[1] - time limit                      */
	/*    characts[2] - number of starts executed       */
	/*    characts[3] - number of solution improvements */
	/*    characts[4] - start no. for the last          */
	/*                  improvement                     */
	/*    characts[5] - lower bound on subgraph's size  */
	/*    characts[6] - upper bound on subgraph's size  */
	/*    characts[7] - subgraph's size                 */
} Results;

void ITS(char *, char *, int, int, double, long, long, Results *);

int main(int argc, char **argv) {
	long start = clock();

	FILE *out;
	Results *pres;
	char in_file_name[80];
	char out_file_name[80];
	char summary_file_name[80];
	int i, j;
	int b1, b2;
	int count;
	int sind;
	long time_limit;
	long iterations_coef;
	double seeds[11] = { 0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000,
			9000, 10000 };
	char out_file[80];
	char numbs[31][2] = { { '0' }, { '1' }, { '2' }, { '3' }, { '4' }, { '5' },
			{ '6' }, { '7' }, { '8' }, { '9' }, { '1', '0' } };
	double av_value = 0., av_time = 0.;

	if (argc <= 3) {
		printf("  specify data and output files");
		exit(1);
	}
	strcpy(in_file_name, argv[1]);
	strcpy(out_file_name, argv[2]);
	strcpy(summary_file_name, argv[3]);
	pres = (Results *) calloc(1, sizeof(Results));
	if ((out = fopen(summary_file_name, "w")) == NULL) {
		printf("  fopen failed for output  %s", summary_file_name);
		exit(1);
	}
	sind = strlen(argv[2]) - 1;

	iterations_coef = 1000;
	time_limit = (long) 1;
	b1 = 30;
	b2 = 30;
	count = 10;
	int pos_ponto;
	for (i = 1; i <= count; i++) {

		for (j = 0; j < 80; j++) {

			if (out_file_name[j]=='.') {
				pos_ponto = j;
				break;
			}
			//out_file[j] = out_file_name[j];

		}
		strncpy(out_file, out_file_name, pos_ponto);
		out_file[pos_ponto] = '0';
		out_file[pos_ponto+1] = i + '0';
		out_file[pos_ponto+2] =  '.';
		//strncpy(out_file, &out_file_name[pos_ponto], 3);
		if (i==1) strcat(out_file, "txt");
		//out_file[]=0;
/*		if (i < 10) {
			for (j = 79; j > sind + 1; j--)
				out_file[j] = out_file[j - 1];
			out_file[sind + 1] = numbs[i][0];
		} else {
			for (j = 79; j > sind + 2; j--)
				out_file[j] = out_file[j - 2];
			out_file[sind + 1] = numbs[i][0];
			out_file[sind + 2] = numbs[i][1];
		}*/
		ITS(in_file_name, out_file, b1, b2, seeds[i], iterations_coef,
				time_limit, pres);
		fprintf(out, "    %11.3lf       %8ld\n", pres->value,
				pres->time_to_opt);
		av_value += pres->value;
		av_time += pres->time_to_opt;
	}
	av_value /= count;
	av_time /= count;
	fprintf(out, "     %11.3lf    %11.3lf\n", av_value, av_time);


	long end = clock();
	long elapsed_time = (long) (end - start) / CLK_TCK;
	fprintf(out, "elapsed time:  %ld\n", elapsed_time);
	fclose(out);
	return 0;
}

