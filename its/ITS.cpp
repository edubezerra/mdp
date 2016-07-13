/* Program: iterated tabu search (ITS) for the maximum diversity 
 problem:
 maximize sum of c_ij*x_i*x_j
 subject to b_1 <= sum of x_i <= b_2,
 where c_ij are any real numbers (also negative).
 Author: Gintaras Palubeckis
 Date: 2006-03-09
 Language: C
 Multistarting mechanism: a solution for each restart is created from
 the current (rather good) solution S (vertex set) by performing
 one of the following operations a specified number of times:
 remove vertex from S (provided |S|>b_1), add vertex to S
 (provided |S|<b_2), exchange two vertices (provided |S|=b_1 or |S|=b_2),
 one being in S and another out of S. A vertex or two vertices as agents
 of these operations are picked out randomly from those candidates
 which, after the move or exchange has been made, would give largest
 values of the objective function.
 Some of the input data are supplied through the parameters and the rest
 through the (input) file. An instance of the problem in this file is
 represented by a list containing all (even zero) coefficients c_ij.
 Internally, the program uses a matrix with entries of 'double' type
 to store these coefficients. The program terminates when a specified
 time limit is reached.
 Parameters:
 - input file name;
 - output file name;
 - b_1;
 - b_2;
 - seed for random number generator;
 - coefficient used to set (by multiplying it by the number of vertices)
 the number of iterations for each tabu search run;
 - time limit on a run of the program (in seconds);
 - a pointer to structure 'Results' for writing down the solution found
 and some characteristics. The structure is defined in file 'ITS.h'.
 It is possible to have the last parameter null: in this case
 information is directed only to the output file.
 Examples of invocation:
 either
 Results *pres;
 char in_file_name[80],out_file_name[80];
 double seed=1000;
 pres=(Results *)calloc(1,sizeof(Results));
 strcpy(in_file_name,"D:\\data\\Type2.1.txt");
 strcpy(out_file_name,"D:\\temp\\Type2.1.res");
 ITS(in_file_name,out_file_name,50,50,seed,1000,20,pres);
 or just
 char in_file_name[80],out_file_name[80];
 double seed=1000;
 strcpy(in_file_name,"D:\\data\\Type2.1.txt");
 strcpy(out_file_name,"D:\\temp\\Type2.1.res");
 ITS(in_file_name,out_file_name,50,50,seed,1000,20,NULL);
 Input file contains:
 - the size n of the instance (the number of vertices of the graph);
 - for each pair i, j, i=1,...,n-1, j=i+1,...,n, the triplet:
 i-1, j-1, c_ij.
 Example of the input file:
 500
 0 1 704.48
 0 2 778.79
 0 3 498.3
 0 4 161.76
 .......
 .......
 */

#include <memory.h>
//#include "process.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "ITS.h"

double random(double *seed, double coef) {
	double rd, rf;
	rd = 16807 * (*seed);
	rf = floor(rd / coef);
	*seed = rd - rf * coef;
	return (*seed / (coef + 1));
}

long take_time(int *time_values, clock_t start) {
	int i;
	int hours, mins;
	long longsecs;
	float elapsed_sec;
	clock_t end;
	end = clock();
	elapsed_sec = (end - start) / CLK_TCK;
	longsecs = elapsed_sec;
	for (i = 1; i <= 4; i++)
		time_values[i] = 0;
	hours = (int) (longsecs / 3600);
	if (hours > 0) /* more than an hour  */
	{
		time_values[1] = hours;
		longsecs -= hours * 3600;
	}
	mins = (int) (longsecs / 60);
	if (mins > 0) /* more than a minute */
	{
		time_values[2] = mins;
		longsecs -= mins * 60;
	}
	time_values[3] = (int) longsecs;
	time_values[4] = elapsed_sec * 1000 - (long) elapsed_sec * 1000;
	return (long) elapsed_sec;
}

double random_start(int size, int b1, int b2, int *cl_size, double coef,
		double *seed, double **pweight, Solution *psol) {
	int i, j, r;
	double sol_value = 0;

	i = random(seed, coef) * (b2 - b1 + 1);
	*cl_size = b1 + i;
	for (i = 1; i <= size; i++) {
		(psol + i)->t = i;
		(psol + i)->sol = 0;
	}
	for (i = 1; i <= *cl_size; i++) {
		r = random(seed, coef) * (size - i + 1);
		r += i;
		(psol + (psol + r)->t)->sol = 1;
		(psol + r)->t = (psol + i)->t;
	}
	for (i = 1; i < size; i++) {
		if ((psol + i)->sol != 1)
			continue;
		for (j = i + 1; j <= size; j++)
			if ((psol + j)->sol == 1)
				sol_value += weight(i,j);
	}
	return sol_value;
}

double steepest_ascent(FILE *out, int size, int b1, int b2, int *cl_size,
		double coef, double *seed, double **pweight, Solution *psol) {
	int i, j, k, q;
	int m0 = 0, m1 = 0;
	int count;
	int cs = 0, ind, val;
	double sol_value = 0, f_value;
	double f = 0.;
	double max_impr, impr0, impr1;
	double db, dr_numb;

	*cl_size = 0;
	if (b1 == b2)
		q = b1;
	else
		q = (b1 + b2) / 2;
	for (i = 1; i <= size; i++) {
		(psol + i)->sol = -1;
		(psol + i)->sf = 0;
		(psol + i)->s1 = 0;
	}
	for (i = 1; i < size; i++)
		for (j = i + 1; j <= size; j++) {
			f += (weight(i,j));
			(psol + i)->sf += (weight(i,j));
			(psol + j)->sf += (weight(i,j));
		}
	f *= (q * q);
	for (k = 1; k <= size; k++) {
		max_impr = NEG_LARGE_LONG;
		ind = -1;
		for (i = 1; i <= size; i++) {
			if ((psol + i)->sol != -1)
				continue;
			db = ((double) q) * (psol + i)->sf
					+ ((double) size) * (psol + i)->s1;
			impr0 = -q * db;
			if (impr0 > max_impr) {
				max_impr = impr0;
				ind = i;
				val = 0;
				count = 1;
			} else if (impr0 == max_impr) {
				count++;
				dr_numb = random(seed, coef);
				if (dr_numb <= 1. / (double) count) {
					ind = i;
					val = 0;
				}
			}
			impr1 = (size - q) * db;
			if (impr1 > max_impr) {
				max_impr = impr1;
				ind = i;
				val = 1;
				count = 1;
			} else if (impr1 == max_impr) {
				count++;
				dr_numb = random(seed, coef);
				if (dr_numb <= 1. / (double) count) {
					ind = i;
					val = 1;
				}
			}
		}
		f += max_impr;
		(psol + ind)->sol = val;
		if (val == 0)
			m0++;
		else {
			m1++;
			(*cl_size)++;
		}
		for (i = 1; i <= size; i++) {
			if ((psol + i)->sol != -1 || weight(i,ind) == 0)
				continue;
			if (val == 1)
				(psol + i)->s1 += weight(i,ind);
			(psol + i)->sf -= weight(i,ind);
		}
		if (m1 == b2) {
			cs = 1;
			break;
		}
		if ((size - m0) == b1) {
			cs = -1;
			break;
		}
	}
	if (cs > 0) {
		for (k = 1; k <= size; k++) {
			if ((psol + k)->sol != -1)
				continue;
			db = ((double) q) * (psol + k)->sf
					+ ((double) size) * (psol + k)->s1;
			impr0 = -q * db;
			f += impr0;
			(psol + k)->sol = 0;
			for (i = k + 1; i <= size; i++) {
				if ((psol + i)->sol != -1 || weight(i,k) == 0)
					continue;
				(psol + i)->sf -= (weight(i,k));
			}
		}
	} else if (cs < 0) {
		for (k = 1; k <= size; k++) {
			if ((psol + k)->sol != -1)
				continue;
			db = ((double) q) * (psol + k)->sf
					+ ((double) size) * (psol + k)->s1;
			impr1 = (size - q) * db;
			f += impr1;
			(psol + k)->sol = 1;
			(*cl_size)++;
			for (i = k + 1; i <= size; i++) {
				if ((psol + i)->sol != -1 || weight(i,k) == 0)
					continue;
				(psol + i)->s1 += (weight(i,k));
				(psol + i)->sf -= (weight(i,k));
			}
		}
	}
	f /= size;
	f /= size;
	f_value = f;
	for (i = 1; i < size; i++) {
		if ((psol + i)->sol != 1)
			continue;
		for (j = i + 1; j <= size; j++)
			if ((psol + j)->sol == 1)
				sol_value += weight(i,j);
	}
	if (f_value < sol_value - 0.00001 || f_value > sol_value + 0.00001) {
		fprintf(out, "!!! discrepancy in solution values: %11.3lf   %11.3lf\n",
				f_value, sol_value);
		exit(1);
	}
	return sol_value;
}

double get_solution(int size, int b1, int b2, int perturb_count,
		int cand_list_size, double init_value, int *cl_size, double coef,
		double *seed, double **pweight, Solution *psol) {
	int i, j, k, m;
	int it = 0;
	int ind, ind1, ind2, minind, cand_count;
	double sol_value;
	double minval, del;

	sol_value = init_value;
	for (i = 1; i <= size; i++)
		(psol + i)->val = 0;

	while (it < perturb_count) {
		cand_count = 0;
		minval = POS_LARGE_INT;
		if (*cl_size < b2)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->val > 0 || (psol + k)->sol == 1)
					continue;
				if (cand_count < cand_list_size) {
					cand_count++;
					(psol + cand_count)->cand1 = k;
					(psol + cand_count)->cand2 = -1;
					(psol + cand_count)->d = (psol + k)->cl;
					if ((psol + k)->cl < minval) {
						minval = (psol + k)->cl;
						minind = cand_count;
					}
				} else if ((psol + k)->cl > minval) {
					(psol + minind)->cand1 = k;
					(psol + minind)->cand2 = -1;
					(psol + minind)->d = (psol + k)->cl;
					minval = (psol + 1)->d;
					minind = 1;
					for (j = 2; j <= cand_count; j++)
						if ((psol + j)->d < minval) {
							minval = (psol + j)->d;
							minind = j;
						}
				}
			}
		if (*cl_size > b1)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->val > 0 || (psol + k)->sol == 0)
					continue;
				if (cand_count < cand_list_size) {
					cand_count++;
					(psol + cand_count)->cand1 = k;
					(psol + cand_count)->cand2 = -1;
					(psol + cand_count)->d = -(psol + k)->cl;
					if (-(psol + k)->cl < minval) {
						minval = -(psol + k)->cl;
						minind = cand_count;
					}
				} else if (-(psol + k)->cl > minval) {
					(psol + minind)->cand1 = k;
					(psol + minind)->cand2 = -1;
					(psol + minind)->d = -(psol + k)->cl;
					minval = (psol + 1)->d;
					minind = 1;
					for (j = 2; j <= cand_count; j++)
						if ((psol + j)->d < minval) {
							minval = (psol + j)->d;
							minind = j;
						}
				}
			}
		if ((*cl_size == b1 || *cl_size == b2))
			for (k = 1; k <= size; k++) {
				if ((psol + k)->val > 0 || (psol + k)->sol == 0)
					continue;
				for (m = 1; m <= size; m++) {
					if ((psol + m)->val > 0 || (psol + m)->sol == 1)
						continue;
					del = (psol + m)->cl - (psol + k)->cl - weight(k,m);
					if (cand_count < cand_list_size) {
						cand_count++;
						(psol + cand_count)->cand1 = k;
						(psol + cand_count)->cand2 = m;
						(psol + cand_count)->d = del;
						if (del < minval) {
							minval = del;
							minind = cand_count;
						}
					} else if (del > minval) {
						(psol + minind)->cand1 = k;
						(psol + minind)->cand2 = m;
						(psol + minind)->d = del;
						minval = (psol + 1)->d;
						minind = 1;
						for (j = 2; j <= cand_count; j++)
							if ((psol + j)->d < minval) {
								minval = (psol + j)->d;
								minind = j;
							}
					}
				}
			}
		ind = random(seed, coef) * cand_count + 1;
		ind1 = (psol + ind)->cand1;
		ind2 = (psol + ind)->cand2;
		if (ind2 == -1) {
			if ((psol + ind1)->sol == 1) {
				for (j = 1; j <= size; j++)
					if (j != ind1)
						(psol + j)->cl -= (weight(j,ind1));
				(*cl_size)--;
			} else {
				for (j = 1; j <= size; j++)
					if (j != ind1)
						(psol + j)->cl += (weight(j,ind1));
				(*cl_size)++;
			}
			(psol + ind1)->sol = 1 - (psol + ind1)->sol;
			(psol + ind1)->val = 1;
			it++;
		} else {
			for (j = 1; j <= size; j++) {
				if (j != ind1)
					(psol + j)->cl -= (weight(j,ind1));
				if (j != ind2)
					(psol + j)->cl += (weight(j,ind2));
			}
			(psol + ind1)->sol = 0;
			(psol + ind2)->sol = 1;
			(psol + ind1)->val = (psol + ind2)->val = 1;
			it += 2;
		}
		sol_value += (psol + ind)->d;
	}
	return sol_value;
}

double local_search(int size, int b1, int b2, int *cl_size, long *it_count,
		double **pweight, Solution *psol) {
	int j, k, m;
	int repeat = 1;
	double del;
	double value_change = 0;

	while (repeat > 0) {
		repeat = 0;
		if (*cl_size < b2)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->sol == 1)
					continue;
				(*it_count)++;
				if ((psol + k)->cl <= 0.00001)
					continue;
				repeat = 1;
				((psol + 2)->performance)++;
				(psol + k)->sol = 1;
				value_change += (psol + k)->cl;
				for (j = 1; j <= size; j++)
					if (j != k)
						(psol + j)->cl += (weight(j,k));
				(*cl_size)++;
				if (*cl_size >= b2)
					break;
			}
		if (*cl_size > b1)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->sol == 0)
					continue;
				(*it_count)++;
				if ((psol + k)->cl >= -0.00001)
					continue;
				repeat = 1;
				((psol + 2)->performance)++;
				(psol + k)->sol = 0;
				value_change -= (psol + k)->cl;
				for (j = 1; j <= size; j++)
					if (j != k)
						(psol + j)->cl -= (weight(j,k));
				(*cl_size)--;
				if (*cl_size <= b1)
					break;
			}
		if ((*cl_size == b1 || *cl_size == b2))
			for (k = 1; k <= size; k++) {
				if ((psol + k)->sol == 0)
					continue;
				for (m = 1; m <= size; m++) {
					if ((psol + m)->sol == 1)
						continue;
					(*it_count)++;
					del = (psol + m)->cl - (psol + k)->cl - weight(k,m);
					if (del <= 0.00001)
						continue;
					repeat = 1;
					((psol + 2)->performance)++;
					(psol + k)->sol = 0;
					(psol + m)->sol = 1;
					value_change += del;
					for (j = 1; j <= size; j++) {
						if (j != k)
							(psol + j)->cl -= (weight(j,k));
						if (j != m)
							(psol + j)->cl += (weight(j,m));
					}
					break;
				}
			}
	}
	return value_change;
}

double tabu_search(int size, int b1, int b2, int keep_tabu_time1,
		int keep_tabu_time2, int start, long time_limit, long it_bound,
		double sol_value, int *cl_size, int *vert1, int *vert2, int *stop_cond,
		int *time_values_opt, double *best_value, clock_t start_time,
		double **pweight, int **ptabu, Solution *psol) {
	int i, j, k, k1, k2, m;
	int ind1, ind2, imp;
	int tl_ln = 0;
	long it = 0;
	long elapsed_time;
	double best_improvement;
	double del;
	clock_t end;

	for (i = 1; i <= size; i++) {
		(psol + i)->t = 0;
		(psol + i)->cl = 0.;
	}
	for (i = 1; i < size; i++)
		for (j = i + 1; j <= size; j++) {
			tabu(i,j) = 0;
			tabu(j,i) = 0;
			if ((psol + i)->sol == 1)
				(psol + j)->cl += (weight(i,j));
			if ((psol + j)->sol == 1)
				(psol + i)->cl += (weight(i,j));
		}
	while (it < it_bound) {
		ind1 = ind2 = -1;
		imp = 0;
		best_improvement = NEG_LARGE_LONG;
		if (*cl_size < b2)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->t > 0 || (psol + k)->sol == 1)
					continue;
				it++;
				if (sol_value + (psol + k)->cl > *best_value + 0.00001) {
					best_improvement = (psol + k)->cl;
					ind1 = k;
					imp = 1;
					break;
				}
				if ((psol + k)->cl > best_improvement) {
					best_improvement = (psol + k)->cl;
					ind1 = k;
				}
			}
		if (*cl_size > b1 && imp == 0)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->t > 0 || (psol + k)->sol == 0)
					continue;
				it++;
				if (sol_value - (psol + k)->cl > *best_value + 0.00001) {
					best_improvement = -(psol + k)->cl;
					ind1 = k;
					imp = 2;
					break;
				}
				if (-(psol + k)->cl > best_improvement) {
					best_improvement = -(psol + k)->cl;
					ind1 = k;
				}
			}
		if ((*cl_size == b1 || *cl_size == b2) && imp == 0)
			for (k = 1; k <= size; k++) {
				if ((psol + k)->sol == 0)
					continue;
				for (m = 1; m <= size; m++) {
					if ((psol + m)->sol == 1 || tabu(k,m) > 0)
						continue;
					it++;
					del = (psol + m)->cl - (psol + k)->cl - weight(k,m);
					if (sol_value + del > *best_value + 0.00001) {
						best_improvement = del;
						ind1 = k;
						ind2 = m;
						imp = 3;
						break;
					}
					if (del > best_improvement) {
						best_improvement = del;
						ind1 = k;
						ind2 = m;
					}
				}
			}
		if (ind2 == -1) {
			if ((psol + ind1)->sol == 1) {
				for (j = 1; j <= size; j++)
					if (j != ind1)
						(psol + j)->cl -= (weight(j,ind1));
				(*cl_size)--;
			} else {
				for (j = 1; j <= size; j++)
					if (j != ind1)
						(psol + j)->cl += (weight(j,ind1));
				(*cl_size)++;
			}
			(psol + ind1)->sol = 1 - (psol + ind1)->sol;
			sol_value += best_improvement;
		} else {
			for (j = 1; j <= size; j++) {
				if (j != ind1)
					(psol + j)->cl -= (weight(j,ind1));
				if (j != ind2)
					(psol + j)->cl += (weight(j,ind2));
			}
			(psol + ind1)->sol = 0;
			(psol + ind2)->sol = 1;
			sol_value += best_improvement;
		}
		if (imp > 0) {
			sol_value += local_search(size, b1, b2, cl_size, &it, pweight,
					psol);
			for (i = 1; i <= size; i++)
				(psol + i)->best_sol = (psol + i)->sol;
			*best_value = sol_value;
			((psol + 1)->performance)++;
			(psol + 3)->performance = start;
			take_time(time_values_opt, start_time);
		}
		for (i = 1; i <= size; i++)
			if ((psol + i)->t > 0)
				((psol + i)->t)--;
		m = -1;
		for (i = 1; i <= tl_ln; i++) {
			k1 = *(vert1 + i);
			k2 = *(vert2 + i);
			(tabu(k1,k2))--;
			(tabu(k2,k1))--;
			if (tabu(k1,k2) == 0)
				m = i;
		}
		if (m > 0) {
			for (i = m; i < tl_ln; i++) {
				*(vert1 + i) = *(vert1 + i + 1);
				*(vert2 + i) = *(vert2 + i + 1);
			}
			tl_ln--;
		}
		if (ind2 == -1)
			(psol + ind1)->t = keep_tabu_time1;
		else {
			tl_ln++;
			*(vert1 + tl_ln) = ind1;
			*(vert2 + tl_ln) = ind2;
			tabu(ind1,ind2) = keep_tabu_time2;
			tabu(ind2,ind1) = keep_tabu_time2;
		}
		end = clock();
		elapsed_time = (long) (end - start_time) / CLK_TCK;
		if (elapsed_time >= time_limit) {
			*stop_cond = 1;
			break;
		}
	}
	return sol_value;
}

double ITS_internal(FILE *out, int size, int b1, int b2, long time_limit,
		int keep_tabu_time1, int keep_tabu_time2, int perturb_count,
		int min_perturb_count, int cand_list_size, long it_bound, int *vert1,
		int *vert2, int *time_values_opt, double *seed1, clock_t start,
		double **pweight, int **ptabu, Solution *psol) {
	int i;
	int st = 1;
	int cl_size;
	int stop_cond = 0;
	double sol_value, best_value;
	double seed2, seed3, coef;

	coef = 2048;
	coef *= 1024;
	coef *= 1024;
	coef -= 1;
	seed2 = 2 * (*seed1);
	seed3 = 3 * (*seed1);
	(psol + 1)->performance = 0;
	(psol + 2)->performance = 0;

	if (size > 200)
		sol_value = random_start(size, b1, b2, &cl_size, coef, seed1, pweight,
				psol);
	else
		sol_value = steepest_ascent(out, size, b1, b2, &cl_size, coef, seed1,
				pweight, psol);
//fprintf(out,"*****Steepest ascent*****   sol_value=%8ld  cl_size=%4d\n",
//sol_value,cl_size);
	best_value = sol_value;
	for (i = 1; i <= size; i++)
		(psol + i)->best_sol = (psol + i)->sol;
	sol_value = tabu_search(size, b1, b2, keep_tabu_time1, keep_tabu_time2, st,
			time_limit, it_bound, sol_value, &cl_size, vert1, vert2, &stop_cond,
			time_values_opt, &best_value, start, pweight, ptabu, psol);
	if (perturb_count > b1)
		perturb_count = b1;
	while (stop_cond == 0) {
		st++;
		if (perturb_count <= min_perturb_count)
			i = perturb_count;
		else {
			i = random(&seed3, coef) * (perturb_count - min_perturb_count + 1);
			i += min_perturb_count;
		}
		sol_value = get_solution(size, b1, b2, i, cand_list_size, sol_value,
				&cl_size, coef, &seed2, pweight, psol);
		sol_value = tabu_search(size, b1, b2, keep_tabu_time1, keep_tabu_time2,
				st, time_limit, it_bound, sol_value, &cl_size, vert1, vert2,
				&stop_cond, time_values_opt, &best_value, start, pweight, ptabu,
				psol);
	}
	(psol + 5)->performance = st;
	return best_value;
}

void ITS(char *in_file_name, char *out_file_name, int b1, int b2, double seed,
		long iterations_coef, long time_limit, Results *pres) {
	FILE *out, *in;
	double **pweight;
	int **ptabu;
	Solution *psol;
	int *vert1;
	int *vert2;

	int i, j;
	int size, cl_size = 0;
	int keep_tabu_time1;
	int perturb_count;
	int time_values[5], time_values_opt[5];
	long lo, e_count;
	long it_bound;
	long time_in_seconds;
	double w;
	double value, value_from_sol = 0;
	clock_t start;

	if ((in = fopen(in_file_name, "r")) == NULL) {
		printf("  fopen failed for input");
		exit(1);
	}
	fscanf(in, "%d", &size);
	if ((out = fopen(out_file_name, "w")) == NULL) {
		printf("  fopen failed for output  %s", out_file_name);
		exit(1);
	}
	ALMF(pweight, size+1)
	for (i = 0; i <= size; i++)
		ALF(*(pweight+i), size+1)
	ALM(ptabu, size+1)
	for (i = 0; i <= size; i++)
		ALI(*(ptabu+i), size+1)
	ALI(vert1, TABU_TIME2+1)
	ALI(vert2, TABU_TIME2+1)
	ALS(psol, Solution, size+1)
	keep_tabu_time1 = TABU_TIME1;
	i = size / TABU_COEF;
	if (i < keep_tabu_time1)
		keep_tabu_time1 = i;
	it_bound = ITERATIONS_FIXED_BOUND;
	lo = ((long) size) * iterations_coef;
	if (it_bound < lo)
		it_bound = lo;
	perturb_count = size * PER_COEF;
	e_count = ((long) size) * (size - 1) / 2;
	for (lo = 1; lo <= e_count; lo++) {
		fscanf(in, "%d %d %lf", &i, &j, &w);
		i++;
		j++;
		weight(i,j) = w;
		weight(j,i) = w;
	}
	start = clock();
	value = ITS_internal(out, size, b1, b2, time_limit, keep_tabu_time1,
			TABU_TIME2, perturb_count, MIN_PER_COUNT, LIST_SIZE, it_bound,
			vert1, vert2, time_values_opt, &seed, start, pweight, ptabu, psol);
	time_in_seconds = take_time(time_values, start);
	for (i = 1; i <= size; i++) {
		if ((psol + i)->best_sol != 1)
			continue;
		cl_size++;
		for (j = i + 1; j <= size; j++)
			if ((psol + j)->best_sol == 1)
				value_from_sol += weight(i,j);
	}
	fprintf(out, "   graph order                    = %5d\n", size);
	fprintf(out, "   lower bound                    = %5d\n", b1);
	fprintf(out, "   upper bound                    = %5d\n", b2);
	fprintf(out, "   time limit                     = %5ld\n", time_limit);
	fprintf(out, "   number of iterations per start = %10ld\n", it_bound);
	fprintf(out, "   number of starts executed      = %3d\n",
			(psol + 5)->performance);
	fprintf(out, "   number of improvements         = %3d\n",
			(psol + 1)->performance);
	fprintf(out, "   last improvement at start no.  = %2d\n",
			(psol + 3)->performance);
	if (value < value_from_sol - 0.00001 || value > value_from_sol + 0.00001)
		fprintf(out,
				"!!! some discrepancy in solution values: %11.3lf   %11.3lf\n",
				value, value_from_sol);
	else
		fprintf(out, "   solution value                 = %11.3lf  %11.3lf\n",
				value_from_sol, value);
	fprintf(out, "   subgraph size                  = %5d\n", cl_size);
	lo = 3600 * (long) time_values_opt[1] + 60 * time_values_opt[2]
			+ time_values_opt[3];
	fprintf(out, "   time to solution: %d : %d : %d.%3d  (=%4ld seconds)\n",
			time_values_opt[1], time_values_opt[2], time_values_opt[3],
			time_values_opt[4], lo);
	fprintf(out, "   total time: %d : %d : %d.%3d  (=%4ld seconds)\n",
			time_values[1], time_values[2], time_values[3], time_values[4],
			time_in_seconds);
	fprintf(out, "\n");
	if (pres != NULL) {
		ALI(pres->sol, size+1)
		for (i = 1; i <= size; i++)
			sol(i) = (psol + i)->best_sol;
		pres->value = value;
		pres->time_to_opt = lo;
		pres->total_time = time_in_seconds;
		pres->characts[0] = size;
		pres->characts[1] = time_limit;
		pres->characts[2] = (psol + 5)->performance;
		pres->characts[3] = (psol + 1)->performance;
		pres->characts[4] = (psol + 3)->performance;
		pres->characts[5] = b1;
		pres->characts[6] = b2;
		pres->characts[7] = cl_size;
	}
	fclose(out);
	fclose(in);
}
