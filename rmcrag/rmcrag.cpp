#include <functional>
#include <queue>
#include <set>
#include <vector>
#include <iostream>
#include <limits>
#include <math.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctime>

using namespace std;

template<typename T> void print_queue(T& q1) {
    T q = q1;
    while(!q.empty()) {
        std::cout << q.top() << " ";
        q.pop();
    }
    std::cout << '\n';
}

void print_set(std::set<int>& a_set) {
	for (std::set<int>::iterator itr = a_set.begin(); itr != a_set.end(); ++itr) {
		std::cout << *itr << ", " << endl;
	}
}

void printClusteringsLabels(std::ostream& os, const std::set<int>& clusterings)
{
    for (int const& label : clusterings)
    {
        os << label << ' ';
    }
}

class Element {

private:
public:
	std::set<int> clusterings;
	double cost;

	class DereferenceCompareElement : public std::binary_function<Element*, Element*, bool>
	{
		public:
		bool operator()(const Element* lhs, const Element* rhs) const
		{
			return lhs->cost < rhs->cost;
		}
	};

	friend ostream& operator<<(std::ostream& os, const Element* dt) {
	    	os << "(labels of clusterings: ";
	    	printClusteringsLabels(os, dt->clusterings);
	    	os << ", cost: " << dt->cost << ")";
		return os;
   	}

	~Element(void) {
//	    cout << "Object is being deleted" << endl;
	}

	void appendTop(std::set<int> initialClusterings, std::set<int> allClusterings, const double **pweight, int k) {
    		std::priority_queue<Element*, std::vector<Element*>, DereferenceCompareElement> q;
	
		for (std::set<int>::iterator it = allClusterings.begin(); it != allClusterings.end(); ++it) {
			const bool is_in = initialClusterings.find(*it) != initialClusterings.end();
			if(!is_in) {
				Element* e = new Element();
				for (std::set<int>::iterator itr = initialClusterings.begin(); itr != initialClusterings.end(); ++itr) {
					e->clusterings.insert(*itr);
				}
				e->clusterings.insert(*it);
				e->cost = getCost(e->clusterings, pweight, k);
				q.push(e);
			}
		}

		this->clusterings = q.top()->clusterings;
		this->cost = getCost(this->clusterings, pweight, k);
		q.pop();
		while(!q.empty()) {
			Element* e = q.top();
			delete e;
			q.pop();
		}
	}

	double getCost(std::set<int>& clusterings, const double **pweight, int k) {
		double result =  mvnmi(clusterings, pweight) + (k - clusterings.size());
		return result;
	}
	
	double mvnmi(std::set<int>& clusterings, const double **pweight) {
		return ganmi(clusterings, pweight) + gvnmi(clusterings, pweight);
	}

	double ganmi(std::set<int>& clusterings, const double **pweight) {
		int p = clusterings.size();

		if(p > 1) {
		   double factor = 1.0 / binomial_coefficient(p, 2);
		   double sum = 0.0;

		   std::set<int>::iterator it1, it2;

		   for (std::set<int>::iterator it1 = clusterings.begin(); it1 != clusterings.end(); ++it1) {
		     int elem1 = *it1;
		     for (std::set<int>::iterator it2 = clusterings.begin(); it2 != clusterings.end(); ++it2) {
		       int elem2 = *it2;
		       if(elem1 < elem2) {
			 sum = sum + pweight[*it1][*it2];
		       }
		     }
		   }

		   return factor * sum;
		} else {
			cout << "Error: number of clusterings must be greater than one!" << endl;
			exit(1);
		}
	}

	double gvnmi(std::set<int>& clusterings, const double **pweight) {
		int p = clusterings.size();

		if(p > 1) {
			double factor = 1.0 / binomial_coefficient(p, 2);
			double sum = 0.0;

			for (std::set<int>::iterator it1 = clusterings.begin(); it1 != clusterings.end(); ++it1) {
			     int elem1 = *it1;
			     for (std::set<int>::iterator it2 = clusterings.begin(); it2 != clusterings.end(); ++it2) {
			       int elem2 = *it2;
			       if(elem1 < elem2) {
				 sum = sum + pow(ganmi(clusterings, pweight) - pweight[*it1][*it2], 2);
			       }
			     }
			} 
			return factor * sum;
		} else {
			cout << "Error: number of clusterings must be greater than one!" << endl;
			exit(1);
		}
	}

	template <class T = unsigned long>
	T binomial_coefficient(unsigned long n, unsigned long k) {
	    unsigned long i;
	    T b;
	    if (0 == k || n == k) {
	    	return 1;
	    }
	    if (k > n) {
    		return 0;
	    }
	    if (k > (n - k)) {
		    k = n - k;
	    }
	    if (1 == k) {
		    return n;
	    }
	    b = 1;
	    for (i = 1; i <= k; ++i) {
    		b *= (n - (k - i));
	    	if (b < 0) return -1; /* Overflow */
	    	b /= i;
	    }
	    return b;
	}
};

Element* rmcrag(std::set<int> clusterings, const double **pweight, int p, unsigned int k);

Element* rmcrag(const std::set<int> clusterings, const double **pweight, int p, unsigned int k) {

	std::priority_queue<Element*, std::vector<Element*>, Element::DereferenceCompareElement> Q;

	for(int i=0; i < p; i++) {
		Element* e = new Element();
		std::set<int> originalClusterings;
		originalClusterings.insert(i);
		e->appendTop(originalClusterings, clusterings, pweight, k);
		Q.push(e);
	}
    
	//print_queue(Q);
    
	Element* anotherE = Q.top();
	Q.pop();
    
	while(anotherE->clusterings.size() < k) {
        	anotherE->appendTop(anotherE->clusterings, clusterings, pweight, k);
        	Q.push(anotherE);
       		anotherE = Q.top();
        	Q.pop();
	}

	return anotherE;
}

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

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

double** load_matrix(const char *in_file_name, unsigned int* p) {
   	FILE  *in;
   	FILE  *out;
	double **pweight;
	int i, j,lo,k;
	double w;
	long e_count;
	unsigned int size = 0;

   	if ((in = fopen(in_file_name, "r")) == NULL) {
		printf("  fopen failed for input\n");
		exit(1);
	}
	fscanf(in, "%d %d", &size,&i);

	*p = size;

	e_count = ((long) size) * (size - 1) / 2;
	
	ALMF(pweight, size)

	for (k = 0; k < size; k++)
		ALF(*(pweight+k), size)
	
	for (lo = 0; lo < e_count; lo++) {
		fscanf(in, "%d %d %lf\n", &i, &j, &w);
		weight(i,j) = w;
		weight(j,i) = w;
	}
	return pweight;
}

int main(int argc, char **argv) {

	time_t tstart, tend; 
	tstart = time(0);

	char* in_file_name;

	unsigned int p;
	unsigned int k;

	int c;

	opterr = 0;

	while ((c = getopt (argc, argv, "k:f:")) != -1) {
	    switch (c)
	      {
	      case 'k':
		k = atoi(optarg);
		break;
	      case 'p':
		p = atoi(optarg);
		break;
	      case 'f':
		in_file_name = optarg;
		break;
	      case '?':
		if (optopt == 'f')
		  fprintf (stderr, "Option -%c requires a file name as argument.\n", optopt);
		else if (optopt == 'k')
		  fprintf (stderr, "Option -%c requires an positive integer as argument.\n", optopt);
		else if (isprint (optopt))
		  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		else
		  fprintf (stderr,
		           "Unknown option character `\\x%x'.\n",
		           optopt);
		return 1;
	      default:
		abort ();
	      }
	}

	double **pweight = load_matrix(in_file_name, &p);

	printf ("k = %d, p = %d, filename = %s\n", k, p, in_file_name);

	if(k >= p) {
		cout << "Invalid values for parameters k and p: " << k << ", " << p << endl;
		exit(1);
	}

	std::set<int> clusterings;

	for (int i=0; i<p; i++) {
		clusterings.insert(i);
	}

	Element* elementWithTopClusterings = rmcrag(clusterings, (const double**) pweight, p, k);

	std::cout << "Top k clusterings: ";
	std::set<int>::iterator it;
	for (it = elementWithTopClusterings->clusterings.begin(); it != elementWithTopClusterings->clusterings.end(); ++it) {
		std::cout << *it << " ";
	}
	std::cout << std::endl;

	tend = time(0); 

	cout << "RMCRAG took "<< difftime(tend, tstart) << " second(s)."<< endl;

	cout << "Solution cost: " << elementWithTopClusterings->cost << endl;

	/** Print cost of solution **/
	double value_from_sol = 0.0;
	for (itA = elementWithTopClusterings->clusterings.begin(); itA != elementWithTopClusterings->clusterings.end(); ++itA) {
		for (itB = elementWithTopClusterings->clusterings.begin(); itB != elementWithTopClusterings->clusterings.end(); ++itB) {
			if(*itA > *itB) {
				value_from_sol += pweight[*itA][*itB];
			}
		}
	}
	cout << "Cost of solution: " << value_from_sol << endl;
	
	return 0;
}

