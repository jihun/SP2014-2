#include <stdio.h>
#include <string>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MAX_ITER 20
#define ALPHA 0.3
// #define EPSILON 0.01

using namespace std;

FILE *BIOGRID;
FILE *prior;
FILE *answer;
FILE *roc_curve;
FILE *fp;

typedef string Gene;

typedef struct {
	Gene start;
	Gene end;
	double w;
} Edge;

bool cmp(const Edge &a, const Edge &b) { 
	return (a.start.compare(b.start) < 0);
} 

class cmp_by_f
{
public:
	bool operator() (const pair<Gene, double> &a, const pair<Gene, double> &b) const {
		return a.second < b.second;
	}	
};

class GeneNetworkSimulator 
{	
	unordered_set<Gene> genes;
	vector<Edge> graph;
	unordered_map<Gene, int> idx;
	unordered_map<Gene, double> Y;
	unordered_map<Gene, double> F;

	char gene1[1000], gene2[1000];
	int prior_count;

	set<Gene> ans, training_ans, validation_ans;
	priority_queue<pair<Gene, double>, vector<pair<Gene,double>>, cmp_by_f> pq;
public:

	void init_graph() 
	{
		double w;
		int count = 0;

		graph.clear();
		genes.clear();
	
		// input edges from BIOGRID
		BIOGRID =  fopen("BIOGRID.txt", "rt");
		while (fscanf(BIOGRID, "%s\t%s\t%lf\n", gene1, gene2, &w) != EOF) {
			Edge e1, e2;
			e1.start = gene1; e1.end = gene2; e1.w = w;
			e2.start = gene2; e2.end = gene1; e2.w = w;
			graph.push_back(e1);
			graph.push_back(e2);

			/*
			count++;
			if (count > 1000)
				break;
			*/
		}
		// build the graph from BIOGRID
		sort(graph.begin(), graph.end(), cmp);
		for (int i=0; i<graph.size(); i++) {
			if (idx.find(graph[i].start) == idx.end()) {
				idx[graph[i].start] = i;
				genes.insert(graph[i].start);
			}
		}
		
		prior = fopen("prior_knowledge.txt", "rt");
	
		// input prior knowledge
		prior_count = 0;
		while (fscanf(prior, "%s\t", gene1) != EOF) {
			fgets(gene2, 1000, prior);
			prior_count++;
		}
		//printf("%d\n",prior_count);
	}
	
	void init_answer()
	{
		answer = fopen("answer.txt", "rt");
		while (fscanf(answer, "%s\t%s\n", gene1, gene2) != EOF) {
			ans.insert(gene1);
		}
	}

	void init_prior_knowledge(bool cross_validation, double limit = 0.0)
	{
		int count = 0;

		if (cross_validation) {
			for (auto iter=training_ans.begin(); iter!=training_ans.end(); iter++) {
				Y[*iter] = 10;
			}
		}
		else {
			fseek(prior, 0L, SEEK_SET);
			while (fscanf(prior, "%s\t", gene1) != EOF) {
				fgets(gene2, 1000, prior);
				double v = (double)(prior_count-count) / prior_count;
				Y[gene1] = v;
				count++;

				if ((double)count/prior_count > limit) {
					break;
				}
			}
		}
	}

	void process()
	{
		int iter;
		int i, j;
		double neighbor_sum, f_sum;

		F.clear();
		
		for (iter=0; iter<MAX_ITER; iter++) {
			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				neighbor_sum = 0;
				for (i=idx[*gene_it]; i<graph.size() && graph[i].start==*gene_it; i++) {
				//	printf("%s %s\n", graph[i].start.c_str(), graph[i].end.c_str());
					neighbor_sum += (F[graph[i].end] * graph[i].w);
				}
				F[*gene_it] = ALPHA * neighbor_sum + (1 - ALPHA) * Y[*gene_it];
			}

			f_sum = 0;
			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				f_sum += F[*gene_it];
			}

			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				F[*gene_it] = F[*gene_it] / f_sum;
			}
		}
	}

	int output(bool cross_validation, int limit)
	{
		int match = 0;

		fp =  fopen("output.txt", "wt");
		
		pq = priority_queue<pair<Gene, double>, vector<pair<Gene,double>>, cmp_by_f>();
		for (auto iter=F.begin(); iter!=F.end(); iter++) {
			pq.push(*iter);
		}

		for (int i=0; i<limit; i++) {
			if (pq.empty()) {
				break;
			}
			string candidate = pq.top().first;

			if (cross_validation) {
				if (validation_ans.find(candidate) != validation_ans.end()) {
					match++;
				}
			}
			else {
				if (ans.find(candidate) != ans.end()) {
					//printf("%lf %s\n", pq.top().second, candidate.c_str());
					match++;
				}
			}

			pq.pop();
		}
		/*
		printf("%d\n", genes.size());
		printf("%d\n", match);
		*/
		return match;
	}

	void ROC_curve_test(double gap) 
	{
		int i, out;
		double limit;

		roc_curve = fopen("roc_curve.txt", "wt");

		init_graph();
		init_answer();

		for (limit = 0; limit <= 1.0; limit += gap) {
			init_prior_knowledge(false, limit);
			process();
			out = output(false, 500);
			fprintf(roc_curve, "%lf %d\n", limit, out);
		}
	}

	void cross_validation(double test_set_ratio)
	{
		int count;
		int out;

		init_graph();
		init_answer();

		count = 0;
		for (auto iter_ans=ans.begin(); iter_ans!=ans.end(); iter_ans++) {
			if ((double)count / ans.size() > test_set_ratio) {
				validation_ans.insert(*iter_ans);
			}
			else {
				training_ans.insert(*iter_ans);
				count++;
			}
		}
		
		init_prior_knowledge(true);
		process();

		out = output(true, 1000);
		printf("%d\n", out);
		printf("%d\n", validation_ans.size());
		
	}
};

int main()
{
	GeneNetworkSimulator gns;
	
	//gns.ROC_curve_test(0.05);
	gns.cross_validation(0.2);
	return 0;
}
