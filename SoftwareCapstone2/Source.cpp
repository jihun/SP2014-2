#include <stdio.h>
#include <string>
#include <queue>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#define MAX_ITER 50
#define EPSILON 0.01
#define ALPHA 0.35

using namespace std;

FILE *BIOGRID = fopen("BIOGRID.txt", "rt");
FILE *prior = fopen("prior_knowledge.txt", "rt");
FILE *answer = fopen("answer.txt", "rt");
FILE *fp = fopen("output.txt", "wt");

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

unordered_set<Gene> genes;
vector<Edge> graph;
unordered_map<Gene, int> idx;
unordered_map<Gene, double> Y;
unordered_map<Gene, double> F;

set<Gene> ans;
priority_queue<pair<Gene, double>, vector<pair<Gene,double>>, cmp_by_f> pq;

void init() 
{
	char gene1[1000], gene2[1000];
	double w;
	int count;
	
	// input edges from BIOGRID
	while (fscanf(BIOGRID, "%s\t%s\t%lf\n", gene1, gene2, &w) != EOF) {
		Edge e1, e2;
		e1.start = gene1; e1.end = gene2; e1.w = w;
		e2.start = gene2; e2.end = gene1; e2.w = w;
		graph.push_back(e1);
		graph.push_back(e2);
	}
	sort(graph.begin(), graph.end(), cmp);
	for (int i=0; i<graph.size(); i++) {
		if (idx.find(graph[i].start) == idx.end()) {
			idx[graph[i].start] = i;
			genes.insert(graph[i].start);
		}
	}

	count = 0;
	while (fscanf(prior, "%s\t", gene1) != EOF) {
		fgets(gene2, 1000, prior);
		count++;
	}

	prior = fopen("prior_knowledge.txt", "rt");
	while (fscanf(prior, "%s\t", gene1) != EOF) {
		fgets(gene2, 1000, prior);

		Y[gene1] = count--;
	}

	while (fscanf(answer, "%s\t%s\n", gene1, gene2) != EOF) {
		ans.insert(gene1);
	}
}

void process()
{
	int iter;
	int i, j;
	
	for (iter=0; iter<MAX_ITER; iter++) {
		for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
			for (i=idx[*gene_it]; i<graph.size() && graph[i].start==*gene_it; i++) {
				F[graph[i].end] = 
					ALPHA * F[graph[i].start] * graph[i].w + 
					(1 - ALPHA) * Y[graph[i].end];
			}
		}
	}
}

void output()
{
	int match = 0;

	for (auto iter=F.begin(); iter!=F.end(); iter++) {
		pq.push(*iter);
	}

	for (int i=0; i<1000; i++) {
		string candidate = pq.top().first;
		
		if (ans.find(candidate) != ans.end()) {
			match++;
		}
		pq.pop();
		if (pq.empty()) {
			break;
		}
	}
	printf("%d\n", genes.size());
	printf("%d\n", match);
}

int main()
{
	init();
	process();
	output();

	return 0;
}
