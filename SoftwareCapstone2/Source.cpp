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

	set<Gene> ans, training_ans, validation_ans;
	priority_queue<pair<Gene, double>, vector<pair<Gene,double>>, cmp_by_f> pq;
public:

  void init_graph(bool entrez) {
    graph.clear();
    genes.clear();
    idx.clear();

    double w;
    FILE *fp;
    
    if (entrez) fp = fopen("BIOGRID-Entrez.txt", "rt");
    else fp = fopen("BIOGRID.txt", "rt");

    while (fscanf(fp, "%s\t%s\t%lf\n", gene1, gene2, &w) != EOF) {
      Edge e1, e2;
      e1.start = gene1; e1.end = gene2; e1.w = w;
      e2.start = gene2; e2.end = gene1; e2.w = w;

      graph.push_back(e1);
      graph.push_back(e2);
    }

    fclose(fp);

    sort(graph.begin(), graph.end(), cmp);
    for (int i = 0; i<graph.size(); i++) {
      if (idx.find(graph[i].start) == idx.end()) {
        idx[graph[i].start] = i;
        genes.insert(graph[i].start);
      }
    }
  }
	
	void init_answer(bool entrez)
	{
    FILE *fp;

    ans.clear();

    if (entrez) {
      fp = fopen("answer-Entrez.txt", "rt");
      while (fscanf(fp, "%s", gene1) != EOF) {
        ans.insert(gene1);
      }
    } else {
      fp = fopen("answer.txt", "rt");
      while (fscanf(fp, "%s", gene1) != EOF) {
        ans.insert(gene1);
        fgets(gene2, 1000, fp);
      }
    }

    fclose(fp);
	}

  void init_prior_knowledge_from_answer_set() {
    for (auto iter = training_ans.begin(); iter != training_ans.end(); iter++) {
      Y[*iter] = 1;
    }
  }

  void init_prior_knowledge_from_miRNA(bool entrez) {
    FILE *fp;

    if (entrez) fp = fopen("prior_knowledge_miRNA-Entrez.txt", "rt");
    else fp = fopen("prior_knowledge_miRNA.txt", "rt");

    while (fscanf(fp, "%s\n", gene1) != EOF) {
      Y[gene1] = 1;
    }

    fclose(fp);
  }

  void init_prior_knowledge_from_publication(bool entrez, double limit = 0.0) {
    FILE *fp;
    int prior_count = 0, rank = 0;

    if (entrez) {
      fp = fopen("prior_knowledge-Entrez.txt", "rt");
      while (fscanf(fp, "%s\n", gene1) != EOF) {
        ++prior_count;
      }
    } else {
      fp = fopen("prior_knowledge.txt", "rt");
      while (fscanf(fp, "%s", gene1) != EOF) {
        fgets(gene2, 1000, fp);
        ++prior_count;
      }
    }

    fseek(fp, 0, SEEK_SET);
    while (fscanf(fp, "%s", gene1) != EOF) {
      if (!entrez) {
        fgets(gene2, 1000, fp);
      }
      double v = (double) (prior_count - rank) / prior_count;
      for (int i = 0, len = strlen(gene1); i < len; ++i) {
        gene1[i] = toupper(gene1[i]);
      }
      Y[gene1] = v;
      ++rank;

      if ((double) rank / prior_count > limit) {
        break;
      }
    }

    fclose(fp);
  }

	void process()
	{
		int iter;
		int i, j;
		double neighbor_sum, f_sum;
    unordered_map<Gene, double> F2;

		F.clear();
		
		for (iter=0; iter<MAX_ITER; iter++) {
			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				neighbor_sum = 0;
				for (i=idx[*gene_it]; i<graph.size() && graph[i].start==*gene_it; i++) {
					neighbor_sum += (F[graph[i].end] * graph[i].w);
				}
				F2[*gene_it] = ALPHA * neighbor_sum + (1 - ALPHA) * Y[*gene_it];
			}

			f_sum = 0;
			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				f_sum += F2[*gene_it];
			}

			for (auto gene_it=genes.begin(); gene_it!=genes.end(); gene_it++) {
				F[*gene_it] = F2[*gene_it] / f_sum;
			}
		}

    pq = priority_queue<pair<Gene, double>, vector<pair<Gene, double>>, cmp_by_f>();
    for (auto iter = F.begin(); iter != F.end(); iter++) {
      pq.push(*iter);
    }
	}

	void match_cross_validation(int limit)
	{
    int matchValidation = 0, matchTraining = 0;

		FILE *fp =  fopen("output.txt", "wt");

    for (int i = 0; !pq.empty(); ++i) {
      string candidate = pq.top().first;
      double score = pq.top().second;

      bool inValidation = validation_ans.find(candidate) != validation_ans.end();
      bool inTraining = training_ans.find(candidate) != training_ans.end();
      bool inAnswer = inValidation || inTraining;

      if (inValidation && i < limit) ++matchValidation;
      if (inTraining && i < limit) ++matchTraining;

      fprintf(fp, "%s\t%.9lf\t%d\n", candidate.c_str(), score, inValidation);

      pq.pop();
    }

    fclose(fp);

    printf("inTraining = %d\tinValidation = %d\n", matchTraining, matchValidation);
	}

  int match_answer_set(int limit) {
    int match = 0;

    FILE *fp = fopen("output.txt", "wt");

    for (int i = 0; !pq.empty(); ++i) {
      string candidate = pq.top().first;
      double score = pq.top().second;

      bool inAnswer = ans.find(candidate) != ans.end();

      if (inAnswer && i < limit) ++match;

      fprintf(fp, "%s\t%.9lf\t%d\n", candidate.c_str(), score, inAnswer);

      pq.pop();
    }

    fclose(fp);

    printf("inAnswer = %d\n", match);

    return match;
  }

	void cross_validation(bool entrez, double training_set_ratio)
	{
		int count;
		int out;

		init_graph(entrez);
		init_answer(entrez);

    training_ans.clear();
    validation_ans.clear();

		count = 0;
		for (auto iter_ans=ans.begin(); iter_ans!=ans.end(); iter_ans++) {
			if ((double)count / ans.size() > training_set_ratio) {
				validation_ans.insert(*iter_ans);
			}
			else {
				training_ans.insert(*iter_ans);
				count++;
			}
		}
		
    //init_prior_knowledge_from_publication(true, 0.2); // if uncomment this line, prior knowledge function Y is initialized with data from publication
    //init_prior_knowledge_from_miRNA(true); // if uncomment this line, prior knowledge function Y is initialized with data from miRNA-Disease-DNA link 
		init_prior_knowledge_from_answer_set();
		process();

		match_cross_validation(200);
		printf("validation set size = %d\n", validation_ans.size());

    printf("\n");
	}
};

int main()
{
	GeneNetworkSimulator gns;
	
  //for (int i = 1; i <= 10; ++i) {
	  gns.cross_validation(true, 0.2);
  //}

  /*gns.init_graph(true);
  gns.init_answer(true);
  //gns.init_prior_knowledge_from_miRNA(true);
  gns.init_prior_knowledge_from_publication(true, 0.2);
  gns.process();
  gns.match_answer_set(100);*/
	return 0;
}
