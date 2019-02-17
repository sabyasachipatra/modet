// Main Function
#include <iostream>
using namespace std;
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Graph.h"
#include "GraphIsomor.h"
#include "Utility.h"
#include "Timer.h"
#include "SubgraphTree.h"
#include "ExpansionTree.h"


#define MAX_BUF 256 // Maximum string buffer size

void parse_cmdline(int argc, char **argv);
void check_inputs();
void initialize();
void prepare_graph();
void compute_real();
void compute_results();
void run_modet();
int compare_results(const void *a, const void *b);

// Variable declarations
char graph_file[MAX_BUF];
char method[MAX_BUF];
int motif_size;
int rand_number;
int num_exchanges;
int num_tries;
int threshold;
Graph *g;
SubgraphTree realSG;
ExpansionTree *realET;
ExpansionTree *randET;

int main(int argc, char **argv) {
	cout << "This program is developed by Sabyasachi Patra, IIIT Bhubaneswar, India" << endl;
	Utility::output_msg("This program is developed by Sabyasachi Patra, IIIT Bhubaneswar, India");
	initialize();
	parse_cmdline(argc, argv);
	check_inputs();
	prepare_graph();
	compute_real();
	compute_results();
	GraphIsomor::finishNauty();
	return 0;
}

// Initialize everything
void initialize() {
	motif_size=-1;
	rand_number=-1;
	num_exchanges = 1;
	num_tries = 10;
	threshold=0;
	strcpy(method, "none");
	strcpy(graph_file, "none");
	FILE *fp;
	fp = fopen("Output.txt","w");
	fprintf (fp, "OUTPUT");
	fprintf (fp, "\n");
	fclose(fp);
	fp = fopen("Error.txt","w");
	fprintf (fp, "ERROR");
	fprintf (fp, "\n");
	fclose(fp);
	GraphIsomor::initNauty(motif_size);
}

// Parse all command line arguments
void parse_cmdline(int argc, char **argv) {
    for (int i=1; i<argc; i++) {
		// Graph file
		if (!strcmp("-g",argv[i]) || !strcmp("--graph",argv[i])) {
			strcpy(graph_file, argv[++i]);
		}
		// Size of motifs to consider
		else if (!strcmp("-s",argv[i]) || !strcmp("--size",argv[i])) {
			motif_size = atoi(argv[++i]);
		}
		// Method for set of subgraphs
		else if (!strcmp("-m",argv[i]) || !strcmp("--method",argv[i])) {
			strcpy(method, argv[++i]);
		}
		// Method for set number o random graphs
		else if (!strcmp("-r",argv[i]) || !strcmp("--random",argv[i])) {
			rand_number = atoi(argv[++i]);
		}
		// Number of exchanges per edge
		else if (!strcmp("-e",argv[i]) || !strcmp("--exchanges",argv[i])) {
			num_exchanges = atoi(argv[++i]);
		}
		// Number of tries per node
		else if (!strcmp("-t",argv[i]) || !strcmp("--tries",argv[i])) {
			num_tries = atoi(argv[++i]);
		}
		// Threshold value
		else if (!strcmp("-th",argv[i]) || !strcmp("--threshold",argv[i])) {
			threshold = atoi(argv[++i]);
		}
	}
}

void check_inputs() {
	if (strcmp(method, "modet") != 0 && strcmp(method, "mdet") != 0) {
		cout << "invalid method" << endl;
		Utility::error_msg("invalid method");
		exit(1);
	}
	if (strcmp(method, "modet") == 0) {
		if (motif_size < 3 || motif_size > 10) {
			cout << "invalid motf size" << endl;
			Utility::error_msg("invalid motf size");
			exit(1);
		}
	}
	if (strcmp(method, "mdet") == 0) {
		if (motif_size < 3 || motif_size > 15) {
			cout << "invalid motf size" << endl;
			Utility::error_msg("invalid motf size");
			exit(1);
		}
	}
	if (strcmp(graph_file, "none") == 0) {
		cout << "no input graph file" << endl;
		Utility::error_msg("no input graph file");
		exit(1);
	}
}

// Prepare the real graph for computation
void prepare_graph() {
	g = new Graph();
	// Read the graph file
	Graph::readGraphFile(g, graph_file);
	// sort and create neighbours array
	g->sortNeighbours();
	g->makeNeighboursArray();
	// Print chosen parameters
	printf("motif size: %d\n", motif_size);
	printf("graph file: %s\n", graph_file);
	printf("%d nodes, %d edges\n", g->numberNodes(), g->numberEdges());
	threshold = int(g->numberNodes()*threshold/100.0);
	cout << "threshold = " << threshold << endl;
}

// Count subgraphs on real network
void compute_real() {
	// Print method name
	if (strcmp(method, "modet") == 0)
		printf("Method: MODET on real network\n");
	else if (strcmp(method, "mdet") == 0)
		printf("Method: MDET on real network\n");
	else
		printf("Invalid method\n");

	// Compute frequency
	printf("\nCounting subgraph frequency on 'REAL NETWORK'\n");
	Timer::start(0);
	if (strcmp(method, "modet") == 0)
		run_modet();
	//else if (strcmp(method, "mdet") == 0)
	//	run_mdet();
	Timer::stop(0);
	printf("%d subgraphs, ", realSG.countSubgraphs());
	printf("%f embeddings\n", realSG.countEmbeddings());
	printf("Time elapsed: %.6fs\n", Timer::elapsed(0));
}

// Run MODET method on graph 'g' and store results on SubgraphTree 'realSG'
void run_modet() {
	Timer::start(0);
	SubgraphTree isomorSG;
	realET = new ExpansionTree();
	realET->create(motif_size, &isomorSG);
	Timer::stop(0);
	printf("Creation time: %.2f\n", Timer::elapsed(0));
	Timer::start(0);
	realET->census(motif_size, g, &realSG);
	//realET->printEmbeddings(motif_size);
	Timer::stop(0);
	printf("census time: %.2f\n", Timer::elapsed(0));
}

// Compute random networks and result
void compute_results() {
	int i, j;
	map<string, int>::iterator ii;
	// Create map and init results
	map<string, int> realMap;
	int resultSize;
	realET->inducedFreq(g,motif_size,threshold,&realSG);  
	realSG.populateMap(&realMap, motif_size);
	ResultType result[realMap.size()];	
	for (ii=realMap.begin(), i=0; ii!=realMap.end(); ii++) {
		if (ii->second >= threshold) {		
			result[i].s = strdup((ii->first).c_str());
			result[i].freq = ii->second;
			result[i].z_score = 0;
			result[i].avgf_rand = 0;
			result[i].devf_rand = 0;
			i++;
		}
	}
	resultSize=i;
	
	// Do we have random networks to compute?
	if (rand_number > 0) {
		map<string, int> randMap[rand_number];
		// Generate all random networks
		printf("Computing random networks: ");
		for (i=0; i<rand_number; i++) {
			cout << "random count = " << i+1 << endl;    
			// Create new random network from previous one
			int nswap=0;
			Graph::randomGraph(g, num_exchanges, num_tries,&nswap);
			//cout << "after:" << nswap << " number of edge swappings\n";			
			g->sortNeighbours();
			g->makeNeighboursArray();
			SubgraphTree randSG;
			SubgraphTree isomorSG;
			randET = new ExpansionTree();
			randET->create(motif_size, &isomorSG);
			randET->census(motif_size, g, &randSG);
			randET->inducedFreq(g,motif_size,threshold,&randSG);  
			delete randET;
			randSG.populateMap(&randMap[i], motif_size);
		}
		// Compute significance
		for (i=0; i<resultSize; i++) {
			// Average frequency
			cout << result[i].freq << " : ";
			double avg = 0;
			for (j=0; j<rand_number; j++) {
				avg += randMap[j][result[i].s];
				cout << randMap[j][result[i].s] << "  ";
			}
			cout << endl;
			avg /= rand_number;
			// Standard deviation
			double dev=0;
			for (j=0; j<rand_number; j++)
				dev += (randMap[j][result[i].s]-avg)*(randMap[j][result[i].s]-avg)/double(rand_number-1);
			dev = sqrt(dev);
			double zscore; 
			if (dev != 0)
				zscore = (result[i].freq - avg)/dev;
			else
				zscore = 0;
			//cout << result[i].freq << ", " << avg << ", " << dev << ", " << zscore << endl; 
			result[i].avgf_rand = avg;
			result[i].devf_rand = dev;
			result[i].z_score = zscore;
		}
	}
	// Sort results
	qsort(result, resultSize, sizeof(ResultType), compare_results);
	FILE *fp;
	fp = fopen("result.txt", "w");
	for (i=0; i<resultSize; i++) {
		if (rand_number > 0) {
			fprintf(fp,"%s : %d : %f\n",result[i].s,result[i].freq,result[i].z_score);
			cout << result[i].s << " : " << result[i].freq << " : " << result[i].z_score << endl;
		}
		else {
			fprintf(fp,"%s : %d\n",result[i].s,result[i].freq);
			cout << result[i].s << " : " << result[i].freq << endl;
		}
	}
 	fprintf(fp,"\nnumber of motifs = %d\n",resultSize);
 	cout << "\nnumber of motifs = " << resultSize << endl;
 	fclose(fp);
}

// Compare two different motif results (for sorting)
int compare_results(const void *a, const void *b) {
  ResultType *r1 = (ResultType *)a;
  ResultType *r2 = (ResultType *)b;

  if (r1->z_score > r2->z_score) return -1;
  if (r1->z_score < r2->z_score) return +1;
  if (r1->freq > r2->freq) return -1;
  if (r1->freq < r2->freq) return +1;
  return 0;
}
