/*Markov Update in Graph: Input: A graph (first line is [r,c] followed by adjacency matrix*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define TRUE 1
#define FALSE 0
#define CONST_E 2.71828

int max(int a, int b){
	return (a >= b)?a:b; 
}

void priorProbOfGraph(int **graph, int N, double *P_graph){
	int i;
	for(i=0;i<N;i++){
		P_graph[i] = 1.0/N; //uniform prior
	}	
}

void updateBelief(int **graph, int N, double *P_graph, double *ClEval, int v, double threshold, int choice) {
	//Before BFS, all the nodes are unvisited
	int *visited = (int *)calloc(N,sizeof(int));
	if(!visited){
		printf("\nError: Memory leak for 'visited'.\n");
		exit(EXIT_FAILURE);
	}
	
	// queue implementation
	int front = 0;
	int rear = 0;
	int *Q = (int *)calloc(N,sizeof(int)); //this is the queue
	if(!Q){
		printf("\nError: Memory leak for Queue.\n");
		exit(EXIT_FAILURE);
	}
	
	//start visit with vertex v
	visited[v] = TRUE; // make vertex v visited
	int numBoxOpened = 0;
	// enqueue vertex v
	Q[rear] = v; // insert at rear
	rear += 1; // increment rear

	while(rear != front){
		// dequeue
		int u = Q[front];
		numBoxOpened += 1;
		//check whether the box has been found or not
		if(ClEval[u] >= threshold){
			FILE *fout = fopen("foundClause.txt","w");
			if(!fout){
				printf("\nError fout: File can't be opened for writing.\n");
				exit(EXIT_FAILURE);
			}
			//printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
			printf("\n\n[SUCCESS] Written in file as: \n[Box, Accuracy, NumBoxOpened]:[%d %lf %d]\n",u,ClEval[u],numBoxOpened);
			fprintf(fout,"%d %lf %d\n",u,ClEval[u],numBoxOpened);
			break;
		}
		
		//printf("\n%d %f", u, P_graph[u]);
		front += 1;

		// check adjacent nodes from u
		int i = 0, k = 0;
		//P_graph[u] = 0; //set its Prior to be updated
		//printf("\n%d(%lf) --> [ ",u,P_graph[u]);
		for (i = 0; i < N; i++) {
			// update the graph prob --> posterior
			if(graph[u][i]){
				switch(choice){
					case 1://adsorption concept
						P_graph[u] = P_graph[u] + P_graph[i];
						k += 1;
						break;
					case 2:
						P_graph[u] = pow(CONST_E, -(threshold - ClEval[u]));
						P_graph[i] = pow(CONST_E, -(threshold - ClEval[i]));
						break;
					case 3:
						P_graph[u] = P_graph[u] * P_graph[i];
						break;
					default:
						printf("\nThe value for \'choice\' is inadequate.\n");
						return;
				}
				//printf("%d(%lf) ",i,P_graph[i]);
			}
		
			// if there is adjacent vertex enqueue it
			if (!visited[i] && graph[u][i]) {
				Q[rear] = i;
				rear += 1;
				visited[i] = TRUE;
			}
		}
		//printf("] [k = %d]",k);
		if(choice == 1){
			if(!k){
				P_graph[u] = P_graph[u];
			}
			else{
				if(P_graph[u]/k > 1.0f){
					P_graph[u] = 1;
				}
				else{
					P_graph[u] = P_graph[u]/k;
				}
			}
		}
		//printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
	}
	printf("\n");
	free(Q);
	return;
}

int main(int *argc, char **argv) {

	//open the adjacency matrix file
	FILE *finp = fopen(argv[1],"r");
	if(!finp){
		printf("\nError finp: File read error.\n");
		exit(EXIT_FAILURE);
	}
	
	int N,row,col,i,j;
	fscanf(finp,"%d%d",&row,&col);
	N = row; //consider row and col as dummy variable

	//create the graph from the file
	int **graph;
	graph = (int **)calloc(N,sizeof(int *));
	if(!graph){
		printf("\nError: Memory leak for 'graph'.\n");
		exit(EXIT_FAILURE);
	}
	for(i=0;i<N;i++){
		graph[i] = (int *)calloc(N,sizeof(int));
		if(!graph[i]){
			printf("\nError: Memory leak for 'graph[i]'.\n");
			exit(EXIT_FAILURE);
		}
		for(j=0;j<N;j++){
			fscanf(finp,"%d",&graph[i][j]);
		}
	}
	fclose(finp);
	
	//make the graph undirected
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			graph[i][j]  = max(graph[i][j],graph[j][i]);
		}
	}
	
	//load the clause performance data
	FILE *fcleval = fopen(argv[2],"r");
	if(!fcleval){
		printf("\nError finp: File read error.\n");
		exit(EXIT_FAILURE);
	}
	double *ClEval = (double *)calloc(N,sizeof(double));
	if(!ClEval){
		printf("\nError: Memory leak for 'ClEval'.\n");
		exit(EXIT_FAILURE);
	}
	int clid,pos,neg,clauselen;//some will be unused
    double accuracy;
    i = 0;
    while (!feof(fcleval)) {
        fscanf(fcleval,"%d%d%d%d%lf",&clid,&pos,&neg,&clauselen,&accuracy);
        ClEval[clid-1] = accuracy;
        i += 1;
    }
    fclose(fcleval);
	
	//Prior probability of graph
	double *P_graph = (double *)calloc(N,sizeof(double));
	if(!P_graph){
		printf("\nError: Memory leak for 'P_graph'.\n");
		exit(EXIT_FAILURE);
	}
	
	// run bfs from 0th vertex (most general clause)
	priorProbOfGraph(graph, N, P_graph); //initialize Prior
	srand(time(NULL));
	
	//int startVertex = rand() % N;
	//convert arguments to adequate datatype
	int startVertex = atoi(argv[3]) - 1; //C standard indexing 1:0, 2:1...
	double threshold = atof(argv[4]);
	int choice = atoi(argv[5]);
	updateBelief(graph, N, P_graph, ClEval, startVertex, threshold, choice); //run updateBelief
	printf("\ndone.\n");
	
	free(P_graph);
	for(i=0;i<N;i++){
		free(graph[i]);
	}
	free(graph);
		
	return 0;
}
