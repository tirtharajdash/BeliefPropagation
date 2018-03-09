/*
	Belief Update in Graph: Input: A graph (first line is [r,c] followed by adjacency matrix
	Version 2
	Programmed by: Tirtharaj Dash, PhD Scholar, BITS Pilani, Goa Campus
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define TRUE 1
#define FALSE 0
#define CONST_E 2.71828
#define INF 99999999

struct stat st = {0};

//Global variables
int numBoxOpened = 0; //Number of boxes opened till now
int printFlag = FALSE; //Do you want to outcome of opening each box during search

//return max of two numbers
int func_max(int a, int b){
	return (a >= b)?a:b; 
}

//return max value of an array
double func_maxArr(double *arr, int N){
	double max = arr[0];
	int i;
	for(i=1;i<N;i++){
		if(max < arr[i])
			max = arr[i];
	}
	return max;
}

//return min value of an array
double func_minArr(double *arr, int N){
	double min = arr[0];
	int i;
	for(i=1;i<N;i++){
		if(min > arr[i])
			min = arr[i];
	}
	return min;
}

//return the index of smallest element in a double array 
int func_indexOfSmallestElemInArr(double *arr, int N){
	int s_index = 0, i;
	double minElem = arr[0];
	for(i=1;i<N;i++){
		if(minElem > arr[i]){
			minElem = arr[i];
			s_index = i;
		}
	}
	return s_index;
}

void priorProbOfGraph(int **graph, int N, double *P_graph){
	int i;
	for(i=0;i<N;i++){
		P_graph[i] = 1.0/N; //uniform prior
	}	
}

/*Random visit from one clause to another*/
int func_RandVisit(int **graph, int N, double *P_graph, double *ClEval, int v, double threshold, int updChoice){
	//initially, all the nodes are unvisited
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
	
	srand(time(NULL));
	
	//start visit with vertex v
	visited[v] = TRUE; // make vertex v visited
	numBoxOpened = 0; //no box has been opened till now
	// enqueue vertex v
	Q[rear] = v; // insert at rear
	rear += 1; // increment rear

	int clauseFound = FALSE;
	while(rear != front){
		// dequeue
		int u = Q[front];
		numBoxOpened += 1;
		//check whether the box has been found or not
		if(ClEval[u] >= threshold){
			if(printFlag)
				printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
			clauseFound = TRUE;
			return u;
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
				switch(updChoice){
					case 1: //adsorption concept
						P_graph[u] = P_graph[u] + P_graph[i];
						k += 1;
						break;
					case 2: //update belief based on closeness to threshold
						P_graph[u] = pow(CONST_E, -(threshold - ClEval[u]));
						P_graph[i] = pow(CONST_E, -(threshold - ClEval[i]));
						break;
					case 3: //Bayesian update P(i) = P(i) * \Pi_{j \in Adj(i)}{P(j|i)}
						P_graph[u] = P_graph[u] * P_graph[i];
						break;
				}
				//printf("%d(%lf) ",i,P_graph[i]);
			}
		}
		//printf("] [k = %d]",k);
		if(updChoice == 1){
			if(!k){
				P_graph[u] = P_graph[u];
			}
			else{
				P_graph[u] = P_graph[u]/k;
			}
		}
		if(printFlag)
			printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
		
		// generate a random unvisited node for next visit
		int cntGen = 0, randNext;
		while(cntGen < N){
			randNext = rand() % N;
			if(!visited[randNext]){
				Q[rear] = randNext;
				rear += 1;
				visited[randNext] = TRUE;
				break;
			}
			cntGen += 1;
		}
	}
	printf("\n");
	free(Q);
	if(!clauseFound)
		return (N+1);
}

/*Random visit from one clause to another by drawing a sample from updated distribution*/
int func_RandDrawnVisit(int **graph, int N, double *P_graph, double *ClEval, int v, double threshold, int updChoice){
	//initially, all the nodes are unvisited
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
	
	srand(time(NULL));
	
	//start visit with vertex v
	visited[v] = TRUE; // make vertex v visited
	numBoxOpened = 0; //no box has been opened till now
	// enqueue vertex v
	Q[rear] = v; // insert at rear
	rear += 1; // increment rear

	int clauseFound = FALSE;
	while(rear != front){
		// dequeue
		int u = Q[front];
		numBoxOpened += 1;
		//check whether the box has been found or not
		if(ClEval[u] >= threshold){
			if(printFlag)
				printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
			clauseFound = TRUE;
			return u;
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
				switch(updChoice){
					case 1: //adsorption concept
						P_graph[u] = P_graph[u] + P_graph[i];
						k += 1;
						break;
					case 2: //update belief based on closeness to threshold
						P_graph[u] = pow(CONST_E, -(threshold - ClEval[u]));
						P_graph[i] = pow(CONST_E, -(threshold - ClEval[i]));
						break;
					case 3: //Bayesian update P(i) = P(i) * \Pi_{j \in Adj(i)}{P(j|i)}
						P_graph[u] = P_graph[u] * P_graph[i];
						break;
				}
				//printf("%d(%lf) ",i,P_graph[i]);
			}
		}
		//printf("] [k = %d]",k);
		if(updChoice == 1){
			if(!k){
				P_graph[u] = P_graph[u];
			}
			else{
				P_graph[u] = P_graph[u]/k;
			}
		}
		if(printFlag)
			printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);

		//Generate next visit clause: our intention is that we want it to be within the upper bound of N
		int cntGen = 0;
		//maximum number of tosses is N
		while(cntGen < N){		
			int qCnt = 0;
			//generate a random number to be used for sampling from updated distribution within the max-min range
			double randProb = func_minArr(P_graph, N) + ((double)rand()/(double)RAND_MAX)*(func_maxArr(P_graph, N) - func_minArr(P_graph, N));
			int *qIndex = (int *)calloc(N, sizeof(int));
			if(!qIndex){
				printf("\nError: Memory Leak in func_randDrawnVisit(): qIndex[].\n");
				exit(EXIT_FAILURE);
			}
			double *qProb = (double *)calloc(N, sizeof(double));
			if(!qProb){
				printf("\nError: Memory Leak in func_randDrawnVisit(): qProb[].\n");
				exit(EXIT_FAILURE);
			}
			for(i=0;i<N;i++){
				if(P_graph[i] >= randProb){
					qProb[qCnt] = P_graph[i];
					qIndex[qCnt] = i; //qualified indices, qCnt: qualified count
					qCnt += 1;
				}
			}
			//printf("qCnt = %d\n",qCnt);
			if(!qCnt){
				continue; //generate new random number for sampling
			}
			
			int tempNextNode = func_indexOfSmallestElemInArr(qProb, qCnt);
			if(!visited[qIndex[tempNextNode]]){
				Q[rear] = qIndex[tempNextNode];
				rear += 1;
				visited[qIndex[tempNextNode]] = TRUE;
				break;
			}

			cntGen += 1;
		}
		if(cntGen == N){
			clauseFound = FALSE;
			return (N+1);
		}
	}
	printf("\n");
	free(Q);
	if(!clauseFound)
		return (N+1);
}

/*Visiting to a random neighbor of till found clauses until the goal clause: uses RandomizedQueue*/
int func_QueueBasedVisit(int **graph, int N, double *P_graph, double *ClEval, int v, double threshold, int updChoice) {
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
	numBoxOpened = 0; //no box has been opened till now
	// enqueue vertex v
	Q[rear] = v; // insert at rear
	rear += 1; // increment rear
	srand(time(NULL));

	int clauseFound = FALSE;
	while(rear != front){
		// dequeue
		int u = Q[front];
		numBoxOpened += 1;
		//check whether the box has been found or not
		if(ClEval[u] >= threshold){
			if(printFlag)
				printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
			clauseFound = TRUE;
			return u;
		}
		
		//printf("\n%d %f", u, P_graph[u]);
		front += 1;

		// check adjacent nodes from u
		int i = 0, k = 0, neighCnt = 0, *neighNodes; //neighCnt: Neighbor count
		//P_graph[u] = 0; //set its Prior to be updated
		//printf("\n%d(%lf) --> [ ",u,P_graph[u]);
		neighNodes = (int *)calloc(N, sizeof(int));
		if(!neighNodes){
			printf("\nError: Memory leak in func_QueueBasedVisit: neighNodes[].\n");
			exit(EXIT_FAILURE);
		}
		for (i = 0; i < N; i++) {
			// update the graph prob --> posterior
			if(graph[u][i]){
				neighNodes[neighCnt] = i;
				neighCnt += 1;
				switch(updChoice){
					case 1: //adsorption concept
						P_graph[u] = P_graph[u] + P_graph[i];
						k += 1;
						break;
					case 2: //update belief based on closeness to threshold
						P_graph[u] = pow(CONST_E, -(threshold - ClEval[u]));
						P_graph[i] = pow(CONST_E, -(threshold - ClEval[i]));
						break;
					case 3: //Bayesian update P(i) = P(i) * \Pi_{j \in Adj(i)}{P(j|i)}
						P_graph[u] = P_graph[u] * P_graph[i];
						break;
				}
				//printf("%d(%lf) ",i,P_graph[i]);
			}
		}
		//printf("] [k = %d]",k);
		if(updChoice == 1){
			if(!k){
				P_graph[u] = P_graph[u];
			}
			else{
				P_graph[u] = P_graph[u]/k;
			}
		}
		if(printFlag)
			printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);

		//suffle the neighboring nodes randomly for random visits
		/*for(i=0;i<neighCnt;i++){
			printf("%d ", neighNodes[i]);
		}
		printf("\n");*/
		srand ( time(NULL) );
		for(i = neighCnt-1; i > 0; i--)
		{
			int j = rand() % (i+1);
			int temp = neighNodes[i];
			neighNodes[i] = neighNodes[j];
			neighNodes[j] = temp;
		}
		/*for(i=0;i<neighCnt;i++){
			printf("%d ", neighNodes[i]);
		}
		break;*/
		for(i=0;i<neighCnt;i++){
     		// if there is unvisited adjacent vertex enqueue it
			if (!visited[neighNodes[i]] && graph[u][neighNodes[i]]) {
				Q[rear] = neighNodes[i];
				rear += 1;
				visited[neighNodes[i]] = TRUE;
			}
		}
		
		free(neighNodes);
	}
	printf("\n");
	free(Q);
	if(!clauseFound)
		return (N+1);
}

/*Visiting neighbor clause with higher belief: uses a Priority Queue*/
int func_PriorityQueueBasedVisit(int **graph, int N, double *P_graph, double *ClEval, int v, double threshold, int updChoice) {
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
	numBoxOpened = 0; //no box has been opened till now
	// enqueue vertex v
	Q[rear] = v; // insert at rear
	rear += 1; // increment rear
	srand(time(NULL));
	
	int clauseFound = FALSE;
	while(rear != front){
		// dequeue
		int u = Q[front];
		numBoxOpened += 1;
		//check whether the box has been found or not
		if(ClEval[u] >= threshold){
			if(printFlag)
				printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);
			clauseFound = TRUE;
			return u;
		}
		
		//printf("\n%d %f", u, P_graph[u]);
		front += 1;

		// check adjacent nodes from u
		int i = 0, k = 0, neighCnt = 0, *neighNodes; //neighCnt: Neighbor count
		double *P_neighNodes; //probability of neighboring nodes
		//P_graph[u] = 0; //set its Prior to be updated
		//printf("\n%d(%lf) --> [ ",u,P_graph[u]);
		neighNodes = (int *)calloc(N, sizeof(int));
		if(!neighNodes){
			printf("\nError: Memory leak in func_PriorityQueueBasedVisit: neighNodes[].\n");
			exit(EXIT_FAILURE);
		}
		P_neighNodes = (double *)calloc(N, sizeof(double));
		if(!neighNodes){
			printf("\nError: Memory leak in func_PriorityQueueBasedVisit: neighNodes[].\n");
			exit(EXIT_FAILURE);
		}
		for (i = 0; i < N; i++) {
			// update the graph prob --> posterior
			if(graph[u][i]){
				neighNodes[neighCnt] = i;
				switch(updChoice){
					case 1: //adsorption concept
						P_graph[u] = P_graph[u] + P_graph[i];
						k += 1;
						break;
					case 2: //update belief based on closeness to threshold
						P_graph[u] = pow(CONST_E, -(threshold - ClEval[u]));
						P_graph[i] = pow(CONST_E, -(threshold - ClEval[i]));
						break;
					case 3: //Bayesian update P(i) = P(i) * \Pi_{j \in Adj(i)}{P(j|i)}
						P_graph[u] = P_graph[u] * P_graph[i];
						break;
				}
				P_neighNodes[neighCnt] = P_graph[i];
				neighCnt += 1;
				//printf("%d(%lf) ",i,P_graph[i]);
			}
		}
		//printf("] [k = %d]",k);
		if(updChoice == 1){
			if(!k){
				P_graph[u] = P_graph[u];
			}
			else{
				P_graph[u] = P_graph[u]/k;
			}
		}
		if(printFlag)
			printf("%d %lf %lf %d\n",u,P_graph[u],ClEval[u],numBoxOpened);

		//Priority based visit to one of the neighbor of the presently dequeued node
		/*for(i=0;i<neighCnt;i++){
			printf("(%lf, %d) ", P_neighNodes[i],neighNodes[i]);
		}
		printf("\n");*/
		int j;
		int swapped = FALSE;
		for (i = 0; i < neighCnt-1; i++){
			swapped = FALSE;
			for (j = 0; j < neighCnt-i-1; j++){
				if (P_neighNodes[j] < P_neighNodes[j+1]){
					double temp = P_neighNodes[j];
					P_neighNodes[j] = P_neighNodes[j+1];
					P_neighNodes[j+1] = temp;
					//swap index
					int tempIndex = neighNodes[j];
					neighNodes[j] = neighNodes[j+1];
					neighNodes[j+1] = tempIndex;
					
					swapped = TRUE;
				}
			}
			if(!swapped)
				break;
		}
				
		/*for(i=0;i<neighCnt;i++){
			printf("(%lf, %d) ", P_neighNodes[i],neighNodes[i]);
		}
		printf("\n\n");*/
		
		for(i=0;i<neighCnt;i++){
     		// if there is adjacent vertex enqueue it
			if (!visited[neighNodes[i]] && graph[u][neighNodes[i]]) {
				Q[rear] = neighNodes[i];
				rear += 1;
				visited[neighNodes[i]] = TRUE;
			}
		}
		free(neighNodes);
	}
	printf("\n");
	free(Q);
	if(!clauseFound)
		return (N+1);
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
			//avoid self-loop
			if(i!=j){
				graph[i][j]  = func_max(graph[i][j],graph[j][i]);
			}
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
		
	//for bunch of simulation
	int numRandVertex = atoi(argv[3]);
	double threshold = atof(argv[4]); //threshold for box quality
	
	//create files for storing results fij:: f: file, i: algo 1-4, j: choice 1--3 
	char *fname = (char *)calloc(20,sizeof(char));
	if(!fname){
		printf("\nError: Memory leak in main() for 'fname'.\n");
		exit(EXIT_FAILURE);
	}
	
	//check if the results directory exists, else create one
	if (stat("./results", &st) == -1) {
    	mkdir("./results", 0700);
	}
	char *dirname = (char *)calloc(20,sizeof(char));
	if(!dirname){
		printf("\nError: Memory leak in main() for 'dirname'.\n");
		exit(EXIT_FAILURE);
	}
	sprintf(dirname,"./results/%d",numRandVertex);
	//create a subdirectory
	if (stat(dirname, &st) == -1) {
    	mkdir(dirname, 0700);
	}
		
	//srand(time(NULL));
	int k=0, u;
	for(k=0;k<numRandVertex;k++){
		int startVertex = rand() % N;
		for(i=1;i<=4;i++){
			switch(i){
				case 1:
					for(j=1;j<=3;j++){
						switch(j){
							case 1:
								sprintf(fname,"%s/rand_ads.dat", dirname);
								FILE *f11 = fopen(fname,"a+");
								u = func_RandVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f11,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f11,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 2:
								sprintf(fname,"%s/rand_exp.dat", dirname);
								FILE *f12 = fopen(fname,"a+");
								u = func_RandVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f12,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f12,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 3:
								sprintf(fname,"%s/rand_bayes.dat", dirname);
								FILE *f13 = fopen(fname,"a+");
								u = func_RandVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f13,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f13,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
						}
					}
					break;
				case 2:
					for(j=1;j<=3;j++){
						switch(j){
							case 1:
								sprintf(fname,"%s/distr_ads.dat", dirname);
								FILE *f21 = fopen(fname,"a+");
								u = func_RandDrawnVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f21,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f21,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 2:
								sprintf(fname,"%s/distr_exp.dat", dirname);
								FILE *f22 = fopen(fname,"a+");
								u = func_RandDrawnVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f22,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f22,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 3:
								sprintf(fname,"%s/distr_bayes.dat", dirname);
								FILE *f23 = fopen(fname,"a+");
								u = func_RandDrawnVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f23,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f23,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
						}
					}
					break;
				case 3:
					for(j=1;j<=3;j++){
						switch(j){
							case 1:
								sprintf(fname,"%s/Q_ads.dat", dirname);
								FILE *f31 = fopen(fname,"a+");
								u = func_QueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f31,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f31,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 2:
								sprintf(fname,"%s/Q_exp.dat", dirname);
								FILE *f32 = fopen(fname,"a+");
								u = func_QueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f32,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f32,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 3:
								sprintf(fname,"%s/Q_bayes.dat", dirname);
								FILE *f33 = fopen(fname,"a+");
								u = func_QueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f33,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f33,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
						}
					}
					break;
				case 4:
					for(j=1;j<=3;j++){
						switch(j){
							case 1:
								sprintf(fname,"%s/PQ_ads.dat", dirname);
								FILE *f41 = fopen(fname,"a+");
								u = func_PriorityQueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f41,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f41,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 2:
								sprintf(fname,"%s/PQ_exp.dat", dirname);
								FILE *f42 = fopen(fname,"a+");
								u = func_PriorityQueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f42,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f42,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
							case 3:
								sprintf(fname,"%s/PQ_bayes.dat", dirname);
								FILE *f43 = fopen(fname,"a+");
								u = func_PriorityQueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, j);
								if(u == (N+1)){
									printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
									fprintf(f43,"%d %lf %d %d\n", u, 0.0, INF, 0);
								}
								else{
									printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
									fprintf(f43,"%d %lf %d %d\n", u, ClEval[u], numBoxOpened, 1);
								}
								break;
						}
					}
					break;
			}
			printf("k = %d, i = %d, j = %d\n", k, i, j);
		}
	}
	
	/*
	//For arg based simulation
	int startVertex = atoi(argv[3]); //C standard indexing 1:0, 2:1...
	double threshold = atof(argv[4]); //threshold for box quality
	int beliefUpdChoice = atoi(argv[5]); //belief update choice 1: Adsorption, 2: Exponential, 3: Bayesian
	
	printFlag = FALSE;
	
	//int u = func_RandVisit(graph, N, P_graph, ClEval, startVertex, threshold, beliefUpdChoice);
	//int u = func_RandDrawnVisit(graph, N, P_graph, ClEval, startVertex, threshold, beliefUpdChoice);
	int u = func_QueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, beliefUpdChoice);
	//int u = func_PriorityQueueBasedVisit(graph, N, P_graph, ClEval, startVertex, threshold, beliefUpdChoice);
	
	if(u == (N+1))
		printf("[FAILURE]\t[Box, Accuracy, NumBoxOpened, 'Tail']:[%d %lf %d]\n", u, 0.0, INF);
	else
		printf("[SUCCESS]\t[Box, Accuracy, NumBoxOpened, 'Head']:[%d %lf %d]\n", u, ClEval[u], numBoxOpened);
	*/

	printf("\n[SEARCH COMPLETE]\n");
	
	free(P_graph);
	for(i=0;i<N;i++){
		free(graph[i]);
	}
	free(graph);
		
	return 0;
}
