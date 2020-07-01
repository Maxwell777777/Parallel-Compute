#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <cstddef>

//#define NODES_NUM 1024000
#define MAX_EDGE_NUM_PER_NODE 10
#define ITERATION 100

using namespace std;
/* pagerank */

double max(long double x, long double y){
    return x>y?x:y;
}

double min(long double x, long double y){
    return x>y?y:x;
}

void print_rank(long double *rank, int NODES_NUM){
    printf("Rank value: ");
    for (int i = 0; i < NODES_NUM;i++){
        printf("%.10Lf ",rank[i]);
    }
    printf("\n");
}


int rand_edge_num(int *indegree_per_node, int nodes_num, int max_edge_num_per_node){
    int edge_num = 0;
    int temp_edge_num;
    for (int i = 0; i < nodes_num; i++){
        srand(i);
        temp_edge_num = (int)(rand()%(max_edge_num_per_node)+1);
        indegree_per_node[i] = temp_edge_num;
        edge_num += temp_edge_num;
    }
    return edge_num;
}


void rand_graph(int *edge_start, int *edge_end, int *indegree_per_node, int *outdegree_per_node, int node_num){
    int current_edge_num = 0;
    int temp_start;
    printf("Edges: ");
    for (int end = 0; end < node_num; end++){
        //rank[i] = 1.0/node_num;
        int *startnodes = new int[indegree_per_node[end]+1]; // startnodes[0] = self
        startnodes[0] = end;
        for(int j = 1; j <= indegree_per_node[end]; j++){
            edge_end[current_edge_num] = end;
            srand(current_edge_num + 1);
            temp_start = (int)(rand() % node_num);
            // avoid shart == end  or repeat
            for (int k = 0; k < j+1; k++){
                if (startnodes[k]==temp_start){
                    temp_start = (int)((temp_start + 1) % node_num);
                    k = 0;
                }
            }
            startnodes[j] = temp_start;
            edge_start[current_edge_num] = temp_start;
            outdegree_per_node[temp_start] += 1;
            printf("(%d,%d) ",temp_start,end);
            current_edge_num++;
        }
    }
    printf("\n");
}



void init_transition(double *transiton, int *outdegree_per_node, int *start_edge, int edge_num){
    // update trainsition
    for (int i = 0; i < edge_num; i++){
        transiton[i] = 1.0/outdegree_per_node[start_edge[i]];
        }
}


void pagerank(long double *nextrank_value, long double *currentrank_value, double *transition,int *start_edge,int *end_edge, int edge_num, int iterator, int thread_count, int NODES_NUM){
    int row;
    int col;
    int iter;
    int idx;
    for (iter = 0; iter < iterator; iter++){
        double tempsum = 0;
#pragma omp parallel for num_threads(thread_count)\
shared(nextrank_value,currentrank_value,transition) \
private(row,col,idx)
        for (idx = 0;idx < edge_num; idx++){
            row = start_edge[idx];
            col = end_edge[idx];
            nextrank_value[col] += currentrank_value[row] * transition[idx];
        }

        // avoid too big or small value
#pragma omp parallel for num_threads(thread_count)\
shared(nextrank_value,tempsum) \
private(row)
        for (row = 0; row < NODES_NUM; row++) {
            nextrank_value[row] = max(1.0 / NODES_NUM / NODES_NUM, nextrank_value[row]);
            nextrank_value[row] = min(1.0 - 1.0 / NODES_NUM / NODES_NUM, nextrank_value[row]);
            tempsum += nextrank_value[row];
        }

        // update and normlize currentrank_value
#pragma omp parallel for num_threads(thread_count)\
shared(nextrank_value,currentrank_value,tempsum) \
private(row)
        for (row = 0; row < NODES_NUM; row++) {
            nextrank_value[row] = nextrank_value[row]/tempsum;
            currentrank_value[row] = nextrank_value[row];
            nextrank_value[row] = 0;
        }

        // print
        //printf("iter %d: ",iter+1);
        //print_rank(currentrank_value);
    }
}


int main(int argc,char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: ./PageRank thread_number node_num\n");
        exit(1);
    }
    int thread_count;
    thread_count=strtol(argv[1],NULL,10);
    int NODES_NUM = atoi(argv[2]);

    // init
    int *inedge_number_per_node = new int[NODES_NUM];
    int *outedge_number_per_node = new int[NODES_NUM];

    int edge_number;
    // decide random edge numbers
    edge_number = rand_edge_num(inedge_number_per_node, NODES_NUM, MAX_EDGE_NUM_PER_NODE);

    printf("Total %d nodes!\n",NODES_NUM);
    printf("Total %d edges!\n",edge_number);

    // init variables
    const std::size_t N = NODES_NUM*10;

    int *edge_start = new int[N];
    int *edge_end = new int[N];

    double *transition = new double[N];

    long double *rank = new long double[NODES_NUM];
    long double *next_rank = new long double[NODES_NUM];


    for (int i = 0; i < NODES_NUM; i++){
        next_rank[i] = 0;
        rank[i] = 1.0/NODES_NUM;
        outedge_number_per_node[i] = 0;
    }
    // print
    printf("Init ", NODES_NUM);
    print_rank(rank, NODES_NUM);

    // random graph
    rand_graph(edge_start, edge_end, inedge_number_per_node, outedge_number_per_node, NODES_NUM);
    // init transition
    for (int i = 0; i < (int)N; i++){
        transition[i] = 0;
    }
    init_transition(transition,outedge_number_per_node,edge_start,edge_number);

    // pagerank
    pagerank(next_rank,rank,transition,edge_start,edge_end,edge_number,ITERATION,thread_count, NODES_NUM);
    // print
    printf("After %d iteration ",ITERATION);
    print_rank(rank, NODES_NUM);

    // free
    delete[] edge_start;
    delete[] edge_end;
    delete[] transition;
    delete[] rank;
    delete[] next_rank;
    delete[] inedge_number_per_node;
    delete[] outedge_number_per_node;

}
