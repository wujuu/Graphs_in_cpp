// Graphs_in_cpp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <climits>
using namespace std;
const int N = 10, E = 40, K = 50;


//Array
void Swap(int* &A, int i, int j) {
	int tmp = A[j];
	A[j] = A[i];
	A[i] = tmp;
}
int Get_min(int* A, bool* &B, int n) {
	int min = A[0], id_min = 0;
	for (int i = 1; i < n; i++) {
		if (!B[i] && A[i] < min) {
			min = A[i];
			id_min = i;
		}
	}
	B[id_min] = false;
	return id_min;
}

//Vertex
struct Vertex {
	int id;
	int w;
	Vertex* next;
};
Vertex* Init_vertex(int id, int w) {
	Vertex* x = new Vertex;
	x->id = id;
	x->w = w;
	x->next = NULL;
	return x;
}
void Push_vertex(Vertex* &first, Vertex* x) {
	x->next = first;
	first = x;
}

//Queue
struct Queue {
	Vertex* head, *tail;
};
Queue* Init_queue() {
	Queue* Q = new Queue;
	Q->head = Q->tail = NULL;
	return Q;
}
bool Empty_queue(Queue* &Q) {
	if (Q->head == NULL) return true;
	else return false;
}
void Push_queue(Queue* &Q, int id) {
	Vertex* x = Init_vertex(id,0);

	if (Empty_queue(Q)) Q->head  = x;
	else Q->tail->next = x;

	Q->tail = x;
}
int Pop_queue(Queue* &Q) {
	if (Empty_queue(Q)) return -1;
	else {
		Vertex* tmp = Q->head; int x = tmp->id;
		Q->head = Q->head->next;
		delete tmp;
		return x;
	}
}
void Print_queue(Queue* Q) {
	for (Vertex* tmp = Q->head; tmp != NULL; tmp = tmp->next) cout << tmp->id <<  " ";
	cout << endl;
}

//Heaps

struct heap {
	Vertex** A;
	int heapsize;
};
void Heap_swap(Vertex** &A, int i, int j) {
	Vertex* tmp = A[i];
	A[i] = A[j];
	A[j] = tmp;
}
heap* Init_heap(int n) {
	heap* H = new heap;
	H->A = new Vertex*[n];
	H->heapsize = 0;
	return H;
}
int Parent(int i) {
	return (i - 1) / 2;
}
int Left(int i) {
	return 2 * i + 1;
}
int Right(int i) {
	return 2 * i + 2;
}
void Heapify(heap* H, int i) {
	int id_min = i;
	if (Left(i) <= H->heapsize - 1 && H->A[Left(i)]->w < H->A[id_min]->w) id_min = Left(i);
	if (Right(i) <= H->heapsize - 1 && H->A[Right(i)]->w < H->A[id_min]->w) id_min = Right(i);
	if (id_min != i) {
		Heap_swap(H->A, i, id_min);
		Heapify(H, id_min);
	}
}
void Decr_key(heap* H, int i, int key) {
	H->A[i]->w = key;
	while (i > 0 && H->A[Parent(i)]->w < H->A[i]->w) {
		Heap_swap(H->A, Parent(i), i);
		i = Parent(i);
	}
}
void Push_heap(heap* H, Vertex* x) {
	H->heapsize++;
	H->A[H->heapsize - 1] = x;
	Decr_key(H, H->heapsize - 1, x->w);
}
int Pop_heap(heap* H) {
	if (H->heapsize - 1 > 0) {
		int max = H->A[0]->w;
		H->A[0] = H->A[H->heapsize - 1];
		H->heapsize--;
		Heapify(H, 0);
		return max;
	}
}

//Lasy zbiorow rozlaczynch
struct node {
	int key, rank;
	node *parent;
};
node* Init_node(int key) {
	node *x;
	x->key = key;
	x->parent = x;
	x->rank = 0;
	return x;
}
node* Find_parent(node* x) {
	if (x->parent == x) return x;
	else x->parent = Find_parent(x->parent);
}
void Union(node* x, node* y) {
	x = Find_parent(x);
	y = Find_parent(y);
	if (x->rank <= y->rank) {
		y->parent = x;
		if (x->rank == y->rank) x->parent++;
	}
	else if (x->rank < y->rank) x->parent = y;
}

//Edge
struct Edge {
	int x, y, w;
};
Edge* Init_edge(int x, int y, int w) {
	Edge* e = new Edge;
	e->x = x;
	e->y = y;
	e->w = w;
	return e;
}

//Quicksort
int Partition(Edge** &A, int start, int end) {
	int i = start - 1, j = end + 1, x = A[start]->w;
	while (true) {
		do i++; while (A[i]->w < x);
		do j--; while (A[j]->w > x);
		if (i < j) {
			int tmp = A[i]->w;
			A[i]->w= A[j]->w;
			A[j]->w = tmp;
		}
		else return j;
	}
}
void Quick_sort(Edge** &A, int end, int start = 0) {
	if (start < end) {
		int middle = Partition(A, end, start);
		Quick_sort(A, middle, start);
		Quick_sort(A, end, middle + 1);
	}
}



//Graphs
struct Graph {
	Vertex** V;
	Edge** E;
	int n;
	int e;
};
Graph* Init_graph(int n, int e, int k) {
	Graph* G = new Graph;
	G->n = n;
	G->e = e;

	G->V = new Vertex*[n]; 

	for (int i = 0; i < n; i++)
		G->V[i] = NULL;

	G->E = new Edge*[e];

	for (int i = 0; i < e; i++) {
		int x = rand() % n, y = rand() % n, w = rand() % k;
		Push_vertex(G->V[x], Init_vertex(y,w));
		G->E[i] = Init_edge(x, y, w);
	}

	return G;
}
int** Transform_graph(Graph* G) {
	int** D = new int*[G->n];
	for (int i = 0; i < G->n; i++)
		D[i] = new int[G->n];

	for (int i = 0; i < G->n; i++)
		for (int j = 0; j < G->n; j++)
			D[i][j] = INT_MAX;

	for (int i = 0; i < G->n; i++)
		for (Vertex* tmp = G->V[i]; tmp != NULL; tmp = tmp->next)
			D[i][tmp->id] = tmp->w;
		
	return D;
}
void Print_graph(Graph* G) {
	for (int i = 0; i < G->n; i++) {
		cout << i << ": ";
		for (Vertex* tmp = G->V[i]; tmp != NULL; tmp = tmp->next) cout << tmp->id << " ";
		cout << endl;
	}
}
void BFS(Graph* G, int* &D, int* &P, int s) {
	//Preparing
	bool* B = new bool[G->n];

	for (int i = 0; i < G->n; i++) {
		D[i] = INT_MAX;
		P[i] = -1;
		B[i] = false;
	}

	D[s] = 0;
	B[s] = true;

	Queue* Q = Init_queue();
	Push_queue(Q, s);

	//BFS
	while (!Empty_queue(Q)) {
		int u = Pop_queue(Q);

		for (Vertex* tmp = G->V[u]; tmp != NULL; tmp = tmp->next) {

			if (!B[tmp->id]) {
				D[tmp->id] = D[u] + 1;
				P[tmp->id] = u;
				B[tmp->id] = true;
				Push_queue(Q, tmp->id);
			}
		}
	}

	

	delete[] B;
}
void DFS_visit(Graph* G, int* &P, int* &I, int * &O, bool* &B, int u, int &time){
	time++;
	I[u] = time;

	B[u] = true;

	for (Vertex* tmp = G->V[u]; tmp != NULL; tmp = tmp->next) 
		if (!B[tmp->id]) {
			P[tmp->id] = u;
			DFS_visit(G, P, I, O, B, tmp->id, time);
		}

	time++;
	O[u] = time;
}
void DFS(Graph* G, int* &P, int* &I, int * &O) {
	//Preparing
	bool* B = new bool[G->n];

	for (int i = 0; i < G->n; i++) {
		P[i] = -1;
		B[i] = false;
	}

	int time = 0;

	//DFS
	for (int i = 0; i < G->n; i++)
		if (!B[i]) DFS_visit(G, P, I, O, B, i, time);

	delete[] B;
}
void Prim(Graph* G, int* &P, int s=0) {
	//Preparing
	bool* B = new bool[G->n];
	int* D = new int[G->n];
	heap* Q = Init_heap(G->n);

	for (int i = 0; i < G->n; i++) {
		D[i] = INT_MAX;
		P[i] = -1;
		B[i] = false;
		Push_heap(Q,Init_vertex(i,D[i]));
	}

	D[s] = 0;
	Decr_key(H,)


	//Prim
	for (int i = 0; i < G->n; i++) {
		int u = Get_min(D, B, G->n);

		for (Vertex* tmp = G->V[u]; tmp != NULL; tmp = tmp->next) 
			if (!B[tmp->id] && tmp->w < D[tmp->id]) {
				P[tmp->id] = u;
				D[tmp->id] = tmp->w;
			}

	}

	delete[] B;
	delete[] D;
}
void Kruskal(Graph* G, int* &P) {
	node** F = new node*[G->n];

	for (int i = 0; i < G->n; i++) {
		P[i] = -1;
		F[i] = Init_node(i);

	}

	Quick_sort(G->E, G->e);

	for (int i = 0; i < G->e; i++) {
		if (Find_parent(F[G->E[i]->x]) != Find_parent(F[G->E[i]->y])) {
			P[G->E[i]->y] = P[G->E[i]->x];
			Union(F[G->E[i]->x], F[G->E[i]->y]);
		}
	}

	delete[] F;
}
bool Relax(int* &D, int *&P,  int start, int end, int weight) {
	if (D[end] > D[start] + weight) {
		D[end] = D[start] + weight;
		P[end] = start;
		return true;
	}
	
	return false;
}
void Dijikstra(Graph* G, int* &D, int* &P, int s) {
	//Preparing
	bool* B = new bool[G->n];

	for (int i = 0; i < G->n; i++) {
		D[i] = INT_MAX;
		P[i] = -1;
		B[i] = false;
	}

	D[s] = 0;

	//Dijikstra
	for (int i = 0; i < N; i++) {

		int u = Get_min(D, B, G->n);

		for (Vertex* tmp = G->V[u]; tmp != NULL; tmp = tmp->next)
			if(!B[tmp->id]) 
				Relax(D, P, u, tmp->id, tmp->w);

	}

	delete[] B;
}
bool Bellmann_Ford(Graph* G, int* &D, int* &P, int s) {
	//Preparing
	for (int i = 0; i < G->n; i++) {
		D[i] = INT_MAX;
		P[i] = -1;
	}

	D[s] = 0;

	for (int i = 0; i < G->n-1; i++) {
		bool test = false;

		for (int j = 0; j < G->n; i++)
			for (Vertex* tmp = G->V[j]; tmp != NULL; tmp = tmp->next)
				test=Relax(D, P, j, tmp->id, tmp->w);

		if (test == false) return true;
	}

	for (int i = 0; i < G->n; i++)
		for (Vertex* tmp = G->V[i]; tmp != NULL; tmp = tmp->next)
			if (Relax(D, P, i, tmp->id, tmp->w)) return false;

	delete[] B;

	return true;
}
int** Floyd_Warshall(Graph* G) {
	int** D = Transform_graph(G);

	for (int k = 0; k < G->n; k++) 
		for (int u = 0; u < G->n; u++) 
			for (int v = 0; v < G->n; v++) 
				if (D[u][v] > D[u][k] + D[k][v]) 
					D[u][v] = D[u][k] + D[k][v];

	return D;
}



int main() {
	srand(time(NULL));
	Graph* G = Init_graph(N, E, K);
	Print_graph(G);

    return 0;
}

