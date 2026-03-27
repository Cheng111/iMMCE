/* Enumerate clique by Bron-Kerbosch Algorithms 
 * Author: Cheng Chen Yun Zhang
 * Date: 2022
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"
#include "graph.h"
#include "bk.h"
//#include "PartiteSize.h"


/* Global Variables */
int LB, UB, PART;
int VERSION;
int PRINT;

vid_t select_tomita_partite_pivot(Graph *G, vid_t *old, int ne, int ce)
{
  vid_t pivot;
  vid_t candidate;
  int best_count;
  int count;
  int i, j;

  pivot = old[ne];
  best_count = -1;

  for (i = ne; i < ce; i++) {
	candidate = old[i];
	count = 0;
	for (j = ne; j < ce; j++) {
	  if (edge_exists(G, candidate, old[j]) || same_category(G, candidate, old[j])) {
		count++;
	  }
	}
	if (count > best_count) {
	  best_count = count;
	  pivot = candidate;
	}
  }

  return pivot;
}

/* ------------------------------------------------------------------- *
 * Function: clique_out()    *                                     *
 * Output the clique or k-partite clique that meets the requirements.  *
 * ------------------------------------------------------------------- */
void clique_out(FILE *fp, Graph *G, vid_t *clique, int len)
{
  int i;

  for (i = 0; i < len-1; i++)
	{
		fprintf(fp, "%s\t", G->_label[clique[i]]);
	}
  fprintf(fp, "%s\n", G->_label[clique[i]]);
  return;
}


/* ------------------------------------------------------------- *
 * Function: clique_profile_out()                                *
 Report the number of cliques or k-partite cliques categorized by their length.
 * ------------------------------------------------------------- */
void clique_profile_out(FILE *fp, u64 *nclique, Graph *G)
{
  unsigned int n = G->_num_vertices;
  u64 sum=0;
  int max=0;
  int i;
  fprintf(fp, "Size\tNumber\n");
  for (i = LB; i <= n; i++) {
	if (nclique[i]) {
	  fprintf(fp, "%d\t%lu\n", i, nclique[i]);
	  sum += nclique[i];
	  max = i;
	}
  }
  fprintf(fp, "\n");
  fprintf(fp, "No. of vertices : %d\n", n);
  fprintf(fp, "No. of edges    : %d\n", G->_num_edges);
  fprintf(fp, "No. of cliques  : %lu\n", sum);
  fprintf(fp, "Max clique size : %d\n", max);
}

/* ------------------------------------------------------------- *
 * Function: clique_find_v1()                                    *
 *   Bron-Kerbosch version 1 : numerical order                   *
 *   Recursive function to find cliques                          *
 * ------------------------------------------------------------- */
void clique_find_v1(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t u;
  int i, j;

  while (ne < ce) 
  {

#ifdef DEBUG
  for (i = 0; i < ne; i++) printf(" %s", G->_label[old[i]]);
  printf("\t|");
  for (i = 0; i < lc; i++) printf(" %s", G->_label[clique[i]]);
  printf("\t|");
  for (i = ne; i < ce; i++) printf(" %s", G->_label[old[i]]);
  printf("\n");
#endif

	u = old[ne];

	/* Set new cand and not */
	memset(new, -1, ce*sizeof(vid_t));
    new_ne = 0;
	// X intersection N(u)
	for (j = 0; j < ne; j++)
	  if (edge_exists(G, u, old[j])) new[new_ne++] = old[j];
	new_ce = new_ne;
	// p intersection N(u)
	for (j = ne+1; j < ce; j++)
	  if (edge_exists(G, u, old[j])) new[new_ce++] = old[j];
	
	/* Output clique or extend */
	clique[lc] = u;
	if (new_ce == 0 && lc+1 >= LB) { //P U X = empty
	  nclique[lc+1]++;
	  if (PRINT) clique_out(fp, G, clique, lc+1);
	}
	else { 
	  if (new_ne < new_ce)
	    clique_find_v1(fp, nclique, G, clique, new, lc+1, new_ne, new_ce);
	}

	
	/* Move u to not */
	ne++;

	/* Bound condition: Stop if any not is a neighbor of all candidates */ 
    for (i = 0; i < ne; i++) {
	  for (j = ne; j < ce; j++) {
	    if (!edge_exists(G, old[i], old[j])) break;
	  }
	  if (j == ce) return;
	}

  }

}


/* ------------------------------------------------------------- *
 * Function: clique_find_v2()                                    *
 *   Bron-Kerbosch version 2                                     *
 *   Recursive function to find cliques                          *
 * ------------------------------------------------------------- */
void clique_find_v2(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t fixp=0, p, u;
  int s=0, pos=0, nod, minnod, count;
  int i, j, k;

  if (lc + (ce - ne) < LB) return;

#ifdef DEBUG
  for (i = 0; i < ne; i++) printf(" %s", G->_label[old[i]]);
  printf("\t|");
  for (i = 0; i < lc; i++) printf(" %s", G->_label[clique[i]]);
  printf("\t|");
  for (i = ne; i < ce; i++) printf(" %s", G->_label[old[i]]);
  printf("\n");
#endif

  /* Choose a vertex, fixp, in old (both not and cand) that
	 has lowest number of non-adjacent vertices in old cand */
  minnod = ce + 1;
  nod = 0;
  for (i = 0; i < ce; i++) {
	count = 0;
	p = old[i];
	for (j = ne; j < ce; j++) {
	  if (!edge_exists(G, p, old[j])) {
		count++;
		pos = j;
	  }
	}
	if (count < minnod) {
	  fixp = p;
	  minnod = count;
	  if (i < ne) { s = pos; }    // if p in not
	  else { s = i; nod = 1; }    // if p in cand
	}
  }
  
  /* Recursively extend clique */
  for (k = minnod+nod; k > 0; k--) {

	/* Swap this candidate to be the next one */
	p = old[s];
	old[s] = old[ne];
	old[ne] = p;

	u = old[ne];

	/* Set new cand and not */
	memset(new, -1, ce*sizeof(vid_t));
    new_ne = 0;
	for (j = 0; j < ne; j++)
	  if (edge_exists(G, u, old[j])) new[new_ne++] = old[j];
	new_ce = new_ne;
	for (j = ne+1; j < ce; j++)
	  if (edge_exists(G, u, old[j])) new[new_ce++] = old[j];
	
	/* Output clique or extend */
	clique[lc] = u;
	if (lc+1 <= UB) {
	  if (new_ce == 0 && lc+1 >= LB) {
	    nclique[lc+1]++;
	    if (PRINT) clique_out(fp, G, clique, lc+1);
	  }
	  else if (new_ne < new_ce) {
	    clique_find_v2(fp, nclique, G, clique, new, lc+1, new_ne, new_ce);
	  }
	}
	
	/* Move u to not */
	ne++;

	/* Bound condition: Stop if fixp is a neighbor of all candidates */ 
    if (k > 1) {
	  for (s = ne; s < ce; s++) {
	    if (!edge_exists(G, fixp, old[s])) break;
	  }
	  if (s == ce) return;
	}

  }

  return;
}


/* ------------------------------------------------------------- *
 * Function: maxclique_find()                                    *
 *   Bron-Kerbosch version 2                                     *
 *   Recursive function to find one of the maximum cliques       *
 * ------------------------------------------------------------- */
void maxclique_find(vid_t *maxclique, int *maxclique_size, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t fixp=0, p, u;
  int s=0, pos=0, nod, minnod, count;
  int i, j, k;

#ifdef DEBUG
  for (i = 0; i < ne; i++) printf(" %d", old[i]);
  printf("\t|");
  for (i = 0; i < lc; i++) printf(" %d", clique[i]);
  printf("\t|");
  for (i = ne; i < ce; i++) printf(" %d", old[i]);
  printf("\n");
#endif

  /* Choose a vertex, fixp, in old (both not and cand) that
	 has lowest number of non-adjacent vertices in old cand */
  minnod = ce + 1;
  nod = 0;
  for (i = 0; i < ce; i++) {
	count = 0;
	p = old[i];
	for (j = ne; j < ce; j++) {
	  if (!edge_exists(G, p, old[j])) {
		count++;
		pos = j;
	  }
	}
	if (count < minnod) {
	  fixp = p;
	  minnod = count;
	  if (i < ne) { s = pos; }    // if p in not
	  else { s = i; nod = 1; }    // if p in cand
	}
  }
  
  /* Recursively extend clique */
  for (k = minnod+nod; k > 0; k--) {

	/* Swap this candidate to be the next one */
	p = old[s];
	old[s] = old[ne];
	old[ne] = p;

	u = old[ne];

	/* Set new cand and not */
	memset(new, -1, ce*sizeof(vid_t));
    new_ne = 0;
	for (j = 0; j < ne; j++)
	  if (edge_exists(G, u, old[j])) new[new_ne++] = old[j];
	new_ce = new_ne;
	for (j = ne+1; j < ce; j++)
	  if (edge_exists(G, u, old[j])) new[new_ce++] = old[j];
	
	/* Record clique or extend */
	clique[lc] = u;
	if (new_ce == 0 && lc+1 >= *maxclique_size) {
	  *maxclique_size = lc + 1;
	  memcpy(maxclique, clique, (*maxclique_size)*sizeof(vid_t));
	  if (PRINT) {
		printf("max size %d\n", lc+1);
		clique_out(stdout, G, clique, lc+1);
	  }
	}
	else if (new_ce-new_ne+lc+1 >= *maxclique_size) {
	  maxclique_find(maxclique, maxclique_size, G, clique, new, lc+1, new_ne, new_ce);
	}
	
	/* Move u to not */
	ne++;

	/* Bound condition: Stop if fixp is a neighbor of all candidates */ 
    if (k > 1) {
	  for (s = ne; s < ce; s++) {
	    if (!edge_exists(G, fixp, old[s])) break;
	  }
	  if (s == ce) return;
	}

  }

  return;
}


/* ------------------------------------------------------------- *
 * Function: clique_enumerate()                                  
 *   For mpiclique version 1, 2
 * ------------------------------------------------------------- */

void clique_enumerate(FILE *fp, u64 *nclique, Graph *G, vid_t *cand, int lcand)
{
  unsigned int n = num_vertices(G);
  vid_t clique[n];
  vid_t vertices[n];
  int ne, ce;
  vid_t u, i, j;
#ifdef DEBUG 
  double utime;
#endif
  
  memset(clique, -1, n*sizeof(vid_t));

  for (i = 0; i < lcand; i++) {

#ifdef DEBUG
    utime = get_cur_time();
#endif
	
	u = cand[i];

	/* Prepare cand and not array */
	ne = 0;
	for (j = 0; j < u; j++)
	  if (edge_exists(G, j, u)) vertices[ne++] = j;
	ce = ne;
	for (j = u+1; j < n; j++)
	  if (edge_exists(G, u, j)) vertices[ce++] = j;

	/* Recursively find cliques containing u if enough candidates */
	if (ce - ne >= LB - 1) {
	  clique[0] = u;
	  if (VERSION == 1)
        clique_find_v1(fp, nclique, G, clique, vertices, 1, ne, ce);
	  else if (VERSION == 2)
        clique_find_v2(fp, nclique, G, clique, vertices, 1, ne, ce);
    }

#ifdef DEBUG
	utime = get_cur_time() - utime;
	printf("task %4d : %d subtasks, %f seconds\n", u, ce-ne, utime);
#endif
  }
  
  return;
}

/* ------------------------------------------------------------- *
 * Function: clique_find_v4()                                    *
 *   Use the Bron–Kerbosch algorithm with pivoting to enumerate  *
 *   maximal k-partite cliques.                                  *
 *   Terminate when the candidate set from one partite is empty  *
 *   and no vertex from this partite is present in the current   *
 *   clique.                                                     *
 * ------------------------------------------------------------- */
 /* vid_t *clique: current candidate partite clique C.
 	vid_t *old: X (vertices cannot extend C) cascade M (vertices may extend the C)
	ne: size of X
	ce: size of X + M
	nclique: used for hist graph for clique length
	lc: clique length - 1 */
void clique_find_v4(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes, int * psizes)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t u, pivot;
  int i, j, pnei;//neibor of pivot vertex should be in P
  int new_psizes[G->Pnum];
  int upid, jpid;
  int parclique = 0;
  //Terminate when the candidate set of a partite is empty and no vertex from that partite is included in the candidate k-partite clique.
  for(i = 0; i < G->Pnum; i++)
  {
	  if(csizes[i] == 0 && psizes[i] == 0)
	  {
		return;}
  }
  pivot = select_tomita_partite_pivot(G, old, ne, ce);
  pnei = ce -1;
  while (ne < ce) 
  {
	u = old[ne];
	upid = G->_category[u];
	if (u != pivot)
	{
        //In this recursion, skip processing vertices that are either adjacent to the pivot or belong to the same partite set as the pivot vertex.
		if (edge_exists(G, u, pivot) || same_category(G, u, pivot))
		{
			old[ne] = old[pnei];
			old[pnei] = u;
			pnei--;

			continue;
		}
	}
	if(pnei < ne)
	{break;}
	/* Set new cand and not */
	memset(new, -1, ce*sizeof(vid_t));
	memset(new_psizes, 0, G->Pnum * sizeof(int));
    new_ne = 0;
	// X intersection (N(u) U P(u))
	for (j = 0; j < ne; j++)
	{
	  if (edge_exists(G, u, old[j]) || same_category(G, u, old[j])) 
	  {
		  jpid = G->_category[old[j]];
		  new[new_ne++] = old[j];
	  }
	}
	new_ce = new_ne;
	// M intersection (N(u) U P(u))
	for (j = ne+1; j < ce; j++)
	{
	  if (edge_exists(G, u, old[j]) || same_category(G, u, old[j]))
	  {
		  jpid = G->_category[old[j]];
		  new_psizes[jpid]++;
		  new[new_ce++] = old[j];//add N(u) nad P(u)(same partite set) to X
	  }
	}
	
	/* Output clique or extend */
	upid = G->_category[u];
	clique[lc] = u;
	csizes[upid]++;
	if(new_ce == 0) 
	{
	  parclique = 0;
	  for(i = 0; i < G->Pnum; i++)
	  {
		  if(csizes[i] < G->lbs[i])
		  {
			  parclique = 1;
			  break;
		  }
	  }
	  if(parclique == 0)
	  {
		  nclique[lc+1]++;
	  	if (PRINT) clique_out(fp, G, clique, lc+1);
	  }
	}
	else { 
	  if (new_ne < new_ce)
	    clique_find_v4(fp, nclique, G, clique, new, lc+1, new_ne, new_ce, csizes, new_psizes);
	}
	csizes[upid]--;

	
	/* Move u to X */
	ne++;

	/* Bound condition: Stop if any not is a neighbor of all candidates */ 
    for (i = 0; i < ne; i++) {
	  for (j = ne; j < ce; j++) {
	    if (!edge_exists(G, old[i], old[j])) break;
	  }
	  if (j == ce) return;
	}


  }

}

void clique_find_v5(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes)
{
	int i, j, k;
	GP * gp;
	vid_t u, v;
	for(k = 0; k < G->Pnum; k++ )
	{

		gp = G->gps + k;
		for(i = 0; i < gp->size; i++)
		{
			for(j = i + 1; j < gp->size; j++)
			{
				u = old[gp->vertices[i]];
				v = old[gp->vertices[j]];
				add_edge(G, u, v);
			}
		}
	}
	clique_find_v5_sub(fp, nclique, G, clique, old, lc, ne, ce, csizes);
}

void clique_find_v5_sub(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t fixp=0, p, u;
  int s=0, pos=0, nod, minnod, count;
  int i, j, k;
  int upid;
  int parclique;

  if (lc + (ce - ne) < LB) return;

#ifdef DEBUG
  for (i = 0; i < ne; i++) printf(" %s", G->_label[old[i]]);
  printf("\t|");
  for (i = 0; i < lc; i++) printf(" %s", G->_label[clique[i]]);
  printf("\t|");
  for (i = ne; i < ce; i++) printf(" %s", G->_label[old[i]]);
  printf("\n");
#endif

  /* Choose a vertex, fixp, in old (both not and cand) that
	 has lowest number of non-adjacent vertices in old cand */
  minnod = ce + 1;
  nod = 0;
  for (i = 0; i < ce; i++) {
	count = 0;
	p = old[i];
	for (j = ne; j < ce; j++) {
	  if (!edge_exists(G, p, old[j])) {
		count++;
		pos = j;
	  }
	}
	if (count < minnod) {
	  fixp = p;
	  minnod = count;
	  if (i < ne) { s = pos; }    // if p in not
	  else { s = i; nod = 1; }    // if p in cand
	}
  }
  
  /* Recursively extend clique */
  for (k = minnod+nod; k > 0; k--) {

	/* Swap this candidate to be the next one */
	p = old[s];
	old[s] = old[ne];
	old[ne] = p;

	u = old[ne];

	/* Set new cand and not */
	memset(new, -1, ce*sizeof(vid_t));
    new_ne = 0;
	for (j = 0; j < ne; j++)
	  if (edge_exists(G, u, old[j])) new[new_ne++] = old[j];
	new_ce = new_ne;
	for (j = ne+1; j < ce; j++)
	  if (edge_exists(G, u, old[j])) new[new_ce++] = old[j];
	
	/* Output clique or extend */
	clique[lc] = u;
	upid = G->_category[u];
	csizes[upid]++;
	if (lc+1 <= UB) {
	  if (new_ce == 0 && lc+1 >= LB) 
	  {
			parclique = 0;
		for(i = 0; i < G->Pnum; i++)
		{
			if(csizes[i] < G->lbs[i])
			{
				parclique = 1;
				break;
			}
		}
		if(PRINT && (parclique == 0))
		{
			nclique[lc+1]++;
			clique_out(fp, G, clique, lc+1);
		}
	  }
	  else if (new_ne < new_ce) 
	  {
	    clique_find_v5_sub(fp, nclique, G, clique, new, lc+1, new_ne, new_ce, csizes);
	  }
	}
	csizes[upid]--;
	
	/* Move u to not */
	ne++;

	/* Bound condition: Stop if fixp is a neighbor of all candidates */ 
    if (k > 1) {
	  for (s = ne; s < ce; s++) {
	    if (!edge_exists(G, fixp, old[s])) break;
	  }
	  if (s == ce) return;
	}

  }

  return;
}

/* ------------------------------------------------------------- 
 * Function: clique_find_v6()                                    
 *   MPCE to find k-partite cliques                
 *  1. Terminate when the candidate set from one partite is empty
 *   and no vertex from this partite is present in the current   
 *   clique.
 *  2. The first k vertices in the candidate clique set are from 
 * 	 the k different partite sets.                                                     
 * ------------------------------------------------------------- */
 /* vid_t *clique(R): current clique
 	vid_t *old(P): X cascade M
	ne: size of X
	ce: size of X + M
	nclique: used for hist graph for clique length
	lc: clique length - 1 */

void clique_find_v6(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes, int * psizes)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t u, pivot;
  int i, j, pnei;//neibor of pivot vertex should be in P
  int ne_copy = ne;
  int flag_color;
  int new_psizes[G->Pnum];
  int upid, jpid;
  int parclique = 0;
  //Terminate when the candidate set of a partite is empty and no vertex from that partite is included in the candidate k-partite clique.
  for(i = 0; i < G->Pnum; i++)
  {
	  if(csizes[i] == 0 && psizes[i] == 0)
	  {
		return;}
  }
  flag_color = ce - 1;
  if(lc < G->Pnum)
  {
	ne_copy = ne;
	while(ne_copy < ce)
	{
		u = old[ne_copy];
		//The first k vertices in the candidate clique set are from the k different partite sets.
		if (csizes[G->_category[u]] > 0)
		{

			old[ne_copy] = old[flag_color];
			old[flag_color] = u;
			flag_color--;
			if(flag_color < ne_copy)
			{
				break;
				}
		}
		else
		{
			ne_copy++;
			if(flag_color < ne_copy)
			{break;}
		}
	}
	pnei = flag_color;
  }
  else
  {
	pnei = ce -1;
	pivot = select_tomita_partite_pivot(G, old, ne, ce);
	ne_copy = ne;
	while(ne_copy <= pnei)
	{
		u = old[ne_copy];
		if (u != pivot)
		{
			//In this recursion, skip processing vertices that are either adjacent to the pivot or belong to the same partite set as the pivot vertex.
			if (edge_exists(G, u, pivot) || same_category(G, u, pivot))
			{
				old[ne_copy] = old[pnei];
				old[pnei] = u;
				pnei--;
				if(pnei < ne_copy)
				{break;}

				continue;
			}
		}
		ne_copy++;
		if(pnei < ne_copy)
		{break;}
	}
  }
  while (ne <= pnei) 
  {
	u = old[ne];
	upid = G->_category[u];
	memset(new, -1, ce*sizeof(vid_t));
	memset(new_psizes, 0, G->Pnum * sizeof(int));
    new_ne = 0;
	// X intersection (N(u) U P(u))
	for (j = 0; j < ne; j++)
	{
	  if (edge_exists(G, u, old[j]) || same_category(G, u, old[j])) 
	  {
		  jpid = G->_category[old[j]];
		  new[new_ne++] = old[j];
	  }
	}
	new_ce = new_ne;
	// M intersection (N(u) U P(u))
	for (j = ne+1; j < ce; j++)
	{
	  if (edge_exists(G, u, old[j]) || same_category(G, u, old[j]))
	  {
		  jpid = G->_category[old[j]];
		  new_psizes[jpid]++;
		  new[new_ce++] = old[j];//add N(u) and P(u)(same partite set) to X
	  }
	}
	
	/* Output clique or extend */
	upid = G->_category[u];
	clique[lc] = u;
	csizes[upid]++;
	if(new_ce == 0) 
	{
	  parclique = 0;
	  for(i = 0; i < G->Pnum; i++)
	  {
		  if(csizes[i] < G->lbs[i])
		  {
			  parclique = 1;
			  break;
		  }
	  }
	  if(parclique == 0)
	  {
		  nclique[lc+1]++;
	  	if (PRINT) clique_out(fp, G, clique, lc+1);
	  }
	}
	else { 
	  if (new_ne < new_ce)
	    clique_find_v6(fp, nclique, G, clique, new, lc+1, new_ne, new_ce, csizes, new_psizes);
	}
	csizes[upid]--;

	
	/* Move u to not */
	ne++;

	/* Bound condition: Stop if any not is a neighbor of all candidates */ 
    for (i = 0; i < ne; i++) {
	  for (j = ne; j < ce; j++) {
		if (!(edge_exists(G, old[i], old[j]) || same_category(G, old[i], old[j]))) break;
	  }
	  if (j == ce) return;
	}


  }

}
