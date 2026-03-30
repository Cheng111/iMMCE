#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utility.h"
#include "graph.h"
#include "bk.h"
#include "ReadKG.h"

/* Global Variables */
extern int VERSION;
extern int PRINT;

static vid_t select_tomita_pivot(Graph *G, vid_t *old, int ne, int ce)
{
  vid_t pivot;
  vid_t candidate;
  int best_count;
  int count;
  int i, j;

  pivot = old[ne];
  best_count = -1;

  for (i = 0; i < ce; i++) {
	candidate = old[i];
	count = 0;
	for (j = ne; j < ce; j++) {
	  if (edge_exists(G, candidate, old[j])) count++;
	}
	if (count > best_count) {
	  best_count = count;
	  pivot = candidate;
	}
  }

  return pivot;
}

void clique_find_v7(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes)
{
	int i, j, k;
	GP * gp;
	vid_t u, v;
	FILE *fout = fopen("filename.txt", "w");
	fclose(fout);
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
	clique_find_v7_sub(fp, nclique, G, clique, old, lc, ne, ce, csizes);
}

void clique_find_v7_sub(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t u, pivot;
  int i, j, pnei;//neibor of pivot vertex should be in P
  int upid;//, jpid;
  int parclique = 0;

  pivot = select_tomita_pivot(G, old, ne, ce);
  pnei = ce -1;
  while (ne < ce) 
  {
	u = old[ne];
	upid = G->_category[u];
	if (u != pivot)
	{
		if(edge_exists(G, u, pivot))
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
    new_ne = 0;
	for (j = 0; j < ne; j++)
	{
	  if (edge_exists(G, u, old[j])) 
	  {
		  new[new_ne++] = old[j];
	  }
	}
	new_ce = new_ne;
	for (j = ne+1; j < ce; j++)
	{
	  if (edge_exists(G, u, old[j]))
	  {
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
	  {
	    clique_find_v7_sub(fp, nclique, G, clique, new, lc+1, new_ne, new_ce, csizes);
	  }
	}
	csizes[upid]--;

	
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
 * Function: clique_find_v8()                                    *
 *   Bron-Kerbosch version 8 : BK with pivot on kpartite graph.  *
 * ------------------------------------------------------------- */
 /* vid_t *clique(R): current clique
 	vid_t *old(P): X cascade P
	ne: size of X
	ce: size of X + P
	nclique: used for hist graph for clique length
	lc: clique length - 1 */
void clique_find_v8(FILE *fp, u64 *nclique, Graph *G, \
		vid_t *clique, vid_t *old, int lc, int ne, int ce, int * csizes, int * psizes)
{
  vid_t new[ce];
  int new_ne, new_ce;
  vid_t u, pivot;
  int i, j, pnei;//neibor of pivot vertex should be in P
  int new_psizes[G->Pnum];
  int upid, jpid;
  int parclique = 0;
  pivot = select_tomita_partite_pivot(G, old, ne, ce);
  pnei = ce -1;
  while (ne < ce) 
  {
	u = old[ne];
	upid = G->_category[u];
	if (u != pivot)
	{
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
	for (j = 0; j < ne; j++)
	{
	  if (edge_exists(G, u, old[j]) || same_category(G, u, old[j])) 
	  {
		  jpid = G->_category[old[j]];
		  new[new_ne++] = old[j];
	  }
	}
	new_ce = new_ne;
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
	    clique_find_v8(fp, nclique, G, clique, new, lc+1, new_ne, new_ce, csizes, new_psizes);
	}
	csizes[upid]--;

	
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
