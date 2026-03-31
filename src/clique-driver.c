/* Driver file for Clique Enumerator by BK
 *
 * Author: Cheng Chen: cchen67@vols.utk.edu; Yun Zhang: yzhang@cs.utk.edu
 * Last modified: 2022
 *
 * Changes: Add the upper bound of clique size to the BK
 *
 * Copyright 2005-2006
 * Department of Computer Science, University of Tennessee, Knoxville
 *
 */

#include <string.h>
#include "graph.h"
#include "bk.h"
#include "ReadKG.h"

#define MAX_PATH_LENGTH 4096 // On Linux

/* Global variables */
extern int LB,UB;  // lower bound and upper bound of clique size
extern int VERSION;
extern int PRINT;
extern int PART;
extern int Spart;
FILE *fp;
char *outfn = NULL;
char *infn;
char * confile;
int lb;


void print_options(void)
{
  fprintf(stderr, "\n Options: \n");
  fprintf(stderr, "  -p             print out clique list \n");
  fprintf(stderr, "                 note: no print out if not specify -p\n");
  fprintf(stderr, "                       print to stdout if not specify -o\n");
  fprintf(stderr, "  -o <filename>  filename to store cliques if choose to print out\n");
  fprintf(stderr, "  -l <value>     least number of vertices in a clique <default = 3>\n");
  fprintf(stderr, "  -u <value>     most number of vertices in a clique <default = graph size>\n");
  fprintf(stderr, "  -v [1~8]     algorithm version <default = 2>\n");
  fprintf(stderr, "                 1 - bron kerbosch version 1 (numerical order)\n");
  fprintf(stderr, "                 2 - bron kerbosch version 2 (improved version)\n");
  fprintf(stderr, "                 3 - modified bron kerbosch to find maximum clique\n");
  fprintf(stderr, "                 4 - maximal partite cliques enumration (BK with piivot)\n");
  fprintf(stderr, "                 5 - maximal partite cliques enumration (MMCE + BK)\n");
  fprintf(stderr, "                 6 - maximal partite cliques enumration (MPCE)\n");
  fprintf(stderr, "                 7 - maximal partite cliques enumration (MMCE + BK with pivot)\n");
  fprintf(stderr, "                 8 - maximal partite cliques enumration (MMCE variant + BK with pivot)\n");
  fprintf(stderr, "\n");
}

void argument_parse(int argc, char **argv)
{
  int i;

  if (argc < 2) {
    fprintf(stderr, "Usage: %s Graph <options>\n", argv[0]);
    print_options();
	exit(1);
  }

  LB = 3; UB = -1;
  VERSION = 2;
  PRINT = 0;
  PART = 0;



  for (i = 2; i < argc; i++) {
	if (!strcmp(argv[i], "-o")) {
	  outfn = strdup(argv[++i]);
	}
	if (!strcmp(argv[i], "-l")) {
	  LB = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-u")) {
	  UB = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-v")) {
	  VERSION = atoi(argv[++i]);
	}
	if (!strcmp(argv[i], "-p")) {
	  PRINT = 1;
	}
  if (!strcmp(argv[i], "-klb")) {
    PART = 1;
	  lb = atoi(argv[++i]);
	}
  }

  infn = strdup(argv[1]);
  if ((fp = fopen(argv[1], "r")) == NULL) {
    fprintf(stderr, "Can't open input file %s\n", argv[1]);
    exit(-1);
  }

  return;
}

int cmp_gp(const void *gp1, const void *gp2) {
    const GP *a = (const GP *)gp1;
    const GP *b = (const GP *)gp2;
    if (a->size > b->size) {return 1;}
    else if (a->size < b->size) {return -1;}
    else {return 0;}
}

void order_vertex(Graph *G, int psizes[], vid_t * vertices, GP * gp)
{
  int i, j, k;
  int c;
  for(i = 0; i < G->Pnum; i++)
  {
    gp[i].vertices = (int *) malloc(psizes[i] * sizeof(int));
    gp[i].index = 0;
  }
  for(i = 0; i < G->_num_vertices; i++)
  {
    c = G->_category[i];
    (gp[c].vertices)[gp[c].index] = i;
    gp[c].index = gp[c].index + 1;
  }
  if((VERSION == 5) | (VERSION == 7))
  {return;}
  qsort(gp, G->Pnum, sizeof(GP), cmp_gp);
  k = 0;
  for(i = 0; i < G->Pnum; i++)
  {
    for(j = 0; j < gp[i].index; j++)
    {
      vertices[k] = gp[i].vertices[j];
      k++;
    }
  }
}

void maximal_clique(char * confile, Graph *G)
{
  FILE *fp1=stdout, *fp2;
  char fname[MAX_PATH_LENGTH];
  double utime;
  unsigned int n = num_vertices(G);
  vid_t clique[n];
  vid_t vertices[n];
  u64 nclique[n+1];
  int i, pid;

  if (outfn != NULL) {
    sprintf(fname, "%s.clique", outfn);
    fp1 = fopen(fname, "w");
  }

  if (strlen(infn) + strlen(".profile") >= sizeof(fname)) {
      // Truncate infn to fit the buffer
      strncpy(fname, infn, sizeof(fname) - strlen(".profile") - 1);
      strcat(fname, ".profile");
  } else {
      snprintf(fname, sizeof(fname), "%s.profile", infn);
  }


  fp2 = fopen(fname, "w");

  utime = get_cur_time();
  memset(nclique, 0, (n+1)*sizeof(u64));
  memset(clique, -1, n*sizeof(vid_t));
  for (i = 0; i < n; i++) vertices[i] = i;
  if (VERSION == 1)
	{clique_find_v1(fp1, nclique, G, clique, vertices, 0, 0, n);}
  else if (VERSION == 2)
	{clique_find_v2(fp1, nclique, G, clique, vertices, 0, 0, n);}
  else
  {

    int psizes[G->Pnum];//each partite has how much nodes
    int csizes[G->Pnum];
    memset(psizes, 0, G->Pnum*sizeof(int));
    memset(csizes, 0, G->Pnum*sizeof(int));
    for (i = 0; i < n; i++)
    {
      pid = G->_category[i];
      psizes[pid]++;
    }
    G->lbs = (int *) malloc(G->Pnum * sizeof(int));
    for(i = 0; i < G->Pnum; i++)
    {
      G->lbs[i] = lb;}

    GP * gp = calloc(G->Pnum, sizeof(GP));
    G->gps = gp;
    for(i = 0; i < G->Pnum; i++)
    {
      gp[i].color = i;
      gp[i].size = psizes[i];
      gp[i].vertices = NULL;
      gp[i].index = 0;
    }
    if(VERSION == 5)
    {
      order_vertex(G, psizes, vertices, gp);
      clique_find_v5(fp1, nclique, G, clique, vertices, 0, 0, n, csizes);
    }
    else if (VERSION == 4)
    {
      clique_find_v4(fp1, nclique, G, clique, vertices, 0, 0, n, csizes, psizes);
    }
    else if(VERSION == 6)
    {
      clique_find_v6(fp1, nclique, G, clique, vertices, 0, 0, n, csizes, psizes);
    }
    else if(VERSION == 7)
    {
      order_vertex(G, psizes, vertices, gp);
      clique_find_v7(fp1, nclique, G, clique, vertices, 0, 0, n, csizes);
    }
    else if (VERSION == 8)
    {
      clique_find_v8(fp1, nclique, G, clique, vertices, 0, 0, n, csizes, psizes);
    }
    else
    {
      fprintf(stderr, "Wrong version number! The version number should be between 1 ~ 8\n");
      exit(1);
    }
  }
  utime = get_cur_time() - utime;


  clique_profile_out(fp2, nclique, G);
  fprintf(fp2, "Time (seconds)  : %.6f\n", utime);

  if (outfn != NULL) {
    free(outfn);
    fclose(fp1);
  }
  fclose(fp2);

}

void maximum_clique(Graph *G)
{
  FILE *fp=stdout;
  unsigned int n = num_vertices(G);
  vid_t clique[n];
  vid_t vertices[n];
  vid_t maxclique[n];
  int maxclique_size=0;
  vid_t i;

  if (outfn != NULL) fp = fopen(outfn, "w");

  memset(clique, -1, n*sizeof(vid_t));
  for (i = 0; i < n; i++) vertices[i] = i;
  maxclique_find(maxclique, &maxclique_size, G, clique, vertices, 0, 0, n);

  if (PRINT) clique_out(fp, G, maxclique, maxclique_size);
  fprintf(fp, "Maximum Clique Size : %d\n", maxclique_size);

  if (outfn != NULL) {
    free(outfn);
    fclose(fp);
  }
  return;
}


int main(int argc, char  **argv)
{
  Graph *G;

  argument_parse(argc, argv);

  G = graph_edgelist_in(fp, PART);
  fclose(fp);
  if (UB == -1) UB = num_vertices(G);

  if (VERSION == 3) maximum_clique(G);
  else maximal_clique(confile, G);

  graph_free(G);

  exit(0);
}
