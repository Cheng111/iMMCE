#define _GNU_SOURCE

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <search.h>
#include "graph.h"

Graph * ReadKG(FILE *fp)
{
    unsigned int N = 0, E = 0, C = 0;
    char * line = NULL;
    char * pch;
    size_t len;
    ssize_t read;
    char * word1;
    char * word2;
    char * color;
    char * delim = " \t";
    ENTRY item;
    ENTRY *found_item;
    ENTRY citem;
    ENTRY *found_citem;
    struct hsearch_data *htab = NULL;
    struct hsearch_data *htabc = NULL;
    int u, v = -1, cid = -1;
    int k=0, edges=0, cid_use = 0;
    Graph * G = NULL;
    int *id;
    int *cids = NULL;
    int *ids = NULL;
    
    while ((read = getline(&line, &len, fp)) != -1) {
        while(line[strlen(line) - 1] == '\n' || line[strlen(line) - 1] == '\x0d'){
            line[strlen(line) - 1] = 0;
            if(strlen(line) == 0)
            {break;}
        }
        if(strlen(line) == 0)
        {continue;}
        pch = NULL;
        pch = strtok(line, delim);
        if(*pch == 'c')
        {continue;}
        else if(*pch == 'p')
        {
            char * edges = "edges";
            pch = strtok (NULL, delim);
            if(strcmp(pch, edges) != 0)
            {
                fprintf(stderr, "Bad file format : p edges N E C incorrect\n");
	            exit(-1);
            }
            pch = strtok (NULL, delim);
            N = atoi(pch);
            pch = strtok (NULL, delim);
            E = atoi(pch);
            pch = strtok (NULL, delim);
            C = atoi(pch);
            G = graph_make(N);
            htab=calloc(1,sizeof(struct hsearch_data));
            hcreate_r(N + 100,htab);
            htabc=calloc(1,sizeof(struct hsearch_data));
            hcreate_r(C + 100,htabc);

            ids = (int *) malloc(N * sizeof(int));
            cids = (int *) malloc(C * sizeof(int));
            G->Pnum = C;
            G->psizes = (int *)malloc(G->Pnum * sizeof(int));
            G->_categoryname = (char **)malloc(G->Pnum * sizeof(char *));
            memset(G->psizes, 0, G->Pnum * sizeof(int));
            memset(G->_categoryname, 0, G->Pnum * sizeof(char *));
        }
        else if(*pch == 'e')
        {
            pch = strtok (NULL, delim);
            word1 = pch;
            pch = strtok (NULL, delim);
            word2 = pch;
            item.key = word1;
            if (hsearch_r(item, FIND, &found_item, htab) != 0) {
	            id = (int *) (found_item->data);
	            u = *id;
	        }
	        else {
	            u = k;
	            G->_label[k++] = strdup(word1);
	            item.key = G->_label[u];
	            ids[u] = u;
	            item.data = (void *) (ids+u);
                hsearch_r(item, ENTER, &found_item, htab);
	        }
            item.key = word2;
            if (hsearch_r(item, FIND, &found_item, htab) != 0) {
	            id = (int *) (found_item->data);
	            v = *id;
	        }
	        else {
	            v = k;
	            G->_label[k++] = strdup(word2);
	            item.key = G->_label[v];
	            ids[v] = v;
	            item.data = (void *) (ids+v);
                hsearch_r(item, ENTER, &found_item, htab);
	        }
            if (k > N) {
	            fprintf(stderr, "Bad file format : too many labels %d %d\n", k, N);
	            exit(-1);
	        }
	        add_edge(G, u, v);
	        edges++;
        }
        else if(*pch == 'n')
        {
            pch = strtok (NULL, delim);
            word1 = pch;
            pch = strtok (NULL, delim);
            color = pch;
            item.key = word1;
            if (hsearch_r(item, FIND, &found_item, htab) != 0) {
	            id = (int *) (found_item->data);
	            v = *id;
            }
            else {
	            v = k;
	            G->_label[k++] = strdup(word1);
	            item.key = G->_label[v];
	            ids[v] = v;
	            item.data = (void *) (ids+v);
                hsearch_r(item, ENTER, &found_item, htab);
	        }
            citem.key = color;
            if (hsearch_r(citem, FIND, &found_citem, htabc) != 0) {
	            id = (int *) (found_citem->data);
	            cid = *id;
            }
            else{
                G->_categoryname[cid_use] = strdup(color);
                cid = cid_use;
                cid_use++;
                cids[cid] = cid;
                citem.key = G->_categoryname[cid];
                citem.data = (void *) (cids + cid);
                hsearch_r(citem, ENTER, &found_citem, htabc);
            }
            G->_category[v] = cid;
            G->psizes[cid]++;
        }
        else{
            fprintf(stderr, "unknown command %c\n", *pch);
        }
    }
    if (edges != E) { 
	    fprintf(stderr, "edgelist_in : # of edges incorrect\n");
	    fprintf(stderr, "edgelist_in : %d vertices, %d edges\n", k, edges);
    }
    if (k != N) {
	    fprintf(stderr, "edgelist_in : # of vertices incorrect\n");
	    fprintf(stderr, "edgelist_in : %d vertices, %d edges\n", k, edges);
	    G->_num_vertices = k;
	    G->_num_active_vertices = k;
    }
    graph_build_same_category(G);

    if (htab != NULL) {
        hdestroy_r(htab);
        free(htab);
    }
    if (htabc != NULL) {
        hdestroy_r(htabc);
        free(htabc);
    }
    if (ids != NULL) free(ids);
    if (cids != NULL) free(cids);
    if (line != NULL) free(line);
    return G;
}

void printGraph(Graph * G, FILE * fout)
{
    int i;
    if(fout != NULL)
    {AdjMatrix_out(fout, G);}
    for(i = 0; i < G->_num_vertices; i++)
    {
        printf("node %d %s color %d %s\n", i, G->_label[i], G->_category[i], G->_categoryname[G->_category[i]]);
    }
}

void GetConfig(char * confile, Graph * G)
{
    char * line = NULL;
    size_t len;
    char * pch;
    ssize_t read;
    char * delim = " \t";
    int i;
    const char * tlb = "lb";
    FILE * fconf;
    if ((fconf = fopen(confile, "r")) == NULL) 
    {
        fprintf(stderr, "Can't open configure file %s\n", confile);
        exit(-1);
    }
    int lb;
    G->lbs = (int *) malloc(G->Pnum * sizeof(int));
    while ((read = getline(&line, &len, fconf)) != -1) {
        while(line[strlen(line) - 1] == '\n' || line[strlen(line) - 1] == '\x0d'){
            line[strlen(line) - 1] = 0;
            if(strlen(line) == 0)
            {break;}
        }
        if(strlen(line) == 0)
        {continue;}
        pch = NULL;
        pch = strtok(line, delim);
        if(strcmp(pch, tlb) == 0){
            pch = strtok (NULL, delim);
            lb = atoi(pch);
            for(i = 0; i < G->Pnum; i++)
            {G->lbs[i] = lb;}
        }
    }
    fclose(fconf);
}
