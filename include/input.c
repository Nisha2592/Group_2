// input.c

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ljmd.h"

char restfile[BLEN], trajfile[BLEN], ergfile[BLEN];
int nprint;

int get_a_line(FILE *fp, char *buf) {
    char tmp[BLEN], *ptr;

    if (fgets(tmp, BLEN, fp)) {
        int i;
        ptr = strchr(tmp, '#');
        if (ptr) *ptr = '\0';
        i = strlen(tmp); --i;
        while (isspace(tmp[i])) {
            tmp[i] = '\0';
            --i;
        }
        ptr = tmp;
        while (isspace(*ptr)) { ++ptr; }
        i = strlen(ptr);
        strcpy(buf, tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

void reading(char *ll, mdsys_t *sys, char *restfile, char *trajfile, char *ergfile, int *nprint) {
    if (get_a_line(stdin, ll)) return;
    sys->natoms = atoi(ll);

    if (get_a_line(stdin, ll)) return;
    sys->mass = atof(ll);

    if (get_a_line(stdin, ll)) return;
    sys->epsilon = atof(ll);

    if (get_a_line(stdin, ll)) return;
    sys->sigma = atof(ll);

    if (get_a_line(stdin, ll)) return;
    sys->rcut = atof(ll);

    if (get_a_line(stdin, ll)) return;
    sys->box = atof(ll);

    if (get_a_line(stdin, restfile)) return;
    if (get_a_line(stdin, trajfile)) return;
    if (get_a_line(stdin, ergfile)) return;

    if (get_a_line(stdin, ll)) return;
    sys->nsteps = atoi(ll);

    if (get_a_line(stdin, ll)) return;
    sys->dt = atof(ll);

    if (get_a_line(stdin, ll)) return;
    *nprint = atoi(ll);
}