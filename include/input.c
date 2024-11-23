// input.c

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ljmd.h"

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
