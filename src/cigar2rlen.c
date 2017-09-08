#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void cigar2rlen(char **cs, int *n, int *result)
{

  const char KEEP[] = "MDNP=XISH";

  int x;
  for (x = 0; x < *n; x++) {
    int csLen = strlen(cs[x]);
    int start = 0;
    int i = 0;
    int len = 0;
    while (start <= csLen - 1) {
      i = strcspn(&cs[x][start], KEEP);
      char flag[1];
      memcpy(flag, &cs[x][start + i], 1);
      if (strpbrk(flag, "ISH") == NULL) {
        char subLen[i];
        memcpy(subLen, &cs[x][start], i);
        len += atoi(subLen);
      }
      start += i + 1;
    }
    result[x] = len;
  }

}

