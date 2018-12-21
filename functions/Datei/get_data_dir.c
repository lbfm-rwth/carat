#include <stdlib.h>
#include <stdio.h>

void get_data_dir(char *result, const char *str) {
  char *dir;
  dir = getenv("CARAT_DIR");
  if (dir)
    sprintf(result, "%s/%s", dir, str);
  else
    sprintf(result, "%s/%s", TOPDIR, str);
}
