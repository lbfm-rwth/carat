#include <stdlib.h>
#include <stdio.h>

const char *get_data_dir(void) {
  char *dir = getenv("CARAT_DIR");
  return dir ? dir : TOPDIR;
}
