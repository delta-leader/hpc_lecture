#include <time.h>

double get_time() {
  return clock() / CLOCKS_PER_SEC;
}
