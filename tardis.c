#include "tardis.h"

int main(int argc, char **argv){
  parameters *params;
  int ret;
  
  init_params(&params);

  ret = parseCommandLine(argc, argv, params);
  if (ret == 0)
    exit(1);
  
  print_params(params);
}
