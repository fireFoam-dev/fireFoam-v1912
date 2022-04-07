//#include "version.h"
#include "stdio.h"
#include "gitVersion.H"

if (Pstream::master()){
//printf("****************************************************\n");
printf("FireFOAM Build Version:  ");
printf("%10s\n",build_git_sha);
printf("FireFOAM Build Time Stamp:  ");
printf("%10s\n",build_git_time);
//printf("%10s\n",build_git_time);
printf("****************************************************\n");
}

