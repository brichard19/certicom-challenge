#ifndef _MPI_HELPER_H
#define _MPI_HELPER_H

#include <mpi.h>

#define MPI_CALL(condition)                                                                        \
  {                                                                                                \
    int err = condition;                                                                           \
    if(err != MPI_SUCCESS) {                                                                       \
      std::cout << "MPI error " << err << " " << __LINE__ << std::endl;                            \
      MPI_Finalize();                                                                              \
      exit(1);                                                                                     \
    }                                                                                              \
  }

#endif