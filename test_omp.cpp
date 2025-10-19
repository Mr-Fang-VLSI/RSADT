#include <omp.h>
#include <iostream>
int main(){ std::cout<<"Threads="<<omp_get_max_threads()<<std::endl; }
