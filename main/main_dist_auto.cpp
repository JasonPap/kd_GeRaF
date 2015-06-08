
#include <mpi.h>
#include <iostream>
#include "../source/Distributed_auto_random_kd_forest.h"

using namespace std;

int main(void)
{	
	size_t N = 11, D = 2;
	Distributed_auto_random_kd_forest<int>* dforest;
	cout << "About to build dist forest" << endl;
	dforest = new Distributed_auto_random_kd_forest<int>(N, D, "testinputfile", 0.0);
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "Dforest created, ready to rock" << endl;
	return 0;
}