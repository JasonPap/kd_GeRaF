
#include <mpi.h>
#include <iostream>
#include "../source/Distributed_auto_random_kd_forest.h"

using namespace std;

int main(void)
{	
	size_t N = 11, D = 10;
	Distributed_auto_random_kd_forest<float>* dforest;
	dforest = new Distributed_auto_random_kd_forest<float>(N, D, "./test_files/data.txt", 0.0);
	MPI_Barrier(MPI_COMM_WORLD);
	cout << "Dforest created, ready to rock" << endl;
	std::vector< std::vector<std::pair<float, int> > > res;
	dforest->perform_queries(9, "./test_files/query.txt", 1, 0.0, res, false);
	delete(dforest);
	return 0;
}