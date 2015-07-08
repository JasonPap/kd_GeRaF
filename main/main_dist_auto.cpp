
#include <mpi.h>
#include <iostream>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "../source/Distributed_auto_random_kd_forest.h"

using namespace std;

int main(void)
{	
	size_t N = 100000, D = 2000, Q = 10000;

    //double start_time = MPI_Wtime();
	MPI_Init(NULL, NULL); 
	int p = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();
	
	Params mypars;
	mypars.points_per_leaf = 1024;
	mypars.trees_no = 128;
	mypars.t = 128;
	mypars.max_leaf_check = 64;
	mypars.rotate_option = No;
	mypars.shuffle_enable = false;
	mypars.sample_size = N/p;

  	clock_t start_time = clock();

	Distributed_auto_random_kd_forest<float>* dforest;
	dforest = new Distributed_auto_random_kd_forest<float>(N, D, "./test_files/data.txt", 0.0, &mypars);
	MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0)
		std::cout<<"Build in "<<(clock() - start_time)/CLOCKS_PER_SEC<<" s."<<std::endl;
	clock_t mid_c = clock();

	std::vector< std::vector<std::pair<float, int> > > res;
	dforest->perform_queries(Q, "./test_files/query.txt", 1, 0.0, res);
	
	if(rank==0)
		std::cout<<"search took: "<<(clock() - mid_c)/CLOCKS_PER_SEC<<" s."<<std::endl;   

	delete(dforest);
	return 0;
}