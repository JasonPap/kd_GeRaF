
#include <mpi.h>
#include <iostream>
  //vv
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "../source/Distributed_auto_random_kd_forest.h"

using namespace std;

int main(void)
{	
	size_t N = 1000000, D = 200;

    //double start_time = MPI_Wtime();
	MPI_Init(NULL, NULL); 
	int p = MPI::COMM_WORLD.Get_size();
	
	Params mypars;
	mypars.points_per_leaf = 64;
	mypars.trees_no = 64;
	mypars.t = 64;
	mypars.max_leaf_check = 512;
	mypars.rotate_option = No;
	mypars.shuffle_enable = false;
	mypars.sample_size = N/p;

	//vv
  	clock_t start_time = clock();

	Distributed_auto_random_kd_forest<float>* dforest;
	dforest = new Distributed_auto_random_kd_forest<float>(N, D, "./test_files/data.txt", 0.0, &mypars);
	MPI_Barrier(MPI_COMM_WORLD);

    double mid_time = MPI_Wtime();
      //vv
	std::cout<<"Build in "<<(clock() - start_time)/CLOCKS_PER_SEC<<" s."<<std::endl;
	clock_t mid_c = clock();

	std::vector< std::vector<std::pair<float, int> > > res;
	dforest->perform_queries(50000, "./test_files/query.txt", 1, 0.0, res, false);
	
	std::cout<<(clock() - mid_c)/CLOCKS_PER_SEC<<" s."<<std::endl;

	MPI_Barrier(MPI_COMM_WORLD);
    double end_time = MPI_Wtime();
    int rank = MPI::COMM_WORLD.Get_rank();
    if (rank == 0)
    {
    	//std::cout<<"Build time: "<<mid_time - start_time << " s."<<std::endl;
    	std::cout<<"Search time: "<<end_time - mid_time << " s."<<std::endl;
    }

	delete(dforest);
	return 0;
}