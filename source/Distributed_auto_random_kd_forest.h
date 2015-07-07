/**
 @file Distributed_auto_random_kd_forest.h
 */

#ifndef DISTRIBUTED_AUTO_RANDOM_KD_TREE_H
#define DISTRIBUTED_AUTO_RANDOM_KD_TREE_H

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctype.h>
#include <mpi.h>
#include <unistd.h>
 //vv
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

#include "Auto_random_kd_forest.h"

template <typename T>
class Distributed_auto_random_kd_forest
{
private:
	/* 
	 * \brief Initializes the MPI framework 
	 * @return	-	the number of MPI processes 
	 */ 
	int initialize_MPI(){ 
	    MPI_Init(NULL, NULL); 
		int p = MPI::COMM_WORLD.Get_size();
	    return p; 
	}

public:
	/**
	 *\brief Contructor for the distributed kd forest. Takes
	 * all the parameters needed to build a auto_random_kd_forest.
	 * @param N           - size of data set
	 * @param D           - dimension of data set
	 * @param datafile    - file that contains the data set
	 * @param epsilon   	- factor of accuracy. Set to 0 for maximum accuracy
	 * @param parameters  - struct of parameters. By default, parameters are auto configured.
	 * @param file_option - how to parse the `datafile`. Default value is 0.
	 */

	Distributed_auto_random_kd_forest(size_t& N, size_t& D, const std::string& datafile,
	const double epsilon, Params* parameters = 0, const int file_option = 0){
		//proc_num = initialize_MPI();
		int proc_num = MPI::COMM_WORLD.Get_size();
		rank = MPI::COMM_WORLD.Get_rank();

		if (file_option == 0) {
	      	 local_datafile = datafile;  
		     local_datafile.append(std::to_string(rank));
	    } else if (file_option == 1) {
	      	std::cout << "File option 1 not supported by distributed kd-forest" << std::endl;
			MPI_Finalize();
			exit(1);
		}

		size_t n = N/proc_num + 1;
		
		RKDf = new Auto_random_kd_forest<T>(n, D, local_datafile, epsilon, parameters, file_option);
	}
	~Distributed_auto_random_kd_forest()
	{
		MPI_Barrier(MPI_COMM_WORLD);
		delete(RKDf);
		MPI_Finalize();
		std::cout<<"dkdf destroyed"<<std::endl;
	}

	void perform_queries(const size_t Q, const std::string& queryfile, const int k, const double epsilon, 
		std::vector<std::vector<std::pair<float, int> > >& results, const bool printPoints){

		clock_t t = clock();
		std::vector<std::vector<std::pair<float, int> > > local_results;
		RKDf->perform_queries(Q, queryfile, k, epsilon, local_results);
		MPI_Barrier(MPI_COMM_WORLD);

		if(rank == 0)
		{
			std::cout<<"Queries done, reducing results."<<std::endl;
			std::cout<<(clock() - t)/CLOCKS_PER_SEC<<" s."<<std::endl;
		}

		if ( k == 1){
			// //then MPI_REDUCE can be used for best efficiency
			// for(int i = 0; i < Q; i++){
			// 	float local_d = local_results[i][0].first;
			// 	float global_d = 0;
			// 	MPI_Reduce((void*)&local_d, &global_d, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);

			// 	// if(rank == 0)
			// 	// {
			// 	// 	std::cout<<"Query: "<<i<<" min distance = "<<global_d<<std::endl;
			// 	// }
			// 	if (printPoints){
			// 		sleep(1);
			// 		MPI_Bcast((void*)&global_d, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
			// 		if(global_d == local_d){
			// 			RKDf->print(local_results[i][0].second);
			// 		}
			// 	}

			// }

			float* local_distances = new float[Q];
			for(int i = 0; i < Q; i++){
				local_distances[i] = (float)local_results[i][0].first;
				//std::cout<<local_distances[i]<<std::endl;
			}
			MPI_Request request;
			MPI_Isend((void*)local_distances, Q, MPI_FLOAT, 0, rank, MPI_COMM_WORLD, &request);

			

			if(rank == 0){
				float** distances = new float*[proc_num];
				MPI_Request* requests = new MPI_Request[proc_num];
				for(int i = 0; i < proc_num; i++){
					distances[i] = new float[Q];
					MPI_Irecv(distances[i], Q, MPI_FLOAT, i, i, MPI_COMM_WORLD, &(requests[i]));
				}

				MPI_Status* statuses = new MPI_Status[proc_num];
				int waitAllRetVal = MPI_Waitall(proc_num, requests, statuses);
				if ( waitAllRetVal != MPI_SUCCESS ){
					std::cout<<"Problem on Receiving data"<<std::endl;
				}
				float min;
				for(int q = 0; q < Q; q++){
					min = distances[0][q];
					for(int i = 1; i < proc_num; i++){
						if(distances[i][q] < min){
							min = distances[i][q];
						}
					}
					//print min distance
					//std::cout<<"min distance for query"<<q<<" = "<<min<<std::endl;
				}
				delete[] statuses;
				delete[] requests;
				for(int i=0; i < proc_num; i++)
					delete[] distances[i];
				delete[] distances;
			}

			//MPI_Status s;
			//MPI_Wait(&request, &s);
			// sleep(rank);
			// std::cout<<rank<<" sent\n";
			//delete[] local_distances;
			// MPI_Barrier(MPI_COMM_WORLD);
			// sleep(rank);
			// std::cout<<"My rank is : "<<rank<<std::endl;

			// float* all_dist = NULL;
			// if(rank == 0){
			// 	std::cout<<"About to allocate a big chunk of memory "<< Q*proc_num<<std::endl;
			// 	all_dist = new float[Q*proc_num];
			// 	all_dist[0] = 4.3;
			// 	if (all_dist == NULL)
			// 		std::cout<<"Error allocation memory"<<std::endl;
			// std::cout<<all_dist[0]<<std::endl;
			// }
			// float* all_dist = new float[Q*proc_num];
			// MPI_Barrier(MPI_COMM_WORLD);
			// MPI_Gather((void*)local_distances, Q, MPI_FLOAT, (void*)all_dist, Q*proc_num, MPI_FLOAT, 0, MPI_COMM_WORLD);
			// std::cout<<"AAAAAAAAAAAAAAA\n";
			// // delete(local_distances);
			// if(rank == 0){
			// 	float min;
			// 	for(int i = 0; i < Q; ++i){
			// 		min = all_dist[i];
			// 		for(int j = 1; j < proc_num; ++j){
			// 			if(all_dist[i+j*Q] < min){
			// 				min = all_dist[i+j*Q];
			// 			}
			// 		}
			// 		std::cout<<min<<std::endl;
			// 	}
			// 	// sleep(10);
			// 	// delete(all_dist);
				
			// }
			// std::cout<<"BBBBBBBBBBBB\n";
			// sleep(111);

		}
		else{ 	// k > 1 Results must all be gathered by the master process
				// this functionality is not implemented, it can't be done with an
				// efficient way with the reduce function.
			std::cout<<"queries with k > 1 not supported by distributed kd-GeRaf"<<std::endl;
			return;
		}
	}



private:
	/**
	 * The kd-forest of each process
	 */
	Auto_random_kd_forest<T>* RKDf;
	/**
	 * Number of MPI processes on the global communicator
	 */
	int proc_num;
	/**
	 * Process MPI rank
	 */
	int rank;
	/**
	 * Dimension of the data set
	 */
	int D;
	/**
	 * Size of the part of the data set of the local kd-forest
	 */
	int N;
	/**
	 * Name of the local file with the points needed for the local
	 * kd-forest
	 */
	std::string local_datafile;

};


 #endif //DISTRIBUTED_AUTO_RANDOM_KD_TREE_H
