/**
 @file Distributed_auto_random_kd_forest.h
 */

#ifndef DISTRIBUTED_AUTO_RANDOM_KD_TREE_H
#define DISTRIBUTED_AUTO_RANDOM_KD_TREE_H

#include <mpi.h>
#include "MPI_functions.h" 
#include "Auto_random_kd_forest.h"

template <typename T>
class Distributed_auto_random_kd_forest
{
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
		proc_num = initialize_MPI();
		rank = MPI::COMM_WORLD.Get_rank();

		if(rank == 0){
			//master process read the file and distributes points to slaves
			if (file_option == 0) {
		      	split_file(datafile, rank, proc_num, 100);
		    } else if (file_option == 1) {
		      	std::cout << "File option 1 not supported by distributed kd-forest" << std::endl;
    			Finalize();
    			exit(1);
    		}
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
	std::string local_datafile

	
	/* 
	 * \brief Initializes the MPI framework 
	 * @return	-	the number of MPI processes 
	 */ 
	int initialize_MPI(){ 
	    MPI::Init(NULL, NULL); 
		int p = MPI::COMM_WORLD.Get_size();
	    return p; 
	}

	void split_file(const std::string& datafile, const int rank, const int size, const int overlap) {
	    MPI_File in, out;
	    //every process opens the datafile for input
	    ierr = MPI::File::Open(MPI::COMM_WORLD, datafile, MPI_MODE_RDONLY, MPI_INFO_NULL, &in);
	    if (ierr) {
	        if (rank == 0) 
	        	std::err << "Couldn't open file " << datafile << std::endl;
	        Finalize();
	        exit(1);
	    }

	    local_datafile = datafile;
	    local_datafile.append(std::to_string(rank));
	    char* cfilename = new char[local_datafile.length() + 1];
	    strcpy(cfilename, local_datafile.c_str());
	    printf("cfilename = %s\n", cfilename);
	    //each process opens its own output file
	    ierr = MPI::File::Open(MPI::COMM_SELF, cfilename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &out);
	    if (ierr) {
	        if (rank == 0) 
	        	std::err << "Couldn't open output file " << std::endl;
	        Finalize();
	        exit(1);
	    }
	    delete[] cfilename;

	    MPI_Offset localsize;
	    MPI_Offset globalstart;
	    int mysize;
	    char *chunk;

	    /* read in relevant chunk of file into "chunk",
	     * which starts at location in the file globalstart
	     * and has size mysize 
	     */
	    {
	        MPI_Offset globalend;
	        MPI_Offset filesize;

	        /* figure out who reads what */
	        MPI_File_get_size(*in, &filesize);
	        filesize--;  /* get rid of text file eof */
	        mysize = filesize/size;
	        globalstart = rank * mysize;
	        globalend   = globalstart + mysize - 1;
	        if (rank == size-1) globalend = filesize-1;

	        /* add overlap to the end of everyone's chunk except last proc... */
	        if (rank != size-1)
	            globalend += overlap;

	        mysize =  globalend - globalstart + 1;

	        /* allocate memory */
	        chunk = malloc( (mysize + 1)*sizeof(char));
	        if (chunk == NULL){
	        	std::err << "Could not allocate memory for new file, exiting..." << std::endl;
	        	Finalize();
	        	exit(1);
	        }

	        /* everyone reads in their part */
	        MPI_File_read_at_all(*in, globalstart, chunk, mysize, MPI_CHAR, MPI_STATUS_IGNORE);
	        chunk[mysize] = '\0';
	    }

	    /*
	     * everyone calculate what their start and end *really* are by going 
	     * from the first newline after start to the first newline after the
	     * overlap region starts (eg, after end - overlap + 1)
	     */
	    int locstart=0, locend=mysize-1;
	    if (rank != 0) {
	        while(chunk[locstart] != '\n') locstart++;
	        locstart++;
	    }
	    if (rank != size-1) {
	        locend-=overlap;
	        while(chunk[locend] != '\n') locend++;
	    }
	    mysize = locend-locstart+1;

	    /* output the processed file */
	    MPI_File_write_at_all(*out, (MPI_Offset)(globalstart+(MPI_Offset)locstart), &(chunk[locstart]), mysize, MPI_CHAR, MPI_STATUS_IGNORE);
	    MPI_File_close(&in);
    	MPI_File_close(&out);

	    return;
}


}


 #endif //DISTRIBUTED_AUTO_RANDOM_KD_TREE_H
