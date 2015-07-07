/** \example main_par.cpp
 * This is an example of how to use parallel building of the forest
 * and then search it efficiently.
 */

#include <string>
  //vv
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */
#include "../source/Auto_random_kd_forest.h"

int main(int argc, char *argv[]) {
  size_t N = 1000000, D = 200, Q = 50000;
  int k = 1;
  double epsilon = 0.0;
  std::string datafile = "test_files/data.txt", queryfile = "test_files/query.txt";
  std::vector< std::vector<std::pair<float, int> > > res;

  
  Params mypars;
  mypars.points_per_leaf = 64;
  mypars.trees_no = 128;
  mypars.t = 64;
  mypars.max_leaf_check = 1024;
  mypars.rotate_option = No;
  mypars.shuffle_enable = false;
  mypars.sample_size = N;

    //vv
  clock_t start_time = clock();

  Auto_random_kd_forest<float>  RKDf(N, D, datafile, epsilon, &mypars);

    //vv
  std::cout<<"Build in "<<(clock() - start_time)/CLOCKS_PER_SEC<<" s."<<std::endl;
  clock_t mid_time = clock();

  RKDf.perform_queries(Q, queryfile, k, epsilon, res);

    //vv
  std::cout<<(clock() - mid_time)/CLOCKS_PER_SEC<<" s."<<std::endl;

  std::cout << "\nRESULTS\n";
  for (std::vector<std::pair<float, int> >::const_iterator it = res[0]
      .begin(); it != res[0].end(); ++it)
    std::cout << it->first << ", index = " << it->second << std::endl;

  std::cout << "main_par successfully terminated" << std::endl;
  return 0;
}
