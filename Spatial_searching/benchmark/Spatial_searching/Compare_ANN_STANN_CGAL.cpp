#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Timer.h>
#include <boost/progress.hpp>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

#include <ANN/ANN.h>
#include <ANN/ANNperf.h>

#include<sfcnn.hpp>



typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Random_points_in_cube_3<Point_3> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;

typedef CGAL::Search_traits_3<Kernel> TreeTraits;
typedef CGAL::Sliding_midpoint<TreeTraits> Splitter;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> OK_search;

template <class Output_iterator>
int read_points(std::string filename,Output_iterator out)
{
  std::ifstream file(filename.c_str());
  int nb_pts;
  Point_3 p;
  file >> nb_pts;
  for (int i=0;i<nb_pts;++i){
    file >> p;
    *out++=p;
  }
  return nb_pts;
}

int main(int argc,char** argv)
{

  if (argc!=3 && argc!=4){
    std::cerr << "Usage: main nb_input_point nb_neighbors_needed\n";
    std::cerr << "or     main Input.xyz Queries.xyz nb_neighbors_needed\n";
    exit(EXIT_FAILURE);
  }

  int nb_input_point,k,nb_queries;
  std::vector<Point_3> points;
  std::vector<Point_3> queries;

  if (argc==3){
    nb_input_point=atoi(argv[1]);
    k=atoi(argv[2]);
    nb_queries=1;
    points.reserve(nb_input_point);
    Random_points_iterator rpit(1.0);
    points.insert(points.begin(),N_Random_points_iterator(rpit,0), N_Random_points_iterator(nb_input_point));
    queries.push_back(Point_3(0,0,0));
  }
  else{
    k=atoi(argv[3]);
    nb_input_point=read_points(argv[1],std::back_inserter(points));
    nb_queries=read_points(argv[2],std::back_inserter(queries));
  }

  std::cout << "Looking for " << k << " nearest neighbors amongst " << nb_input_point << " points, with "<< nb_queries <<" query points.\n";

  std::random_shuffle(points.begin(),points.end());

  double STANN_time=0,ANN_time=0,OK_time=0;

  CGAL::Timer time;
//Building trees
//--STANN
  time.start();
  sfcnn<Point_3, 3, double> NN(&points[0], nb_input_point);
  time.stop();
  double STANN_tree=time.time();
  time.reset();
//--ANN
  ANNpointArray    dataPts;     // data points
  ANNpoint         queryPt;     // query point
  ANNidxArray      nnIdx;       // near neighbor indices
  ANNdistArray     dists;       // near neighbor distances
  ANNkd_tree*      kdTree;      // search structure
  dataPts = annAllocPts(nb_input_point, 3);      // allocate data points


  queryPt = annAllocPt(3);          // allocate query point
  nnIdx = new ANNidx[k];            // allocate near neigh indices
  dists = new ANNdist[k];            // allocate near neighbor dists
  //set points
  for (int i=0;i<nb_input_point;++i){
     dataPts[i][0]=points[i][0];
     dataPts[i][1]=points[i][1];
     dataPts[i][2]=points[i][2];
  }
  time.start();
  kdTree = new ANNkd_tree(dataPts,nb_input_point,3);
  time.stop();
  double ANN_tree=time.time();
  time.reset();
//--OK_search
  time.start();
  Splitter splitter(10); //bucket size can be changed here
  OK_search::Tree ok_tree(points.begin(),points.end(),splitter);
  ok_tree.build();
  time.stop();
  double OK_tree=time.time();
  time.reset();

boost::progress_display show_progress( nb_queries );

//running NN algorithms
  for (std::vector<Point_3>::const_iterator it=queries.begin();it!=queries.end();++it)
  {
    const Point_3& query=*it;

  //STANN
    CGAL::Timer time;
    std::vector<long unsigned int> answer;
    time.start();
    NN.ksearch(query, k, answer);
    time.stop();
    STANN_time+=time.time();
    time.reset();
  //ANN
    queryPt[0]=query[0];
    queryPt[1]=query[1];
    queryPt[2]=query[2];

    time.start();
    kdTree->annkSearch(            // search
                    queryPt,            // query point
                    k,                // number of near neighbors
                    nnIdx,              // nearest neighbors (returned)
                    dists,              // distance (returned)
                    0);              // error bound

    time.stop();
    ANN_time+=time.time();
    time.reset();

    //~ ANNkdStats stats;
    //~ kdTree->getStats(stats);
    //~ std::cout << "====ANN stats ====\n";
    //~ std::cout << "dimension of space " << stats.dim << "\n";
    //~ std::cout << "no. of points " << stats.n_pts << "\n";
    //~ std::cout << "bucket size " << stats.bkt_size << "\n";
    //~ std::cout << "no. of leaves (including trivial) " << stats.n_lf << "\n";
    //~ std::cout << "no. of trivial leaves (no points) " << stats.n_tl << "\n";
    //~ std::cout << "no. of splitting nodes " << stats.n_spl << "\n";
    //~ std::cout << "no. of shrinking nodes (for bd-trees) " << stats.n_shr << "\n";
    //~ std::cout << "depth of tree " << stats.depth << "\n";
    //~ std::cout << "sum of leaf aspect ratios " << stats.sum_ar << "\n";
    //~ std::cout << "average leaf aspect ratio " << stats.avg_ar << "\n";
    //~ std::cout << "==================\n";

    //~ std::cout << "====ANN tree  ====\n";
    //~ kdTree->Print(ANNtrue,std::cout);
    //~ std::cout << "==================\n";

  //Ortho-k-NN
    std::vector<std::pair<Point_3,double> > ok_result;
    ok_result.reserve(k);
    time.start();
    OK_search ok_search(ok_tree, query,k);
    time.stop();
    std::copy(ok_search.begin(),ok_search.end(),std::back_inserter(ok_result));
    OK_time+=time.time();

    //~ std::cout << std::endl; ok_tree.statistics(std::cout); std::cout << std::endl;

    //~ std::cout << "====CGAL tree ====\n";
    //~ ok_tree.print();
    //~ std::cout << "==================\n";

    for (int i = 0; i < k; i++) {      //check results are the same
      if ( nnIdx[i]!=answer[i] ){
        std::cerr << "STANN and ANN produced different results\n";
        exit(EXIT_FAILURE);
      }
      if ( ok_result[i].first!=points[answer[i]] ){
        std::cerr << "CGAL::OK_search and STANN (ANN) produced different results\n";
        exit(EXIT_FAILURE);
      }
    }
    ++show_progress;
  }

  std::cout << "Time to build trees\n";
  std::cout << "STANN ANN OK_search\n";
  std::cout << STANN_tree << " " << ANN_tree << " " << OK_tree << "\n";
  std::cout << "Time spent for all queries\n";
  std::cout << "STANN ANN OK_search\n";
  std::cout << STANN_time << " " << ANN_time << " " << OK_time << "\n";

  delete [] nnIdx;              // clean things up
  delete [] dists;
  delete kdTree;
  annClose();
}
