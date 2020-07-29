#include <CGAL/Simple_cartesian.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/point_generators_3.h>
#include <vector>
#include <CGAL/Timer.h>
#include <cstdlib>
#include <CGAL/Memory_sizer.h>


typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Search_traits_3<K> TreeTraits;



typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits> OK_search;
typedef CGAL::K_neighbor_search<TreeTraits> K_search;
typedef CGAL::Incremental_neighbor_search<TreeTraits> I_search;
typedef CGAL::Orthogonal_incremental_neighbor_search<TreeTraits> OI_search;

typedef CGAL::Random_points_in_cube_3<Point> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;


typedef OK_search::Point_d test1;
typedef OK_search::FT test2;
typedef OK_search::Point_with_transformed_distance test3;
typedef OK_search::Query_item test4;
typedef OK_search::iterator test5;
typedef OK_search::Tree test6;


void one_run(int NB_INPUT_POINTS,int VALUE_OF_K,int NB_NEIGH){
  if(VALUE_OF_K < NB_NEIGH || NB_NEIGH > NB_INPUT_POINTS){
    std::cout << "must have VALUE_OF_K >= NB_NEIGH and NB_NEIGH <= NB_INPUT_POINTS) " << std::endl;
    exit(EXIT_FAILURE);
  }

  CGAL::Timer time;
  CGAL::Memory_sizer m;
  std::size_t ok_mem,k_mem,oi_mem,i_mem;
  double ok_time,k_time,oi_time,i_time;

  std::cout << "Nb input points: " <<NB_INPUT_POINTS << std::endl;
  std::cout << "Value of K: " <<VALUE_OF_K << std::endl;
  std::cout << "Nb neighbor requested: " <<NB_NEIGH << std::endl;

  Point ok_pt(0,0,0),k_pt(0,0,0),oi_pt(0,0,0),i_pt(0,0,0);


  time.start();
  // generator for random data points in the cube ( (-1,-1,-1), (1,1,1) )
  std::cout << "Generating points... ";
  Random_points_iterator rpit( 1.0);
  std::vector<Point> points; points.reserve(NB_INPUT_POINTS);
  points.insert(points.begin(),N_Random_points_iterator(rpit,0), N_Random_points_iterator(NB_INPUT_POINTS));
  Point query(0,0,0);
  time.stop(); std::cout << time.time() << " "; time.reset();
  std::size_t mem_points=m.virtual_size()+m.resident_size();
  std::cout << mem_points << std::endl;

  {
    //Orthogonal K search
    std::cout << "OK search... "; time.start();
    OK_search::Tree ok_tree(points.begin(),points.end());
    OK_search ok_search(ok_tree, query,VALUE_OF_K);
    OK_search::iterator ok_it=ok_search.begin(), ok_end=ok_search.end();
    for (int j=0; (j < NB_NEIGH-1)&&(ok_it!=ok_end); ++j,++ok_it) {}
    time.stop(); ok_time=time.time(); std::cout << ok_time << " "; time.reset();
    if (NB_INPUT_POINTS!=0){
      assert(ok_it!=ok_end);
      ok_pt=(*ok_it).first;
    }
    else
      assert(ok_it==ok_end);

    ok_mem=(m.virtual_size()+m.resident_size())-mem_points;
    std::cout << ok_mem << std::endl;
  }


  {
    //K search
    std::cout << "K search... "; time.start();
    K_search::Tree k_tree(points.begin(),points.end());
    K_search k_search(k_tree, query,VALUE_OF_K);
    K_search::iterator k_it=k_search.begin(), k_end=k_search.end();
    for (int j=0; (j < NB_NEIGH-1)&&(k_it!=k_end); ++j,++k_it) {}
    time.stop(); k_time=time.time(); std::cout << k_time << " "; time.reset();
    if (NB_INPUT_POINTS!=0){
      assert(k_it!=k_end);
      k_pt=(*k_it).first;
    }
    else
      assert(k_it==k_end);
    k_mem=(m.virtual_size()+m.resident_size())-mem_points;
    std::cout << k_mem << std::endl;
  }

  {
    //Orthogonal incremental search
    std::cout << "OI search... "; time.start();
    OI_search::Tree oi_tree(points.begin(),points.end());
    OI_search oi_search(oi_tree, query);
    OI_search::iterator oi_it=oi_search.begin(), oi_end=oi_search.end();
    for (int j=0; (j < NB_NEIGH-1)&&(oi_it!=oi_end); ++j,++oi_it) {}
    time.stop(); oi_time=time.time(); std::cout << oi_time << " "; time.reset();
    if (NB_INPUT_POINTS!=0){
      assert(oi_it!=oi_end);
      oi_pt=(*oi_it).first;
    }
    else
      assert(oi_it==oi_end);
    oi_mem=(m.virtual_size()+m.resident_size())-mem_points;
    std::cout << oi_mem << std::endl;
  }

  {
    //incremental search
    std::cout << "I search... "; time.start();
    I_search::Tree i_tree(points.begin(),points.end());
    I_search i_search(i_tree, query);
    I_search::iterator i_it=i_search.begin(), i_end=i_search.end();
    for (int j=0; (j < NB_NEIGH-1)&&(i_it!=i_end); ++j,++i_it) {}
    time.stop(); i_time=time.time(); std::cout << i_time << " "; time.reset();
    if (NB_INPUT_POINTS!=0){
      assert(i_it!=i_end);
      i_pt=(*i_it).first;
    }
    else
      assert(i_it==i_end);
    i_mem=(m.virtual_size()+m.resident_size())-mem_points;
    std::cout << i_mem << std::endl;
  }

  if ( ok_pt != k_pt ) std::cout << "K different\n";
  if ( ok_pt != oi_pt ) std::cout << "OI different\n";
  if ( ok_pt != i_pt ) std::cout << "I different\n";

  assert (ok_pt == k_pt);
  assert (ok_pt == oi_pt);
  assert (ok_pt == i_pt);

  std::cerr << NB_INPUT_POINTS << " " << VALUE_OF_K << " " << NB_NEIGH << " ";
  std::cerr << ok_time  << " "  << ok_mem  << " ";
  std::cerr << k_time   << " "  << k_mem   << " ";
  std::cerr << oi_time  << " "  << oi_mem  << " ";
  std::cerr << i_time   << " "  << i_mem   << "\n";
}


int main(int argc,char** argv) {

  if(argc !=4){
    std::cout << "./Compare_methods NB_INPUT_POINTS VALUE_OF_K NB_NEIGH" << std::endl;
    int NB_INPUT_POINTS =100000;
    int VALUE_OF_K=10000;
    int NB_NEIGH = 5000;
    std::cout << "Normal run  .....\n";
    one_run(NB_INPUT_POINTS,VALUE_OF_K,NB_NEIGH);
    NB_INPUT_POINTS =100000;
    VALUE_OF_K=1000000;
    NB_NEIGH = 5000;
    std::cout << "Run with k > nb input points  .....\n";
    one_run(NB_INPUT_POINTS,VALUE_OF_K,NB_NEIGH);
    NB_INPUT_POINTS =0;
    VALUE_OF_K=1000000;
    NB_NEIGH = 0;
    std::cout << "Run with empty tree  .....\n";
    one_run(NB_INPUT_POINTS,VALUE_OF_K,NB_NEIGH);
  }
  else
  {
    int NB_INPUT_POINTS = atoi(argv[1]);
    int VALUE_OF_K = atoi(argv[2]);
    int NB_NEIGH = atoi(argv[3]);
    one_run(NB_INPUT_POINTS,VALUE_OF_K,NB_NEIGH);
  }

  return 0;
}





