#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/random.hpp>


#include <cassert>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point          Point;

typedef CGAL::Creator_uniform_2<double, Point>             Creator;
typedef CGAL::Random_points_in_square_2<Point, Creator> Point_generator;


int main(int argc,char** argv )
{
  int n_segments=100000;
  if (argc==2) n_segments=atoi(argv[1]);
  
  
  CGAL::Random rand(0);
  std::vector<Point> point_set;
  point_set.reserve(2*n_segments);
  
  CGAL::cpp11::copy_n(Point_generator(1,rand), 2*n_segments,std::back_inserter(point_set));

  std::cout << point_set.size()/2  << " segments" << std::endl;
  
  CDT cdt;

  CGAL::Timer time;
  time.start();
  cdt.insert( point_set.begin(),point_set.end() );
  time.stop();
  
  std::cout << "Inserting points in " << time.time() << std::endl;
  time.reset();
  
  std::vector<CDT::Vertex_handle> vertex_handles;
  vertex_handles.reserve(2*n_segments);
  
  for (CDT::Finite_vertices_iterator vit=cdt.finite_vertices_begin(),
                                     vit_end=cdt.finite_vertices_end();vit!=vit_end;++vit)
  {
    vertex_handles.push_back(vit);
  }
  
  boost::rand48 random;
  boost::random_number_generator<boost::rand48> rng(random);
  std::random_shuffle(vertex_handles.begin(),vertex_handles.end(),rng);
  
  time.start();
  for (int i=0;i<n_segments;++i){
    cdt.insert_constraint( vertex_handles[2*i],vertex_handles[2*i+1] );
  }
  time.stop();
  
  std::cout << "Adding constraints in " << time.time() << std::endl;
  
  std::cout << cdt.number_of_vertices()-static_cast<unsigned>(2*n_segments) << " intersection points\n";
  
  return 0;
}
