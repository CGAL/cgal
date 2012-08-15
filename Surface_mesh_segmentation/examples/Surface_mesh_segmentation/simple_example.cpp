#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_polyhedron_triangle_primitive<K,Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

int main(int argc,char** argv)
{
  if (argc==1){
    std::cerr << "Please provide a file a name to an off file" << std::endl;
    return EXIT_FAILURE;
  }
  std::ifstream input(argv[1]);
  if (!input){
    std::cerr << "Could not open " << argv[1] << std::endl;
    return EXIT_FAILURE;
  }
  //create an empty polyhedron
  Polyhedron polyhedron;
  //read polyhedron from off file
  input >> polyhedron;
  
  //create an aabb tree with the polyhedron facets
  Tree aabb_tree(polyhedron.facets_begin(),polyhedron.facets_end());
  
  //create a ray form source and vector
  K::Ray_3 ray(K::Point_3(0,0,0),K::Vector_3(1,0,0));
  
  typedef std::list< Tree::Object_and_primitive_id > List;
  List all_intersections;
  
  //compute all facets intersected and intersection
  aabb_tree.all_intersections ( ray, std::back_inserter(all_intersections) );
  
  for (List::iterator  it=all_intersections.begin();it!=all_intersections.end();++it)
  {
    //check the ray intersect facets in point
    const K::Point_3* point = CGAL::object_cast<K::Point_3>( &(it->first) );
    if ( point ){
      std::cout << *point << std::endl;
    }
  }
}