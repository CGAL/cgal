#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/corefinement_operations.h>
#include <CGAL/iterator.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polyhedron_corefinement<Polyhedron> Corefinement;
typedef std::vector<Kernel::Point_3> Polyline;



int main(int argc,char** argv)
{
  if (argc!=3){
    std::cerr << "Usage: " << argv[0] << "file1.off file2.off" << std::endl;
    return 1;
  }

  //read input polyhedra
  Polyhedron A,B;
  std::ifstream input;
  
  input.open(argv[1]);
  if (!input) {std::cerr << "Error cannot read " << argv[1] << std::endl; return 1;}
  input >> A;
  input.close();

  input.open(argv[2]);
  if (!input) {std::cerr << "Error cannot read " << argv[2] << std::endl; return 1;}
  input >> B;
  input.close();
  
  //check validity of polyhedra
  if (!A.is_pure_triangle() || !B.is_pure_triangle()){
    std::cerr << "Inputs polyhedra must be triangulated." << std::endl;
    return 1;
  }

  if (!A.size_of_vertices() || !B.size_of_vertices()){
    std::cerr << "Inputs polyhedra must not be empty." << std::endl;
    return 1;
  }

  if (!A.is_valid() || !B.is_valid()){
    std::cerr << "Inputs polyhedra must be valid." << std::endl;
    return 1;
  }
  
  //define the vector to contain result.
  std::vector<std::pair <Polyhedron*,int> > result;
  result.reserve(2);
  
  //Define the functor performing corefinement and boolean operations
  Corefinement coref;
  CGAL::Emptyset_iterator polyline_output;  //if you are interested by the intersection polyline,
  //                                          declare std::list<Polyline> polylines;  and use 
  //                                          std::back_inserter(polylines) instead of polyline_output
  
  //run the corefinement. The sum of tags indicates we are interested in the union and the intersection of
  // the input polyhedra.
  coref(A,B,polyline_output,std::back_inserter(result),Corefinement::P_minus_Q_tag);
  
  std::ofstream output;
 
  //write result into off files.
  output.open("res.off");
  output << *( result[0].first );
  output.close();

  //delete result polyhedra when no longer needed
  delete result[0].first;

  return 0;
}

