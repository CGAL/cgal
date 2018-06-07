#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <Eigen/Sparse>
#include <Eigen/Dense>
typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;
typedef CGAL::Surface_mesh<Point>                            Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Heat_method_3::Heat_method_3<Mesh,Kernel> Heat_method;
typedef CGAL::Heat_method_3::Heat_method_Eigen_traits_3::SparseMatrix SparseMatrix;

int main()
{
  Mesh sm;
  std::ifstream in("data/pyramid0.off");
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  //source set tests
  Heat_method hm(sm);
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  std::set<vertex_descriptor> source_copy = hm.get_sources();
  assert(*(source_copy.begin()) == source);;
  assert(*(hm.sources_begin()) == source);
  assert(hm.remove_source(source));
  assert((hm.get_sources()).empty());
  assert(hm.add_source(*(vertices(sm).first)));
  assert(*(hm.sources_end()) == source);
  assert(*(hm.sources_begin()) == source);
  assert(hm.add_source(*(std::next(vertices(sm).first,3))));
  assert(source != *(std::next(vertices(sm).first,3)));
  assert(*(hm.sources_begin()) == source);
  assert(!(hm.add_source(source)));
  assert(hm.remove_source(*(std::next(vertices(sm).first,3))));
  //cotan matrix tests
  const SparseMatrix& M = hm.get_mass_matrix();
  const SparseMatrix& c = hm.get_cotan_matrix();
  double sum = 0;
  for(int k = 0; k<c.outerSize(); ++k)
  {
    for(SparseMatrix::InnerIterator it(c,k); it; ++it)
    {
      sum +=it.value();
    }
  }
  //Every row should sum up to 0
  assert(sum == 0);

  for(int k = 0; k<M.outerSize(); ++k)
  {
    for(SparseMatrix::InnerIterator it(M,k); it; ++it)
    {
      sum +=it.value();
    }
  }
  //total Area matrix should be equal to the sum of all faces on the mesh
  //have to allow for the error because of rounding issues: Andreas might be able to help with this?
  assert((sum-1.866025)<=0.000005);

  double time_step = hm.get_time_step();
  double length_sum = hm.summation_of_edges();
  //there are 6 edges in pyramid
  double time_step_computed = (1./6)*length_sum;
  assert(time_step_computed ==time_step);

  //mass matrix tests

  const SparseMatrix& K = hm.get_kronecker_delta();
  assert(K.nonZeros()==1);


  std::cout<<"SUCCESS";
  return 0;
}
