#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Dynamic_property_map.h>
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
//typedef CGAL::Polyhedron_3<Kernel> Mesh;

typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
typedef boost::property_map<Mesh, Vertex_distance_tag >::type Vertex_distance_map;
 
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef CGAL::Heat_method_3::Heat_method_3<Mesh,Kernel,Vertex_distance_map> Heat_method;
typedef CGAL::Heat_method_3::Heat_method_Eigen_traits_3::SparseMatrix SparseMatrix;

void source_set_tests(Heat_method hm, Mesh sm)
{
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  const std::set<vertex_descriptor>& source_copy = hm.get_sources();
  assert(*(source_copy.begin()) == source);;
  assert(*(hm.sources_begin()) == source);
  assert(hm.remove_source(source));
  assert((hm.get_sources()).empty());
  assert(hm.add_source(*(vertices(sm).first)));
  assert(*(hm.sources_begin()) == source);
  assert(hm.add_source(*(std::next(vertices(sm).first,3))));
  assert(source != *(std::next(vertices(sm).first,3)));
  assert(*(hm.sources_begin()) == source);
  assert(!(hm.add_source(source)));
  assert(hm.remove_source(*(std::next(vertices(sm).first,3))));
}

void cotan_matrix_test(SparseMatrix c)
{
  double sum = 0;
  for(int k = 0; k<c.outerSize(); ++k)
  {
    for(SparseMatrix::InnerIterator it(c,k); it; ++it)
    {
      sum +=it.value();
    }
  }
  //Every row should sum up to 0, allow for slight error for large meshes
  assert(sum < 0.000000001);
}

void mass_matrix_test(SparseMatrix M)
{
  double sum = 0;
  for(int k = 0; k<M.outerSize(); ++k)
    {
      for(SparseMatrix::InnerIterator it(M,k); it; ++it)
      {
        sum +=it.value();
      }
    }
    //total Area matrix should be equal to the sum of all faces on the mesh
    //have to allow for the error because of rounding issues: Andreas might be able to help with this?
    //this will only work for the pyramid mesh
    assert((sum-1.866025)<=0.000005);
}

void check_for_zero(Eigen::VectorXd u)
{
  for(int c_i = 0; c_i<4; c_i++)
  {
      assert(u(c_i,0)<0.00001);
  }
}

void check_for_unit(Eigen::MatrixXd X, int dimension)
{
  for(int k = 0; k<dimension; k++)
  {
    double sum = CGAL::sqrt(X(k,0)*X(k,0) + X(k,1)*X(k,1) + X(k,2)*X(k,2));
    assert((sum-1)<0.00001);
  }
}

int main()
{
  Mesh sm;
  Vertex_distance_map vertex_distance_map = get(Vertex_distance_tag(),sm);
  
  
  std::ifstream in("data/pyramid0.off");
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  //source set tests
  Heat_method hm(sm, vertex_distance_map);
  source_set_tests(hm,sm);
  //cotan matrix tests
  const SparseMatrix& M = hm.get_mass_matrix();
  const SparseMatrix& c = hm.get_cotan_matrix();
  cotan_matrix_test(c);
  mass_matrix_test(M);


  double time_step = hm.get_time_step();
  double length_sum = hm.summation_of_edges();
  //there are 6 edges in pyramid
  double time_step_computed = (1./6)*length_sum;
  assert(time_step_computed ==time_step);


  const SparseMatrix& K = hm.get_kronecker_delta();
  // AF: I commented the assert as I commented in build()
  assert(K.nonZeros()==1);
  Eigen::VectorXd solved_u = hm.solve_cotan_laplace(M,c,K,time_step,4);
  Eigen::VectorXd check_u = ((M+time_step*c)*solved_u)-K;
  check_for_zero(check_u);
  Eigen::MatrixXd X = hm.compute_unit_gradient(solved_u);
  check_for_unit(X,3);

  SparseMatrix XD = hm.compute_divergence(X,4);

  Eigen::VectorXd solved_dist = hm.solve_phi(c, XD,4);

  std::cout<<"PHASE 1 DONE \n";
  std::cout<<"PHASE 2 DONE \n";
  std::cout<<"PHASE 3 DONE \n";

  Mesh sm2;
  Vertex_distance_map vertex_distance_map2 = get(Vertex_distance_tag(),sm2);
  
  std::ifstream llets("data/sphere.off");
  llets>>sm2;
  if(!llets|| num_vertices(sm2) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  Heat_method hm2(sm2, vertex_distance_map2);
  //Eigen::VectorXd solved_dist_sphere = hm2.get_distances();
  const SparseMatrix& M2 = hm2.get_mass_matrix();
  const SparseMatrix& c2 = hm2.get_cotan_matrix();
  cotan_matrix_test(c2);
  //mass_matrix_test(M2);
  const SparseMatrix& K2 = hm2.get_kronecker_delta();
  // AF: I commented the assert as I commented in build()
  assert(K2.nonZeros()==1);
  double time_step_2 = hm2.get_time_step();

  Eigen::VectorXd solved_u2 = hm2.solve_cotan_laplace(M2,c2,K2,time_step_2,43562);
  Eigen::VectorXd check_u2 = ((M2+time_step_2*c2)*solved_u2)-K2;
  check_for_zero(check_u2);
  Eigen::MatrixXd X2 = hm2.compute_unit_gradient(solved_u2);
  check_for_unit(X2, 87120);
  SparseMatrix XD2 = hm.compute_divergence(X2,43562);
  Eigen::VectorXd solved_dist2 = hm.solve_phi(c2, XD2,43562);
  //verified a few of the actual values against the estimated values, avg. error was 0.0001
  //In future, want to check performance against other solver


  return 0;
}
