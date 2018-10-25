#define CGAL_TESTSUITE

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
typedef Kernel::Point_2                                      Point_2;
typedef CGAL::Surface_mesh<Point>                            Surface_mesh;
//typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;

typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
typedef boost::property_map<Surface_mesh, Vertex_distance_tag >::type Vertex_distance_map;

typedef CGAL::Heat_method_3::Heat_method_3<Surface_mesh,Kernel, CGAL::Tag_true> Heat_method;
typedef CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<CGAL::Eigen_sparse_matrix<double>::EigenType > > Solver_traits;
typedef Solver_traits::Matrix SparseMatrix;


struct Heat_method_3_private_tests {
  
#if 0
void source_set_tests(Heat_method hm, const Idt& sm)
{
  vertex_descriptor source = *(vertices(sm).first);
  hm.add_source(source);
  assert(*(hm.sources_begin()) == source);
  assert(hm.remove_source(source));
  assert(hm.add_source(*(vertices(sm).first)));
  assert(*(hm.sources_begin()) == source);
  assert(hm.add_source(*(std::next(vertices(sm).first,3))));
  assert(source != *(std::next(vertices(sm).first,3)));
  assert(*(hm.sources_begin()) == source);
  assert(!(hm.add_source(source)));
  assert(hm.remove_source(*(std::next(vertices(sm).first,3))));

}

#endif

#if 1

void cotan_matrix_test(const SparseMatrix& c)
{
  double sum = 0;
  for(int k = 0; k<c.eigen_object().outerSize(); ++k)
  {
    for(SparseMatrix::EigenType::InnerIterator it(c.eigen_object(),k); it; ++it)
    {
      sum +=it.value();
    }
  }
  //Every row should sum up to 0, allow for slight error for large meshes
  std::cout<<"sum is: "<< sum << "\n";
  assert(sum < 1e-3);
}

void mass_matrix_test(const SparseMatrix& M)
{
  double sum = 0;
  for(int k = 0; k<M.eigen_object().outerSize(); ++k)
    {
      for(SparseMatrix::EigenType::InnerIterator it(M.eigen_object(),k); it; ++it)
      {
        sum +=it.value();
      }
    }
    //total Area matrix should be equal to the sum of all faces on the mesh
    //have to allow for the error because of rounding issues: Andreas might be able to help with this?
    //this will only work for the pyramid mesh
    assert((sum-1.866025)<=0.000005);
}

void check_for_zero(const Eigen::VectorXd& u)
{
  for(int c_i = 0; c_i<4; c_i++)
  {
      assert(u(c_i,0)<0.00001);
  }
}

void check_for_unit(const Eigen::MatrixXd& X, int dimension)
{
  for(int k = 0; k<dimension; k++)
  {
    double sum = CGAL::sqrt(X(k,0)*X(k,0) + X(k,1)*X(k,1) + X(k,2)*X(k,2));
    assert((sum-1)<0.00001);
  }
}
#if 0
void check_no_update(const Idt& sm, const Vertex_distance_map& original, const Vertex_distance_map& updated)
{
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm))
  {
    assert(get(original, vd) == get(updated,vd));
  }
}
#endif

#endif


int main()
{
  Surface_mesh sm;

  std::ifstream in("../data/disk.off");   //(argc>1)?argv[1]:"data/pyramid1.off");
  in >> sm;
  if(!in || num_vertices(sm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }
  Vertex_distance_map vdm = get(Vertex_distance_tag(),sm);

  //source set tests
  Heat_method hm(sm);

  hm.add_source(* vertices(sm).first);

  hm.fill_distance_map(vdm);

  BOOST_FOREACH(boost::graph_traits<Surface_mesh>::vertex_descriptor vd, vertices(sm)){
    std::cout << get(vdm,vd) << std::endl;
  }

  // source_set_tests(hm,idt);

#if 1

  //cotan matrix tests
  //const SparseMatrix& M = hm.mass_matrix();
  //std::cout<<"and M is: "<< Eigen::MatrixXd(M) << "\n";
  const SparseMatrix& c = hm.cotan_matrix();
  cotan_matrix_test(c);
  //mass_matrix_test(M);

  //double time_step = hm.time_step();
  //double length_sum = hm.summation_of_edges();
  //there are 6 edges in pyramid
  //double time_step_computed = (1./6)*length_sum;
  ///assert(time_step_computed*time_step_computed ==time_step);


  const SparseMatrix& K = hm.kronecker_delta();
  // AF: I commented the assert as I commented in build()
  assert(K.eigen_object().nonZeros()==1);
//  Eigen::VectorXd solved_u = hm.solve_cotan_laplace(M,c,K,time_step,19768);
//  Eigen::VectorXd check_u = ((M+time_step*c)*solved_u)-K;
//  check_for_zero(check_u);
//  Eigen::MatrixXd X = hm.compute_unit_gradient(solved_u);
//  check_for_unit(X,3);
//  SparseMatrix XD = hm.compute_divergence(X,19768);
//  Eigen::VectorXd solved_dist = hm.solve_phi(c, XD, 19768);

#endif

  return 0;
}

};

int main()
{
  Heat_method_3_private_tests tests;
  return tests.main();
}
