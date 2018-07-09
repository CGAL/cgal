#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Heat_method_3/Heat_method_3.h>
#include <CGAL/Heat_method_3/Intrinsic_Delaunay_Triangulation_3.h>
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
typedef CGAL::Surface_mesh<Point>                            BaseMesh;
typedef CGAL::dynamic_halfedge_property_t<Point_2> Halfedge_coordinate_tag;
typedef boost::property_map<BaseMesh, Halfedge_coordinate_tag >::type Halfedge_coordinate_map;

typedef CGAL::Intrinsic_Delaunay_Triangulation_3::Intrinsic_Delaunay_Triangulation_3<BaseMesh,Kernel, Halfedge_coordinate_map> Mesh;


typedef CGAL::dynamic_vertex_property_t<double> Vertex_distance_tag;
typedef boost::property_map<Mesh, Vertex_distance_tag >::type Vertex_distance_map;


typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

typedef CGAL::Heat_method_3::Heat_method_3<Mesh,Kernel,Vertex_distance_map> Heat_method;
typedef CGAL::Heat_method_3::Heat_method_Eigen_traits_3::SparseMatrix SparseMatrix;


#if 0
void source_set_tests(Heat_method hm, const Mesh& sm)
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

void cotan_matrix_test(const SparseMatrix& c)
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
  std::cout<<"sum is: "<< sum << "\n";
  assert(sum < 0.000000001);
}

void mass_matrix_test(const SparseMatrix& M)
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

void check_no_update(const Mesh& sm, const Vertex_distance_map& original, const Vertex_distance_map& updated)
{
  BOOST_FOREACH(vertex_descriptor vd, vertices(sm))
  {
    assert(get(original, vd) == get(updated,vd));
  }
}

#endif


int main()
{
  BaseMesh bm;
  Halfedge_coordinate_map hcm;
  Mesh sm(bm,hcm);

  Vertex_distance_map vertex_distance_map = get(Vertex_distance_tag(),sm);
  bool idf = false;

  std::ifstream in("data/pyramid0.off");
  in >> bm;
  if(!in || num_vertices(bm) == 0) {
    std::cerr << "Problem loading the input data" << std::endl;
    return 1;
  }

  put(vertex_distance_map, * vertices(sm).first, 1.0);
  put(vertex_distance_map, * halfedges(sm).first, 1.0);
  

  //source set tests
  Heat_method hm(sm, vertex_distance_map);
#if 0
  source_set_tests(hm,sm);
  //cotan matrix tests
  const SparseMatrix& M = hm.mass_matrix();
  //std::cout<<"and M is: "<< Eigen::MatrixXd(M) << "\n";
  const SparseMatrix& c = hm.cotan_matrix();
  cotan_matrix_test(c);
  mass_matrix_test(M);


  double time_step = hm.time_step();
  double length_sum = hm.summation_of_edges();
  //there are 6 edges in pyramid
  double time_step_computed = (1./6)*length_sum;
  assert(time_step_computed ==time_step);


  const SparseMatrix& K = hm.kronecker_delta();
  // AF: I commented the assert as I commented in build()
  assert(K.nonZeros()==1);
  Eigen::VectorXd solved_u = hm.solve_cotan_laplace(M,c,K,time_step,4);
  Eigen::VectorXd check_u = ((M+time_step*c)*solved_u)-K;
  check_for_zero(check_u);
  Eigen::MatrixXd X = hm.compute_unit_gradient(solved_u);
  check_for_unit(X,3);

  SparseMatrix XD = hm.compute_divergence(X,4);

  Eigen::VectorXd solved_dist = hm.solve_phi(c, XD,4);

#endif
  
  return 0;
}
