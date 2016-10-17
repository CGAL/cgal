#ifndef CGAL_ARAP_PARAMETERIZER_3_H
#define CGAL_ARAP_PARAMETERIZER_3_H

#include <CGAL/assertions.h>
#include <CGAL/circulator.h>
#include <CGAL/Timer.h>

#include <CGAL/OpenNL/linear_solver.h>
#include <CGAL/Eigen_solver_traits.h>

#include <CGAL/Parameterizer_traits_3.h>
#include <CGAL/Two_vertices_parameterizer_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/parameterize.h>
#include <CGAL/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/LSCM_parameterizer_3.h>

#include <CGAL/Algebraic_kernel_d_2.h>
#include <CGAL/Conic_misc.h> // @tmp used for solving conic equations

#include <iostream>
#include <fstream>
#include <set>
#include <string>
#include <vector>

/// \file ARAP_parameterizer_3.h

// @todo Determine the proper name of this file
// @todo Have an initial parameterization (LSCM or MVC) that depends on the number
//       of components
// @todo Use the algeabric kernel to solve the cubic equation
// @todo Add the post processing step
// @todo Add distortion measures
// @todo Make it work with a surface mesh type
// @todo Use a boost array for the roots
// @todo Handle the case cot(0) with a local parameterization aligned with the axes
//       (this produces C2=0 which is problematic to compute a & b)

// @todo Add support for the OpenNL solver?
// @todo The two systems A Xu = Bu and A Xv = BV could be merged in one system
//       using complex numbers?

namespace CGAL {

// ------------------------------------------------------------------------------------
// Declaration
// ------------------------------------------------------------------------------------

template
<
  class TriangleMesh, ///< a model of `FaceGraph`
  class BorderParameterizer_3
    = Two_vertices_parameterizer_3<TriangleMesh>,
    ///< Strategy to parameterize the surface border.
    ///< The minimum is to parameterize two vertices.
  class SparseLinearAlgebraTraits_d
    = Eigen_solver_traits< > // defaults to Eigen::BICGSTAB with Eigen_sparse_matrix
    ///< Traits class to solve a sparse linear system.
>
class ARAP_parameterizer_3
  : public Parameterizer_traits_3<TriangleMesh>
{
// Private types
private:
  // Superclass
  typedef Parameterizer_traits_3<TriangleMesh>        Base;

// Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef typename Base::Error_code                   Error_code;
  /// @endcond

  /// Export BorderParameterizer_3 template parameter.
  typedef BorderParameterizer_3                       Border_param;

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator        face_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator      vertex_iterator;

  typedef CGAL::Vertex_around_target_circulator<TriangleMesh>    vertex_around_target_circulator;
  typedef CGAL::Halfedge_around_target_circulator<TriangleMesh>  halfedge_around_target_circulator;
  typedef CGAL::Halfedge_around_face_circulator<TriangleMesh>  halfedge_around_face_circulator;

  typedef boost::unordered_set<vertex_descriptor>       Vertex_set;
  typedef std::vector<face_descriptor>                  Faces_vector;

  // Mesh_Adaptor_3 subtypes:
  typedef Parameterizer_traits_3<TriangleMesh>  Traits;
  typedef typename Traits::NT                   NT;
  typedef typename Traits::Point_2              Point_2;
  typedef typename Traits::Point_3              Point_3;
  typedef typename Traits::Vector_2             Vector_2;
  typedef typename Traits::Vector_3             Vector_3;

  // SparseLinearAlgebraTraits_d subtypes:
  typedef SparseLinearAlgebraTraits_d                  Sparse_LA;

  typedef typename Sparse_LA::Vector                   Vector;
  typedef typename Sparse_LA::Matrix                   Matrix;

  // Memory maps
    // Each triangle is associated a linear transformation matrix
  typedef std::pair<NT, NT>                                         Lt_matrix;
  typedef CGAL::Unique_hash_map<face_descriptor, Lt_matrix>         Lt_hash_map;
  typedef boost::associative_property_map<Lt_hash_map>              Lt_map;

    // Each angle (uniquely determined by the opposite half edge) has a cotangent
  typedef CGAL::Unique_hash_map<halfedge_descriptor, NT>            Cot_hm;
  typedef boost::associative_property_map<Cot_hm>                   Cot_map;

    // Each face has a local 2D isometric parametrisation
  typedef std::pair<int, int>                                       Local_indices;
  typedef CGAL::Unique_hash_map<halfedge_descriptor, Local_indices> Lp_hm;
  typedef boost::associative_property_map<Lp_hm>                    Lp_map;
  typedef std::vector<Point_2>                                      Local_points;

  // Angle computations
  using Base::cotangent;

// Private fields
private:
  /// Object that maps (at least two) border vertices onto a 2D space
  Border_param m_borderParameterizer;

  /// Traits object to solve a sparse linear system
  Sparse_LA m_linearAlgebra;

  /// Controlling parameter
  const NT m_lambda;
  const NT m_lambda_tolerance;

  /// Energy minimization parameters
  const unsigned int m_iterations;
  const NT m_tolerance;

// Private accessors
private:
  /// Get the object that maps the surface's border onto a 2D space.
  Border_param& get_border_parameterizer() { return m_borderParameterizer; }

  /// Get the sparse linear algebra (traits object to access the linear system).
  Sparse_LA& get_linear_algebra_traits() { return m_linearAlgebra; }

// Private utilities
private:
  template <typename VertexUVmap>
  void output_uvmap(const std::string filename,
                    const TriangleMesh& mesh,
                    const Faces_vector& faces,
                    const VertexUVmap uvmap) const
  {
    std::ofstream out(filename.c_str());
    BOOST_FOREACH(face_descriptor fd, faces){
      halfedge_descriptor hd = halfedge(fd, mesh);
      out << "4 " << uvmap[target(hd, mesh)] << " 0 ";
      hd = next(hd, mesh);
      BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)){
        out << uvmap[vd] << " 0 ";
      }
      out << std::endl;
    }
  }

  template <typename VertexUVmap>
  void output_uvmap(const std::string filename,
                    const unsigned int iter,
                    const TriangleMesh& mesh,
                    const Faces_vector& faces,
                    const VertexUVmap uvmap) const
  {
    std::ostringstream outss;
    outss << filename << iter << ".txt" << std::ends;
    output_uvmap(outss.str(), mesh, faces, uvmap);
  }

// Private operations
private:
  void initialize_containers(const TriangleMesh& mesh,
                             halfedge_descriptor bhd,
                             Vertex_set& vertices,
                             Faces_vector& faces) const
  {
    internal::ContainersFiller<TriangleMesh> fc(mesh, faces, vertices);
    CGAL::Polygon_mesh_processing::connected_component(face(opposite(bhd, mesh), mesh),
                                                       mesh,
                                                       boost::make_function_output_iterator(fc));
    CGAL_postcondition(vertices.size() == num_vertices(mesh) &&
                       faces.size() == num_faces(mesh));
  }

  template <typename VertexUVmap>
  Error_code compute_initial_uv_map(TriangleMesh& mesh,
                                    halfedge_descriptor bhd,
                                    VertexUVmap uvmap) const
  {
    // Get the initial uv map from a LSCM run
    typedef CGAL::LSCM_parameterizer_3<TriangleMesh,
                                       Border_param>         LSCM_parameterizer;
//    typedef CGAL::Mean_value_coordinates_parameterizer_3<TriangleMesh,
//                                       Border_param>         MVC_parameterizer;

    Error_code status = CGAL::parameterize(mesh, LSCM_parameterizer(), bhd, uvmap);
    return status;
  }

  template <typename VertexUVmap, typename VertexParameterizedMap>
  Error_code parameterize_border(const TriangleMesh& mesh,
                                 const Vertex_set& vertices,
                                 halfedge_descriptor bhd,
                                 VertexUVmap uvmap,
                                 VertexParameterizedMap vpmap)
  {
    // Compute (u,v) for (at least two) border vertices and mark them as "parameterized"
    Error_code status = Base::OK;

#define FIXED_VERTICES_IN_PARAMETERIZATION_SPACE
#ifdef FIXED_VERTICES_IN_PARAMETERIZATION_SPACE
    // Find the two farthest vertices in the initial parameterization

    // @fixme brute force algorithm; use Convex hull + rotating caliphers instead
    NT max_dist = (std::numeric_limits<NT>::min)();
    vertex_descriptor vd1_max, vd2_max;

    BOOST_FOREACH(vertex_descriptor vd1, vertices){
      BOOST_FOREACH(vertex_descriptor vd2, vertices){
        Point_2 uv1 = get(uvmap, vd1);
        Point_2 uv2 = get(uvmap, vd2);

        NT sq_d = CGAL::squared_distance(uv1, uv2);
        if(sq_d > max_dist)
        {
          max_dist = sq_d;
          vd1_max = vd1;
          vd2_max = vd2;
        }
      }
    }

    // the value in uvmap is already set
    put(vpmap, vd1_max, true);
    put(vpmap, vd2_max, true);
#else
    // This fixes two vertices that are far in the original geometry. We lose
    // the UV information of the initial param
    status = get_border_parameterizer().parameterize_border(mesh, bhd, uvmap, vpmap);
#endif
    return status;
  }

  void compute_cotangent_angle(const TriangleMesh& mesh,
                               halfedge_descriptor hd,
                               vertex_descriptor vi,
                               vertex_descriptor vj,
                               vertex_descriptor vk,
                               Cot_map ctmap) const
  {
    typedef typename boost::property_map<TriangleMesh,
                                         boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(vertex_point, mesh);

    Point_3 position_vi = get(ppmap, vi);
    Point_3 position_vj = get(ppmap, vj);
    Point_3 position_vk = get(ppmap, vk);

    double cot = cotangent(position_vi, position_vj, position_vk);
    put(ctmap, hd, cot);
  }

  // Filling ctmap with the cotangent of the angles opposite of halfedges
  Error_code compute_cotangent_angles(const TriangleMesh& mesh,
                                      const Faces_vector& faces,
                                      Cot_map ctmap) const
  {
    BOOST_FOREACH(face_descriptor fd, faces){
      halfedge_descriptor hd = halfedge(fd, mesh), hdb = hd;

      vertex_descriptor vi = target(hd, mesh);
      hd = next(hd, mesh);
      vertex_descriptor vj = target(hd, mesh);
      hd = next(hd, mesh);
      vertex_descriptor vk = target(hd, mesh);
      hd = next(hd, mesh);

      if(hd != hdb){ // make sure that it is a triangular face
        return Base::ERROR_NON_TRIANGULAR_MESH;
      }

      compute_cotangent_angle(mesh, hd, vk, vj, vi , ctmap); // angle at v_j
      compute_cotangent_angle(mesh, next(hd, mesh), vi, vk, vj , ctmap); // angle at v_k
      compute_cotangent_angle(mesh, prev(hd, mesh), vj, vi, vk , ctmap); // angle at v_i
    }

    return Base::OK;
  }

  /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
  NT compute_w_ij(const TriangleMesh& mesh,
                  halfedge_descriptor hd,
                  const Cot_map ctmap) const
  {
    // Note that halfedge and vertex circulators move clockwise!!

    // coefficient corresponding to the angle at vk if vk is the vertex before vj
    // while circulating around vi
    NT c_k = get(ctmap, opposite(hd, mesh));

    // coefficient corresponding to the angle at vl if vl is the vertex after vj
    // while circulating around vi
    NT c_l = get(ctmap, hd);

    double weight = c_k + c_l;
    return weight;
  }

  /// Compute the line i of matrix A
  /// - call compute_w_ij() to compute the A coefficient w_ij for each neighbor v_j.
  /// - compute w_ii = - sum of w_ijs.
  ///
  /// \pre Vertices must be indexed.
  /// \pre Vertex i musn't be already parameterized.
  /// \pre Line i of A must contain only zeros.
  template <typename VertexIndexMap>
  Error_code fill_linear_system_matrix(Matrix& A,
                                       const TriangleMesh& mesh,
                                       vertex_descriptor vertex,
                                       const Cot_map ct_map,
                                       VertexIndexMap vimap) const
  {
    int i = get(vimap, vertex);

    // Circulate over the vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    halfedge_around_target_circulator hc(vertex, mesh), end = hc;
    CGAL_For_all(hc, end){
      halfedge_descriptor hd = *hc;
      CGAL_assertion(target(hd, mesh) == vertex);

      NT w_ij = -1.0 * compute_w_ij(mesh, hd, ct_map);

      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      vertex_descriptor v_j = source(hd, mesh);
      int j = get(vimap, v_j);

      // Set w_ij in matrix
      A.set_coef(i, j, w_ij, true /*new*/);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return Base::ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i, i, w_ii, true /*new*/);
    return Base::OK;
  }

  /// Initialize the (constant) matrix A in the linear system "A*X = B",
  /// after (at least two) border vertices parameterization.
  template <typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code initialize_matrix_A(const TriangleMesh& mesh,
                                 const Vertex_set& vertices,
                                 const Cot_map ctmap,
                                 VertexIndexMap vimap,
                                 VertexParameterizedMap vpmap,
                                 Matrix& A) const
  {
    Error_code status = Base::OK;

    // compute A
    unsigned int count = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices){
      if(!get(vpmap, vd)){ // not yet parameterized
        // Compute the line i of the matrix A
        status = fill_linear_system_matrix(A, mesh, vd, ctmap, vimap);
        if(status != Base::OK)
          return status;
      } else { // fixed vertices
        int index = get(vimap, vd);
        std::cout << index << " is fixed in A " << std::endl;
        A.set_coef(index, index, 1, true /*new*/);
        ++count;
      }
    }
    return status;
  }

  /// Solves the equation a3 x^3 + a2 x^2 + a1 x + a0 = 0.
  int solve_cubic_equation(const NT a3, const NT a2, const NT a1, const NT a0,
                           std::vector<NT>& roots) const
  {
    CGAL_precondition(roots.empty());
    NT r1 = 0, r2 = 0, r3 = 0; // roots of the cubic equation

    int root_n = CGAL::solve_cubic(a3, a2, a1, a0, r1, r2, r3);
    CGAL_postcondition(root_n > 0);

    roots.push_back(r1);
    if(root_n > 1)
      roots.push_back(r2);
    if(root_n > 2)
      roots.push_back(r3);

    std::cout << "coeffs: " << a3 << " " << a2 << " " << a1 << " " << a0 << std::endl;
    std::cout << root_n << " roots: " << r1 << " " << r2 << " " << r3 << std::endl;

    return roots.size();
  }

  /// Solves the equation a3 x^3 + a2 x^2 + a1 x + a0 = 0 using CGAL's algeabric kernel.
  int solve_cubic_equation_with_AK(const NT a3, const NT a2,
                                   const NT a1, const NT a0,
                                   std::vector<NT>& roots) const
  {
    CGAL_precondition(roots.empty());

    typedef CGAL::GMP_arithmetic_kernel                       AK;
    typedef CGAL::Algebraic_kernel_d_1<AK::Rational>          Algebraic_kernel_d_1;
    typedef typename Algebraic_kernel_d_1::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_d_1::Algebraic_real_1   Algebraic_real_1;
    typedef typename Algebraic_kernel_d_1::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_d_1::Coefficient        Coefficient;

    typedef CGAL::Polynomial_traits_d<Polynomial_1>           Polynomial_traits_1;

    typedef Algebraic_kernel_d_1::Solve_1                     Solve_1;

    Algebraic_kernel_d_1 ak_1;
    const Solve_1 solve_1 = ak_1.solve_1_object();
    typename Polynomial_traits_1::Construct_polynomial construct_polynomial_1;
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_x[1]
      = {std::make_pair(CGAL::Exponent_vector(1,0),Coefficient(1))};
    Polynomial_1 x=construct_polynomial_1(coeffs_x, coeffs_x+1);

    AK::Rational a3q(a3);
    AK::Rational a2q(a2);
    AK::Rational a1q(a1);
    AK::Rational a0q(a0);

    Polynomial_1 pol = a3q*CGAL::ipower(x,3) + a2q*CGAL::ipower(x,2)
                       + a1q*x + a0q;

    std::vector<std::pair<Algebraic_real_1, Multiplicity_type> > ak_roots;
    solve_1(pol, std::back_inserter(ak_roots));

    for(std::size_t i=0; i<ak_roots.size(); i++)
    {
      roots.push_back(ak_roots[i].first.to_double());
      std::cout << CGAL::to_double(ak_roots[i].first) << std::endl;
    }

    return roots.size();
  }

  /// Solve the bivariate system
  /// { C1 * a + 2 * lambda * a ( a^2 + b^2 - 1 ) = C2
  /// { C1 * b + 2 * lambda * b ( a^2 + b^2 - 1 ) = C3
  /// using CGAL's algeabric kernel.
  int solve_bivariate_system(const NT C1, const NT C2, const NT C3,
                             std::vector<NT>& a_roots, std::vector<NT>& b_roots) const
  {
    typedef CGAL::GMP_arithmetic_kernel                       AK;
    typedef CGAL::Algebraic_kernel_d_2<AK::Rational>          Algebraic_kernel_d_2;
    typedef typename Algebraic_kernel_d_2::Polynomial_2       Polynomial_2;
    typedef typename Algebraic_kernel_d_2::Algebraic_real_2   Algebraic_real_2;
    typedef typename Algebraic_kernel_d_2::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_d_2::Coefficient        Coefficient;

    typedef CGAL::Polynomial_traits_d<Polynomial_2>           Polynomial_traits_2;

    typedef Algebraic_kernel_d_2::Solve_2             Solve_2;
    typedef Algebraic_kernel_d_2::Is_coprime_2        Is_coprime_2;
    typedef Algebraic_kernel_d_2::Make_coprime_2      Make_coprime_2;
    typedef Algebraic_kernel_d_2::Is_square_free_2    Is_square_free_2;

    Algebraic_kernel_d_2 ak_2;
    const Is_coprime_2 is_coprime_2 = ak_2.is_coprime_2_object();
    const Solve_2 solve_2 = ak_2.solve_2_object();
    const Is_square_free_2 is_square_free_2 = ak_2.is_square_free_2_object();
    const Make_coprime_2 make_coprime_2 = ak_2.make_coprime_2_object();

    typename Polynomial_traits_2::Construct_polynomial construct_polynomial_2;
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_x[1]
      = {std::make_pair(CGAL::Exponent_vector(1,0),Coefficient(1))};
    Polynomial_2 x=construct_polynomial_2(coeffs_x,coeffs_x+1);
    std::pair<CGAL::Exponent_vector, Coefficient> coeffs_y[1]
      = {std::make_pair(CGAL::Exponent_vector(0,1),Coefficient(1))};
    Polynomial_2 y=construct_polynomial_2(coeffs_y,coeffs_y+1);

    AK::Rational C1q(C1);
    AK::Rational C2q(C2);
    AK::Rational C3q(C3);

    Polynomial_2 pol1 = C1q * x + 2 * m_lambda * x * ( x*x + y*y - 1) - C2q;
    Polynomial_2 pol2 = C1q * y + 2 * m_lambda * y * ( x*x + y*y - 1) - C3q;

    std::cout << "init" << std::endl;
    std::cout << "pol1: " << pol1 << std::endl << "pol2: " << pol2 << std::endl;

    CGAL_precondition(is_square_free_2(pol1));
    CGAL_precondition(is_square_free_2(pol2));
    if(!is_coprime_2(pol1, pol2)){
      std::cout << "not coprime" << std::endl;
      CGAL_assertion(false); // @todo handle that case

      Polynomial_2 g, q1, q2;
      make_coprime_2(pol1, pol2, g, q1, q2);
      std::cout << "g: " << g << std::endl;
      pol1 = q1;
      pol2 = q2;
    }

    std::vector<std::pair<Algebraic_real_2, Multiplicity_type> > roots;
    solve_2(pol1, pol2, std::back_inserter(roots));
    std::cout << "solved with " << roots.size() << " roots." << std::endl;

    for(std::size_t i=0; i<roots.size(); i++)
    {
      a_roots.push_back(roots[i].first.to_double().first);
      b_roots.push_back(roots[i].first.to_double().second);
      std::cout << CGAL::to_double(roots[i].first.to_double().first) << std::endl;
      std::cout << CGAL::to_double(roots[i].first.to_double().second) << std::endl;
    }

    return a_roots.size();
  }

  /// Compute the root that gives the lowest face energy.
  template <typename VertexUVmap>
  std::size_t compute_root_with_lowest_energy(const TriangleMesh& mesh,
                                              face_descriptor fd,
                                              const Cot_map ctmap,
                                              const Local_points& lp,
                                              const Lp_map lpmap,
                                              const VertexUVmap uvmap,
                                              const NT C2_denom, const NT C3,
                                              const std::vector<NT>& roots) const
  {
    CGAL_precondition(!roots.empty());
    NT E_min = (std::numeric_limits<NT>::max)();
    std::size_t index_arg = 0.;
    for(std::size_t i=0; i<roots.size(); ++i)
    {
      const NT a = roots[i];
      const NT b = C3 * C2_denom * a;
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap, a, b);
      if(Ef < E_min){
        E_min = Ef;
        index_arg = i;
      }
    }
    return index_arg;
  }

  /// Compute the root that gives the lowest face energy.
  template <typename VertexUVmap>
  std::size_t compute_root_with_lowest_energy(const TriangleMesh& mesh,
                                              face_descriptor fd,
                                              const Cot_map ctmap,
                                              const Local_points& lp,
                                              const Lp_map lpmap,
                                              const VertexUVmap uvmap,
                                              const std::vector<NT>& a_roots,
                                              const std::vector<NT>& b_roots) const
  {
    CGAL_precondition(!a_roots.empty() && a_roots.size() == b_roots.size());
    NT E_min = (std::numeric_limits<NT>::max)();
    std::size_t index_arg = -1;
    for(std::size_t i=0; i<a_roots.size(); ++i)
    {
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap,
                                          a_roots[i], b_roots[i]);
      if(Ef < E_min){
        E_min = Ef;
        index_arg = i;
      }
    }
    return index_arg;
  }

  /// Compute the optimal values of the linear transformation matrices Lt.
  template <typename VertexUVmap>
  Error_code compute_optimal_Lt_matrices(const TriangleMesh& mesh,
                                         const Faces_vector& faces,
                                         const Cot_map ctmap,
                                         const Local_points& lp,
                                         const Lp_map lpmap,
                                         const VertexUVmap uvmap,
                                         Lt_map ltmap) const
  {
    Error_code status = Base::OK;

    BOOST_FOREACH(face_descriptor fd, faces){
      // Compute the coefficients C1, C2, C3
      NT C1 = 0., C2 = 0., C3 = 0.;

      halfedge_around_face_circulator hc(halfedge(fd, mesh), mesh), end(hc);
      CGAL_For_all(hc, end){
        halfedge_descriptor hd = *hc;
        NT c = get(ctmap, hd);

        // UV positions
        Point_2 uvpi = get(uvmap, source(hd, mesh));
        Point_2 uvpj = get(uvmap, target(hd, mesh));
        NT diff_x = uvpi.x() - uvpj.x();
        NT diff_y = uvpi.y() - uvpj.y();

        // local positions (in the isometric 2D param)
        Local_indices li = get(lpmap, hd);
        Point_2 ppi = lp[ li.first ];
        Point_2 ppj = lp[ li.second ];
        NT p_diff_x = ppi.x() - ppj.x();
        NT p_diff_y = ppi.y() - ppj.y();

        std::cout << "c: " << c << std::endl;
        std::cout << "diff: " << diff_x << " " << diff_y << std::endl;
        std::cout << "pdiff: " << p_diff_x << " " << p_diff_y << std::endl;

        std::cout << "ADD1: " << c * ( p_diff_x*p_diff_x + p_diff_y*p_diff_y ) << std::endl;
        std::cout << "ADD2: " << c * ( diff_x*p_diff_x + diff_y*p_diff_y ) << std::endl;
        std::cout << "ADD3: " << c * ( diff_x*p_diff_y - diff_y*p_diff_x ) << std::endl;

        C1 += c * ( p_diff_x*p_diff_x + p_diff_y*p_diff_y );
        C2 += c * ( diff_x*p_diff_x + diff_y*p_diff_y );
        C3 += c * ( diff_x*p_diff_y - diff_y*p_diff_x );
      }

      std::cout << "C1: " << C1 << " , C2: " <<  C2 << " , C3: " << C3 << std::endl;

      // Compute a and b
      NT a = 0., b = 0.;

      if(m_lambda == 0.){ // ASAP
        CGAL_precondition(C1 != 0.);
        a = C2 / C1;
        b = C3 / C1;
      }
      else if( std::abs(C1) < m_lambda_tolerance * m_lambda &&
               std::abs(C2) < m_lambda_tolerance * m_lambda ){ // ARAP
        // If lambda is large compared to C1 and C2, the cubic equation that
        // determines a and b can be simplified to a simple quadric equation

        CGAL_precondition(C2*C2 + C3*C3 != 0.);
        NT denom = 1. / CGAL::sqrt(C2*C2 + C3*C3);
        a = C2 * denom;
        b = C3 * denom;
      }
      else { // general case
#define SOLVE_CUBIC_EQUATION
#ifdef SOLVE_CUBIC_EQUATION
        CGAL_precondition(C2 != 0.);
        NT C2_denom = 1. / C2;
        NT a3_coeff = 2. * m_lambda * (C2 * C2 + C3 * C3) * C2_denom * C2_denom;

        std::vector<NT> roots;
        solve_cubic_equation(a3_coeff, 0., (C1 - 2. * m_lambda), -C2, roots);

        // The function above is correct up to 10^{-14}, but below can be used
        // if more precision is needed (should never be the case)
//        std::vector<NT> roots_with_AK;
//        solve_cubic_equation_with_AK(a3_coeff, 0., (C1 - 2. * m_lambda), -C2, roots_with_AK);

        std::size_t ind = compute_root_with_lowest_energy(mesh, fd,
                                                          ctmap, lp, lpmap, uvmap,
                                                          C2_denom, C3, roots);
        a = roots[ind];
        b = C3 * C2_denom * a;
#else // solve the bivariate system
        std::vector<NT> a_roots;
        std::vector<NT> b_roots;
        solve_bivariate_system(C1, C2, C3, a_roots, b_roots);

        std::size_t ind = compute_root_with_lowest_energy(mesh, fd,
                                                          ctmap, lp, lpmap, uvmap,
                                                          a_roots, b_roots);
        a = a_roots[ind];
        b = b_roots[ind];
#endif

        std::cout << "ab: " << a << " " << b << std::endl;
      }

      // Update the map faces --> optimal Lt matrices
      Lt_matrix ltm = std::make_pair(a, b);
      put(ltmap, fd, ltm);
    }

    return status;
  }

  /// Computes the coordinates of the vertices p0, p1, p2
  /// in a local 2D orthonormal basis of the triangle's plane.
  void project_triangle(const Point_3& p0, const Point_3& p1, const Point_3& p2,  // in
                        Point_2& z0, Point_2& z1, Point_2& z2) const              // out
  {
    Vector_3 X = p1 - p0;
    NT X_norm = std::sqrt(X * X);
    if(X_norm != 0.0)
      X = X / X_norm;

    Vector_3 Z = CGAL::cross_product(X, p2 - p0);
    NT Z_norm = std::sqrt(Z * Z);
    if(Z_norm != 0.0)
      Z = Z / Z_norm;

    Vector_3 Y = CGAL::cross_product(Z, X);

    NT x0 = 0;
    NT y0 = 0;
    NT x1 = std::sqrt( (p1 - p0)*(p1 - p0) );
    NT y1 = 0;
    NT x2 = (p2 - p0) * X;
    NT y2 = (p2 - p0) * Y;

    z0 = Point_2(x0, y0);
    z1 = Point_2(x1, y1);
    z2 = Point_2(x2, y2);
  }

  /// Compute the local parameterization (2D) of a face and store it in memory.
  void project_face(const TriangleMesh& mesh,
                    vertex_descriptor vi,
                    vertex_descriptor vj,
                    vertex_descriptor vk,
                    Local_points& lp) const
  {
    typedef typename boost::property_map<TriangleMesh,
                                         boost::vertex_point_t>::const_type PPmap;
    PPmap ppmap = get(vertex_point, mesh);

    Point_3 position_vi = get(ppmap, vi);
    Point_3 position_vj = get(ppmap, vj);
    Point_3 position_vk = get(ppmap, vk);

    Point_2 pvi, pvj, pvk;
    project_triangle(position_vi, position_vj, position_vk,
                     pvi, pvj, pvk);

    // Add the newly computed 2D points to the vector of local positions
    lp.push_back(pvi);
    lp.push_back(pvj);
    lp.push_back(pvk);
  }

  /// Utility for fill_linear_system_rhs():
  /// Compute the local isometric parametrization (2D) of the faces of the mesh.
  Error_code compute_local_parametrisation(const TriangleMesh& mesh,
                                           const Faces_vector& faces,
                                           Local_points& lp,
                                           Lp_map lpmap) const
  {
    int global_index = 0;

    BOOST_FOREACH(face_descriptor fd, faces){
      halfedge_descriptor hd = halfedge(fd, mesh), hdb = hd;

      vertex_descriptor vi = target(hd, mesh); // hd is k -- > i
      put(lpmap, hd, std::make_pair(global_index + 2, global_index));

      hd = next(hd, mesh);
      vertex_descriptor vj = target(hd, mesh); // hd is i -- > j
      put(lpmap, hd, std::make_pair(global_index, global_index + 1));

      hd = next(hd, mesh);
      vertex_descriptor vk = target(hd, mesh); // hd is j -- > k
      put(lpmap, hd, std::make_pair(global_index + 1, global_index + 2));

      hd = next(hd, mesh);
      if(hd != hdb){ // to make sure that it is a triangular face
        return Base::ERROR_NON_TRIANGULAR_MESH;
      }

      project_face(mesh, vi, vj, vk, lp);
      global_index += 3;
    }

    if(lp.size() != 3 * faces.size())
      return Base::ERROR_NON_TRIANGULAR_MESH;

    return Base::OK;
  }

  /// Compute the coefficient b_ij = (i, j) of the matrix B, for j neighbor vertex of i.
  void compute_b_ij(const TriangleMesh& mesh,
                    halfedge_descriptor hd,
                    const Cot_map ctmap,
                    const Local_points& lp,
                    const Lp_map lpmap,
                    const Lt_map ltmap,
                    NT& x, NT& y) const // x for Bu and y for Bv
  {
    CGAL_precondition(x == 0.0 && y == 0.0);

    // Note that :
    // -- Circulators move in a clockwise manner
    // -- The halfedge hd points to vi

    // -- Face IJK with vk before vj while circulating around vi
    halfedge_descriptor hd_opp = opposite(hd, mesh); // hd_opp points to vj
    face_descriptor fd_k = face(hd_opp, mesh);

    // Get the matrix L_t corresponding to the face ijk
    Lt_matrix ltm_k = get(ltmap, fd_k);
    NT a_k = ltm_k.first;
    NT b_k = ltm_k.second;

    // Get the local parametrisation in the triangle ijk
    Local_indices loc_indices_k = get(lpmap, hd_opp);
    const Point_2& pvi_k = lp[ loc_indices_k.first ];
    const Point_2& pvj_k = lp[ loc_indices_k.second ];

    // x_i - x_j in the local parameterization of the triangle ijk
    NT diff_k_x = pvi_k.x() - pvj_k.x();
    NT diff_k_y = pvi_k.y() - pvj_k.y();

    // get the cotangent angle at vk
    NT ct_k = get(ctmap, hd_opp);

    // cot_k * Lt_k * (xi - xj)_k
    x = ct_k * (  a_k * diff_k_x + b_k * diff_k_y );
    y = ct_k * ( -b_k * diff_k_x + a_k * diff_k_y );

    // -- Face IJL with vl after vj while circulating around vi
    face_descriptor fd_l = face(hd, mesh);

    // Get the matrix L_t corresponding to the face
    Lt_matrix ltm_l = get(ltmap, fd_l);
    NT a_l = ltm_l.first;
    NT b_l = ltm_l.second;

    // Get the local parametrisation in the triangle ijl
    Local_indices loc_indices_l = get(lpmap, hd);
    const Point_2& pvi_l = lp[ loc_indices_l.second ]; // because hd points to vi
    const Point_2& pvj_l = lp[ loc_indices_l.first ];

    // x_i - x_j in the local parameterization of the triangle ijl
    NT diff_l_x = pvi_l.x() - pvj_l.x();
    NT diff_l_y = pvi_l.y() - pvj_l.y();

    // get the cotangent angle at vl
    NT ct_l = get(ctmap, hd);

    // cot_l * Lt_l * (xi - xj)_l
    x += ct_l * (  a_l * diff_l_x + b_l * diff_l_y );
    y += ct_l * ( -b_l * diff_l_x + a_l * diff_l_y );

    std::cout << "akbk: " << a_k << " " << b_k << std::endl;
    std::cout << "albl: " << a_l << " " << b_l << std::endl;
    std::cout << "xy: " << x << " " << y << std::endl;
  }

  /// Compute the line i of vectors Bu and Bv
  /// - call compute_b_ij() for each neighbor v_j to compute the B coefficient b_i
  ///
  /// \pre Vertices must be indexed.
  /// \pre Vertex i musn't be already parameterized.
  /// \pre Lines i of Bu and Bv must be zero.
  template <typename VertexIndexMap>
  Error_code fill_linear_system_rhs(const TriangleMesh& mesh,
                                    vertex_descriptor vertex,
                                    const Cot_map ctmap,
                                    const Local_points& lp,
                                    const Lp_map lpmap,
                                    const Lt_map ltmap,
                                    VertexIndexMap vimap,
                                    Vector& Bu, Vector& Bv) const
  {
    int i = get(vimap, vertex);

    // Circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT bu_i = 0, bv_i = 0;
    int vertexIndex = 0;

    halfedge_around_target_circulator hc(vertex, mesh), end = hc;
    CGAL_For_all(hc, end){
      halfedge_descriptor hd = *hc;
      CGAL_assertion(target(hd, mesh) == vertex);

      NT x = 0., y = 0.;
      compute_b_ij(mesh, hd, ctmap, lp, lpmap, ltmap, x, y);

      bu_i += x;
      bv_i += y;
      vertexIndex++;
    }

    // Set the entries for the right hand side vectors
    Bu.set(i, bu_i);
    Bv.set(i, bv_i);

    if (vertexIndex < 2)
      return Base::ERROR_NON_TRIANGULAR_MESH;

    return Base::OK;
  }

  /// Compute the entries of the right hand side of the linear system.
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code compute_rhs(const TriangleMesh& mesh,
                         const Vertex_set& vertices,
                         const Cot_map ctmap,
                         const Local_points& lp,
                         const Lp_map lpmap,
                         const Lt_map ltmap,
                         VertexUVmap uvmap,
                         VertexIndexMap vimap,
                         VertexParameterizedMap vpmap,
                         Vector& Bu, Vector& Bv) const
  {
    // Initialize the right hand side B in the linear system "A*X = B"
    Error_code status = Base::OK;

    unsigned int count = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices){
      if(!get(vpmap, vd)){ // not yet parameterized
        // Compute the lines i of the vectors Bu and Bv
        status = fill_linear_system_rhs(mesh, vd, ctmap, lp, lpmap,
                                                  ltmap, vimap, Bu, Bv);
        if(status != Base::OK)
          return status;
      } else { // fixed vertices
        int index = get(vimap, vd);
        Point_2 uv = get(uvmap, vd);
        Bu.set(index, uv.x());
        Bv.set(index, uv.y());
        ++count;
        std::cout << index << " is fixed in B: ";
        std::cout << uv.x() << " " << uv.y() << std::endl;
      }
    }
    return status;
  }

  /// Copy the data from two vectors to the UVmap.
  template <typename VertexUVmap, typename VertexIndexMap>
  void assign_solution(const Vector& Xu, const Vector& Xv,
                       const Vertex_set& vertices,
                       VertexUVmap uvmap,
                       VertexIndexMap vimap) const
  {
    BOOST_FOREACH(vertex_descriptor vd, vertices){
      int index = get(vimap, vd);
      NT u = Xu(index);
      NT v = Xv(index);
      put(uvmap, vd, Point_2(u, v));
    }
  }

  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code update_solution(const TriangleMesh& mesh,
                             const Vertex_set& vertices,
                             const Cot_map ctmap,
                             const Local_points& lp,
                             const Lp_map lpmap,
                             const Lt_map ltmap,
                             VertexUVmap uvmap,
                             VertexIndexMap vimap,
                             VertexParameterizedMap vpmap,
                             const Matrix& A)
  {
    Error_code status = Base::OK;

    // Create two sparse linear systems "A*Xu = Bu" and "A*Xv = Bv"
    int nbVertices = num_vertices(mesh);
    Vector Xu(nbVertices), Xv(nbVertices), Bu(nbVertices), Bv(nbVertices);

    // Fill the right hand side vectors
    status = compute_rhs(mesh, vertices, ctmap, lp, lpmap, ltmap,
                                         uvmap, vimap, vpmap, Bu, Bv);
    if (status != Base::OK)
      return status;

    std::cout << "A, B: " << std::endl << A.eigen_object() << std::endl << Bu << std::endl << Bv << std::endl;

    // Solve "A*Xu = Bu". On success, solution is (1/Du) * Xu.
    // Solve "A*Xv = Bv". On success, solution is (1/Dv) * Xv.
    NT Du, Dv;
    if(!get_linear_algebra_traits().linear_solver(A, Bu, Xu, Du) ||
       !get_linear_algebra_traits().linear_solver(A, Bv, Xv, Dv))
    {
      std::cout << "Could not solve linear system" << std::endl;
      status = Base::ERROR_CANNOT_SOLVE_LINEAR_SYSTEM;
      return status;
    }

    // WARNING: this package does not support homogeneous coordinates!
    CGAL_assertion(Du == 1.0);
    CGAL_assertion(Dv == 1.0);

    CGAL_postcondition_code
    (
      // make sure that the constrained vertices have not been moved
      BOOST_FOREACH(vertex_descriptor vd, vertices){
        if(get(vpmap, vd)){
          int index = get(vimap, vd);
          std::cout << "at: " << index << " " << Xu[index] << " " << Xv[index] << std::endl;
          CGAL_postcondition(std::abs(Xu[index] - Bu[index] ) < 1e-10);
          CGAL_postcondition(std::abs(Xv[index] - Bv[index] ) < 1e-10);
        }
      }
    )

    assign_solution(Xu, Xv, vertices, uvmap, vimap);
    return status;
  }


  /// Compute the current energy of a face, given a linear transformation matrix.
  template <typename VertexUVmap>
  NT compute_current_face_energy(const TriangleMesh& mesh,
                                 face_descriptor fd,
                                 const Cot_map ctmap,
                                 const Local_points& lp,
                                 const Lp_map lpmap,
                                 const VertexUVmap uvmap,
                                 const NT a, const NT b) const
  {
    NT Ef = 0.;

    halfedge_around_face_circulator hc(halfedge(fd, mesh), mesh), end(hc);
    CGAL_For_all(hc, end){
      halfedge_descriptor hd = *hc;
      NT cot = get(ctmap, hd);
      NT nabla_x = 0., nabla_y = 0.;

      // UV positions
      Point_2 pi = get(uvmap, source(hd, mesh));
      Point_2 pj = get(uvmap, target(hd, mesh));
      NT diff_x = pi.x() - pj.x();
      NT diff_y = pi.y() - pj.y();

      // local positions (in the 2D param)
      Local_indices li = get(lpmap, hd);
      Point_2 ppi = lp[ li.first ];
      Point_2 ppj = lp[ li.second ];
      NT p_diff_x = ppi.x() - ppj.x();
      NT p_diff_y = ppi.y() - ppj.y();

      nabla_x = diff_x - (  a * p_diff_x + b * p_diff_y );
      nabla_y = diff_y - ( -b * p_diff_x + a * p_diff_y );

      NT sq_nabla_norm = nabla_x * nabla_x + nabla_y * nabla_y;
      Ef += cot * sq_nabla_norm;
    }

    NT s = a * a + b * b - 1;
    Ef += m_lambda * s * s;
    return Ef;
  }

  /// Compute the current energy of a face.
  template <typename VertexUVmap>
  NT compute_current_face_energy(const TriangleMesh& mesh,
                                 face_descriptor fd,
                                 const Cot_map ctmap,
                                 const Local_points& lp,
                                 const Lp_map lpmap,
                                 const Lt_map ltmap,
                                 const VertexUVmap uvmap) const
  {
    NT Ef = 0.;

    Lt_matrix ltm = get(ltmap, fd); // the (current) optimal linear transformation
    NT a = ltm.first;
    NT b = ltm.second;

    return compute_current_face_energy(mesh, fd, ctmap, lp, lpmap, uvmap, a, b);
  }

  /// Compute the current energy of the parameterization.
  template <typename VertexUVmap>
  NT compute_current_energy(const TriangleMesh& mesh,
                            const Faces_vector& faces,
                            const Cot_map ctmap,
                            const Local_points& lp,
                            const Lp_map lpmap,
                            const Lt_map ltmap,
                            const VertexUVmap uvmap) const
  {
    NT E = 0.;

    BOOST_FOREACH(face_descriptor fd, faces){
      NT Ef = compute_current_face_energy(mesh, fd, ctmap, lp, lpmap,
                                          ltmap, uvmap);
      E += Ef;
    }

    E *= 0.5;
    return E;
  }

  // Post processing
  Error_code post_process() const
  {
    // convex virtual boundary algorithm of Karni et al.[2005]

    return Base::OK;
  }

// Public operations
public:
  template <typename VertexUVmap,
            typename VertexIndexMap,
            typename VertexParameterizedMap>
  Error_code parameterize(TriangleMesh& mesh,
                          halfedge_descriptor bhd,
                          VertexUVmap uvmap,
                          VertexIndexMap vimap,
                          VertexParameterizedMap vpmap)
  {
    Error_code status = Base::OK;

    int nbVertices = num_vertices(mesh);
    Matrix A(nbVertices, nbVertices); // the constant matrix using in the linear system A*X = B

    // vertices and faces containers
    Faces_vector faces;
    Vertex_set vertices;
    initialize_containers(mesh, bhd, vertices, faces);

    // linear transformation matrices L_t
    // Only need to store 2 indices since the complete matrix is {{a,b},{-b,a}}
    Lt_hash_map lt_hm;
    Lt_map ltmap(lt_hm); // will be filled in 'compute_optimal_Lt_matrices()'

    // Compute the initial parameterization of the mesh
    status = compute_initial_uv_map(mesh, bhd, uvmap);
    output_uvmap("ARAP_initial_param.txt", mesh, faces, uvmap);
    if(status != Base::OK)
      return status;

    // Fix two vertices on the border
    if(true || m_lambda == 0)
    {
      status = parameterize_border(mesh, vertices, bhd, uvmap, vpmap);
      if(status != Base::OK)
        return status;
    }

    // Compute all cotangent angles
    Cot_hm cthm;
    Cot_map ctmap(cthm);
    status = compute_cotangent_angles(mesh, faces, ctmap);
    if(status != Base::OK)
      return status;

    // Compute all local 2D parametrisation
    Lp_hm lphm;
    Lp_map lpmap(lphm);
    Local_points lp;
    status = compute_local_parametrisation(mesh, faces, lp, lpmap);
    if(status != Base::OK)
      return status;

    // The matrix A is constant and can be initialized outside of the loop
    status = initialize_matrix_A(mesh, vertices, ctmap, vimap, vpmap, A);
    if(status != Base::OK)
      return status;

    NT energy_this = compute_current_energy(mesh, faces, ctmap, lp, lpmap,
                                            ltmap, uvmap);
    NT energy_last;
    std::cout << "Initial energy: " << energy_this << std::endl;

    // main loop
    for(unsigned int ite=1; ite<=m_iterations; ++ite)
    {
      compute_optimal_Lt_matrices(mesh, faces, ctmap, lp, lpmap, uvmap, ltmap);
      status = update_solution(mesh, vertices, ctmap, lp, lpmap, ltmap,
                                               uvmap, vimap, vpmap, A);

      // Output the current situation
      output_uvmap("ARAP_iteration_", ite, mesh, faces, uvmap);
      energy_last = energy_this;
      energy_this = compute_current_energy(mesh, faces, ctmap, lp, lpmap,
                                                        ltmap, uvmap);
      std::cout << "Energy at iteration " << ite << " : " << energy_this << std::endl;
      CGAL_warning(energy_this >= 0);

      if(status != Base::OK)
        return status;

      // energy based termination
      if(m_tolerance > 0.0 && ite <= m_iterations) // if tolerance <= 0 then don't compute energy
      {                                                 // also no need compute energy if this iteration is the last iteration
        double energy_diff = std::abs((energy_last - energy_this) / energy_this);
        if(energy_diff < m_tolerance){
          std::cout << "Minimization process over after: "
                    << ite + 1 << " iterations. "
                    << "Energy diff: " << energy_diff << std::endl;
          break;
        }
      }
    }

    // Use post processing to handle flipped elements
    status = post_process();

    return status;
  }

public:
  /// Constructor
  ARAP_parameterizer_3(Border_param border_param = Border_param(),
                       ///< Object that maps the surface's border to 2D space
                       Sparse_LA sparse_la = Sparse_LA(),
                       ///< Traits object to access a sparse linear system
                       NT lambda = 1,
                       unsigned int iterations = 10,
                       NT tolerance = 1e-4)
    :
      m_borderParameterizer(border_param),
      m_linearAlgebra(sparse_la),
      m_lambda(lambda),
      m_lambda_tolerance(1e-10),
      m_iterations(iterations),
      m_tolerance(tolerance)
  { }

  // Default copy constructor and operator=() are fine
};

} // namespace CGAL

#endif // CGAL_ARAP_PARAMETERIZER_3_H
