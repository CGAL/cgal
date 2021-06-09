#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_PROBABILISTIC_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {

// forward-declare template class
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_policies;

template<typename TriangleMesh, typename GeomTraits>
using Cost_property = CGAL::dynamic_vertex_property_t<Eigen::Matrix<typename GeomTraits::FT, 4, 4,
  Eigen::DontAlign>>;

template<typename TriangleMesh, typename GeomTraits>
using Vertex_cost_map = typename boost::property_map<TriangleMesh, CGAL::dynamic_vertex_property_t<
  Eigen::Matrix<typename GeomTraits::FT, 4, 4, Eigen::DontAlign>>>::type;

template<typename TriangleMesh, typename GeomTraits>
using Cost_base = internal::GarlandHeckbert_cost_base<
    Vertex_cost_map<TriangleMesh, GeomTraits>, 
    GeomTraits, 
    GarlandHeckbert_policies<TriangleMesh, GeomTraits>>;

template<typename TriangleMesh, typename GeomTraits>
using Placement_base = internal::GarlandHeckbert_placement_base<
    Vertex_cost_map<TriangleMesh, GeomTraits>, 
    GeomTraits, 
    GarlandHeckbert_policies<TriangleMesh, GeomTraits>>;

  // derived class implements functions used in the base class that
  // takes the derived class as template argument - see "CRTP"
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_probabilistic_policies : 
  public Placement_base<TriangleMesh, GeomTraits>,
  public Cost_base<TriangleMesh, GeomTraits>
{

  public:

    typedef Cost_base<TriangleMesh, GeomTraits> Get_cost;
    typedef Placement_base<TriangleMesh, GeomTraits> Get_placement;

    // these using directives are needed to choose between the definitions of these types
    // in Cost_base and Placement_base (even though they are the same)
    // TODO alternatives - e.g. rename base class types so they don't clash
    using typename Cost_base<TriangleMesh, GeomTraits>::Mat_4;
    using typename Cost_base<TriangleMesh, GeomTraits>::Col_4;
    using typename Cost_base<TriangleMesh, GeomTraits>::Point_3;
    using typename Cost_base<TriangleMesh, GeomTraits>::Vector_3;

    typedef typename GeomTraits::FT FT;

    // TODO good default values
    GarlandHeckbert_probabilistic_policies(
        TriangleMesh& tmesh, 
        FT dm = FT(100),
        Vector_3 mv = Vector_3(0.01, 0.01, 0.01), 
        Vector_3 mn = Vector_3(0.01, 0.01, 0.01), 
        FT sdn = 0.01,
        FT sdp = 0.01) 
      : Get_cost(dm), mean_vec(mv), mean_normal(mn), std_dev_normal(sdn), std_dev_pos(sdp) 
    {
      Vertex_cost_map<TriangleMesh, GeomTraits> vcm_ = get(Cost_property<TriangleMesh, GeomTraits>(),
          tmesh);

      /**
       * initialize the two base class cost matrices (protected members)
       */
      Get_cost::m_cost_matrices = vcm_;
      Get_placement::m_cost_matrices = vcm_;
    }

    Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const 
    {
      Mat_4 X;
      X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

      Col_4 opt_pt;

      opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
      
      return opt_pt;

    }

    Mat_4 construct_quadric_from_normal(const Vector_3& normal, const Point_3& point,
        const GeomTraits& gt) const {
      const auto squared_length = gt.compute_squared_length_3_object();
      const auto dot_product = gt.compute_scalar_product_3_object();
      const auto construct_vec_3 = gt.construct_vector_3_object();

      const FT dot_mnmv = dot_product(mean_normal, mean_vec);

      // Eigen column vector of length 3
      const Col_4 mean_n_col = vector_3_to_col_vec(mean_normal);

      // intermediate values for simplicity
      const FT sdev_n_2 = square(std_dev_normal);
      const FT sdev_p_2 = square(std_dev_pos);
      
      // start by setting values along the diagonal
      Mat_4 mat = sdev_n_2 * Mat_4::Identity();

      // add outer product of the mean normal with itself
      // to the upper left 3x3 block
      mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();

      // set the first 3 values of the last row and the first 
      // 3 values of the last column
      //TODO why do we have to flip this sign as well? Probably linked to
      // our weird cross product order in the three point overloads, 
      const auto b = - vector_3_to_col_vec(dot_mnmv * mean_normal + sdev_n_2 * mean_vec);
      mat.col(3).head(3) = b;
      mat.row(3).head(3) = b.transpose();

      // set the value in the bottom right corner
      mat(3, 3) = CGAL::square(dot_mnmv) 
        + sdev_n_2 * squared_length(mean_vec) 
        + sdev_p_2 * squared_length(mean_normal)
        + 3 * sdev_n_2 + sdev_p_2;

      return mat;
    }

  private:
    Vector_3 mean_vec;
    Vector_3 mean_normal;
    FT std_dev_normal;
    FT std_dev_pos;

    Eigen::Matrix<FT, 3, 1> vector_3_to_col_vec(const Vector_3& v)
    {
      return Eigen::Matrix<FT, 3, 1>(v.x(), v.y(), v.z());
    }
};

} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif
