#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_QUADRICS_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_INTERNAL_GARLAND_HECKBERT_QUADRICS_H

#include <CGAL/license/Surface_mesh_simplification.h>
#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <Eigen/Dense>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {
namespace quadrics 
{
  template<typename Traits>
  using Mat4 = Eigen::Matrix<typename Traits::FT, 4, 4, Eigen::DontAlign>;

  template<typename Traits>
  using Point_3 = typename Traits::Point_3;

  template<typename Traits>
  using Vector_3 = typename Traits::Vector_3;

  template<typename Traits>
  using FT = typename Traits::FT;
  
  template<typename Traits>
  Eigen::Matrix<FT<Traits>, 3, 1> vector_3_to_col_vec(const Vector_3<Traits> & v) {
    return Eigen::Matrix<FT<Traits>, 3, 1>(v.x(), v.y(), v.z());
  }
    
  template<typename Traits>
  Mat4<Traits> construct_probabilistic_plane_quadric(
    const Vector_3<Traits> & mean_normal,
    const Point_3<Traits> & mean_pos,
    FT<Traits> std_dev_normal, 
    FT<Traits> std_dev_pos,
    const Traits traits)
  {
    const auto squared_length = traits.compute_squared_length_3_object();
    const auto dot_product = traits.compute_scalar_product_3_object();
    const auto construct_vec_3 = traits.construct_vector_3_object();
    
    const auto sn2 = CGAL::square(std_dev_normal);
    const auto sp2 = CGAL::square(std_dev_pos);

    const auto mean_vec = construct_vec_3(CGAL::ORIGIN, mean_pos);
    const auto dot_mnmv = dot_product(mean_normal, mean_vec);
    
    // Eigen column vector of length 3
    const auto mean_n_col = vector_3_to_col_vec<Traits>(mean_normal);

    // start by setting values along the diagonal
    Mat4<Traits> mat = sn2 * Mat4<Traits>::Identity();
    
    // add outer product of the mean normal with itself
    // to the upper left 3x3 block
    mat.block(0, 0, 3, 3) += mean_n_col * mean_n_col.transpose();
    
    // set the first 3 values of the last row and the first 
    // 3 values of the last column
    //TODO why do we have to flip this sign as well? Probably linked to
    // our weird cross product order in the three point overloads, 
    // but users will run into problems using this overload
    const auto b = - vector_3_to_col_vec<Traits>(dot_mnmv * mean_normal + sn2 * mean_vec);
    mat.col(3).head(3) = b;
    mat.row(3).head(3) = b.transpose();
    
    // set the value in the bottom right corner
    mat(3, 3) = CGAL::square(dot_mnmv) 
      + sn2 * squared_length(mean_vec) 
      + sp2 * squared_length(mean_normal)
      + 3 * sn2 + sp2;

    return mat;
  }  
  
  template<typename Traits>
  Mat4<Traits> construct_probabilistic_plane_quadric(
      const Point_3<Traits> & p,
      const Point_3<Traits> & q,
      const Point_3<Traits> & r,
      FT<Traits> std_dev_norm,
      FT<Traits> std_dev_pos,
      const Traits traits)
  {
    const auto unit_norm = traits.construct_unit_normal_3_object();
    
    //TODO how does a different sequence of r, p, q affect output? 
    return construct_probabilistic_plane_quadric<Traits>(unit_norm(r, p, q), 
        r, std_dev_norm, std_dev_pos, traits);
  }   
  
  template<typename Traits>
  Mat4<Traits> construct_classical_plane_quadric(
      const Vector_3<Traits> & n,
      const Point_3<Traits> & p,
      const Traits traits) 
  {

    const auto dot_product = traits.compute_scalar_product_3_object();
    const auto construct_vec_3 = traits.construct_vector_3_object();

    //TODO why is a negative sign here?
    const auto d = - dot_product(n, construct_vec_3(CGAL::ORIGIN, p));

    const auto row = Eigen::Matrix<typename Traits::FT, 1, 4>(n.x(), n.y(), n.z(), d);
    return row.transpose() * row;
  }

  /**
   * construction of a classical quadric, taken over literally from 
   * a previous version of GarlandHeckbert_core.h, so can potentially be simplified
   * TODO chec optimizations / simplifications
   */
  template<typename Traits>
  Mat4<Traits> construct_classical_plane_quadric(
      const Point_3<Traits> & p,
      const Point_3<Traits> & q,
      const Point_3<Traits> & r,
      const Traits traits)
  {
    const auto unit_norm = traits.construct_unit_normal_3_object();
    
    //TODO how does a different sequence of r, p, q affect output? 
    return construct_classical_plane_quadric<Traits>(unit_norm(r, p, q), r, traits);
  }
};
} //namespace CGAL
} //namespace Surface_mesh_simplification
} //namespace internal

#endif
