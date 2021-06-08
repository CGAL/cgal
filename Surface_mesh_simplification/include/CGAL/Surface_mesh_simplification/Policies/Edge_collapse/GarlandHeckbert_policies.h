#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>

#include <CGAL/tags.h>

#include <Eigen/Dense>

#include <boost/optional/optional.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {

// forward-declare template class
template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_policies;

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
class GarlandHeckbert_policies : 
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

    GarlandHeckbert_policies() { }
    GarlandHeckbert_policies(TriangleMesh tm) : tm_(tm) { }  

    /**
     * TODO implementations of all functions needed by the base implementation
     */
    Col_4 construct_optimal_point(const Mat_4& aQuadric, const Col_4& p0, const Col_4& p1) const 
    {
      Mat_4 X;
      X << aQuadric.block(0, 0, 3, 4), 0, 0, 0, 1;

      Col_4 opt_pt;

      if(X.determinant() == 0)
      {
        // not invertible
        const Col_4 p1mp0 = std::move(p1 - p0);
        const FT a = (p1mp0.transpose() * aQuadric * p1mp0)(0, 0);
        const FT b = 2 * (p0.transpose() * aQuadric * p1mp0)(0, 0);

        if(a == 0)
        {
          if(b < 0)
            opt_pt = p1;
          else if(b == 0)
            opt_pt = 0.5 * (p0 + p1);
          else
            opt_pt = p0;
        }
        else
        {
          FT ext_t = -b/(2*a);
          if(ext_t < 0 || ext_t > 1 || a < 0)
          {
            // one of endpoints
            FT p0_cost = (p0.transpose() * aQuadric * p0)(0, 0);
            FT p1_cost = (p1.transpose() * aQuadric * p1)(0, 0);
            if(p0_cost > p1_cost)
              opt_pt = p1;
            else
              opt_pt = p0;
          }
          else
          {
            // extremum of the parabola
            opt_pt = p0 + ext_t * (p1 - p0);
          }
        }
      }
      else // invertible
      {
        opt_pt = X.inverse().col(3); // == X.inverse() * (0 0 0 1)
      }
      return opt_pt;
    }

    Mat_4 construct_quadric_from_normal(const Vector_3& normal, const Point_3& point) const {
      // negative dot product between the normal and the position vector
      const auto d = - normal * Vector_3(ORIGIN, point);

      // row vector given by d appended to the normal
      const auto row = Eigen::Matrix<FT, 1, 4>(normal.x(), normal.y(), normal.z(), d);

      // outer product
      return row.transpose() * row;
    }

  private:
    TriangleMesh tm_;
};

} //namespace Surface_mesh_simplification
} //namespace CGAL

#endif
