// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_INTERNAL_GARLANDHECKBERT_POLICIES_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_INTERNAL_GARLANDHECKBERT_POLICIES_BASE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/internal/Common.h>

#include <CGAL/tags.h>

#include <Eigen/Dense>

#include <optional>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template <typename Mat>
Mat combine_matrices(const Mat& a, const Mat& b)
{
  return a + b;
}

template <typename GeomTraits>
struct GarlandHeckbert_matrix_types
{
  typedef typename GeomTraits::FT                                               FT;

  typedef Eigen::Matrix<FT, 3, 3, Eigen::DontAlign>                             Mat_3;
  typedef Eigen::Matrix<FT, 3, 1>                                               Col_3;
  typedef Eigen::Matrix<FT, 4, 4, Eigen::DontAlign>                             Mat_4;
  typedef Eigen::Matrix<FT, 4, 1>                                               Col_4;
  typedef Eigen::Matrix<FT, 1, 4>                                               Row_4;
};

// Storage is initialized by the most-derived class (e.g. GarlandHeckbert_plane_policies)
template <typename QuadricCalculator, typename TriangleMesh, typename GeomTraits>
struct GarlandHeckbert_quadrics_storage
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

  typedef Mat_4                                                                Cost_matrix;
  typedef CGAL::dynamic_vertex_property_t<Cost_matrix>                         Cost_property;
  typedef typename boost::property_map<TriangleMesh, Cost_property>::type      Vertex_cost_map;

  typedef QuadricCalculator                                                    Quadric_calculator;

protected:
  Vertex_cost_map m_cost_matrices;
  Quadric_calculator m_quadric_calculator;

public:
  GarlandHeckbert_quadrics_storage() = delete;

  GarlandHeckbert_quadrics_storage(TriangleMesh& tmesh,
                                   const Quadric_calculator& quadric_calculator)
    : m_quadric_calculator(quadric_calculator)
  {
    m_cost_matrices = get(Cost_property(), tmesh);
  }
};

template <typename QuadricCalculator, typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_cost_and_placement
  : public GarlandHeckbert_quadrics_storage<QuadricCalculator, TriangleMesh, GeomTraits>
{
  typedef QuadricCalculator                                                    Quadric_calculator;
  typedef GarlandHeckbert_quadrics_storage<
            Quadric_calculator, TriangleMesh, GeomTraits>                      Base;

public:
  // Tells the main function of 'Edge_collapse' that these
  // policies must call "initialize" and "update" functions.
  typedef CGAL::Tag_true                                                       Update_tag;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

  typedef typename GeomTraits::FT                                              FT;
  typedef typename GeomTraits::Point_3                                         Point_3;
  typedef typename GeomTraits::Vector_3                                        Vector_3;

  typedef typename Base::Mat_4                                                 Mat_4;
  typedef typename Base::Col_4                                                 Col_4;

private:
  FT discontinuity_multiplier;

public:
  GarlandHeckbert_cost_and_placement(TriangleMesh& tmesh,
                                     const Quadric_calculator& quadric_calculator,
                                     const FT dm = FT(100))
    : Base(tmesh, quadric_calculator), discontinuity_multiplier(dm)
  { }

  decltype(auto) vcm() const { return this->m_cost_matrices; }
  const Quadric_calculator& quadric_calculator() const { return this->m_quadric_calculator; }

public:
  static Col_4 point_to_homogenous_column(const Point_3& point)
  {
    return Col_4 { point.x(), point.y(), point.z(), FT(1) };
  }

  Col_4 construct_optimum(const Mat_4& mat, const Col_4& p0, const Col_4& p1) const
  {
    return quadric_calculator().construct_optimal_point(mat, p0, p1);
  }

  static bool is_discontinuity_edge(const halfedge_descriptor h,
                                    const TriangleMesh& tmesh)
  {
    return is_border_edge(h, tmesh);
  }

public:
  template <typename VertexPointMap>
  Mat_4 construct_quadric(const halfedge_descriptor he,
                          const TriangleMesh& tmesh,
                          const VertexPointMap vpm,
                          const GeomTraits& gt) const
  {
    return quadric_calculator().construct_quadric_from_edge(he, tmesh, vpm, gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric(const face_descriptor f,
                          const TriangleMesh& tmesh,
                          const VertexPointMap vpm,
                          const GeomTraits& gt) const
  {
    return quadric_calculator().construct_quadric_from_face(f, tmesh, vpm, gt);
  }

public:
  // initialize all quadrics
  template <typename VertexPointMap>
  void initialize(const TriangleMesh& tmesh,
                  const VertexPointMap vpm,
                  const GeomTraits& gt) const
  {
    Mat_4 zero_mat = Mat_4::Zero();

    for(vertex_descriptor v : vertices(tmesh))
      put(vcm(), v, zero_mat);

    for(face_descriptor f : faces(tmesh))
    {
      if(f == boost::graph_traits<TriangleMesh>::null_face())
        continue;

      const halfedge_descriptor h = halfedge(f, tmesh);

      // construtct the (4 x 4) matrix representing the plane quadric
      const Mat_4 quadric = construct_quadric(f, tmesh, vpm, gt);

      for(halfedge_descriptor shd : halfedges_around_face(h, tmesh))
      {
        const vertex_descriptor vs = source(shd, tmesh);
        const vertex_descriptor vt = target(shd, tmesh);

        put(vcm(), vs, combine_matrices(get(vcm(), vs), quadric));

        if(!is_discontinuity_edge(shd, tmesh))
          continue;

        const Mat_4 discontinuous_quadric =
          discontinuity_multiplier * construct_quadric(shd, tmesh, vpm, gt);

        put(vcm(), vs, combine_matrices(get(vcm(), vs), discontinuous_quadric));
        put(vcm(), vt, combine_matrices(get(vcm(), vt), discontinuous_quadric));
      }
    }
  }

  template <typename Profile, typename VertexDescriptor>
  void update_after_collapse(const Profile& profile,
                             const VertexDescriptor new_v) const
  {
    put(vcm(), new_v, combine_matrices(get(vcm(), profile.v0()),
                                       get(vcm(), profile.v1())));
  }

public:
  // Cost
  template <typename Profile>
  std::optional<typename Profile::FT>
  operator()(const Profile& profile,
             const std::optional<typename Profile::Point>& placement) const
  {
    typedef std::optional<typename Profile::FT>                              Optional_FT;

    if(!placement)
      return std::optional<typename Profile::FT>();

    CGAL_precondition(!get(vcm(), profile.v0()).isZero(0));
    CGAL_precondition(!get(vcm(), profile.v1()).isZero(0));

    const Mat_4 combined_matrix = combine_matrices(get(vcm(), profile.v0()),
                                                   get(vcm(), profile.v1()));
    const Col_4 pt = point_to_homogenous_column(*placement);
    const Optional_FT cost = (pt.transpose() * combined_matrix * pt)(0, 0);

    return cost;
  }

public:
  // Placement
  template <typename Profile>
  std::optional<typename Profile::Point> operator()(const Profile& profile) const
  {
    CGAL_precondition(!get(vcm(), profile.v0()).isZero(0));
    CGAL_precondition(!get(vcm(), profile.v1()).isZero(0));

    // the combined matrix has already been computed in the evaluation of the cost...
    const Mat_4 combined_matrix = combine_matrices(get(vcm(), profile.v0()),
                                                   get(vcm(), profile.v1()));

    const Col_4 p0 = point_to_homogenous_column(profile.p0());
    const Col_4 p1 = point_to_homogenous_column(profile.p1());

    const Col_4 opt = construct_optimum(combined_matrix, p0, p1);

    std::optional<typename Profile::Point> pt = typename Profile::Point(opt(0) / opt(3),
                                                                          opt(1) / opt(3),
                                                                          opt(2) / opt(3));

    return pt;
  }
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_INTERNAL_GARLANDHECKBERT_POLICIES_BASE_H
