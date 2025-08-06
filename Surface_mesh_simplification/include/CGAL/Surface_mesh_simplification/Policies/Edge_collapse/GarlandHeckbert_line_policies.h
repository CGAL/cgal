// Copyright (c) 2025  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Leo Valque

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_policy_base.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_functions.h>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {

template <typename TriangleMesh, typename GeomTraits>
class Line_quadric_calculator
{
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;
  typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Row_4             Row_4;

public:
  Line_quadric_calculator() { }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_vertex(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
                                      const TriangleMesh& tmesh,
                                      const VertexPointMap point_map,
                                      const GeomTraits& gt) const
  {
    return construct_line_quadric_from_vertex(v, tmesh, point_map, gt);
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return Mat_4::Zero();
  }

  template <typename VertexPointMap>
  Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                    const TriangleMesh& tmesh,
                                    const VertexPointMap point_map,
                                    const GeomTraits& gt) const
  {
    return Mat_4::Zero();
  }


  Col_4 construct_optimal_point(const Mat_4& quadric,
                                const Col_4& p0,
                                const Col_4& p1) const
  {
    return construct_optimal_point_singular<GeomTraits>(quadric, p0, p1);
  }
};

} // namespace internal

template<typename TriangleMesh, typename GeomTraits>
class GarlandHeckbert_line_policies
  : public internal::GarlandHeckbert_cost_and_placement<
             internal::Line_quadric_calculator<TriangleMesh, GeomTraits>, TriangleMesh, GeomTraits>
{
public:
  typedef internal::Line_quadric_calculator<TriangleMesh, GeomTraits>         Quadric_calculator;

private:
  typedef internal::GarlandHeckbert_cost_and_placement<
            Quadric_calculator, TriangleMesh, GeomTraits>                      Base;
  typedef GarlandHeckbert_line_policies<TriangleMesh, GeomTraits>             Self;

public:
  typedef Self                                                                 Get_cost;
  typedef Self                                                                 Get_placement;

  typedef typename GeomTraits::FT                                              FT;

public:
  GarlandHeckbert_line_policies(TriangleMesh& tmesh,
                                 const FT dm = FT(100))
    : Base(tmesh, Quadric_calculator(), dm)
  { }

public:
  const Get_cost& get_cost() const { return *this; }
  const Get_placement& get_placement() const { return *this; }

  using Base::operator();
};

} // namespace Surface_mesh_simplification
} // namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H


// ////////

// namespace CGAL {
// namespace Surface_mesh_simplification {
// namespace internal {

//   // Storage is initialized by the most-derived class (e.g. GarlandHeckbert_plane_policies)
// template <typename TriangleMesh, typename GeomTraits>
// struct GarlandHeckbert_quadrics_storage_custom
// {
//   typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Mat_4             Mat_4;
//   typedef typename GarlandHeckbert_matrix_types<GeomTraits>::Col_4             Col_4;

//   typedef Mat_4                                                                Cost_matrix;
//   typedef CGAL::dynamic_vertex_property_t<Cost_matrix>                         Cost_property;
//   typedef typename boost::property_map<TriangleMesh, Cost_property>::type      Vertex_cost_map;

// protected:
//   Vertex_cost_map m_cost_matrices;

// public:
//   GarlandHeckbert_quadrics_storage_custom() = delete;

//   GarlandHeckbert_quadrics_storage_custom(TriangleMesh& tmesh)
//   {
//     m_cost_matrices = get(Cost_property(), tmesh);
//   }
// };
// }

// template <typename TriangleMesh, typename GeomTraits>
// class GarlandHeckbert_line_policies
//   : public internal::GarlandHeckbert_quadrics_storage_custom<TriangleMesh, GeomTraits>
// {
//   typedef internal::GarlandHeckbert_quadrics_storage_custom<TriangleMesh, GeomTraits>                      Base;

// public:
//   // Tells the main function of 'Edge_collapse' that these
//   // policies must call "initialize" and "update" functions.
//   typedef CGAL::Tag_true                                                       Update_tag;

//   typedef GarlandHeckbert_line_policies<TriangleMesh, GeomTraits>          Self;

//   typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor        vertex_descriptor;
//   typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor      halfedge_descriptor;
//   typedef typename boost::graph_traits<TriangleMesh>::face_descriptor          face_descriptor;

//   typedef typename GeomTraits::FT                                              FT;
//   typedef typename GeomTraits::Point_3                                         Point_3;
//   typedef typename GeomTraits::Vector_3                                        Vector_3;

//   typedef typename Base::Mat_4                                                 Mat_4;
//   typedef typename Base::Col_4                                                 Col_4;

// private:
//   FT discontinuity_multiplier;
//   FT line_factor;

// public:
//   GarlandHeckbert_line_policies(TriangleMesh& tmesh, const FT lf=FT(0.01), const FT dm = FT(100))
//     : Base(tmesh), discontinuity_multiplier(dm), line_factor(lf)
//   { }

//   decltype(auto) vcm() const { return this->m_cost_matrices; }

// public:
//   static Col_4 point_to_homogenous_column(const Point_3& point)
//   {
//     return Col_4 { point.x(), point.y(), point.z(), FT(1) };
//   }

//   Col_4 construct_optimum(const Mat_4& mat, const Col_4& p0, const Col_4& p1) const
//   {
//     // return quadric_calculator().construct_optimal_point(mat, p0, p1);
//     return internal::construct_optimal_point_singular<GeomTraits>(mat, p0, p1);
//   }

//   template <typename VertexPointMap>
//   Mat_4 construct_quadric_from_edge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
//                                     const TriangleMesh& tmesh,
//                                     const VertexPointMap point_map,
//                                     const GeomTraits& gt) const
//   {
//     return internal::construct_classic_plane_quadric_from_edge(he, tmesh, point_map, gt);
//   }

//   template <typename VertexPointMap>
//   Mat_4 construct_quadric_from_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
//                                     const TriangleMesh& tmesh,
//                                     const VertexPointMap point_map,
//                                     const GeomTraits& gt) const
//   {
//     return internal::construct_classic_plane_quadric_from_face(f, tmesh, point_map, gt);
//   }

//   static bool is_discontinuity_edge(const halfedge_descriptor h,
//                                     const TriangleMesh& tmesh)
//   {
//     return is_border_edge(h, tmesh);
//   }

// public:
//   // initialize all quadrics
//   template <typename VertexPointMap>
//   void initialize(const TriangleMesh& tmesh,
//                   const VertexPointMap vpm,
//                   const GeomTraits& gt) const
//   {
//     Mat_4 zero_mat = Mat_4::Zero();

//     for(vertex_descriptor v : vertices(tmesh))
//       put(vcm(), v, line_factor*internal::construct_line_quadric_from_vertex(v, tmesh, vpm, gt));


//     for(face_descriptor f : faces(tmesh))
//     {
//       if(f == boost::graph_traits<TriangleMesh>::null_face())
//         continue;

//       const halfedge_descriptor h = halfedge(f, tmesh);

//       // construtct the (4 x 4) matrix representing the plane quadric
//       const Mat_4 quadric = construct_quadric_from_face(f, tmesh, vpm, gt);

//       for(halfedge_descriptor shd : halfedges_around_face(h, tmesh))
//       {
//         const vertex_descriptor vs = source(shd, tmesh);
//         const vertex_descriptor vt = target(shd, tmesh);

//         put(vcm(), vs, internal::combine_matrices(get(vcm(), vs), quadric));

//         if(!is_discontinuity_edge(shd, tmesh))
//           continue;

//         const Mat_4 discontinuous_quadric =
//           discontinuity_multiplier * construct_quadric_from_edge(shd, tmesh, vpm, gt);

//         put(vcm(), vs, internal::combine_matrices(get(vcm(), vs), discontinuous_quadric));
//         put(vcm(), vt, internal::combine_matrices(get(vcm(), vt), discontinuous_quadric));
//       }
//     }
//   }

//   template <typename Profile, typename VertexDescriptor>
//   void update_after_collapse(const Profile& profile,
//                              const VertexDescriptor new_v) const
//   {
//     put(vcm(), new_v, internal::combine_matrices(get(vcm(), profile.v0()),
//                                        get(vcm(), profile.v1())));
//   }

// public:
//   // Cost
//   template <typename Profile>
//   std::optional<typename Profile::FT>
//   operator()(const Profile& profile,
//              const std::optional<typename Profile::Point>& placement) const
//   {
//     typedef std::optional<typename Profile::FT>                              Optional_FT;

//     if(!placement)
//       return std::optional<typename Profile::FT>();

//     CGAL_precondition(!get(vcm(), profile.v0()).isZero(0));
//     CGAL_precondition(!get(vcm(), profile.v1()).isZero(0));

//     const Mat_4 combined_matrix = internal::combine_matrices(get(vcm(), profile.v0()),
//                                                    get(vcm(), profile.v1()));
//     const Col_4 pt = point_to_homogenous_column(*placement);
//     const Optional_FT cost = (pt.transpose() * combined_matrix * pt)(0, 0);

//     return cost;
//   }

// public:
//   typedef Self                                                                 Get_cost;
//   typedef Self                                                                 Get_placement;

//   // Placement
//   template <typename Profile>
//   std::optional<typename Profile::Point> operator()(const Profile& profile) const
//   {
//     CGAL_precondition(!get(vcm(), profile.v0()).isZero(0));
//     CGAL_precondition(!get(vcm(), profile.v1()).isZero(0));

//     // the combined matrix has already been computed in the evaluation of the cost...
//     const Mat_4 combined_matrix = internal::combine_matrices(get(vcm(), profile.v0()),
//                                                    get(vcm(), profile.v1()));

//     const Col_4 p0 = point_to_homogenous_column(profile.p0());
//     const Col_4 p1 = point_to_homogenous_column(profile.p1());

//     const Col_4 opt = construct_optimum(combined_matrix, p0, p1);

//     std::optional<typename Profile::Point> pt = typename Profile::Point(opt(0) / opt(3),
//                                                                           opt(1) / opt(3),
//                                                                           opt(2) / opt(3));

//     return pt;
//   }

//   const Get_cost& get_cost() const { return *this; }
//   const Get_placement& get_placement() const { return *this; }
// };

// // } // namespace internal
// } // namespace Surface_mesh_simplification
// } // namespace CGAL

// #endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_GARLANDHECKBERT_LINE_POLICIES_H
