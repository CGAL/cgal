// Copyright (c) 2023 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//

#ifndef CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H
#define CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H

#include <CGAL/license/Straight_skeleton_extrusion_2.h>

#ifdef CGAL_SLS_DEBUG_DRAW
  #include <CGAL/draw_straight_skeleton_2.h>
  #include <CGAL/draw_polygon_2.h>
  #include <CGAL/draw_polygon_with_holes_2.h>
  #include <CGAL/draw_triangulation_2.h>
#endif

#include <CGAL/create_weighted_straight_skeleton_2.h>
#include <CGAL/create_weighted_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/create_weighted_offset_polygons_2.h>
#include <CGAL/create_weighted_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>

#include <CGAL/General_polygon_with_holes_2.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#ifdef CGAL_SLS_OUTPUT_FILES
#include <CGAL/IO/polygon_soup_io.h>
#endif

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>

#include <CGAL/Cartesian/function_objects.h>
#include <CGAL/enum.h>
#include <CGAL/IO/helpers.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_2.h>

#include <boost/algorithm/clamp.hpp>
#include <boost/range/value_type.hpp>

#include <algorithm>
#include <iostream>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <optional>
#include <memory>

namespace CGAL {
namespace Straight_skeleton_extrusion {
namespace internal {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Below is to handle vertical slabs.

enum class Slope
{
  UNKNOWN = 0,
  INWARD,
  OUTWARD,
  VERTICAL
};

template <typename FT>
inline constexpr FT default_extrusion_height()
{
  return (std::numeric_limits<double>::max)();
}

// @todo Maybe this postprocessing is not really necessary? Do users really care if the point
// is not perfectly above the input contour edge (it generally cannot be anyway if the kernel is not exact except for some
// specific cases)?
#define CGAL_SLS_SNAP_TO_VERTICAL_SLABS
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS

// The purpose of this visitor is to snap back almost-vertical (see preprocessing_weights()) edges
// to actual vertical slabs.
template <typename HDS, typename GeomTraits>
typename GeomTraits::Point_2
snap_point_to_contour_halfedge_plane(const typename GeomTraits::Point_2& op,
                                     typename HDS::Halfedge_const_handle ch)
{
  using FT = typename GeomTraits::FT;
  using Segment_2 = typename GeomTraits::Segment_2;
  using Line_2 = typename GeomTraits::Line_2;

  using Vertex_const_handle = typename HDS::Vertex_const_handle;

  const Vertex_const_handle sv = ch->opposite()->vertex();
  const Vertex_const_handle tv = ch->vertex();

  if(sv->point().x() == tv->point().x())
  {
    // vertical edge
    // std::cout << "vertical edge, snapping " << op << " to " << sv->point().x() << " "  << op.y() << std::endl;
    return { sv->point().x(), op.y() };
  }
  else if(sv->point().y() == tv->point().y())
  {
    // horizontal edge
    // std::cout << "horizontal edge, snapping " << op << " to " << op.x() << " " << sv->point().y() << std::endl;
    return { op.x(), sv->point().y() };
  }
  else
  {
    // Project orthogonally onto the halfedge
    // @todo should the projection be along the direction of the other offset edge sharing this point?
    Segment_2 s { sv->point(), tv->point() };
    std::optional<Line_2> line = CGAL_SS_i::compute_normalized_line_coeffC2(s);
    CGAL_assertion(bool(line)); // otherwise the skeleton would have failed already

    FT px, py;
    CGAL::line_project_pointC2(line->a(),line->b(),line->c(), op.x(),op.y(), px,py);
    // std::cout << "snapping " << op << " to " << px << " " << py << std::endl;
    return { px, py };
  }
}

template <typename HDS, typename GeomTraits>
void snap_skeleton_vertex(typename HDS::Halfedge_const_handle hds_h,
                          typename HDS::Halfedge_const_handle contour_h,
                          std::map<typename GeomTraits::Point_2,
                                   typename GeomTraits::Point_2>& snapped_positions)
{
  typename HDS::Vertex_const_handle hds_tv = hds_h->vertex();

  // this re-applies snapping towards contour_h even if the point was already snapped towards another contour
  auto insert_result = snapped_positions.emplace(hds_tv->point(), hds_tv->point());
  insert_result.first->second = snap_point_to_contour_halfedge_plane<HDS, GeomTraits>(insert_result.first->second, contour_h);

  // std::cout << "snap_skeleton_vertex(V" << hds_tv->id() << " pt: " << hds_h->vertex()->point() << ")"
  //           << " to " << insert_result.first->second << std::endl;
}

template <typename GeomTraits, typename PointRange>
void apply_snapping(PointRange& points,
                    const std::map<typename GeomTraits::Point_2,
                                   typename GeomTraits::Point_2>& snapped_positions)
{
  using Point_3 = typename GeomTraits::Point_3;

  for(Point_3& p3 : points)
  {
    auto it = snapped_positions.find({ p3.x(), p3.y() });
    if(it != snapped_positions.end())
      p3 = Point_3{it->second.x(), it->second.y(), p3.z()};
  }
}
#endif // CGAL_SLS_SNAP_TO_VERTICAL_SLABS

template <typename StraightSkeleton_2,
          typename OffsetBuilderTraits,
          typename GeomTraits>
class Skeleton_offset_correspondence_builder_visitor
  : public CGAL::Default_polygon_offset_builder_2_visitor<OffsetBuilderTraits, StraightSkeleton_2>
{
  using Base = CGAL::Default_polygon_offset_builder_2_visitor<OffsetBuilderTraits, StraightSkeleton_2>;

  using FT = typename GeomTraits::FT;
  using Point_2 = typename GeomTraits::Point_2;

  using SS_Halfedge_const_handle = typename StraightSkeleton_2::Halfedge_const_handle;

  using HDS = typename StraightSkeleton_2::Base;
  using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;

public:
  Skeleton_offset_correspondence_builder_visitor(const StraightSkeleton_2& ss,
                                                 std::unordered_map<HDS_Halfedge_const_handle, Point_2>& offset_points
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                                               , const FT vertical_weight
                                               , std::map<Point_2, Point_2>& snapped_positions
#endif
                                                )
    : m_ss(ss)
    , m_offset_points(offset_points)
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    , m_vertical_weight(vertical_weight)
    , m_snapped_positions(snapped_positions)
#endif
  { }

public:
  void on_offset_contour_started() const
  {
    // std::cout << "~~ new contour ~~" << std::endl;
  }

  // can't modify the position yet because we need arrange_polygons() to still work properly
  void on_offset_point(const Point_2& op,
                       SS_Halfedge_const_handle hook) const
  {
    CGAL_assertion(hook->is_bisector());

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    // @fixme on paper one could create a polygon thin-enough w.r.t. the max weight value such that
    // there is a skeleton vertex that wants to be snapped to two different sides...
    CGAL_assertion(m_snapped_positions.count(op) == 0);

    HDS_Halfedge_const_handle canonical_hook = (hook < hook->opposite()) ? hook : hook->opposite();

    SS_Halfedge_const_handle contour_h1 = hook->defining_contour_edge();
    CGAL_assertion(contour_h1->opposite()->is_border());
    SS_Halfedge_const_handle contour_h2 = hook->opposite()->defining_contour_edge();
    CGAL_assertion(contour_h2->opposite()->is_border());

    const bool is_h1_vertical = (contour_h1->weight() == m_vertical_weight);
    const bool is_h2_vertical = (contour_h2->weight() == m_vertical_weight);

    // this can happen when the offset is passing through vertices
    m_offset_points[canonical_hook] = op;

    // if both are vertical, it's the common vertex (which has to exist)
    if(is_h1_vertical && is_h2_vertical)
    {
      CGAL_assertion(contour_h1->vertex() == contour_h2->opposite()->vertex() ||
                     contour_h2->vertex() == contour_h1->opposite()->vertex());
      if(contour_h1->vertex() == contour_h2->opposite()->vertex())
        m_snapped_positions[op] = contour_h1->vertex()->point();
      else
        m_snapped_positions[op] = contour_h2->vertex()->point();
    }
    else if(is_h1_vertical)
    {
      m_snapped_positions[op] = snap_point_to_contour_halfedge_plane<HDS, GeomTraits>(op, contour_h1);
    }
    else if(is_h2_vertical)
    {
      m_snapped_positions[op] = snap_point_to_contour_halfedge_plane<HDS, GeomTraits>(op, contour_h2);
    }
#endif
  }

private:
  const StraightSkeleton_2& m_ss;
  std::unordered_map<HDS_Halfedge_const_handle, Point_2>& m_offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
  const FT m_vertical_weight;
  std::map<Point_2, Point_2>& m_snapped_positions;
#endif
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename GeomTraits>
class Extrusion_builder
{
  using Geom_traits = GeomTraits;

  using FT = typename Geom_traits::FT;
  using Point_2 = typename Geom_traits::Point_2;
  using Segment_2 = typename Geom_traits::Segment_2;
  using Line_2 = typename Geom_traits::Line_2;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;

  using Polygon_2 = CGAL::Polygon_2<Geom_traits>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<Geom_traits>;

  using Offset_polygons = std::vector<std::shared_ptr<Polygon_2> >;
  using Offset_polygons_with_holes = std::vector<std::shared_ptr<Polygon_with_holes_2> >;

  using Straight_skeleton_2 = CGAL::Straight_skeleton_2<Geom_traits>;
  using Straight_skeleton_2_ptr = std::shared_ptr<Straight_skeleton_2>;

  using SS_Vertex_const_handle = typename Straight_skeleton_2::Vertex_const_handle;
  using SS_Halfedge_const_handle = typename Straight_skeleton_2::Halfedge_const_handle;

  using HDS = typename Straight_skeleton_2::Base;
  using HDS_Vertex_const_handle = typename HDS::Vertex_const_handle;
  using HDS_Halfedge_const_handle = typename HDS::Halfedge_const_handle;
  using HDS_Face_handle = typename HDS::Face_handle;

  // Standard CDT2 for the horizontal (z constant) faces
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, Geom_traits>;
  using Vbb = CGAL::Triangulation_vertex_base_2<Geom_traits, Vb>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<Geom_traits>;
  using TDS = CGAL::Triangulation_data_structure_2<Vb,Fb>;
  using Itag = CGAL::No_constraint_intersection_requiring_constructions_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<Geom_traits, TDS, Itag>;
  using CDT_Vertex_handle = typename CDT::Vertex_handle;
  using CDT_Face_handle = typename CDT::Face_handle;

  // Projection CDT2 for the lateral faces
  using PK = CGAL::Projection_traits_3<Geom_traits>;
  using PVbb = CGAL::Triangulation_vertex_base_with_info_2<std::size_t, PK>;
  using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
  using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
  using PTDS = CGAL::Triangulation_data_structure_2<PVb,PFb>;
  using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
  using PCDT_Vertex_handle = typename PCDT::Vertex_handle;
  using PCDT_Face_handle = typename PCDT::Face_handle;

  using Offset_builder_traits = CGAL::Polygon_offset_builder_traits_2<Geom_traits>;
  using Visitor = Skeleton_offset_correspondence_builder_visitor<Straight_skeleton_2, Offset_builder_traits, Geom_traits>;
  using Offset_builder = CGAL::Polygon_offset_builder_2<Straight_skeleton_2,
                                                        Offset_builder_traits,
                                                        Polygon_2,
                                                        Visitor>;

private:
  Geom_traits m_gt;

public:
  Extrusion_builder(const Geom_traits& gt) : m_gt(gt) { }

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  template <typename PolygonWithHoles,
            typename PointRange, typename FaceRange>
  void construct_horizontal_faces(const PolygonWithHoles& p,
                                  const FT altitude,
                                  PointRange& points,
                                  FaceRange& faces,
                                  const bool invert_faces = false)
  {
#ifdef CGAL_SLS_DEBUG_DRAW
    CGAL::draw(p);
#endif

    CDT cdt;

    try
    {
      cdt.insert_constraint(p.outer_boundary().begin(), p.outer_boundary().end(), true /*close*/);
      for(auto h_it=p.holes_begin(); h_it!=p.holes_end(); ++h_it)
        cdt.insert_constraint(h_it->begin(), h_it->end(), true /*close*/);
    }
    catch(const typename CDT::Intersection_of_constraints_exception&)
    {
      std::cerr << "Warning: Failed to triangulate horizontal face" << std::endl;
      return;
    }

    std::size_t id = points.size(); // point ID offset (previous faces inserted their points)
    for(CDT_Vertex_handle vh : cdt.finite_vertex_handles())
    {
      points.emplace_back(cdt.point(vh).x(), cdt.point(vh).y(), altitude);
      vh->info() = id++;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // CGAL::draw(cdt);
#endif

    std::unordered_map<CDT_Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<CDT_Face_handle, bool> > in_domain(in_domain_map);

    CGAL::mark_domain_in_triangulation(cdt, in_domain);

    for(CDT_Face_handle f : cdt.finite_face_handles())
    {
      if(!get(in_domain, f))
        continue;

      // invert faces for the z=0 plane (bottom face)
      if(invert_faces)
        faces.push_back({f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info()});
      else
        faces.push_back({f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()});
    }
  }

  template <typename PointRange, typename FaceRange>
  void construct_horizontal_faces(const Offset_polygons_with_holes& p_ptrs,
                                  const FT altitude,
                                  PointRange& points,
                                  FaceRange& faces)
  {
    for(const auto& p_ptr : p_ptrs)
      construct_horizontal_faces(*p_ptr, altitude, points, faces);
  }

  template <typename SLSFacePoints, typename PointRange, typename FaceRange>
  void triangulate_skeleton_face(SLSFacePoints& face_points,
                                 const bool invert_faces,
                                 PointRange& points,
                                 FaceRange& faces)
  {
    CGAL_precondition(face_points.size() >= 3);

    // shift once to ensure that face_points[0] and face_points[1] are at z=0 and thus the normal is correct
    std::rotate(face_points.rbegin(), face_points.rbegin() + 1, face_points.rend());
    CGAL_assertion(face_points[0][2] == 0 && face_points[1][2] == 0);

    const Vector_3 n = CGAL::cross_product(face_points[1] - face_points[0], face_points[2] - face_points[0]);
    PK traits(n);
    PCDT pcdt(traits);

    try
    {
      pcdt.insert_constraint(face_points.begin(), face_points.end(), true /*close*/);
    }
    catch(const typename PCDT::Intersection_of_constraints_exception&)
    {
      std::cerr << "Warning: Failed to triangulate skeleton face" << std::endl;
      return;
    }

    std::size_t id = points.size(); // point ID offset (previous faces inserted their points);
    for(PCDT_Vertex_handle vh : pcdt.finite_vertex_handles())
    {
      points.push_back(pcdt.point(vh));
      vh->info() = id++;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // CGAL::draw(pcdt);
#endif

    std::unordered_map<PCDT_Face_handle, bool> in_domain_map;
    boost::associative_property_map< std::unordered_map<PCDT_Face_handle, bool> > in_domain(in_domain_map);

    CGAL::mark_domain_in_triangulation(pcdt, in_domain);

    for(PCDT_Face_handle f : pcdt.finite_face_handles())
    {
      if(!get(in_domain, f))
        continue;

      // invert faces for exterior skeletons
      if(invert_faces)
        faces.push_back({f->vertex(0)->info(), f->vertex(2)->info(), f->vertex(1)->info()});
      else
        faces.push_back({f->vertex(0)->info(), f->vertex(1)->info(), f->vertex(2)->info()});
    }
  }

  // This version is for default height, so just gather all the full faces
  void construct_lateral_faces(const Straight_skeleton_2& ss,
                               std::vector<Point_3>& points,
                               std::vector<std::vector<std::size_t> >& faces,
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                               const FT& vertical_weight,
                               std::map<Point_2, Point_2>& snapped_positions,
#endif
                               const bool ignore_frame_faces = false,
                               const bool invert_faces = false)
  {
    const HDS& hds = static_cast<const HDS&>(ss);

    std::size_t fc = 0;

    for(const HDS_Face_handle hds_f : CGAL::faces(hds))
    {
      std::vector<Point_3> face_points;

      // If they exist (exterior skeleton), the first four faces of the SLS correspond
      // to the outer frame, and should be ignored.
      if(ignore_frame_faces && fc++ < 4)
        continue;

      HDS_Halfedge_const_handle hds_h = hds_f->halfedge(), done = hds_h;
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
      CGAL_assertion(hds_h == contour_h);
      const bool is_vertical = (contour_h->weight() == vertical_weight);
#endif

      do
      {
        HDS_Vertex_const_handle hds_tv = hds_h->vertex();

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
        // this computes the snapped position but does not change the geometry of the skeleton
        if(is_vertical && !hds_tv->is_contour())
          snap_skeleton_vertex<HDS, Geom_traits>(hds_h, contour_h, snapped_positions);
#endif

        face_points.emplace_back(hds_tv->point().x(), hds_tv->point().y(), hds_tv->time());

        hds_h = hds_h->next();
      }
      while(hds_h != done);

      if(face_points.size() < 3)
      {
        std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
        continue;
      }

      triangulate_skeleton_face(face_points, invert_faces, points, faces);
    }
  }

  void construct_lateral_faces(const Straight_skeleton_2& ss,
                              const Offset_builder& offset_builder,
                              const FT height,
                              std::vector<Point_3>& points,
                              std::vector<std::vector<std::size_t> >& faces,
                              const std::unordered_map<HDS_Halfedge_const_handle, Point_2>& offset_points,
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
                              const FT& vertical_weight,
                              std::map<Point_2, Point_2>& snapped_positions,
#endif
                              const bool ignore_frame_faces = false,
                              const bool invert_faces = false)
  {
    CGAL_precondition(height != default_extrusion_height<FT>());

    const bool extrude_upwards = is_positive(height);
    const FT abs_height = CGAL::abs(height);

    const HDS& hds = static_cast<const HDS&>(ss);

    std::size_t fc = 0;
    for(const HDS_Face_handle hds_f : CGAL::faces(hds))
    {
      // If they exist (exterior skeleton), the first four faces of the SLS correspond
      // to the outer frame, and should be ignored.
      if(ignore_frame_faces && fc++ < 4)
        continue;

      std::vector<Point_3> face_points;

      HDS_Halfedge_const_handle hds_h = hds_f->halfedge(), done = hds_h;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      HDS_Halfedge_const_handle contour_h = hds_h->defining_contour_edge();
      CGAL_assertion(hds_h == contour_h);
      const bool is_vertical = (contour_h->weight() == vertical_weight);
#endif

      do
      {
        HDS_Vertex_const_handle hds_sv = hds_h->opposite()->vertex();
        HDS_Vertex_const_handle hds_tv = hds_h->vertex();

        // Compare_offset_against_event_time compares offset to node->time(),
        // so when the offset is greater or equal than the node->time(), the node is the face
        auto compare_time_to_offset = [&](HDS_Vertex_const_handle node) -> CGAL::Comparison_result
        {
          if(node->is_contour())
            return CGAL::LARGER; // offset > 0 and contour nodes' time is 0
          else
            return offset_builder.Compare_offset_against_event_time(abs_height, node);
        };

        const CGAL::Comparison_result sc = compare_time_to_offset(hds_sv);
        const CGAL::Comparison_result tc = compare_time_to_offset(hds_tv);

        // if the offset is crossing at the source, it will be added when seen as a target
        // from the previous halfedge

        if(sc != tc && sc != CGAL::EQUAL && tc != CGAL::EQUAL)
        {
          // std::cout << "sc != tc" << std::endl;
          CGAL_assertion(sc != CGAL::EQUAL && tc != CGAL::EQUAL);

          HDS_Halfedge_const_handle hds_off_h = hds_h;
          if(hds_h->slope() == CGAL::NEGATIVE) // ensure same geometric point on both sides
            hds_off_h = hds_off_h->opposite();

          // The offset point must already been computed in the offset builder visitor
          auto off_p = offset_points.find(hds_off_h);
          CGAL_assertion(off_p != offset_points.end());

          face_points.emplace_back(off_p->second.x(), off_p->second.y(), height);
        }

        if(tc != CGAL::SMALLER)
        {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
          if(is_vertical && !hds_tv->is_contour())
            snap_skeleton_vertex<HDS, Geom_traits>(hds_h, contour_h, snapped_positions);
#endif

          const Point_2& off_p = hds_tv->point();
          // ->time() could be an approximation of the true time. If we are here, the target's time
          // is smaller than the height, but we still sanitize it in case of numerical errors in ->time().
          const FT time = boost::algorithm::clamp(hds_tv->time(), FT(0), abs_height);
          face_points.emplace_back(off_p.x(), off_p.y(), extrude_upwards ? time : - time);
        }

        hds_h = hds_h->next();
      }
      while(hds_h != done);

      if(face_points.size() < 3)
      {
        std::cerr << "Warning: sm_vs has size 1 or 2: offset crossing face at a single point?" << std::endl;
        continue;
      }

      triangulate_skeleton_face(face_points, invert_faces, points, faces);
    }
  }

public:
  // this is roughly "CGAL::create_interior_weighted_skeleton_and_offset_polygons_with_holes_2()",
  // but we want to know the intermediate straight skeleton to build the lateral faces of the 3D mesh
  template <typename PolygonWithHoles,
            typename PointRange, typename FaceRange>
  bool inward_construction(const PolygonWithHoles& pwh,
                          const std::vector<std::vector<FT> >& speeds,
                          const FT vertical_weight,
                          const FT height,
                          PointRange& points,
                          FaceRange& faces)
  {
    CGAL_precondition(!is_zero(height));

    const FT abs_height = abs(height);

    // bottom face (z=0)
    construct_horizontal_faces(pwh, 0 /*height*/, points, faces, true /*invert faces*/);

    // Avoid recomputing offset points multiple times
    std::unordered_map<HDS_Halfedge_const_handle, Point_2> offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    // This is to deal with vertical slabs: we temporarily give a non-vertical slope to be able
    // to construct the SLS, and we then snap back the position to the vertical planes.
    // Note that points in non-vertical faces are also changed a bit
    std::map<Point_2, Point_2> snapped_positions;
#endif

    Straight_skeleton_2_ptr ss_ptr;

    bool is_default_extrusion_height = (height == default_extrusion_height<FT>());

    // @partial_wsls_pwh interior SLS of weighted polygons with holes can have skeleton faces with holes
    // The current postprocessing is in the function EnforceSimpleConnectedness, but it has not yet
    // been made compatible with partial skeletons.
    if(is_default_extrusion_height || pwh.number_of_holes() != 0)
    {
      ss_ptr = CGAL::create_interior_weighted_straight_skeleton_2(
                 CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                 CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                 pwh.holes_begin(), pwh.holes_end(),
                 std::begin(speeds[0]), std::end(speeds[0]),
                 std::next(std::begin(speeds)), std::end(speeds),
                 m_gt);
    }
    else
    {
      ss_ptr = CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
                 abs_height,
                 CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                 CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                 pwh.holes_begin(), pwh.holes_end(),
                 std::begin(speeds[0]), std::end(speeds[0]),
                 std::next(std::begin(speeds)), std::end(speeds),
                 m_gt);
    }

    if(!ss_ptr)
    {
      std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
      return false;
    }

#ifdef CGAL_SLS_DEBUG_DRAW
    // print_straight_skeleton(*ss_ptr);
    CGAL::draw(*ss_ptr);
#endif

    if(is_default_extrusion_height)
    {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      construct_lateral_faces(*ss_ptr, points, faces, vertical_weight, snapped_positions);
#else
      construct_lateral_faces(*ss_ptr, points, faces);
#endif
    }
    else // height != default_extrusion_height<FT>()
    {
#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      Visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
      Visitor visitor(*ss_ptr, vertical_weight, offset_points);
#endif
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      Offset_polygons raw_output;
      ob.construct_offset_contours(abs_height, std::back_inserter(raw_output));

      Offset_polygons_with_holes output = CGAL::arrange_offset_polygons_2<Polygon_with_holes_2>(raw_output);
      construct_horizontal_faces(output, height, points, faces);

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points, vertical_weight, snapped_positions);
#else
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points);
#endif
    }

    // if the height is negative (extruding downwards), we need to invert all the faces
    if(is_negative(height))
    {
      for(auto& f : faces)
        std::swap(f[0], f[1]);
    }

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    apply_snapping<Geom_traits>(points, snapped_positions);
#endif

    return true;
  }

  template <typename PolygonWithHoles,
            typename PointRange, typename FaceRange>
  bool outward_construction(const PolygonWithHoles& pwh,
                            const std::vector<std::vector<FT> >& speeds,
                            const FT vertical_weight,
                            const FT height,
                            PointRange& points,
                            FaceRange& faces)
  {
    CGAL_precondition(!is_zero(height));
    CGAL_precondition(height != default_extrusion_height<FT>()); // was checked before, this is just a reminder

    const FT abs_height = abs(height);

    // bottom face (z=0)
    construct_horizontal_faces(pwh, 0 /*height*/, points, faces, true /*invert faces*/);

    // Avoid recomputing offset points multiple times
    std::unordered_map<HDS_Halfedge_const_handle, Point_2> offset_points;

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    // This is to deal with vertical slabs: we temporarily give a non-vertical slope to be able
    // to construct the SLS, and we then snap back the position to the vertical planes.
    // Note that points in non-vertical faces are also changed a bit
    std::map<Point_2, Point_2> snapped_positions;
#endif

    Offset_polygons raw_output; // accumulates for both the outer boundary and the holes

    // the exterior of a polygon with holes is the exterior of its outer boundary,
    // and the interior of its inverted holes
    //
    // Start with the outer boundary
    {
      Straight_skeleton_2_ptr ss_ptr = CGAL_SS_i::create_partial_exterior_weighted_straight_skeleton_2(
                                          abs_height,
                                          CGAL_SS_i::vertices_begin(pwh.outer_boundary()),
                                          CGAL_SS_i::vertices_end(pwh.outer_boundary()),
                                          std::begin(speeds[0]), std::end(speeds[0]),
                                          m_gt);

      if(!ss_ptr)
      {
        std::cerr << "Error: encountered an error during outer skeleton construction" << std::endl;
        return false;
      }

#ifdef CGAL_SLS_DEBUG_DRAW
      // print_straight_skeleton(*ss_ptr);
      CGAL::draw(*ss_ptr);
#endif

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      Visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
      Visitor visitor(*ss_ptr, offset_points);
#endif
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      ob.construct_offset_contours(abs_height, std::back_inserter(raw_output));

      // Manually filter the offset of the outer frame
      std::swap(raw_output[0], raw_output.back());
      raw_output.pop_back();

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points,
                              vertical_weight, snapped_positions,
                              true /*ignore frame faces*/, true /*invert faces*/);
#else
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points,
                              true /*ignore frame faces*/, true /*invert faces*/);
#endif
    }

    // now, deal with the holes

    std::size_t hole_id = 0;
    for(auto hit=pwh.holes_begin(); hit!=pwh.holes_end(); ++hit, ++hole_id)
    {
      Polygon_2 hole = *hit; // intentional copy
      hole.reverse_orientation();

      // this is roughly "CGAL::create_exterior_weighted_skeleton_and_offset_polygons_with_holes_2()",
      // but we want to know the intermediate straight skeleton to build the lateral faces of the 3D mesh

      std::vector<Polygon_2> no_holes;
      std::vector<std::vector<FT> > no_speeds;
      Straight_skeleton_2_ptr ss_ptr = CGAL_SS_i::create_partial_interior_weighted_straight_skeleton_2(
                                          abs_height,
                                          CGAL_SS_i::vertices_begin(hole),
                                          CGAL_SS_i::vertices_end(hole),
                                          std::begin(no_holes), std::end(no_holes),
                                          std::begin(speeds[hole_id]), std::end(speeds[hole_id]),
                                          std::begin(no_speeds), std::end(no_speeds),
                                          m_gt);

      if(!ss_ptr)
      {
        std::cerr << "Error: encountered an error during skeleton construction" << std::endl;
        return EXIT_FAILURE;
      }

#ifdef CGAL_SLS_DEBUG_DRAW
      // print_straight_skeleton(*ss_ptr);
      CGAL::draw(*ss_ptr);
#endif

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      Visitor visitor(*ss_ptr, offset_points, vertical_weight, snapped_positions);
#else
      Visitor visitor(*ss_ptr, offset_points);
#endif
      Offset_builder ob(*ss_ptr, Offset_builder_traits(), visitor);
      ob.construct_offset_contours(abs_height, std::back_inserter(raw_output));

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points,
                              vertical_weight, snapped_positions,
                              false /*no outer frame*/, true /*invert faces*/);
#else
      construct_lateral_faces(*ss_ptr, ob, height, points, faces, offset_points,
                              false /*no outer frame*/, true /*invert faces*/);
#endif
    }

    // - the exterior offset of the outer boundary is built by creating an extra frame and turning
    // the outer boundary into a hole. Hence, it needs to be reversed back to proper orientation
    // - the exterior offset of the holes is built by reversing the holes and computing an internal
    // skeleton. Hence, the result also needs to be reversed.
    for(std::shared_ptr<Polygon_2> ptr : raw_output)
      ptr->reverse_orientation();

    Offset_polygons_with_holes output = CGAL::arrange_offset_polygons_2<Polygon_with_holes_2>(raw_output);
    construct_horizontal_faces(output, height, points, faces);

    // if the height is negative (extruding downwards), we need to invert all the faces
    if(is_negative(height))
    {
      for(auto& f : faces)
        std::swap(f[0], f[1]);
    }

#ifdef CGAL_SLS_SNAP_TO_VERTICAL_SLABS
    apply_snapping<Geom_traits>(points, snapped_positions);
#endif

    return true;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

// convert angles (in degrees) to weights, and handle vertical angles
template <typename FT, typename AngleRange>
void convert_angles(AngleRange& angles)
{
  CGAL_precondition(!angles.empty());

  auto angle_to_weight = [](const FT& angle) -> FT
  {
    CGAL_precondition(0 < angle && angle < 180);

    // @todo should this be an epsilon around 90°? As theta goes to 90°, tan(theta) goes to infinity
    // and thus we could get numerical issues (overflows) if the kernel is not exact
    if(angle == 90)
      return 0;
    else
      return std::tan(CGAL::to_double(angle * CGAL_PI / 180));
  };

  for(auto& contour_angles : angles)
    for(FT& angle : contour_angles)
      angle = angle_to_weight(angle);
}

// handle vertical angles (inf speed)
// returns
// - whether the weights are positive or negative (inward / outward)
// - whether the input weights are valid
// - the weight of vertical slabs
template <typename FT, typename WeightRange>
std::tuple<Slope, bool, FT>
preprocess_weights(WeightRange& weights)
{
  CGAL_precondition(!weights.empty());

  Slope slope = Slope::UNKNOWN;

  FT max_value = 0; // non-inf, maximum absolute value
  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      // '0' means a vertical slab, aka 90° angle (see preprocess_angles())
      if(w == 0)
        continue;

      // determine whether weights indicate all inward or all outward
      if(slope == Slope::UNKNOWN)
      {
        // w is neither 0 nor inf here
        slope = (w > 0) ? Slope::INWARD : Slope::OUTWARD;
      }
      else if(slope == Slope::INWARD && w < 0)
      {
        std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }
      else if(slope == Slope::OUTWARD && w > 0)
      {
        std::cerr << "Error: mixing positive and negative weights is not yet supported" << std::endl;
        return {Slope::UNKNOWN, false, FT(-1)};
      }

      // if we are going outwards, it is just an interior skeleton with opposite weights
      w = CGAL::abs(w);
      if(w > max_value)
        max_value = w;
    }
  }

  if(slope == Slope::UNKNOWN)
  {
    std::cerr << "Warning: all edges vertical?" << std::endl;
    slope = Slope::VERTICAL;
  }

  // Take a weight which is a large % of the max value to ensure there's no ambiguity
  //
  // Since the max value might not be very close to 90°, take the max between of the large-% weight
  // and the weight corresponding to an angle of 89.9999999°
  const FT weight_of_89d9999999 = 572957787.3425436; // tan(89.9999999°)
  const FT scaled_max = (std::max)(weight_of_89d9999999, FT(1e3) * max_value);

  for(auto& contour_weights : weights)
  {
    for(FT& w : contour_weights)
    {
      if(w == FT(0))
        w = scaled_max;
    }
  }

  return {slope, true, scaled_max};
}

template <typename PolygonWithHoles,
          typename WeightRange,
          typename PolygonMesh,
          typename NamedParameters>
bool extrude_skeleton(const PolygonWithHoles& pwh,
                      WeightRange& weights,
                      PolygonMesh& out,
                      const NamedParameters& np)
{
  namespace PMP = ::CGAL::Polygon_mesh_processing;

  using Polygon_2 = typename PolygonWithHoles::Polygon_2;
  using Point = typename boost::range_value<Polygon_2>::type;

  using Default_kernel = typename Kernel_traits<Point>::type;
  using Geom_traits = typename internal_np::Lookup_named_param_def<internal_np::geom_traits_t,
                                                                   NamedParameters,
                                                                   Default_kernel>::type;
  using FT = typename Geom_traits::FT;

  using Point_3 = typename Geom_traits::Point_3;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);
  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, CGAL::internal_np::geom_traits));

  const FT height = choose_parameter(get_parameter(np, internal_np::maximum_height), default_extrusion_height<FT>());

  Slope slope;
  bool valid_input;
  FT vertical_weight;
  std::tie(slope, valid_input, vertical_weight) = preprocess_weights<FT>(weights);

  if(!valid_input)
  {
    if(verbose)
      std::cerr << "Error: invalid input weights" << std::endl;
    return false;
  }

  if(verbose)
  {
    switch(slope)
    {
      case Slope::UNKNOWN: std::cout << "Slope is UNKNOWN??" << std::endl; break;
      case Slope::INWARD: std::cout << "Slope is INWARD" << std::endl; break;
      case Slope::OUTWARD: std::cout << "Slope is OUTWARD" << std::endl; break;
      case Slope::VERTICAL: std::cout << "Slope is VERTICAL" << std::endl; break;
    }
  }

  // End of preprocessing, start the actual skeleton computation

  if(slope != Slope::INWARD && height == default_extrusion_height<FT>())
  {
    if(verbose)
      std::cerr << "Error: height must be specified when using an outward (or vertical) slope" << std::endl;
    return false;
  }

  // build a soup, to be converted to a mesh afterwards
  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > faces;
  points.reserve(2 * pwh.outer_boundary().size()); // just a reasonable guess
  faces.reserve(2 * pwh.outer_boundary().size() + 2*pwh.number_of_holes());

  Extrusion_builder<Geom_traits> builder(gt);
  bool res;
  if(slope != Slope::OUTWARD) // INWARD or VERTICAL
    res = builder.inward_construction(pwh, weights, vertical_weight, height, points, faces);
  else
    res = builder.outward_construction(pwh, weights, vertical_weight, height, points, faces);

  if(!res)
    return false;

#ifdef CGAL_SLS_OUTPUT_FILES
  // This soup provides one connected component per edge of the input polygon
  CGAL::IO::write_polygon_soup("extruded_skeleton_soup.off", points, faces, CGAL::parameters::stream_precision(17));
#endif

  // Convert the triangle soup to a triangle mesh

  PMP::merge_duplicate_points_in_polygon_soup(points, faces);
  if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    PMP::orient_polygon_soup(points, faces);
  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));

  PMP::polygon_soup_to_polygon_mesh(points, faces, out);

  CGAL_warning(is_valid_polygon_mesh(out) && is_closed(out));

  return true;
}

} // namespace internal
} // namespace Straight_skeleton_extrusion

// Documented in the Straight_skeleton_2 package
template <typename PolygonWithHoles,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool extrude_skeleton(const PolygonWithHoles& pwh,
                      PolygonMesh& out,
                      const NamedParameters& np = parameters::default_values(),
                      std::enable_if_t<CGAL_SS_i::has_Hole_const_iterator<PolygonWithHoles>::value>* = nullptr)
{
  namespace SSEI = Straight_skeleton_extrusion::internal;

  using Polygon_2 = typename PolygonWithHoles::General_polygon_2;
  using K = typename Kernel_traits<typename boost::range_value<Polygon_2>::type>::type;
  using FT = typename K::FT;
  using Default_speed_type = std::vector<std::vector<FT> >;

  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter_reference;

  const bool has_weights = !(is_default_parameter<NamedParameters, internal_np::weights_param_t>::value);
  const bool has_angles = !(is_default_parameter<NamedParameters, internal_np::angles_param_t>::value);

  Default_speed_type def_speeds; // never used, but needed for compilation

  if(has_weights)
  {
    // not using "Angles" here is on purpose, I want to copy the range but the NP is ref only
    auto weights = choose_parameter(get_parameter_reference(np, internal_np::weights_param), def_speeds);
    return SSEI::extrude_skeleton(pwh, weights, out, np);
  }
  else if(has_angles)
  {
    // not using "Weights" here is on purpose, I want to copy the range but the NP is ref only
    auto angles = choose_parameter(get_parameter_reference(np, internal_np::angles_param), def_speeds);
    SSEI::convert_angles<FT>(angles);
    return SSEI::extrude_skeleton(pwh, angles, out, np);
  }
  else // neither angles nor weights were provided
  {
    std::vector<std::vector<FT> > uniform_weights;
    uniform_weights.reserve(pwh.number_of_holes() + 1);

    uniform_weights.push_back(std::vector<FT>(pwh.outer_boundary().size(), FT(1)));
    for(const auto& hole : pwh.holes())
      uniform_weights.push_back(std::vector<FT>(hole.size(), FT(1)));

    return SSEI::extrude_skeleton(pwh, uniform_weights, out, np);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename Polygon, // without holes
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool extrude_skeleton(const Polygon& p,
                      PolygonMesh& out,
                      const NamedParameters& np = parameters::default_values(),
                      std::enable_if_t<!CGAL_SS_i::has_Hole_const_iterator<Polygon>::value>* = nullptr)
{
  using Polygon_with_holes = CGAL::General_polygon_with_holes_2<Polygon>;

  return extrude_skeleton(Polygon_with_holes(p), out, np);
}

} // namespace CGAL

#endif // CGAL_SLS_EXTRUSION_EXTRUDE_SKELETON_H
