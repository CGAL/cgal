// Copyright (c) 2019-2023 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Pierre Alliez
//                 Cedric Portaneri,
//                 Mael Rouxel-Labbé
//                 Andreas Fabri
//                 Michael Hemmer
//
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_2_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_2_H

#ifdef CGAL_AW2_DEBUG_PP
 #ifndef CGAL_AW2_DEBUG
  #define CGAL_AW2_DEBUG
  #define CGAL_AW2_DEBUG_INITIALIZATION
  #define CGAL_AW2_DEBUG_STEINER_COMPUTATION
  #define CGAL_AW2_DEBUG_QUEUE
  #define CGAL_AW2_DEBUG_EDGE_STATUS
  #define CGAL_AW2_DEBUG_MANIFOLDNESS
 #endif
#endif

// @bug
// - restore hints in closest point call
// - point_set_wrap broken with unordered queue
// - mixed_input broken (ordered or unordered queue)
// - fix return (*this)(p, bb, bound, Boolean_tag<internal::Has_static_filters<GeomTraits>::value>());

// @todo
// - Add initialize_with_cavities
// - arrange polygons if we have no purged cavities

// @todo long
// - check what could be factorized with AW3 (e.g. oracle base?)
// - test relaxed sphere marching
// - test expensive assertions in T3
// - document and test new kernel functors
// - document and test DT2 face base classes
// - bench queues and whatnot
// - if(is_neighbor_cc_in_offset) should be removed if sphere marching is used
// - test new traversable criterion (ball fits)

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_triangulation_face_base_2.h>
#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_triangulation_vertex_base_2.h>
#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_2/internal/gate_priority_queue.h>
#include <CGAL/Alpha_wrap_2/internal/geometry_utils.h>
#include <CGAL/Alpha_wrap_2/internal/oracles.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_2.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Default.h>
#include <CGAL/Named_function_parameters.h>
#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
 #include <CGAL/Modifiable_priority_queue.h>
#endif
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // only if non-manifoldness is not treated
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>

#include <algorithm>
#include <array>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

namespace {

namespace AW2i = ::CGAL::Alpha_wraps_2::internal;

} // unnamed namespace

struct Wrapping_default_visitor
{
  Wrapping_default_visitor() { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_begin(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_flood_fill_begin(const AlphaWrapper&) { }

  // Whether the flood filling process should continue
  template <typename Wrapper>
  constexpr bool go_further(const Wrapper&) { return true; }

  template <typename AlphaWrapper, typename Gate>
  void before_edge_treatment(const AlphaWrapper&, const Gate&) { }

  template <typename Wrapper, typename Point>
  void before_Steiner_point_insertion(const Wrapper&, const Point&) { }

  template <typename Wrapper, typename VertexHandle>
  void after_Steiner_point_insertion(const Wrapper&, VertexHandle) { }

  template <typename AlphaWrapper>
  void on_flood_fill_end(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_end(const AlphaWrapper&) { };
};

template <typename Oracle_,
          typename Triangulation_ = CGAL::Default>
class Alpha_wrapper_2
{
  using Oracle = Oracle_;

  // Triangulation
  using Base_GT = typename Oracle::Geom_traits;
  using Default_Gt = Robust_circumcenter_filtered_traits_2<Base_GT>;

  using Default_Vb = Alpha_wrap_triangulation_vertex_base_2<Default_Gt>;
  using Default_Fb = Alpha_wrap_triangulation_face_base_2<Default_Gt>;
  using Default_Fbt = Face_base_with_timestamp<Default_Fb>; // for determinism
  using Default_Tds = CGAL::Triangulation_data_structure_2<Default_Vb, Default_Fbt>;
  using Default_Triangulation = CGAL::Delaunay_triangulation_2<Default_Gt, Default_Tds>;

public:
  using Triangulation = typename Default::Get<Triangulation_, Default_Triangulation>::type;

  // Use the geom traits from the triangulation, and trust the (advanced) user that provided it
  using Geom_traits = typename Triangulation::Geom_traits;

private:
  using Face_handle = typename Triangulation::Face_handle;
  using Face_circulator = typename Triangulation::Face_circulator;
  using Edge = typename Triangulation::Edge;
  using Vertex_handle = typename Triangulation::Vertex_handle;
  using Vertex_circulator = typename Triangulation::Vertex_circulator;
  using Locate_type = typename Triangulation::Locate_type;

  using Gate = internal::Gate<Triangulation>;

  // A sorted queue is a priority queue sorted by circumradius, and is experimentally significantly
  // slower. However, intermediate results are usable: at each point of the algorithm, the wrap
  // has a somewhat uniform mesh element size.
  //
  // An unsorted queue is a LIFO queue, and is experimentally much faster (~35%),
  // but intermediate results are not useful: a LIFO carves deep before than wide,
  // and can thus for example leave scaffolding faces till almost the end of the refinement.
  //
  // @todo bench this in 2D
#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
  using Alpha_PQ = Modifiable_priority_queue<Gate, Less_gate, Gate_ID_PM<Triangulation>, CGAL_BOOST_PAIRING_HEAP>;
#else
  using Alpha_PQ = std::stack<Gate>;
#endif

  using FT = typename Geom_traits::FT;
  using Point_2 = typename Geom_traits::Point_2;
  using Vector_2 = typename Geom_traits::Vector_2;
  using Disk_2 = typename Geom_traits::Disk_2;
  using Iso_rectangle_2 = typename Geom_traits::Iso_rectangle_2;

  using SC = Simple_cartesian<double>;
  using SC_Point_2 = SC::Point_2;
  using SC_Vector_2 = SC::Vector_2;
  using SC_Iso_rectangle_2 = SC::Iso_rectangle_2;
  using SC2GT = Cartesian_converter<SC, Geom_traits>;

  using Seeds = std::vector<Point_2>;

protected:
  Oracle m_oracle;
  SC_Iso_rectangle_2 m_bbox;

  FT m_alpha = FT(-1), m_sq_alpha = FT(-1);
  FT m_offset = FT(-1), m_sq_offset = FT(-1);

  Seeds m_seeds;

  Triangulation m_tr;

  Alpha_PQ m_queue;

public:
  Alpha_wrapper_2()
#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
      // '4096' is an arbitrary, not-too-small value for the largest ID in queue initialization
    : m_queue(4096)
#endif
  {
    // Due to the Steiner point computation being a dichotomy, the algorithm is inherently inexact
    // and passing exact kernels is explicitly disabled to ensure no misunderstanding.
    static_assert(std::is_floating_point<FT>::value);
  }

  Alpha_wrapper_2(const Oracle& oracle)
    :
      m_oracle(oracle),
      m_tr(Geom_traits(oracle.geom_traits()))
#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
      , m_queue(4096)
#endif
  {
    static_assert(std::is_floating_point<FT>::value);
  }

public:
  const Geom_traits& geom_traits() const { return m_tr.geom_traits(); }
  Oracle& oracle() { return m_oracle; }
  const Oracle& oracle() const { return m_oracle; }
  Triangulation& triangulation() { return m_tr; }
  const Triangulation& triangulation() const { return m_tr; }
  const Alpha_PQ& queue() const { return m_queue; }

  double default_alpha() const
  {
    const Bbox_2 bbox = m_oracle.bbox();
    const double diag_length = std::sqrt(square(bbox.xmax() - bbox.xmin()) +
                                         square(bbox.ymax() - bbox.ymin()));

    return diag_length / 20.;
  }

private:
  const Point_2& circumcenter(const Face_handle f) const
  {
    // We only cross an infinite edge once, so this isn't going to be recomputed many times
    if(m_tr.is_infinite(f))
    {
      const int inf_index = f->index(m_tr.infinite_vertex());
      f->set_circumcenter(
            geom_traits().construct_circumcenter_2_object()(m_tr.point(f, Triangulation::ccw(inf_index)),
                                                            m_tr.point(f, Triangulation::cw(inf_index))));
    }

    return f->circumcenter(geom_traits());
  }

public:
  template <typename OutputPolygons,
            typename InputNamedParameters = parameters::Default_named_parameters,
            typename OutputNamedParameters = parameters::Default_named_parameters>
  void operator()(const double alpha, // = default_alpha()
                  const double offset, // = alpha / 30.
                  OutputPolygons& output_polygons,
                  const InputNamedParameters& in_np = parameters::default_values(),
                  const OutputNamedParameters& out_np = parameters::default_values())
  {
    namespace PMP = Polygon_mesh_processing;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    //
    using Visitor = typename internal_np::Lookup_named_param_def<
                      internal_np::visitor_t,
                      InputNamedParameters,
                      Wrapping_default_visitor // default
                    >::reference;

    Wrapping_default_visitor default_visitor;
    Visitor visitor = choose_parameter(get_parameter_reference(in_np, internal_np::visitor), default_visitor);

    // Points used to create initial cavities
    m_seeds = choose_parameter(get_parameter_reference(in_np, internal_np::seed_points), Seeds());

    // Whether or not some faces should be reflagged as "inside" after the refinement+carving loop
    // as ended, as to ensure that the result is not only combinatorially manifold, but also
    // geometrically manifold.
    //
    // -- Warning --
    // Manifoldness postprocessing will be performed even if the wrapping is interrupted (and
    // this option is enabled).
    const bool do_enforce_manifoldness = choose_parameter(get_parameter(in_np, internal_np::do_enforce_manifoldness), true);

    // Whether to keep pockets of "outside" faces that are not connected to the exterior (or to the
    // initial cavities, if used).
    //
    // -- Warning --
    // Pockets of "outside" faces will be purged even if the wrapping is interrupted (and
    // this option is enabled).
    const bool keep_inner_ccs = choose_parameter(get_parameter(in_np, internal_np::keep_inner_connected_components), false);

    // This parameter enables avoiding recomputing the triangulation from scratch when wrapping
    // the same input for multiple values of alpha (and typically the same offset values).
    //
    // -- Warning --
    // If this is enabled, the 2D triangulation will NOT be re-initialized at launch.
    // This means that the triangulation is NOT cleared, even if:
    // - you use an alpha value that is greater than what was used in a previous run; you will
    //   obtain the same result as the last run.
    // - you use a different offset value between runs, you might then get points that are not
    //   on the offset surface corresponding to that corresponding to the latter offset value.
    const bool refining = choose_parameter(get_parameter(in_np, internal_np::refine_triangulation), false);

#ifdef CGAL_AW2_TIMER
    CGAL::Real_timer t;
    t.start();
#endif

    visitor.on_alpha_wrapping_begin(*this);

    if(!initialize(alpha, offset, refining))
      return;

#ifdef CGAL_AW2_TIMER
    t.stop();
    std::cout << "Initialization took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

#ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation("starting_wrap.off");
#endif

    alpha_flood_fill(visitor);

#ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation("flood_filled_wrap.off");
#endif

#ifdef CGAL_AW2_TIMER
    t.stop();
    std::cout << "Flood filling took: " << t.time() << " s." << std::endl;
    t.reset();
    t.start();
#endif

    if(do_enforce_manifoldness)
    {
      make_manifold();

#ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation("manifold_wrap.off");
#endif

#ifdef CGAL_AW2_TIMER
      t.stop();
      std::cout << "Manifoldness post-processing took: " << t.time() << " s." << std::endl;
      t.reset();
      t.start();
#endif
    }

    if(!keep_inner_ccs)
    {
      // We could purge *before* manifold enforcement, but making the mesh manifold is
      // very cheap in most cases, so it is better to keep the code simpler.
      purge_inner_connected_components();

#ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
      dump_triangulation("purged_wrap.off");
#endif
    }

    extract_boundary(output_polygons);

#ifdef CGAL_AW2_TIMER
    t.stop();
    std::cout << "Surface extraction took: " << t.time() << " s." << std::endl;
#endif

#ifdef CGAL_AW2_DEBUG
    std::cout << "Alpha wrap polygons:  " << output_polygons.size() << std::endl;

 #ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation("final_wrap.off");
 #endif
#endif

    visitor.on_alpha_wrapping_end(*this);
  }

  // Convenience overloads
  template <typename OutputPolygons>
  void operator()(const double alpha,
                  OutputPolygons& output_polygons)
  {
    return operator()(alpha, alpha / 30. /*offset*/, output_polygons);
  }

  template <typename OutputPolygons>
  void operator()(OutputPolygons& output_polygons)
  {
    return operator()(default_alpha(), output_polygons);
  }

  // This function is public only because it is used in the tests
  SC_Iso_rectangle_2 construct_bbox(const double offset)
  {
    // Input axis-aligned bounding box
    SC_Iso_rectangle_2 bbox = m_oracle.bbox();
    const SC_Point_2 bbox_centroid = midpoint((bbox.min)(), (bbox.max)());

    // Scale a bit to create the initial points not too close to the input
    double scaling = 1.2;
    CGAL::Aff_transformation_2<SC> scale(SCALING, scaling);
    bbox = SC_Iso_rectangle_2(scale.transform((bbox.min)()), scale.transform((bbox.max)()));

    // Translate bbox back to initial centroid
    const SC_Point_2 bbox_transformed_centroid = midpoint((bbox.min)(), (bbox.max)());
    const SC_Vector_2 diff_centroid = bbox_centroid - bbox_transformed_centroid;
    CGAL::Aff_transformation_2<SC> centroid_translate(TRANSLATION, diff_centroid);
    bbox = SC_Iso_rectangle_2(centroid_translate.transform((bbox.min)()),
                           centroid_translate.transform((bbox.max)()));

    // Add the offset
    SC_Vector_2 offset_ext = std::sqrt(2.) * offset * SC_Vector_2(1, 1);
    CGAL::Aff_transformation_2<SC> translate_m(TRANSLATION, - offset_ext);
    CGAL::Aff_transformation_2<SC> translate_M(TRANSLATION,   offset_ext);
    bbox = SC_Iso_rectangle_2(translate_m.transform((bbox.min)()), translate_M.transform((bbox.max)()));

    return bbox;
  }

private:
  // The distinction between inside boundary and outside boundary is the presence of faces
  // being flagged for manifoldness: inside boundary considers those outside, and outside
  // boundary considers them inside.
  bool is_on_inside_boundary(Face_handle fh, Face_handle nh) const
  {
    return (fh->is_inside() != nh->is_inside()); // one is "inside", the other is not
  }

  bool is_on_outside_boundary(Face_handle fh, Face_handle nh) const
  {
    return (fh->is_outside() != nh->is_outside()); // one is "outside", the other is not
  }

private:
  // Adjust the bbox & insert its corners to construct the starting triangulation
  void insert_bbox_corners()
  {
    m_bbox = construct_bbox(CGAL::to_double(m_offset));

#ifdef CGAL_AW2_DEBUG_INITIALIZATION
    std::cout << "Insert Bbox vertices" << std::endl;
#endif

    // insert in dt the eight corner vertices of the input loose bounding box
    for(int i=0; i<8; ++i)
    {
      const Point_2 bp = SC2GT()(m_bbox.vertex(i));
      Vertex_handle bv = m_tr.insert(bp);
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
      std::cout << "\t" << bp << std::endl;
#endif
      bv->type() = AW2i::Vertex_type::BBOX_VERTEX;
    }
  }

  // Two criteria:
  // - Faces that are intersecting the input are inside
  // - Faces whose circumcenter is in the offset volume are inside: this is because
  // we need to have outside face circumcenters outside of the volume to ensure
  // that the refinement point is separated from the existing point set.
  Face_label cavity_face_label(const Face_handle fh)
  {
    CGAL_precondition(!m_tr.is_infinite(fh));

    const Triangle_with_outside_info<Geom_traits> tr(fh, geom_traits());
    if(m_oracle.do_intersect(tr))
      return Face_label::INSIDE;

    const Point_2& fh_cc = circumcenter(fh);
    typename Geom_traits::Construct_disk_2 disk = geom_traits().construct_disk_2_object();
    const Disk_2 fh_cc_offset_disk = disk(fh_cc, m_sq_offset);
    const bool is_cc_in_offset = m_oracle.do_intersect(fh_cc_offset_disk);

    return is_cc_in_offset ? Face_label::INSIDE : Face_label::OUTSIDE;
  }

  // Create a cavity using seeds rather than starting from the infinity.
  //
  // For each seed, insert the seeds and its neighboring vertices on a regular hexagon.
  // The idea behind a hexagon rather than e.g. a simple triangle is twofold:
  // - Improve the odds of both starting the algorithm, but also of not immediately stopping,
  //   which could happen if a conflict zone included all the initial cavity when using a simple cavity.
  // - Base the cavity diameter on alpha, and allow its faces to intersect the input. If a single
  //   triangle is used, one is forced to make it small-enough such that it does not intersect
  //   the input, and then it could be forced to be smaller than 'alpha' and the algorithm cannot start,
  //   for example if the seed is close to the offset. If we have many triangles, this is far less
  //   likely to happen.
  //
  // For an edge length a of a regular hexagon, the radius is:
  //   r = a
  // Triangles in a regular hexagon are equilateral triangles with side length a and circumradius a/sqrt(3)
  // Since we want faces of the hexagon to be traversable, we want a such that:
  //   a / sqrt(3) > alpha
  // Hence for radius r:
  //   r > alpha * sqrt(3) ~= 1.732 * alpha
  //
  // Another way is to simply make faces incident to the seed always traversable, and then
  // we only have to ensure faces opposite of the seed are traversable (i.e., radius ~= 1.65 * alpha)
  bool initialize_with_cavities()
  {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
    std::cout << "> Dig cavities" << std::endl;
    std::cout << m_seeds.size() << " seed(s)" << std::endl;
#endif

    CGAL_precondition(!m_seeds.empty());

    Iso_rectangle_2 bbox = SC2GT()(m_bbox);

    std::vector<Vertex_handle> seed_vs;
    for(const Point_2& seed_p : m_seeds)
    {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
      std::cout << "Initialize from seed " << seed_p << std::endl;
#endif

      if(bbox.has_on_unbounded_side(seed_p))
      {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is outside the bounding box" << std::endl;
#endif
        continue;
      }

      // get the closest point on the input
      const Point_2 closest_pt = m_oracle.closest_point(seed_p);
      const FT sq_d_to_closest = geom_traits().compute_squared_distance_2_object()(seed_p, closest_pt);

      if(sq_d_to_closest < m_sq_offset)
      {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is in the offset" << std::endl;
#endif
        continue;
      }

      // Mark the seeds and hexagon vertices as "scaffolding" vertices such that the edges
      // incident to these vertices are always traversable regardless of their circumradius.
      // This is done because otherwise some cavities can appear on the mesh: non-traversable edges
      // with two vertices on the offset, and the third being a deeper inside seed / hex_seed.
      // Making them traversable will cause more refinement than "alpha", but they will eventually
      // not appear anymore in the inside/outside boundary and the shape will look smoother.

      Vertex_handle seed_v = m_tr.insert(seed_p);
      seed_v->type() = AW2i::Vertex_type::SEED_VERTEX;
      seed_vs.push_back(seed_v);

      // Regular hexagon vertices
      const Point_2 center = seed_p;
      const FT radius = FT(1.74) * m_alpha; // sqrt(3) ~= 1.732, rounded up
      const FT half_radius = radius / FT(2);
      const FT half_height = radius * CGAL::approximate_sqrt(FT(3)) / FT(2);

      std::array<Point_2, 6> hex_ps =
      {
        Point_2(center.x() + radius, center.y()),
        Point_2(center.x() + half_radius, center.y() + half_height),
        Point_2(center.x() - half_radius, center.y() + half_height),
        Point_2(center.x() - radius, center.y()),
        Point_2(center.x() - half_radius, center.y() - half_height),
        Point_2(center.x() + half_radius, center.y() - half_height)
      };

      for(const Point_2& seed_neighbor_p : hex_ps)
      {
#ifdef CGAL_AW2_DEBUG_PP
        std::cout << seed_neighbor_p << std::endl;
#endif
        if(bbox.has_on_unbounded_side(seed_neighbor_p))
          continue;

        Vertex_handle hex_v = m_tr.insert(seed_neighbor_p, seed_v->face() /*hint*/);
        hex_v->type() = AW2i::Vertex_type::SEED_VERTEX;
      }
    }

    if(seed_vs.empty())
    {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
      std::cerr << "Error: no acceptable seed was provided" << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW2_DEBUG_INITIALIZATION
    std::cout << m_tr.number_of_vertices() - 4 /*bbox*/ << " vertice(s) due to seeds" << std::endl;
#endif

    for(Vertex_handle seed_v : seed_vs)
    {
      Face_circulator f = m_tr.incident_faces(seed_v), done(f);
      do
      {
        if(!m_tr.is_infinite(f))
          f->set_label(cavity_face_label(f));
      }
      while(++f != done);
    }

    // Should be cheap enough to go through the full triangulation as only seeds have been inserted
    for(Face_handle fh : m_tr.all_face_handles())
    {
      if(fh->is_inside())
        continue;

      // When the algorithm starts from a manually dug hole, infinite faces are initialized
      // as "inside" such that they do not appear on the boundary
      CGAL_assertion(!m_tr.is_infinite(fh));

      for(int i=0; i<3; ++i)
        push_edge(std::make_pair(fh, i));
    }

    if(m_queue.empty())
    {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
      std::cerr << "Could not initialize the algorithm with these seeds, and alpha|offset values" << std::endl;
#endif
      return false;
    }

    return true;
  }

  // tag all infinite faces "outside" and all finite faces "inside"
  // init queue with all convex hull edges
  bool initialize_from_infinity()
  {
    for(Face_handle fh : m_tr.all_face_handles())
    {
      if(m_tr.is_infinite(fh))
      {
        fh->set_label(Face_label::OUTSIDE);
        const int inf_index = fh->index(m_tr.infinite_vertex());
        push_edge(std::make_pair(fh, inf_index));
      }
      else
      {
        fh->set_label(Face_label::INSIDE);
      }
    }

    return true;
  }

  void reset_manifold_labels()
  {
    // No erase counter increment, or it will mess up with a possibly non-empty queue.
    for(Face_handle fh : m_tr.all_face_handles())
      if(fh->label() == Face_label::MANIFOLD)
        fh->set_label(Face_label::OUTSIDE);
  }

  // This function is used in the case of resumption of a previous run: m_tr is not cleared,
  // and we fill the queue with the new parameters.
  bool initialize_from_existing_triangulation()
  {
#ifdef CGAL_AW2_DEBUG_INITIALIZATION
    std::cout << "Restart from a DT of " << m_tr.number_of_faces() << " faces" << std::endl;
#endif

    for(Face_handle fh : m_tr.all_face_handles())
    {
      if(fh->is_inside())
        continue;

      for(int i=0; i<3; ++i)
      {
        if(fh->neighbor(i)->is_inside())
          push_edge(std::make_pair(fh, i));

      }
    }

    return true;
  }

public:
  // Since we duplicate points in polygons, manifold or not does not matter much.
  //
  // It could be tempting to walk "inside" to get a single polygon from pinched polygons,
  // but it would not be better: neither walk could get a proper result for e.g. an "8"-shaped
  // wrap and a crescent-with-ends-touching-shaped wrap.
  template <typename OutputPolygons>
  void extract_boundary(OutputPolygons& output_polygons) const
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "> Extract wrap..." << std::endl;
#endif

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

    for(auto eit = m_tr.finite_edges_begin(), eend = m_tr.finite_edges_end(); eit != eend; ++eit)
    {
      Edge e = *eit;
      if(!e.first->is_outside())
        e = m_tr.mirror_edge(e);

      const Face_handle fh = e.first;
      if(fh->tds_data().processed())
        continue;

      const int s = e.second;
      const Face_handle nh = fh->neighbor(s);
      if(!is_on_outside_boundary(fh, nh))
        continue;

      std::cout << "Starting new polygon" << std::endl;

      std::vector<Point_2> polygon;

      // rotate around the edge's source until finding the next boundary edge
      // and do so until we are the first edge again
      const Face_handle start_face = fh;
      const int start_s = s;

      Face_handle curr_face = start_face;
      int curr_s = start_s;
      do {
        // Add the vertex in CCW order when looking from outside
        polygon.push_back(m_tr.point(curr_face, Triangulation::ccw(curr_s)));
        std::cout << "Adding " << polygon.back() << std::endl;

        curr_face->tds_data().mark_processed();

        curr_s = Triangulation::cw(curr_s);
        for (;;)
        {
          Face_handle curr_nh = curr_face->neighbor(curr_s);
          if(is_on_outside_boundary(curr_face, curr_nh))
            break;

            curr_s = Triangulation::cw(m_tr.mirror_index(curr_face, curr_s));
            curr_face = curr_nh;
        }
      }
      while (curr_face != start_face || curr_s != start_s);

      output_polygons.emplace_back(polygon.begin(), polygon.end());
    }

    // Reset the conflict flags
    for(Face_handle fh : m_tr.all_face_handles())
      fh->tds_data().clear();
  }

private:
  bool is_traversable(const Edge& e) const
  {
    return less_squared_radius_of_min_empty_circle(m_sq_alpha, e, m_tr);
  }

  bool compute_steiner_point(const Face_handle fh,
                             const Face_handle neighbor,
                             Point_2& steiner_point) const
  {
    CGAL_precondition(!m_tr.is_infinite(neighbor));

    typename Geom_traits::Construct_disk_2 disk = geom_traits().construct_disk_2_object();
    typename Geom_traits::Construct_vector_2 vector = geom_traits().construct_vector_2_object();
    typename Geom_traits::Construct_translated_point_2 translate = geom_traits().construct_translated_point_2_object();
    typename Geom_traits::Construct_scaled_vector_2 scale = geom_traits().construct_scaled_vector_2_object();

    const Point_2& neighbor_cc = circumcenter(neighbor);
    const Disk_2 neighbor_cc_offset_disk = disk(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_disk);

#ifdef CGAL_AW2_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*fh << ", " << &*neighbor << ")" << std::endl;

    const Point_2& chc = circumcenter(fh);
    std::cout << "FH" << std::endl;
    std::cout << "\t" << fh->vertex(0)->point() << std::endl;
    std::cout << "\t" << fh->vertex(1)->point() << std::endl;
    std::cout << "\t" << fh->vertex(2)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // fh's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_2& fh_cc = circumcenter(fh);)
    CGAL_assertion_code(const Disk_2 fh_cc_offset_disk = disk(fh_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(fh_cc_offset_disk));

    if(is_neighbor_cc_in_offset)
    {
      const Point_2& fh_cc = circumcenter(fh);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(fh_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW2_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Triangle_with_outside_info<Geom_traits> tr(neighbor, geom_traits());
    if(m_oracle.do_intersect(tr))
    {
      // steiner point is the closest point on input from face centroid with offset
      const Point_2 closest_pt = m_oracle.closest_point(neighbor_cc);
#ifdef CGAL_AW2_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring triangle intersecting the input" << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
#endif
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_2 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_2_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW2_DEBUG_STEINER_COMPUTATION
      std::cout << "Direction: " << unit << std::endl;
      std::cout << "Steiner: " << steiner_point << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW2_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

private:
  // A permissive gate is a gate that we can traverse without checking its circumradius
  enum class Edge_status
  {
    IRRELEVANT = 0,
    HAS_INFINITE_NEIGHBOR, // the face incident to the mirrored edge is infinite (permissive)
    SCAFFOLDING, // incident to a SEED or BBOX vertex (permissive)
    TRAVERSABLE
  };

  inline const char* get_status_message(const Edge_status status)
  {
    constexpr std::size_t status_count = 4;

    // Messages corresponding to Error_code list above. Must be kept in sync!
    static const char* message[status_count] = {
      "Irrelevant edge",
      "Edge incident to infinite neighbor",
      "Edge with a bbox/seed vertex",
      "Traversable edge"
    };

    const std::size_t status_id = static_cast<std::size_t>(status);
    if(status_id > status_count || status_id < 0)
      return "Unknown status";
    else
      return message[status_id];
  }

public:
  // @speed some decent time may be spent constructing Edges (pairs) for no reason as it's always
  // just to grab the .first and .second as soon as it's constructed, and not due to API requirements
  Edge_status edge_status(const Edge& e) const
  {
    CGAL_precondition(!m_tr.is_infinite(e));

#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
    std::cout << "edge status: "
              << m_tr.point(e.first, Triangulation::ccw(e.second)) << " "
              << m_tr.point(e.first, Triangulation::cw(e.second)) << std::endl;
#endif

    // skip if neighbor is "outside" or infinite
    const Face_handle fh = e.first;
    const int id = e.second;
    CGAL_precondition(fh->label() == Face_label::INSIDE || fh->label() == Face_label::OUTSIDE);

    if(!fh->is_outside()) // "inside" or "manifold"
    {
#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
      std::cout << "Edge is inside" << std::endl;
#endif
      return Edge_status::IRRELEVANT;
    }

    const Face_handle nh = fh->neighbor(id);
    CGAL_precondition(fh->label() == Face_label::INSIDE || fh->label() == Face_label::OUTSIDE);

    if(m_tr.is_infinite(nh))
    {
#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
      std::cout << "Infinite neighboring face" << std::endl;
#endif
      return Edge_status::HAS_INFINITE_NEIGHBOR;
    }

    if(nh->is_outside())
    {
#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
      std::cout << "Neighbor already outside" << std::endl;
#endif
      return Edge_status::IRRELEVANT;
    }

    // push if edge is connected to scaffolding vertices
    for(int i=0; i<2; ++i)
    {
      const Vertex_handle vh = fh->vertex((id+i+1)%3);
      if(vh->type() == AW2i::Vertex_type::BBOX_VERTEX ||
         vh->type() == AW2i::Vertex_type::SEED_VERTEX)
      {
#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
        std::cout << "Scaffolding edge due to vertex #" << i << std::endl;
#endif
        return Edge_status::SCAFFOLDING;
      }
    }

    // skip if e's min empty circle radius is smaller than alpha
    if(is_traversable(e))
    {
#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
      std::cout << "traversable" << std::endl;
#endif
      return Edge_status::TRAVERSABLE;
    }

#ifdef CGAL_AW2_DEBUG_EDGE_STATUS
    std::cout << "not traversable" << std::endl;
#endif
    return Edge_status::IRRELEVANT;
  }

private:
  bool push_edge(const Edge& e)
  {
    CGAL_precondition(e.first->is_outside());

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
    // skip if e is already in queue
    if(m_queue.contains_with_bounds_check(Gate(e)))
      return false;
#endif

    // @todo could avoid useless edge_status() calls by doing it after the zombie check
    // for the unsorted priority queue, but AFAIR, it doesn't save noticeable time (and that
    // increases the queue size).
    const Edge_status status = edge_status(e);
    if(status == Edge_status::IRRELEVANT)
      return false;

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
    const FT sqr = smallest_squared_radius_2(e, m_tr);
    const bool is_permissive = (status == Edge_status::HAS_INFINITE_NEIGHBOR ||
                                status == Edge_status::SCAFFOLDING);
    m_queue.resize_and_push(Gate(e, sqr, is_permissive));
#else
    m_queue.push(Gate(e, m_tr));
#endif

#ifdef CGAL_AW2_DEBUG_QUEUE
    const Face_handle fh = e.first;
    const int s = e.second;
    const Point_2& p0 = m_tr.point(fh, Triangulation::ccw(s));
    const Point_2& p1 = m_tr.point(fh, Triangulation::cw(s));

    static int gid = 0;
    std::cout << "Queue insertion #" << gid++ << "\n"
              << "  fh = " << &*fh << " (" << m_tr.is_infinite(fh) << ") " << "\n"
              << "\t" << p0 << "\n\t" << p1 << std::endl;
    std::cout << "  Status: " << get_status_message(status) << std::endl;
 #ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
    std::cout << "  SQR: " << sqr << std::endl;
    std::cout << "  Permissiveness: " << is_permissive << std::endl;

    CGAL_assertion(is_permissive || sqr >= m_sq_alpha);
 #endif // CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
#endif // CGAL_AW2_DEBUG_QUEUE

    return true;
  }

private:
  bool initialize(const double alpha,
                  const double offset,
                  const bool refining)
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "> Initialize..." << std::endl;
#endif

    const bool resuming = refining && (alpha == m_alpha) && (offset == m_offset);

#ifdef CGAL_AW2_DEBUG
    std::cout << "\tAlpha: " << alpha << std::endl;
    std::cout << "\tOffset: " << offset << std::endl;
    std::cout << "\tRefining? " << refining << std::endl;
    std::cout << "\tResuming? " << resuming << std::endl;
#endif

    if(!is_positive(alpha) || !is_positive(offset))
    {
#ifdef CGAL_AW2_DEBUG
      std::cerr << "Error: invalid input parameters: " << alpha << " and " << offset << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW2_DEBUG
    if(refining && alpha > m_alpha)
      std::cerr << "Warning: refining with an alpha greater than the last iteration's!" << std::endl;
    if(refining && offset != m_offset)
      std::cerr << "Warning: refining with a different offset value!" << std::endl;
#endif

    m_alpha = FT(alpha);
    m_sq_alpha = square(m_alpha);
    m_offset = FT(offset);
    m_sq_offset = square(m_offset);

    // Resuming means that we do not need to re-initialize the queue: either we have finished
    // and there is nothing to do, or the interruption was due to a user callback in the visitor,
    // and we can resume with the current queue
    if(resuming)
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Resuming with a queue of size: " << m_queue.size() << std::endl;
#endif

      reset_manifold_labels();
      return true;
    }

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
    m_queue.clear();
#else
    m_queue = { };
#endif

    if(refining)
    {
      // If we are reusing the triangulation, change the label of the extra elements
      // that we have added to ensure a manifold result back to external ("manifold" -> "outside")
      reset_manifold_labels();

      return initialize_from_existing_triangulation();
    }
    else
    {
      m_tr.clear();

      insert_bbox_corners();

      if(m_seeds.empty())
        return initialize_from_infinity();
      else
        return initialize_with_cavities();
    }
  }

  template <typename Visitor>
  bool alpha_flood_fill(Visitor& visitor)
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "> Flood fill..." << std::endl;
#endif

    visitor.on_flood_fill_begin(*this);

    // Explore all finite faces that are reachable from one of the initial outside faces.
    while(!m_queue.empty())
    {
      if(!visitor.go_further(*this))
        return false;

#ifdef CGAL_AW2_DEBUG_QUEUE_PP
      check_queue_sanity();
#endif

      // const& to something that will be popped, but safe as `fh` && `id` are extracted before the pop
      const Gate& gate = m_queue.top();

#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
      if(gate.is_zombie())
      {
        m_queue.pop();
        continue;
      }
#endif

      const Edge& e = gate.edge();
      const Face_handle fh = e.first;
      const int s = e.second;
      const Face_handle nh = fh->neighbor(s);

#ifdef CGAL_AW2_DEBUG_QUEUE
      static int eid = 0;
      std::cout << m_tr.number_of_vertices() << " DT vertices" << std::endl;
      std::cout << m_queue.size() << " edges in the queue" << std::endl;
      std::cout << "Edge " << eid++ << " ["
                << m_tr.point(fh, Triangulation::ccw(s)) << " "
                << m_tr.point(fh, Triangulation::cw(s)) << "]" << std::endl;
      std::cout << "fh = " << &*fh << " (inf: " << m_tr.is_infinite(fh) << "; label: " << fh->label() <<  "), nh = "
                           << &*nh << " (inf: " << m_tr.is_infinite(nh) << "; label: " << nh->label() <<  ")" << "\n";

# ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
      std::cout << "Priority: " << gate.priority() << " (sq alpha: " << m_sq_alpha << ")" << std::endl;
      std::cout << "Permissiveness: " << gate.is_permissive_edge() << std::endl;
# endif
#endif

#ifdef CGAL_AW2_DEBUG_DUMP_EVERY_STEP
      static int i = 0;
      std::string step_name = "results/steps/step_" + std::to_string(static_cast<int>(i++)) + ".off";
      dump_triangulation(step_name);
#endif

      CGAL_precondition(!m_tr.is_infinite(e));
      CGAL_precondition(fh->is_outside());
      CGAL_precondition(nh->label() == Face_label::INSIDE || nh->label() == Face_label::OUTSIDE);

      visitor.before_edge_treatment(*this, gate);

      m_queue.pop();

      if(m_tr.is_infinite(nh))
      {
        nh->set_label(Face_label::OUTSIDE);
#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif
        continue;
      }

      Point_2 steiner_point;
      if(compute_steiner_point(fh, nh, steiner_point))
      {
//        std::cout << CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset)
//                  << " vs " << 1e-2 * m_offset << std::endl;
        CGAL_assertion(CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset) <= 1e-2 * m_offset);

        // locate faces that are going to be destroyed and remove their edges from the queue
        int li = 0;
        Locate_type lt;
        const Face_handle conflict_face = m_tr.locate(steiner_point, lt, li, nh);
        CGAL_assertion(lt != Triangulation::VERTEX);

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
        // Using small vectors like in Triangulation_2 does not bring any runtime improvement
        std::vector<Edge> boundary_edges;
        std::vector<Face_handle> conflict_zone;
        boundary_edges.reserve(16);
        conflict_zone.reserve(16);

        m_tr.get_conflicts_and_boundary(steiner_point,
                                        std::back_inserter(conflict_zone),
                                        std::back_inserter(boundary_edges),
                                        conflict_face);

        // Purge the queue of edges that will be deleted/modified by the Steiner point insertion,
        // and which might have been gates
        for(const Face_handle& cch : conflict_zone)
        {
          for(int i=0; i<3; ++i)
          {
            const Edge cf = std::make_pair(cch, i);
            if(m_queue.contains_with_bounds_check(Gate(cf)))
              m_queue.erase(Gate(cf));
          }
        }

        for(const Edge& e : boundary_edges)
        {
          const Edge me = m_tr.mirror_edge(e); // boundary edges have incident faces in the CZ
          if(m_queue.contains_with_bounds_check(Gate(me)))
            m_queue.erase(Gate(me));
        }
#endif

        visitor.before_Steiner_point_insertion(*this, steiner_point);

        // Actual insertion of the Steiner point
        // We could use TDS functions to avoid recomputing the conflict zone, but in practice
        // it does not bring any runtime improvements
        Vertex_handle vh = m_tr.insert(steiner_point, lt, conflict_face, li);
        vh->type() = AW2i::Vertex_type:: DEFAULT;

        visitor.after_Steiner_point_insertion(*this, vh);

        Face_circulator new_fh = m_tr.incident_faces(vh), done(new_fh);
        do
        {
          // std::cout << "new face has time stamp " << new_fh->time_stamp() << std::endl;
          new_fh->set_label(m_tr.is_infinite(new_fh) ? Face_label::OUTSIDE : Face_label::INSIDE);
        }
        while(++new_fh != done);

        // Push all new boundary edges to the queue.
        // It is not performed by looking at the edges on the boundary of the conflict zones
        // because we need to handle internal edges, infinite edges, and also more subtle changes
        // such as a new face being marked inside which now creates a boundary
        // with its incident "outside" flagged face.
        do
        {
          for(int i=0; i<3; ++i)
          {
            if(m_tr.is_infinite(new_fh, i))
              continue;

            const Face_handle new_nh = new_fh->neighbor(i);
            if(new_nh->label() == new_fh->label()) // not on a boundary
              continue;

            const Edge boundary_f = std::make_pair(new_fh, i);
            if(new_fh->is_outside())
              push_edge(boundary_f);
            else
              push_edge(m_tr.mirror_edge(boundary_f));
          }
        }
        while(++new_fh != done);
      }
      else // no need for a Steiner point, carve through and continue
      {
        nh->set_label(Face_label::OUTSIDE);
#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
        nh->increment_erase_counter();
#endif

        // for each finite edge of neighbor, push it to the queue
        const int mi = m_tr.mirror_index(fh, s);
        for(int i=1; i<3; ++i)
        {
          const Edge neighbor_f = std::make_pair(nh, (mi+i)%3);
          push_edge(neighbor_f);
        }
      }
    } // while(!queue.empty())

    visitor.on_flood_fill_end(*this);

    // Check that no useful edge has been ignored
    CGAL_postcondition_code(for(auto eit=m_tr.finite_edges_begin(), eend=m_tr.finite_edges_end(); eit!=eend; ++eit) {)
    CGAL_postcondition_code(  Face_handle fh = eit->first; Face_handle nh = eit->first->neighbor(eit->second); )
    CGAL_postcondition_code(  if(fh->label() == nh->label()) continue;)
    CGAL_postcondition_code(  Edge e = *eit;)
    CGAL_postcondition_code(  if(fh->is_inside()) e = m_tr.mirror_edge(e);)
    CGAL_postcondition(       edge_status(e) == Edge_status::IRRELEVANT);
    CGAL_postcondition_code(})

    return true;
  }

  // Any outside face that isn't reachable from infinity is a cavity that can be discarded.
  std::size_t purge_inner_connected_components()
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "> Purge inner islands..." << std::endl;

    std::size_t pre_counter = 0;
    for(Face_handle fh : m_tr.all_face_handles())
      if(fh->is_outside())
        ++pre_counter;
    std::cout << pre_counter << " / " << m_tr.all_face_handles().size() << " (pre purge)" << std::endl;
#endif

    std::size_t label_change_counter = 0;

    std::stack<Face_handle> faces_to_visit;

    if(!m_seeds.empty())
    {
      for(const Point_2& seed : m_seeds)
      {
        Face_handle fh = m_tr.locate(seed);

        if(!fh->is_outside())
        {
          std::cerr << "Warning: face containing seed is not outside?!" << std::endl;
          continue;
        }

        faces_to_visit.push(fh);
      }
    }
    else // typical flooding from outside
    {
      faces_to_visit.push(m_tr.infinite_vertex()->face());
    }

    while(!faces_to_visit.empty())
    {
      Face_handle curr_c = faces_to_visit.top();
      faces_to_visit.pop();

      CGAL_assertion(curr_c->is_outside());

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      for(int j=0; j<3; ++j)
      {
        Face_handle neigh_c = curr_c->neighbor(j);
        if(neigh_c->tds_data().processed() || !neigh_c->is_outside())
          continue;

        faces_to_visit.push(neigh_c);
      }
    }

    for(Face_handle fh : m_tr.all_face_handles())
    {
      if(fh->tds_data().is_clear() && fh->is_outside())
      {
        fh->set_label(Face_label::INSIDE);
#ifndef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
        fh->increment_erase_counter();
#endif
        ++label_change_counter;
      }
    }

    // reset the conflict flags
    for(Face_handle fh : m_tr.all_face_handles())
      fh->tds_data().clear();

#ifdef CGAL_AW2_DEBUG
    std::size_t post_counter = 0;
    for(Face_handle fh : m_tr.all_face_handles())
      if(fh->is_outside())
        ++post_counter;
    std::cout << post_counter << " / " << m_tr.all_face_handles().size() << " (pre purge)" << std::endl;

    std::cout << label_change_counter << " label changes" << std::endl;
#endif

    return label_change_counter;
  }

private:
  bool is_non_manifold(Vertex_handle v) const
  {
    CGAL_precondition(!m_tr.is_infinite(v));

    bool is_non_manifold = false;

    // Flood one inside and outside CC within the face umbrella of the vertex.
    // Process both an inside and an outside CC to also detect edge pinching.
    // If there are still unprocessed afterwards, there is a non-manifoldness issue.
    //
    // Squat the conflict face data to mark visits
    Face_handle inside_start = Face_handle();
    Face_handle outside_start = Face_handle();

    Face_circulator inc_f = m_tr.incident_faces(v), done(inc_f);
    do
    {
      inc_f->tds_data().clear();
      if(inc_f->is_outside())
        outside_start = inc_f;
      else if(inside_start == Face_handle())
        inside_start = inc_f; // can be "inside" or "manifold"
    }
    while(++inc_f != done);

    // fully inside / outside
    if(inside_start == Face_handle() || outside_start == Face_handle())
      return false;

    std::stack<Face_handle> faces_to_visit;
    faces_to_visit.push(inside_start);
    faces_to_visit.push(outside_start);
    while(!faces_to_visit.empty())
    {
      Face_handle curr_c = faces_to_visit.top();
      faces_to_visit.pop();

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      int v_pos = -1;
      CGAL_assertion_code(bool res = ) curr_c->has_vertex(v, v_pos);
      CGAL_assertion(res);

      for(int j=0; j<3; ++j)
      {
        if(j == v_pos)
          continue;

        Face_handle neigh_c = curr_c->neighbor(j);
        CGAL_assertion(neigh_c->has_vertex(v));

        if(neigh_c->tds_data().processed() ||
           is_on_outside_boundary(neigh_c, curr_c)) // do not cross the boundary
        {
          continue;
        }

        faces_to_visit.push(neigh_c);
      }
    }

    do
    {
      if(inc_f->tds_data().is_clear())
        is_non_manifold = true;

      inc_f->tds_data().clear();
    }
    while(++inc_f != done);

    return is_non_manifold;
  }

  bool is_non_manifold(Face_handle f) const
  {
    CGAL_precondition(!m_tr.is_infinite(f));

    for(int i=0; i<3; ++i)
    {
      Vertex_handle v = f->vertex(i);
      if(is_non_manifold(v))
        return true;
    }

    return false;
  }

  bool is_manifold() const
  {
    // Not the best complexity, but it's not important: this function is purely for information
    // Better complexity --> see PMP::non_manifold_vertices + throw
    for(const Vertex_handle v : m_tr.finite_vertex_handles())
      if(is_non_manifold(v))
        return true;

    return false;
  }

  // Remove bbox vertices, if they are not necessary (i.e., no "inside" incident face)
  // This is to try and avoid having long tets with bbox vertices being tagged "inside" as part
  // of the manifold re-tagging
  bool remove_bbox_vertices()
  {
    bool do_remove = true;
    auto vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;

      Face_circulator f = m_tr.incident_faces(v), done(f);
      do
      {
        if(!f->is_outside())
        {
          do_remove = false;
          break;
        }
      }
      while(++f != done);

      if(!do_remove)
        break;
    }

    std::cout << "Removing bbox vertices? " << do_remove << std::endl;
    if(!do_remove)
      return false;

    vit = m_tr.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;
      m_tr.remove(v);
    }

    return true;
  }

public:
  // Not the best complexity, but it's very cheap compared to the rest of the algorithm.
  void make_manifold()
  {
#ifdef CGAL_AW2_DEBUG
    std::cout << "> Make manifold..." << std::endl;

    auto wrap_area = [&]()
    {
      FT a = 0;
      for(Face_handle fh : m_tr.finite_face_handles())
        if(!fh->is_outside())
          a += area(m_tr.point(fh, 0), m_tr.point(fh, 1), m_tr.point(fh, 2));

      return a;
    };

 #ifdef CGAL_AW2_DEBUG_DUMP_INTERMEDIATE_WRAPS
    dump_triangulation("carved_wrap.off");
 #endif

    FT base_a = wrap_area();
    if(!is_positive(base_a))
      std::cerr << "Warning: wrap with non-positive area?" << std::endl;
#endif

    // This ends up more harmful than useful after the priority queue has been introduced since
    // it usually results in a lot of flat faces into the triangulation, which then get added
    // to the mesh during manifoldness fixing.
//    remove_bbox_vertices();

    std::stack<Vertex_handle> non_manifold_vertices; // @todo sort somehow?
    for(Vertex_handle v : m_tr.finite_vertex_handles())
    {
      if(is_non_manifold(v))
      {
#ifdef CGAL_AW2_DEBUG_MANIFOLDNESS_PP
        std::cout << v->point() << " is non-manifold" << std::endl;
#endif
        non_manifold_vertices.push(v);
      }
    }

    // Some lambdas for the comparer
    auto has_scaffolding_vertex = [](Face_handle f) -> bool
    {
      for(int i=0; i<3; ++i)
      {
        if(f->vertex(i)->type() == AW2i::Vertex_type::BBOX_VERTEX ||
           f->vertex(i)->type() == AW2i::Vertex_type::SEED_VERTEX)
        {
          return true;
        }
      }

      return false;
    };

    // This originally seemed like a good idea, but in the end it can have strong cascading issues,
    // whereas some faces with much smaller area could have solved the non-manifoldness.
//    auto is_on_boundary = [](Face_handle f, int i) -> bool
//    {
//      return is_on_outside_boundary(f, f->neighbor(i));
//    };
//
//    auto count_boundary_edges = [&](Face_handle f, Vertex_handle v) -> int
//    {
//      const int v_index_in_f = f->index(v);
//      int boundary_edges = 0;
//      for(int i=0; i<2; ++i) // also take into account the opposite edge?
//      {
//        if(i == v_index_in_f)
//          continue;
//
//        boundary_edges += is_on_boundary(f, i);
//      }
//
//      return boundary_edges;
//    };

    // Experimentally, longest edge works better
//    auto sq_circumradius = [&](Face_handle f) -> FT
//    {
//      const Point_2& cc = circumcenter(f);
//      return geom_traits().compute_squared_distance_2_object()(m_tr.point(f, 0), cc);
//    };

    // The reasoning behind using longest edge rather than area is that we want to avoid
    // spikes (which would have a small area), and can often happen since we do not spend
    // any care on the quality of triangles.
    auto sq_longest_edge = [&](Face_handle f) -> FT
    {
      return (std::max)({ squared_distance(m_tr.point(f, 0), m_tr.point(f, 1)),
                          squared_distance(m_tr.point(f, 1), m_tr.point(f, 2)),
                          squared_distance(m_tr.point(f, 2), m_tr.point(f, 0)) });
    };

#ifdef CGAL_AW2_DEBUG_MANIFOLDNESS_PP
    std::cout << non_manifold_vertices.size() << " initial NMV" << std::endl;
#endif

    while(!non_manifold_vertices.empty())
    {
#ifdef CGAL_AW2_DEBUG_MANIFOLDNESS_PP
      std::cout << non_manifold_vertices.size() << " NMV in queue" << std::endl;
#endif

      Vertex_handle v = non_manifold_vertices.top();
      non_manifold_vertices.pop();

      if(!is_non_manifold(v))
        continue;

      // Prioritize:
      // - faces without bbox vertices
      // - small faces when equal number of boundary edges
      //
      // Note that these are properties that do not depend on the face labels, and so we only need
      // to sort once. However, if a criterion such as the number of incident inside faces were added,
      // we would need to sort after each modification of "inside"/"outside" labels.
      auto comparer = [&](Face_handle l, Face_handle r) -> bool
      {
        CGAL_precondition(!m_tr.is_infinite(l) && !m_tr.is_infinite(r));

        if(has_scaffolding_vertex(l) != has_scaffolding_vertex(r))
          return has_scaffolding_vertex(r);

        return sq_longest_edge(l) < sq_longest_edge(r);
      };

      std::vector<Face_handle> inc_faces;
      inc_faces.reserve(16);
      Face_circulator inc_f = m_tr.incident_faces(v), done(inc_f);
      do
      {
        inc_faces.push_back(inc_f);
      }
      while(++inc_f != done);

      std::vector<Face_handle> finite_outside_inc_faces;
      finite_outside_inc_faces.reserve(16);
      std::copy_if(inc_faces.begin(), inc_faces.end(), std::back_inserter(finite_outside_inc_faces),
                   [&](Face_handle f) -> bool { return !m_tr.is_infinite(f) && f->is_outside(); });

      // 'std::stable_sort' to have determinism without having to write something like:
      //     if(longest_edge(l) == longest_edge(r)) return ...
      // in the comparer. It's almost always a small range, so the extra cost does not matter.
      std::stable_sort(finite_outside_inc_faces.begin(), finite_outside_inc_faces.end(), comparer);

      for(Face_handle inc_f : finite_outside_inc_faces)
      {
        CGAL_assertion(!m_tr.is_infinite(inc_f) && inc_f->is_outside());

        // This is where new material is added
        inc_f->set_label(Face_label::MANIFOLD);

#ifdef CGAL_AW2_DEBUG_DUMP_EVERY_STEP
        static int i = 0;
        std::string step_name = "results/steps_manifold/step" + std::to_string(static_cast<int>(i++)) + ".off";
        dump_triangulation(step_name);
#endif

        // @speed could update the manifold status while tagging
        if(!is_non_manifold(v))
          break;
      }

      CGAL_assertion(!is_non_manifold(v));

      // Check if the new material has not created a non-manifold configuration.
      // @speed this could be done on only the vertices of faces whose labels have changed.
      Vertex_circulator adj_v = m_tr.incident_vertices(v), adj_done(adj_v);
      do
      {
        if(m_tr.is_infinite(adj_v))
          continue;
        if(is_non_manifold(adj_v))
          non_manifold_vertices.push(adj_v);
      }
      while(++adj_v != adj_done);
    }

    CGAL_assertion_code(for(Vertex_handle v : m_tr.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

#ifdef CGAL_AW2_DEBUG
    std::size_t nm_faces_counter = 0;
    for(Face_handle fh : m_tr.all_face_handles())
      if(fh->label() == Face_label::MANIFOLD)
        ++nm_faces_counter;
    std::cout << "Number of added faces: " << nm_faces_counter << std::endl;

    if(!is_zero(base_a))
    {
      const FT manifold_a = wrap_area();
      const FT ratio = manifold_a / base_a;

      std::cout << "Areas post-manifoldness fix:\n"
                << "before: " << base_a << "\n"
                << "after:  " << manifold_a << "\n"
                << "ratio:  " << ratio << std::endl;
      if(ratio > 1.1) // more than 10% extra volume
        std::cerr << "Warning: large increase of volume after manifoldness resolution" << std::endl;
    }
#endif
  }

private:
  void check_queue_sanity()
  {
    std::cout << "\t~~~ Check queue sanity ~~~" << std::endl;

    std::vector<Gate> queue_gates;
    Gate previous_top_gate = m_queue.top();
    while(!m_queue.empty())
    {
      const Gate& current_gate = m_queue.top();
      queue_gates.push_back(current_gate);
      const Edge& current_f = current_gate.edge();
      const Face_handle fh = current_f.first;
      const int id = current_f.second;
      const Point_2& p0 = m_tr.point(fh, Triangulation::ccw(id));
      const Point_2& p1 = m_tr.point(fh, Triangulation::cw(id));

      std::cout << "Edge with VID " << get(Gate_ID_PM<Triangulation>(), current_gate) << "\n";
      std::cout << "\t" << p0 << "\n\t" << p1 << "\n";

#ifdef CGAL_AW2_USE_SORTED_PRIORITY_QUEUE
      std::cout << "  Permissiveness: " << current_gate.is_permissive_edge() << "\n";
      std::cout << "  SQR: " << geom_traits().compute_squared_radius_2_object()(p0, p1) << "\n";
      std::cout << "  Priority " << current_gate.priority() << std::endl;

      if(Less_gate()(current_gate, previous_top_gate))
        std::cerr << "Error: current gate has higher priority than the previous top" << std::endl;

      previous_top_gate = current_gate;
#else
      if(current_gate.is_zombie())
        std::cout << "Gate is a zombie!" << std::endl;
#endif

      m_queue.pop();
    }

    std::cout << "\t~~~ End queue sanity check ~~~" << std::endl;

    // Rebuild
    CGAL_assertion(m_queue.empty());
    for(const auto& g : queue_gates)
      m_queue.push(g); // no need for a resize here since the vector capacity is unchanged
  }

  void dump_triangulation(const std::string filename)
  {
    std::stringstream vertices_ss;
    vertices_ss.precision(17);
    std::stringstream faces_ss;
    faces_ss.precision(17);

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;
    std::size_t nf = 0;

    // First collect all vertices from finite faces
    for(auto fit = m_tr.finite_faces_begin(); fit != m_tr.finite_faces_end(); ++fit)
    {
      for(int i = 0; i < 3; ++i)
      {
        Vertex_handle v = fit->vertex(i);
        if(vertex_to_id.find(v) == vertex_to_id.end())
        {
          vertex_to_id[v] = nv++;
          vertices_ss << m_tr.point(v) << " 0" << "\n";
        }
      }
    }

    // Then output faces with colors
    for(auto fit = m_tr.finite_faces_begin(); fit != m_tr.finite_faces_end(); ++fit)
    {
      faces_ss << "3";
      for(int i = 0; i < 3; ++i)
        faces_ss << " " << vertex_to_id[fit->vertex(i)];

      // Color format: r g b
      if(fit->is_inside())
        faces_ss << " 51 153 51"; // green for inside
      else if(fit->label() == Face_label::MANIFOLD)
        faces_ss << " 246 211 45"; // yellow for manifold
      else
        faces_ss << " 192 28 40"; // red for outside

      faces_ss << "\n";
      ++nf;
    }

    std::ofstream out(filename.c_str());
    out << "COFF\n" << nv << " " << nf << " 0\n";
    out << vertices_ss.str() << faces_ss.str() << std::endl;
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_ALPHA_WRAP_2_H
