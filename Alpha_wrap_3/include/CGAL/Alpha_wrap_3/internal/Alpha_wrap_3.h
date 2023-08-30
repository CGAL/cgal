// Copyright (c) 2019-2022 Google LLC (USA).
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
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H

#ifdef CGAL_AW3_DEBUG_PP
 #ifndef CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG
  #define CGAL_AW3_DEBUG_INITIALIZATION
  #define CGAL_AW3_DEBUG_STEINER_COMPUTATION
  #define CGAL_AW3_DEBUG_QUEUE
  #define CGAL_AW3_DEBUG_FACET_STATUS
  #define CGAL_AW3_DEBUG_MANIFOLDNESS
 #endif
#endif

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/gate_priority_queue.h>
#include <CGAL/Alpha_wrap_3/internal/geometry_utils.h>
#include <CGAL/Alpha_wrap_3/internal/oracles.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Modifiable_priority_queue.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h> // only if non-manifoldness is not treated
#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_with_circumcenter_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Robust_weighted_circumcenter_filtered_traits_3.h>

#include <array>
#include <algorithm>
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
namespace Alpha_wraps_3 {
namespace internal {

template <typename Cb>
class Cell_base_with_timestamp
  : public Cb
{
  std::size_t time_stamp_;

public:
  template <typename... Args>
  Cell_base_with_timestamp(const Args&... args) : Cb(args...), time_stamp_(-1) { }

  Cell_base_with_timestamp(const Cell_base_with_timestamp& other) : Cb(other), time_stamp_(other.time_stamp_) { }

  typedef CGAL::Tag_true Has_timestamp;

  std::size_t time_stamp() const { return time_stamp_; }
  void set_time_stamp(const std::size_t& ts) { time_stamp_ = ts; }

  template <class TDS>
  struct Rebind_TDS
  {
    typedef typename Cb::template Rebind_TDS<TDS>::Other Cb2;
    typedef Cell_base_with_timestamp<Cb2> Other;
  };
};

struct Wrapping_default_visitor
{
  Wrapping_default_visitor() { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_begin(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_flood_fill_begin(const AlphaWrapper&) { }

  template <typename AlphaWrapper, typename Gate>
  void before_facet_treatment(const AlphaWrapper&, const Gate&) { }

  template <typename Wrapper, typename Point>
  void before_Steiner_point_insertion(const Wrapper&, const Point&) { }

  template <typename Wrapper, typename VertexHandle>
  void after_Steiner_point_insertion(const Wrapper&, VertexHandle) { }

  template <typename AlphaWrapper>
  void on_flood_fill_end(const AlphaWrapper&) { }

  template <typename AlphaWrapper>
  void on_alpha_wrapping_end(const AlphaWrapper&) { };
};

template <typename Oracle>
class Alpha_wrap_3
{
  using Base_GT = typename Oracle::Geom_traits;
  using Geom_traits = Robust_circumcenter_filtered_traits_3<Base_GT>;

  using FT = typename Geom_traits::FT;
  using Point_3 = typename Geom_traits::Point_3;
  using Vector_3 = typename Geom_traits::Vector_3;
  using Ball_3 = typename Geom_traits::Ball_3;
  using Iso_cuboid_3 = typename Geom_traits::Iso_cuboid_3;

  using SC = Simple_cartesian<double>;
  using SC_Point_3 = SC::Point_3;
  using SC_Vector_3 = SC::Vector_3;
  using SC_Iso_cuboid_3 = SC::Iso_cuboid_3;
  using SC2GT = Cartesian_converter<SC, Geom_traits>;

  struct Cell_info
  {
    bool is_outside = false;
  };

  enum Vertex_info
  {
    DEFAULT = 0,
    BBOX_VERTEX,
    SEED_VERTEX
  };

  using Vb = Triangulation_vertex_base_3<Geom_traits>;
  using Vbi = Triangulation_vertex_base_with_info_3<Vertex_info, Geom_traits, Vb>;
  using Cbb = Delaunay_triangulation_cell_base_3<Geom_traits>;
  using Cb = Delaunay_triangulation_cell_base_with_circumcenter_3<Geom_traits, Cbb>;
  using Cbi = Triangulation_cell_base_with_info_3<Cell_info, Geom_traits, Cb>;
  using Cbt = Cell_base_with_timestamp<Cbi>;
  using Tds = Triangulation_data_structure_3<Vbi, Cbt>;
  using Dt = Delaunay_triangulation_3<Geom_traits, Tds, Fast_location>;

  using Cell_handle = typename Dt::Cell_handle;
  using Facet = typename Dt::Facet;
  using Vertex_handle = typename Dt::Vertex_handle;
  using Locate_type = typename Dt::Locate_type;

  using Gate = internal::Gate<Dt>;
  using Alpha_PQ = Modifiable_priority_queue<Gate, Less_gate, Gate_ID_PM<Dt>, CGAL_BOOST_PAIRING_HEAP>;

protected:
  const Oracle m_oracle;
  SC_Iso_cuboid_3 m_bbox;

  FT m_alpha, m_sq_alpha;
  FT m_offset, m_sq_offset;

  Dt m_dt;
  Alpha_PQ m_queue;

public:
  // Main constructor
  Alpha_wrap_3(const Oracle& oracle)
    :
      m_oracle(oracle),
      m_dt(Geom_traits(oracle.geom_traits())),
      // used to set up the initial MPQ, use some arbitrary not-too-small value
      m_queue(4096)
  {
    // Due to the Steiner point computation being a dichotomy, the algorithm is inherently inexact
    // and passing exact kernels is explicitly disabled to ensure no misunderstanding.
    static_assert(std::is_floating_point<FT>::value);
  }

public:
  const Geom_traits& geom_traits() const { return m_dt.geom_traits(); }
  Dt& triangulation() { return m_dt; }
  const Dt& triangulation() const { return m_dt; }
  const Alpha_PQ& queue() const { return m_queue; }

  double default_alpha() const
  {
    const Bbox_3 bbox = m_oracle.bbox();
    const double diag_length = std::sqrt(square(bbox.xmax() - bbox.xmin()) +
                                         square(bbox.ymax() - bbox.ymin()) +
                                         square(bbox.zmax() - bbox.zmin()));

    return diag_length / 20.;
  }

private:
  const Point_3& circumcenter(const Cell_handle c) const
  {
    // We only cross an infinite facet once, so this isn't going to be recomputed many times
    if(m_dt.is_infinite(c))
    {
      const int inf_index = c->index(m_dt.infinite_vertex());
      c->set_circumcenter(
            geom_traits().construct_circumcenter_3_object()(m_dt.point(c, (inf_index+1)&3),
                                                            m_dt.point(c, (inf_index+2)&3),
                                                            m_dt.point(c, (inf_index+3)&3)));
    }

    return c->circumcenter(geom_traits());
  }

public:
  template <typename OutputMesh,
            typename InputNamedParameters = parameters::Default_named_parameters,
            typename OutputNamedParameters = parameters::Default_named_parameters>
  void operator()(const double alpha, // = default_alpha()
                  const double offset, // = alpha / 30.
                  OutputMesh& output_mesh,
                  const InputNamedParameters& in_np = parameters::default_values(),
                  const OutputNamedParameters& out_np = parameters::default_values())
  {
    namespace PMP = Polygon_mesh_processing;

    using parameters::get_parameter;
    using parameters::get_parameter_reference;
    using parameters::choose_parameter;

    using OVPM = typename CGAL::GetVertexPointMap<OutputMesh, OutputNamedParameters>::type;
    OVPM ovpm = choose_parameter(get_parameter(out_np, internal_np::vertex_point),
                                 get_property_map(vertex_point, output_mesh));

    typedef typename internal_np::Lookup_named_param_def <
      internal_np::visitor_t,
      InputNamedParameters,
      Wrapping_default_visitor // default
    >::reference                                                                 Visitor;

    Wrapping_default_visitor default_visitor;
    Visitor visitor = choose_parameter(get_parameter_reference(in_np, internal_np::visitor), default_visitor);

    std::vector<Point_3> no_seeds;
    using Seeds = typename internal_np::Lookup_named_param_def<
                    internal_np::seed_points_t, InputNamedParameters, std::vector<Point_3> >::reference;
    Seeds seeds = choose_parameter(get_parameter_reference(in_np, internal_np::seed_points), no_seeds);

    const bool do_enforce_manifoldness = choose_parameter(get_parameter(in_np, internal_np::do_enforce_manifoldness), true);

#ifdef CGAL_AW3_TIMER
    CGAL::Real_timer t;
    t.start();
#endif

    visitor.on_alpha_wrapping_begin(*this);

    if(!initialize(alpha, offset, seeds))
      return;

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
    extract_surface(output_mesh, ovpm, true /*tolerate non manifoldness*/);
    CGAL::IO::write_polygon_mesh("initial_cavities.off", output_mesh,
                                 CGAL::parameters::vertex_point_map(ovpm).stream_precision(17));
#endif

    alpha_flood_fill(visitor);

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Flood filling took: " << t.time() << " s." << std::endl;
#endif

    if(do_enforce_manifoldness)
    {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS
      std::cout << "> Make manifold..." << std::endl;

      extract_surface(output_mesh, ovpm, true /*tolerate non manifoldness*/);

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
      dump_triangulation_faces("intermediate_dt3.off", false /*only_boundary_faces*/);
      IO::write_polygon_mesh("intermediate.off", output_mesh,
                             CGAL::parameters::vertex_point_map(ovpm).stream_precision(17));
#endif

      FT base_vol = 0;
      if(is_closed(output_mesh)) // might not be due to manifoldness
        base_vol = PMP::volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
      else
        std::cerr << "Warning: couldn't compute volume before manifoldness fixes (mesh is not closed)" << std::endl;
#endif

#ifdef CGAL_AW3_TIMER
    t.reset();
    t.start();
#endif

      make_manifold();

#ifdef CGAL_AW3_TIMER
      t.stop();
      std::cout << "\nManifoldness post-processing took: " << t.time() << " s." << std::endl;
#endif

#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS
      if(!is_zero(base_vol))
      {
        extract_surface(output_mesh, ovpm, false /*do not tolerate non-manifoldness*/);

        const FT manifold_vol = PMP::volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
        const FT ratio = manifold_vol / base_vol;

        std::cout << "Volumes post-manifoldness fix:\n"
                  << "before: " << base_vol << "\n"
                  << "after:  " << manifold_vol << "\n"
                  << "ratio:  " << ratio << std::endl;
        if(ratio > 1.1) // more than 10% extra volume
          std::cerr << "Warning: large increase of volume after manifoldness resolution" << std::endl;
      }
#endif
    } // do_enforce_manifoldness

#ifdef CGAL_AW3_TIMER
    t.reset();
    t.start();
#endif

    extract_surface(output_mesh, ovpm, !do_enforce_manifoldness);

#ifdef CGAL_AW3_TIMER
    t.stop();
    std::cout << "Surface extraction took: " << t.time() << " s." << std::endl;
#endif

#ifdef CGAL_AW3_DEBUG
    std::cout << "Alpha wrap vertices:  " << vertices(output_mesh).size() << std::endl;
    std::cout << "Alpha wrap faces:     " << faces(output_mesh).size() << std::endl;

 #ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
    IO::write_polygon_mesh("final.off", output_mesh, CGAL::parameters::stream_precision(17));
    dump_triangulation_faces("final_dt3.off", false /*only_boundary_faces*/);
 #endif
#endif

    visitor.on_alpha_wrapping_end(*this);
  }

  // Convenience overloads
  template <typename OutputMesh>
  void operator()(const double alpha,
                  OutputMesh& output_mesh)
  {
    return operator()(alpha, alpha / 30. /*offset*/, output_mesh);
  }

  template <typename OutputMesh>
  void operator()(OutputMesh& output_mesh)
  {
    return operator()(default_alpha(), output_mesh);
  }

  // This function is public only because it is used in the tests
  SC_Iso_cuboid_3 construct_bbox(const double offset)
  {
    // Input axis-aligned bounding box
    SC_Iso_cuboid_3 bbox = m_oracle.bbox();
    const SC_Point_3 bbox_centroid = midpoint((bbox.min)(), (bbox.max)());

    // Scale a bit to create the initial points not too close to the input
    double scaling = 1.2;
    CGAL::Aff_transformation_3<SC> scale(SCALING, scaling);
    bbox = SC_Iso_cuboid_3(scale.transform((bbox.min)()), scale.transform((bbox.max)()));

    // Translate bbox back to initial centroid
    const SC_Point_3 bbox_transformed_centroid = midpoint((bbox.min)(), (bbox.max)());
    const SC_Vector_3 diff_centroid = bbox_centroid - bbox_transformed_centroid;
    CGAL::Aff_transformation_3<SC> centroid_translate(TRANSLATION, diff_centroid);
    bbox = SC_Iso_cuboid_3(centroid_translate.transform((bbox.min)()),
                           centroid_translate.transform((bbox.max)()));

    // Add the offset
    SC_Vector_3 offset_ext = std::sqrt(3.) * offset * SC_Vector_3(1, 1, 1);
    CGAL::Aff_transformation_3<SC> translate_m(TRANSLATION, - offset_ext);
    CGAL::Aff_transformation_3<SC> translate_M(TRANSLATION,   offset_ext);
    bbox = SC_Iso_cuboid_3(translate_m.transform((bbox.min)()), translate_M.transform((bbox.max)()));

    return bbox;
  }

private:
  // Adjust the bbox & insert its corners to construct the starting triangulation
  void insert_bbox_corners()
  {
    m_bbox = construct_bbox(CGAL::to_double(m_offset));

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "Insert Bbox vertices" << std::endl;
#endif

    // insert in dt the eight corner vertices of the input loose bounding box
    for(int i=0; i<8; ++i)
    {
      const Point_3 bp = SC2GT()(m_bbox.vertex(i));
      Vertex_handle bv = m_dt.insert(bp);
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "\t" << bp << std::endl;
#endif
      bv->info() = BBOX_VERTEX;
    }
  }

  // Two criteria:
  // - Cells that are intersecting the input are inside
  // - Cells whose circumcenter is in the offset volume are inside: this is because
  // we need to have outside cell circumcenters outside of the volume to ensure
  // that the refinement point is separated from the existing point set.
  bool cavity_cell_outside_tag(const Cell_handle ch)
  {
    CGAL_precondition(!m_dt.is_infinite(ch));

    const Tetrahedron_with_outside_info<Geom_traits> tet(ch, geom_traits());
    if(m_oracle.do_intersect(tet))
      return false;

    const Point_3& ch_cc = circumcenter(ch);
    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);
    const bool is_cc_in_offset = m_oracle.do_intersect(ch_cc_offset_ball);

    return !is_cc_in_offset;
  }

  // Create a cavity using seeds rather than starting from the infinity.
  //
  // For each seed, insert the seeds and its neighboring vertices on a regular icosahedron.
  // The idea behind an icosahedron rather than e.g. a simple tetrahedron is twofold:
  // - Improve the odds of both starting the algorithm, but also of not immediately stopping,
  //   which could happen if a conflict zone included all the initial cavity when using a simple cavity.
  // - Base the cavity diameter on alpha, and allow its cells to intersect the input. If a single
  //   tetrahedron is used, one is forced to make it small-enough such that it does not intersect
  //   the input, and then it could be forced to be smaller than 'alpha' and the algorithm cannot start,
  //   for example if the seed is close to the offset. If we have many tetrahedra, this is far less
  //   likely to happen.
  //
  // For an edge length a, the radius of an icosahedron is
  //   r = 0.5 * a * sqrt(phi * sqrt(5)), phi being the golden ratio
  // which yields
  //   r = a * sin(2pi / 5) =~ 0.95105651629515353 * a
  // Faces of an icosahedron are equilateral triangles with size a and circumradius a / sqrt(3)
  // Since we want faces of the icosahedron to be traversable, we want a such that
  //   a / sqrt(3) > alpha
  // Hence r such that
  //   r / (sqrt(3) * sin(2pi/5)) > alpha
  //   r > alpha * sqrt(3) * sin(2pi / 5) ~= 1.6472782070926637 * alpha
  //
  // Furthermore, the triangles between edges of the icosahedron and the center of the icosahedron
  // are not equilateral triangles since a is slightly bigger than r. They are
  // slightly flattened isocele triangles with base 'a' and the circumradius is smaller.
  // The circumradius is
  //   r_iso = r² / (2 * h) = r² / (2 * sqrt(r² - (a / 2)²))
  // Since r = a * sin(2pi / 5)
  //  r_iso = r² / (2 * sqrt(r² - (r / 2 * sin(2pi / 5))²))
  //        = r / (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²))
  // So we want r_iso > alpha, i.e.
  //   r > (2 * sqrt(1 - 1 / (2 * sin(2pi / 5))²)) * alpha ~= 1.7013016167040798 * alpha
  //
  // Another way is to simply make faces incident to the seed always traversable, and then
  // we only have to ensure faces opposite of the seed are traversable (i.e., radius ~= 1.65 * alpha)
  template <typename SeedRange>
  bool initialize_with_cavities(const SeedRange& seeds)
  {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << "> Dig cavities" << std::endl;
    std::cout << seeds.size() << " seed(s)" << std::endl;
#endif

    CGAL_precondition(!seeds.empty());

    // Get a double value approximating the scaling factors
//    std::cout << sqrt(3) * sin(2pi / 5) << std::endl;
//    std::cout << (2. * std::sqrt(1. - 1. / square(2 * std::sin(2 * CGAL_PI / 5)))) << std::endl;

    Iso_cuboid_3 bbox = SC2GT()(m_bbox);

    std::vector<Vertex_handle> seed_vs;
    for(const Point_3& seed_p : seeds)
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cout << "Initialize from seed " << seed_p << std::endl;
#endif

      if(bbox.has_on_unbounded_side(seed_p))
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is outside the bounding box" << std::endl;
#endif
        continue;
      }

      // get the closest point on the input
      const Point_3 closest_pt = m_oracle.closest_point(seed_p);
      const FT sq_d_to_closest = geom_traits().compute_squared_distance_3_object()(seed_p, closest_pt);

      if(sq_d_to_closest < m_sq_offset)
      {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
        std::cerr << "Warning: seed " << seed_p << " is in the offset" << std::endl;
#endif
        continue;
      }

      // Mark the seeds and icosahedron vertices as "artificial vertices" such that the facets
      // incident to these vertices are always traversable regardless of their circumcenter.
      // This is done because otherwise some cavities can appear on the mesh: non-traversable facets
      // with two vertices on the offset, and the third being a deeper inside seed / ico_seed.
      // Making them traversable will cause more refinement than "alpha", but they will eventually
      // not appear anymore in the inside/outside boundary and the surface will look smoother.
      //
      // This problem only appears when the seed and icosahedron vertices are close to the offset surface,
      // which usually happens for large alpha values.

      Vertex_handle seed_v = m_dt.insert(seed_p);
      seed_v->info() = SEED_VERTEX;
      seed_vs.push_back(seed_v);

      // Icosahedron vertices (see also BGL::make_icosahedron())
      const Point_3 center = seed_p;
      const FT radius = 1.65 * m_alpha;
      const FT phi = (FT(1) + approximate_sqrt(FT(5))) / FT(2);
      const FT t = radius / approximate_sqrt(1 + square(phi));
      const FT t_phi = t * phi;

      std::array<Point_3, 12> ico_ps =
      {
        Point_3(center.x(), center.y() + t, center.z() + t_phi),
        Point_3(center.x(), center.y() + t, center.z() - t_phi),
        Point_3(center.x(), center.y() - t, center.z() + t_phi),
        Point_3(center.x(), center.y() - t, center.z() - t_phi),

        Point_3(center.x() + t, center.y() + t_phi, center.z()),
        Point_3(center.x() + t, center.y() - t_phi, center.z()),
        Point_3(center.x() - t, center.y() + t_phi, center.z()),
        Point_3(center.x() - t, center.y() - t_phi, center.z()),

        Point_3(center.x() + t_phi, center.y(), center.z() + t),
        Point_3(center.x() + t_phi, center.y(), center.z() - t),
        Point_3(center.x() - t_phi, center.y(), center.z() + t),
        Point_3(center.x() - t_phi, center.y(), center.z() - t)
      };

      for(const Point_3& seed_neighbor_p : ico_ps)
      {
#ifdef CGAL_AW3_DEBUG_PP
        std::cout << seed_neighbor_p << std::endl;
#endif
        if(bbox.has_on_unbounded_side(seed_neighbor_p))
          continue;

        Vertex_handle ico_v = m_dt.insert(seed_neighbor_p, seed_v /*hint*/);
        ico_v->info() = SEED_VERTEX;
      }
    }

    if(seed_vs.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Error: no acceptable seed was provided" << std::endl;
#endif
      return false;
    }

#ifdef CGAL_AW3_DEBUG_INITIALIZATION
    std::cout << m_dt.number_of_vertices() - 8 /*bbox*/ << " vertice(s) due to seeds" << std::endl;
#endif

    for(Vertex_handle seed_v : seed_vs)
    {
      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_dt.incident_cells(seed_v, std::back_inserter(inc_cells));
      for(Cell_handle ch : inc_cells)
        ch->info().is_outside = cavity_cell_outside_tag(ch);
    }

    // Might as well go through the full triangulation since only seeds should have been inserted
    for(Cell_handle ch : m_dt.all_cell_handles())
    {
      if(!ch->info().is_outside)
        continue;

      // When the algorithm starts from a manually dug hole, infinite cells are tagged "inside"
      CGAL_assertion(!m_dt.is_infinite(ch));

      for(int i=0; i<4; ++i)
        push_facet(std::make_pair(ch, i));
    }

    if(m_queue.empty())
    {
#ifdef CGAL_AW3_DEBUG_INITIALIZATION
      std::cerr << "Could not initialize the algorithm with these seeds, and alpha|offset values" << std::endl;
#endif
      return false;
    }

    return true;
  }

  // tag all infinite cells OUTSIDE and all finite cells INSIDE
  // init queue with all convex hull facets
  bool initialize_from_infinity()
  {
    for(Cell_handle ch : m_dt.all_cell_handles())
    {
      if(m_dt.is_infinite(ch))
      {
        ch->info().is_outside = true;
        const int inf_index = ch->index(m_dt.infinite_vertex());
        push_facet(std::make_pair(ch, inf_index));
      }
      else
      {
        ch->info().is_outside = false;
      }
    }

    return true;
  }

public:
  // Manifoldness is tolerated while debugging and extracting at intermediate states
  // Not the preferred way because it uses 3*nv storage
  template <typename OutputMesh, typename OVPM>
  void extract_possibly_non_manifold_surface(OutputMesh& output_mesh,
                                             OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract possibly non-manifold wrap... ()" << std::endl;
#endif

    clear(output_mesh);

    CGAL_assertion_code(for(auto cit=m_dt.finite_cells_begin(), cend=m_dt.finite_cells_end(); cit!=cend; ++cit))
    CGAL_assertion(cit->tds_data().is_clear());

    for(auto cit=m_dt.finite_cells_begin(), cend=m_dt.finite_cells_end(); cit!=cend; ++cit)
    {
      Cell_handle seed = cit;
      if(seed->info().is_outside || seed->tds_data().processed())
        continue;

      std::queue<Cell_handle> to_visit;
      to_visit.push(seed);

      std::vector<Point_3> points;
      std::vector<std::vector<size_t> > faces;
      std::size_t idx = 0;

      while(!to_visit.empty())
      {
        const Cell_handle cell = to_visit.front();
        CGAL_assertion(!cell->info().is_outside && !m_dt.is_infinite(cell));

        to_visit.pop();

        if(cell->tds_data().processed())
          continue;

        cell->tds_data().mark_processed();

        for(int fid=0; fid<4; ++fid)
        {
          const Cell_handle neighbor = cell->neighbor(fid);
          if(neighbor->info().is_outside)
          {
            // There shouldn't be any artificial vertex on the inside/outside boundary
            // (past initialization)
//            CGAL_assertion(cell->vertex((fid + 1)&3)->info() == DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 2)&3)->info() == DEFAULT);
//            CGAL_assertion(cell->vertex((fid + 3)&3)->info() == DEFAULT);

            points.push_back(m_dt.point(cell, Dt::vertex_triple_index(fid, 0)));
            points.push_back(m_dt.point(cell, Dt::vertex_triple_index(fid, 1)));
            points.push_back(m_dt.point(cell, Dt::vertex_triple_index(fid, 2)));
            faces.push_back({idx, idx + 1, idx + 2});
            idx += 3;
          }
          else
          {
            to_visit.push(neighbor);
          }
        }
      }

      PMP::duplicate_non_manifold_edges_in_polygon_soup(points, faces);

      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));
      PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                        CGAL::parameters::default_values(),
                                        CGAL::parameters::vertex_point_map(ovpm));

      PMP::stitch_borders(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
      CGAL_assertion(is_closed(output_mesh));
    }

    for(auto cit=m_dt.finite_cells_begin(), cend=m_dt.finite_cells_end(); cit!=cend; ++cit)
      cit->tds_data().clear();

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
  }

  template <typename OutputMesh, typename OVPM>
  void extract_manifold_surface(OutputMesh& output_mesh,
                                OVPM ovpm) const
  {
    namespace PMP = Polygon_mesh_processing;

#ifdef CGAL_AW3_DEBUG
    std::cout << "> Extract wrap... ()" << std::endl;
#endif

    CGAL_assertion_code(for(Vertex_handle v : m_dt.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));

    clear(output_mesh);

    // boundary faces to polygon soup
    std::vector<Point_3> points;
    std::vector<std::array<std::size_t, 3> > faces;

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;

    for(auto fit=m_dt.finite_facets_begin(), fend=m_dt.finite_facets_end(); fit!=fend; ++fit)
    {
      Facet f = *fit;
      if(!f.first->info().is_outside)
        f = m_dt.mirror_facet(f);

      const Cell_handle c = f.first;
      const int s = f.second;
      const Cell_handle nh = c->neighbor(s);
      if(c->info().is_outside == nh->info().is_outside)
        continue;

      std::array<std::size_t, 3> ids;
      for(int pos=0; pos<3; ++pos)
      {
        Vertex_handle vh = c->vertex(Dt::vertex_triple_index(s, pos));
        auto insertion_res = vertex_to_id.emplace(vh, nv);
        if(insertion_res.second) // successful insertion, never-seen-before vertex
        {
          points.push_back(m_dt.point(vh));
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      faces.emplace_back(std::array<std::size_t, 3>{ids[0], ids[1], ids[2]});
    }

    if(faces.empty())
      return;

    if(!PMP::is_polygon_soup_a_polygon_mesh(faces))
    {
      CGAL_warning_msg(false, "Could NOT extract mesh...");
      return;
    }

    PMP::polygon_soup_to_polygon_mesh(points, faces, output_mesh,
                                      CGAL::parameters::default_values(),
                                      CGAL::parameters::vertex_point_map(ovpm));

    CGAL_postcondition(!is_empty(output_mesh));
    CGAL_postcondition(is_valid_polygon_mesh(output_mesh));
    CGAL_postcondition(is_closed(output_mesh));

    PMP::orient_to_bound_a_volume(output_mesh, CGAL::parameters::vertex_point_map(ovpm));
  }

  template <typename OutputMesh, typename OVPM>
  void extract_surface(OutputMesh& output_mesh,
                       OVPM ovpm,
                       const bool tolerate_non_manifoldness = false) const
  {
    if(tolerate_non_manifoldness)
      extract_possibly_non_manifold_surface(output_mesh, ovpm);
    else
      extract_manifold_surface(output_mesh, ovpm);
  }

private:
  bool is_traversable(const Facet& f) const
  {
    return less_squared_radius_of_min_empty_sphere(m_sq_alpha, f, m_dt);
  }

  bool compute_steiner_point(const Cell_handle ch,
                             const Cell_handle neighbor,
                             Point_3& steiner_point) const
  {
    CGAL_precondition(!m_dt.is_infinite(neighbor));

    typename Geom_traits::Construct_ball_3 ball = geom_traits().construct_ball_3_object();
    typename Geom_traits::Construct_vector_3 vector = geom_traits().construct_vector_3_object();
    typename Geom_traits::Construct_translated_point_3 translate = geom_traits().construct_translated_point_3_object();
    typename Geom_traits::Construct_scaled_vector_3 scale = geom_traits().construct_scaled_vector_3_object();

    const Point_3& neighbor_cc = circumcenter(neighbor);
    const Ball_3 neighbor_cc_offset_ball = ball(neighbor_cc, m_sq_offset);
    const bool is_neighbor_cc_in_offset = m_oracle.do_intersect(neighbor_cc_offset_ball);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "Compute_steiner_point(" << &*ch << ", " << &*neighbor << ")" << std::endl;

    const Point_3& chc = circumcenter(ch);
    std::cout << "CH" << std::endl;
    std::cout << "\t" << ch->vertex(0)->point() << std::endl;
    std::cout << "\t" << ch->vertex(1)->point() << std::endl;
    std::cout << "\t" << ch->vertex(2)->point() << std::endl;
    std::cout << "\t" << ch->vertex(3)->point() << std::endl;
    std::cout << "cc: " << chc << std::endl;
    std::cout << "CC Distance to input: " << CGAL::squared_distance(chc, m_oracle.closest_point(chc)) << std::endl;

    std::cout << "NCH" << std::endl;
    std::cout << "\t" << neighbor->vertex(0)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(1)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(2)->point() << std::endl;
    std::cout << "\t" << neighbor->vertex(3)->point() << std::endl;
    std::cout << "ncc: " << neighbor_cc << std::endl;
    std::cout << "NC Distance to input: " << CGAL::squared_distance(neighbor_cc, m_oracle.closest_point(neighbor_cc)) << std::endl;
    std::cout << "is_neighbor_cc_in_offset = " << is_neighbor_cc_in_offset << std::endl;

    std::cout << "squared offset " << m_sq_offset << std::endl;
#endif

    // ch's circumcenter should not be within the offset volume
    CGAL_assertion_code(const Point_3& ch_cc = circumcenter(ch);)
    CGAL_assertion_code(const Ball_3 ch_cc_offset_ball = ball(ch_cc, m_sq_offset);)
    CGAL_assertion(!m_oracle.do_intersect(ch_cc_offset_ball));

    if(is_neighbor_cc_in_offset)
    {
      const Point_3& ch_cc = circumcenter(ch);

      // If the voronoi edge intersects the offset, the steiner point is the first intersection
      if(m_oracle.first_intersection(ch_cc, neighbor_cc, steiner_point, m_offset))
      {
#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
        std::cout << "Steiner found through first_intersection(): " << steiner_point << std::endl;
#endif
        return true;
      }
    }

    Tetrahedron_with_outside_info<Geom_traits> tet(neighbor, geom_traits());
    if(m_oracle.do_intersect(tet))
    {
      // steiner point is the closest point on input from cell centroid with offset
      const Point_3 closest_pt = m_oracle.closest_point(neighbor_cc);
      CGAL_assertion(closest_pt != neighbor_cc);

      Vector_3 unit = vector(closest_pt, neighbor_cc);

      // PMP::internal::normalize() requires sqrt()
      unit = scale(unit, m_offset / approximate_sqrt(geom_traits().compute_squared_length_3_object()(unit)));
      steiner_point = translate(closest_pt, unit);

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
      std::cout << "Steiner found through neighboring tet intersecting the input: " << steiner_point << std::endl;
      std::cout << "Closest point: " << closest_pt << std::endl;
      std::cout << "Direction: " << vector(closest_pt, neighbor_cc) << std::endl;
#endif

      return true;
    }

#ifdef CGAL_AW3_DEBUG_STEINER_COMPUTATION
    std::cout << "No Steiner point" << std::endl;
#endif

    return false;
  }

private:
  enum Facet_queue_status
  {
    IRRELEVANT = 0,
    ARTIFICIAL_FACET,
    TRAVERSABLE
  };

  // @speed some decent time may be spent constructing Facet (pairs) for no reason as it's always
  // just to grab the .first and .second as soon as it's constructed, and not due to API requirements
  // e.g. from DT3
  Facet_queue_status facet_status(const Facet& f) const
  {
    CGAL_precondition(!m_dt.is_infinite(f));

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "facet status: "
              << f.first->vertex((f.second + 1)&3)->point() << " "
              << f.first->vertex((f.second + 2)&3)->point() << " "
              << f.first->vertex((f.second + 3)&3)->point() << std::endl;
#endif

    // skip if neighbor is OUTSIDE or infinite
    const Cell_handle ch = f.first;
    const int id = f.second;
    const Cell_handle nh = ch->neighbor(id);
    if(m_dt.is_infinite(nh))
      return TRAVERSABLE;

    if(nh->info().is_outside)
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "Neighbor already outside" << std::endl;
#endif
      return IRRELEVANT;
    }

    // push if facet is connected to artificial vertices
    for(int i=0; i<3; ++i)
    {
      const Vertex_handle vh = ch->vertex(Dt::vertex_triple_index(id, i));
      if(vh->info() == BBOX_VERTEX || vh->info() == SEED_VERTEX)
      {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
        std::cout << "artificial facet due to artificial vertex #" << i << std::endl;
#endif
        return ARTIFICIAL_FACET;
      }
    }

    // skip if f min empty sphere radius is smaller than alpha
    if(is_traversable(f))
    {
#ifdef CGAL_AW3_DEBUG_FACET_STATUS
      std::cout << "traversable" << std::endl;
#endif
      return TRAVERSABLE;
    }

#ifdef CGAL_AW3_DEBUG_FACET_STATUS
    std::cout << "not traversable" << std::endl;
#endif
    return IRRELEVANT;
  }

  bool push_facet(const Facet& f)
  {
    CGAL_precondition(f.first->info().is_outside);

    // skip if f is already in queue
    if(m_queue.contains_with_bounds_check(Gate(f)))
      return false;

    const Facet_queue_status s = facet_status(f);
    if(s == IRRELEVANT)
      return false;

    const Cell_handle ch = f.first;
    const int id = f.second;
    const Point_3& p0 = m_dt.point(ch, (id+1)&3);
    const Point_3& p1 = m_dt.point(ch, (id+2)&3);
    const Point_3& p2 = m_dt.point(ch, (id+3)&3);

    // @todo should prob be the real value we compare to alpha instead of squared_radius
    const FT sqr = geom_traits().compute_squared_radius_3_object()(p0, p1, p2);
    m_queue.resize_and_push(Gate(f, sqr, (s == ARTIFICIAL_FACET)));

    return true;
  }

private:
  template <typename SeedRange>
  bool initialize(const double alpha,
                  const double offset,
                  const SeedRange& seeds)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Initialize..." << std::endl;
    std::cout << "Alpha: " << alpha << std::endl;
    std::cout << "Offset: " << offset << std::endl;
#endif

    if(!is_positive(alpha) || !is_positive(offset))
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Error: invalid input parameters" << std::endl;
#endif
      return false;
    }

    m_alpha = FT(alpha);
    m_sq_alpha = square(m_alpha);
    m_offset = FT(offset);
    m_sq_offset = square(m_offset);

    m_dt.clear();
    m_queue.clear();

    insert_bbox_corners();

    if(seeds.empty())
      return initialize_from_infinity();
    else
      return initialize_with_cavities(seeds);
  }

  template <typename Visitor>
  void alpha_flood_fill(Visitor& visitor)
  {
#ifdef CGAL_AW3_DEBUG
    std::cout << "> Flood fill..." << std::endl;
#endif

    visitor.on_flood_fill_begin(*this);

    // Explore all finite cells that are reachable from one of the initial outside cells.
    while(!m_queue.empty())
    {
#ifdef CGAL_AW3_DEBUG_QUEUE_PP
      check_queue_sanity();
#endif

      // const& to something that will be popped, but safe as `ch` && `id` are extracted before the pop
      const Gate& gate = m_queue.top();
      const Facet& f = gate.facet();
      CGAL_precondition(!m_dt.is_infinite(f));

      const Cell_handle ch = f.first;
      const int id = f.second;
      const Cell_handle neighbor = ch->neighbor(id);

#ifdef CGAL_AW3_DEBUG_QUEUE
      static int fid = 0;
      std::cout << m_dt.number_of_vertices() << " DT vertices" << std::endl;
      std::cout << m_queue.size() << " facets in the queue" << std::endl;
      std::cout << "Face " << fid++ << "\n"
                << "c = " << &*ch << " (" << m_dt.is_infinite(ch) << "), n = " << &*neighbor << " (" << m_dt.is_infinite(neighbor) << ")" << "\n"
                << m_dt.point(ch, (id+1)&3) << "\n" << m_dt.point(ch, (id+2)&3) << "\n" << m_dt.point(ch, (id+3)&3) << std::endl;
      std::cout << "Priority: " << gate.priority() << std::endl;
#endif

      visitor.before_facet_treatment(*this, gate);

      m_queue.pop();

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
      static int i = 0;
      std::string step_name = "results/steps/step_" + std::to_string(static_cast<int>(i)) + ".off";
      dump_triangulation_faces(step_name, true /*only_boundary_faces*/);

      std::string face_name = "results/steps/face_" + std::to_string(static_cast<int>(i++)) + ".xyz";
      std::ofstream face_out(face_name);
      face_out.precision(17);
      face_out << "3\n" << m_dt.point(ch, (id+1)&3) << "\n" << m_dt.point(ch, (id+2)&3) << "\n" << m_dt.point(ch, (id+3)&3) << std::endl;
      face_out.close();
#endif

      if(m_dt.is_infinite(neighbor))
      {
        neighbor->info().is_outside = true;
        continue;
      }

      Point_3 steiner_point;
      if(compute_steiner_point(ch, neighbor, steiner_point))
      {
//        std::cout << CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset)
//                  << " vs " << 1e-2 * m_offset << std::endl;
        CGAL_assertion(CGAL::abs(CGAL::approximate_sqrt(m_oracle.squared_distance(steiner_point)) - m_offset) <= 1e-2 * m_offset);

        // locate cells that are going to be destroyed and remove their facet from the queue
        int li, lj = 0;
        Locate_type lt;
        const Cell_handle conflict_cell = m_dt.locate(steiner_point, lt, li, lj, neighbor);
        CGAL_assertion(lt != Dt::VERTEX);

        std::vector<Facet> boundary_facets;
        std::vector<Cell_handle> conflict_zone;
        boundary_facets.reserve(32);
        conflict_zone.reserve(32);
        m_dt.find_conflicts(steiner_point, conflict_cell,
                            std::back_inserter(boundary_facets),
                            std::back_inserter(conflict_zone));

        // Purge the queue of facets that will be deleted/modified by the Steiner point insertion,
        // and which might have been gates
        for(const Cell_handle& cch : conflict_zone)
        {
          for(int i=0; i<4; ++i)
          {
            const Facet cf = std::make_pair(cch, i);
            if(m_queue.contains_with_bounds_check(Gate(cf)))
              m_queue.erase(Gate(cf));
          }
        }

        for(const Facet& f : boundary_facets)
        {
          const Facet mf = m_dt.mirror_facet(f); // boundary facets have incident cells in the CZ
          if(m_queue.contains_with_bounds_check(Gate(mf)))
            m_queue.erase(Gate(mf));
        }

        visitor.before_Steiner_point_insertion(*this, steiner_point);

        // Actual insertion of the Steiner point
        Vertex_handle vh = m_dt.insert(steiner_point, lt, conflict_cell, li, lj);
        vh->info() = DEFAULT;

        visitor.after_Steiner_point_insertion(*this, vh);

        std::vector<Cell_handle> new_cells;
        new_cells.reserve(32);
        m_dt.incident_cells(vh, std::back_inserter(new_cells));
        for(const Cell_handle& ch : new_cells)
        {
          // std::cout << "new cell has time stamp " << ch->time_stamp() << std::endl;
          ch->info().is_outside = m_dt.is_infinite(ch);
        }

        // Push all new boundary facets to the queue.
        // It is not performed by looking at the facets on the boundary of the conflict zones
        // because we need to handle internal facets, infinite facets, and also more subtle changes
        // such as a new cell being marked inside which now creates a boundary
        // with its incident "outside" flagged cell.
        for(Cell_handle ch : new_cells)
        {
          for(int i=0; i<4; ++i)
          {
            if(m_dt.is_infinite(ch, i))
              continue;

            const Cell_handle nh = ch->neighbor(i);
            if(nh->info().is_outside == ch->info().is_outside) // not on the boundary
              continue;

            const Facet boundary_f = std::make_pair(ch, i);
            if(ch->info().is_outside)
              push_facet(boundary_f);
            else
              push_facet(m_dt.mirror_facet(boundary_f));
          }
        }
      }
      else
      {
        // tag neighbor as OUTSIDE
        neighbor->info().is_outside = true;

        // for each finite facet of neighbor, push it to the queue
        for(int i=0; i<4; ++i)
        {
          const Facet neighbor_f = std::make_pair(neighbor, i);
          push_facet(neighbor_f);
        }
      }
    } // while(!queue.empty())

    visitor.on_flood_fill_end(*this);

    // Check that no useful facet has been ignored
    CGAL_postcondition_code(for(auto fit=m_dt.finite_facets_begin(), fend=m_dt.finite_facets_end(); fit!=fend; ++fit) {)
    CGAL_postcondition_code(  if(fit->first->info().is_outside == fit->first->neighbor(fit->second)->info().is_outside) continue;)
    CGAL_postcondition_code(  Facet f = *fit;)
    CGAL_postcondition_code(  if(!fit->first->info().is_outside) f = m_dt.mirror_facet(f);)
    CGAL_postcondition(       facet_status(f) == IRRELEVANT);
    CGAL_postcondition_code(})
  }

private:
  bool is_non_manifold(Vertex_handle v) const
  {
    CGAL_precondition(!m_dt.is_infinite(v));

    bool is_non_manifold = false;

    std::vector<Cell_handle> inc_cells;
    inc_cells.reserve(64);
    m_dt.incident_cells(v, std::back_inserter(inc_cells));

    // Flood one inside and outside CC.
    // Process both an inside and an outside CC to also detect edge pinching.
    // If there are still unprocessed afterwards, there is a non-manifoldness issue.
    //
    // Squat the conflict cell data to mark visits
    Cell_handle inside_start = Cell_handle();
    Cell_handle outside_start = Cell_handle();

    for(Cell_handle ic : inc_cells)
    {
      ic->tds_data().clear();
      if(ic->info().is_outside)
        outside_start = ic;
      else if(inside_start == Cell_handle())
        inside_start = ic;
    }

    // fully inside / outside
    if(inside_start == Cell_handle() || outside_start == Cell_handle())
      return false;

    std::stack<Cell_handle> cells_to_visit;
    cells_to_visit.push(inside_start);
    cells_to_visit.push(outside_start);
    while(!cells_to_visit.empty())
    {
      Cell_handle curr_c = cells_to_visit.top();
      cells_to_visit.pop();

      if(curr_c->tds_data().processed())
        continue;
      curr_c->tds_data().mark_processed();

      int v_pos = -1;
      CGAL_assertion_code(bool res = ) curr_c->has_vertex(v, v_pos);
      CGAL_assertion(res);

      for(int j=0; j<4; ++j)
      {
        if(j == v_pos)
          continue;

        Cell_handle neigh_c = curr_c->neighbor(j);
        CGAL_assertion(neigh_c->has_vertex(v));

        if(neigh_c->tds_data().processed() ||
           neigh_c->info().is_outside != curr_c->info().is_outside) // do not cross the boundary
          continue;

        cells_to_visit.push(neigh_c);
      }
    }

    for(Cell_handle ic : inc_cells)
    {
      if(ic->tds_data().is_clear()) // <=> unvisited cell
      {
        is_non_manifold = true;
        break;
      }
    }

    // reset the conflict flags
    for(Cell_handle ic : inc_cells)
      ic->tds_data().clear();

    return is_non_manifold;
  }

  bool is_non_manifold(Cell_handle c) const
  {
    CGAL_precondition(!m_dt.is_infinite(c));

    for(int i=0; i<4; ++i)
    {
      Vertex_handle v = c->vertex(i);
      if(is_non_manifold(v))
        return true;
    }

    return false;
  }

  bool is_manifold() const
  {
    // Not the best complexity, but it's not important: this function is purely for information
    // Better complexity --> see PMP::non_manifold_vertices + throw
    for(const Vertex_handle v : m_dt.finite_vertex_handles())
      if(is_non_manifold(v))
        return true;

    return false;
  }

  // Remove bbox vertices, if they are not necessary (i.e., no "inside" incident cell)
  // This is to try and avoid having long tets with bbox vertices being tagged "inside" as part
  // of the manifold re-tagging
  bool remove_bbox_vertices()
  {
    bool do_remove = true;
    auto vit = m_dt.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_dt.finite_incident_cells(v, std::back_inserter(inc_cells));

      for(Cell_handle c : inc_cells)
      {
        if(!c->info().is_outside)
        {
          do_remove = false;
          break;
        }
      }

      if(!do_remove)
        break;
    }

    std::cout << "Removing bbox vertices? " << do_remove << std::endl;
    if(!do_remove)
      return false;

    vit = m_dt.finite_vertices_begin();
    for(std::size_t i=0; i<8; ++i)
    {
      Vertex_handle v = vit++;
      m_dt.remove(v);
    }

    return true;
  }

public:
  // Not the best complexity, but it's very cheap compared to the rest of the algorithm.
  void make_manifold()
  {
    namespace PMP = Polygon_mesh_processing;

    // This seems more harmful than useful after the priority queue has been introduced since
    // it adds a lot of flat cells into the triangulation, which then get added to the mesh
    // during manifoldness fixing.
//    remove_bbox_vertices();

    std::stack<Vertex_handle> non_manifold_vertices; // @todo sort somehow?
    for(Vertex_handle v : m_dt.finite_vertex_handles())
    {
      if(is_non_manifold(v))
        non_manifold_vertices.push(v);
    }

    // Some lambdas for the comparer
    auto has_artificial_vertex = [](Cell_handle c) -> bool
    {
      for(int i=0; i<4; ++i)
        if(c->vertex(i)->info() == BBOX_VERTEX || c->vertex(i)->info() == SEED_VERTEX)
          return true;

      return false;
    };

    auto is_on_boundary = [](Cell_handle c, int i) -> bool
    {
      return (c->info().is_outside != c->neighbor(i)->info().is_outside);
    };

    auto count_boundary_facets = [&](Cell_handle c, Vertex_handle v) -> int
    {
      const int v_index_in_c = c->index(v);
      int boundary_facets = 0;
      for(int i=0; i<3; ++i) // also take into account the opposite facet?
      {
        if(i == v_index_in_c)
          continue;

        boundary_facets += is_on_boundary(c, i);
      }

      return boundary_facets;
    };

    // longest edge works better
//    auto sq_circumradius = [&](Cell_handle c) -> FT
//    {
//      const Point_3& cc = circumcenter(c);
//      return geom_traits().compute_squared_distance_3_object()(m_dt.point(c, 0), cc);
//    };

    auto sq_longest_edge = [&](Cell_handle c) -> FT
    {
      return (std::max)({ squared_distance(m_dt.point(c, 0), m_dt.point(c, 1)),
                          squared_distance(m_dt.point(c, 0), m_dt.point(c, 2)),
                          squared_distance(m_dt.point(c, 0), m_dt.point(c, 3)),
                          squared_distance(m_dt.point(c, 1), m_dt.point(c, 2)),
                          squared_distance(m_dt.point(c, 3), m_dt.point(c, 3)),
                          squared_distance(m_dt.point(c, 2), m_dt.point(c, 3)) });
    };

#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS
    std::cout << non_manifold_vertices.size() << " initial NMV" << std::endl;
#endif

    while(!non_manifold_vertices.empty())
    {
#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS
      std::cout << non_manifold_vertices.size() << " NMV in queue" << std::endl;
#endif

      Vertex_handle v = non_manifold_vertices.top();
      non_manifold_vertices.pop();

#ifdef CGAL_AW3_DEBUG_MANIFOLDNESS
      std::cout << "·";
#endif

      if(!is_non_manifold(v))
        continue;

      // Prioritize:
      // - cells without bbox vertices
      // - cells that already have a large number of boundary facets
      // - small cells when equal number of boundary facets
      // @todo give topmost priority to cells with > 1 non-manifold vertex?
      auto comparer = [&](Cell_handle l, Cell_handle r) -> bool
      {
        if(has_artificial_vertex(l))
          return false;
        if(has_artificial_vertex(r))
          return true;

        const int l_bf_count = count_boundary_facets(l, v);
        const int r_bf_count = count_boundary_facets(r, v);
        if(l_bf_count != r_bf_count)
          return l_bf_count > r_bf_count;

        return sq_longest_edge(l) < sq_longest_edge(r);
      };

      std::vector<Cell_handle> inc_cells;
      inc_cells.reserve(64);
      m_dt.finite_incident_cells(v, std::back_inserter(inc_cells));

#define CGAL_AW3_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
#ifndef CGAL_AW3_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
      std::sort(inc_cells.begin(), inc_cells.end(), comparer); // sort once
#endif

      for(auto cit=inc_cells.begin(), cend=inc_cells.end(); cit!=cend; ++cit)
      {
#ifdef CGAL_AW3_USE_BRUTE_FORCE_MUTABLE_PRIORITY_QUEUE
        // sort at every iteration since the number of boundary facets evolves
        std::sort(cit, cend, comparer);
#endif
        Cell_handle ic = *cit;
        CGAL_assertion(!m_dt.is_infinite(ic));

        // This is where new material is added
        ic->info().is_outside = false;

#ifdef CGAL_AW3_DEBUG_DUMP_EVERY_STEP
        static int i = 0;
        std::string step_name = "results/steps_manifold/step" + std::to_string(static_cast<int>(i++)) + ".off";
        dump_triangulation_faces(step_name, true /*only_boundary_faces*/);
#endif

        // @speed could update the manifold status while tagging
        if(!is_non_manifold(v))
          break;
      }

      CGAL_assertion(!is_non_manifold(v));

      std::vector<Vertex_handle> adj_vertices;
      adj_vertices.reserve(64);
      m_dt.finite_adjacent_vertices(v, std::back_inserter(adj_vertices));

      for(Vertex_handle nv : adj_vertices)
        if(is_non_manifold(nv))
          non_manifold_vertices.push(nv);
    }

    CGAL_assertion_code(for(Vertex_handle v : m_dt.finite_vertex_handles()))
    CGAL_assertion(!is_non_manifold(v));
  }

private:
  void check_queue_sanity()
  {
    std::cout << "Check queue sanity..." << std::endl;
    std::vector<Gate> queue_gates;
    Gate previous_top_gate = m_queue.top();
    while(!m_queue.empty())
    {
      const Gate& current_gate = m_queue.top();
      queue_gates.push_back(current_gate);
      const Facet& current_f = current_gate.facet();
      const Cell_handle ch = current_f.first;
      const int id = current_f.second;
      const Point_3& p0 = m_dt.point(ch, (id+1)&3);
      const Point_3& p1 = m_dt.point(ch, (id+2)&3);
      const Point_3& p2 = m_dt.point(ch, (id+3)&3);
      const FT sqr = geom_traits().compute_squared_radius_3_object()(p0, p1, p2);

      std::cout << "At Facet with VID " << get(Gate_ID_PM<Dt>(), current_gate)  << std::endl;

      if(current_gate.priority() != sqr)
        std::cerr << "Error: facet in queue has wrong priority" << std::endl;

      if(Less_gate()(current_gate, previous_top_gate))
        std::cerr << "Error: current gate has higher priority than the previous top" << std::endl;

      previous_top_gate = current_gate;

      m_queue.pop();
    }
    std::cout << "End sanity check" << std::endl;

    // Rebuild
    CGAL_assertion(m_queue.empty());
    for(const auto& g : queue_gates)
      m_queue.push(g); // no need for a resize here since the vector capacity is unchanged
  }

  void dump_triangulation_faces(const std::string filename,
                                bool only_boundary_faces = false)
  {
    std::stringstream vertices_ss;
    vertices_ss.precision(17);

    std::stringstream facets_ss;
    facets_ss.precision(17);

    std::unordered_map<Vertex_handle, std::size_t> vertex_to_id;
    std::size_t nv = 0;
    std::size_t nf = 0;

    for(auto fit=m_dt.finite_facets_begin(), fend=m_dt.finite_facets_end(); fit!=fend; ++fit)
    {
      Cell_handle c = fit->first;
      int s = fit->second;

      Cell_handle nc = c->neighbor(s);
      if(only_boundary_faces && (c->info().is_outside == nc->info().is_outside))
        continue;

      std::array<std::size_t, 3> ids;
      for(std::size_t pos=0; pos<3; ++pos)
      {
        Vertex_handle v = c->vertex((s+pos+1)&3);
        auto insertion_res = vertex_to_id.emplace(v, nv);
        if(insertion_res.second)
        {
          vertices_ss << m_dt.point(v) << "\n";
          ++nv;
        }

        ids[pos] = insertion_res.first->second;
      }

      facets_ss << "3 " << ids[0] << " " << ids[1] << " " << ids[2] << "\n";
      ++nf;
    }

    std::ofstream out(filename.c_str());
    out << "OFF\n" << nv << " " << nf << " 0\n";
    out << vertices_ss.str() << "\n" << facets_ss.str() << std::endl;
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_ALPHA_WRAP_3_H
