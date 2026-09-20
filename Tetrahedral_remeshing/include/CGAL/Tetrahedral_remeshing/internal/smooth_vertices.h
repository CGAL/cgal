// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_SMOOTH_VERTICES_H
#define CGAL_INTERNAL_SMOOTH_VERTICES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Vector_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/FMLS.h>
#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>

#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#include <CGAL/AABB_segment_primitive_3.h>
#include <CGAL/use.h>

#include <optional>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include <cmath>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
template <typename C3t3, typename SizingFunction, typename CellSelector>
class Vertex_smoothing_context
{
public:
  using Tr                  = typename C3t3::Triangulation;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  using Cell_handle         = typename Tr::Cell_handle;
  using Vertex_handle       = typename Tr::Vertex_handle;
  using Edge                = typename Tr::Edge;
  using Facet               = typename Tr::Facet;

  static constexpr bool is_parallel
    = std::is_convertible_v<typename Tr::Concurrency_tag, CGAL::Parallel_tag>;

  // One finite edge and whether it is in the 1-D complex, so that the parallel
  // collection answers both questions in a single scan.
  struct Classified_edge
  {
    Edge e;
    bool in_complex;
  };

  using Gt       = typename Tr::Geom_traits;
  using Vector_3 = typename Gt::Vector_3;
  using Point_3  = typename Gt::Point_3;
  using FT       = typename Gt::FT;

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Segment_vec = std::vector<typename Gt::Segment_3>;
  using Segment_iter = typename Segment_vec::iterator;
  using Segment_primitive = CGAL::AABB_segment_primitive_3<Gt, Segment_iter>;
  using AABB_segment_traits = CGAL::AABB_traits_3<Gt, Segment_primitive>;
  using AABB_segment_tree = CGAL::AABB_tree<AABB_segment_traits>;

private:
  Triangle_vec m_aabb_triangles;
  Segment_vec m_aabb_segments;

public:
  using FMLS = CGAL::Tetrahedral_remeshing::internal::FMLS<Gt>;
  std::vector<FMLS> subdomain_FMLS;
  std::unordered_map<Surface_patch_index, std::size_t, boost::hash<Surface_patch_index>> subdomain_FMLS_indices;

public:
  AABB_triangle_tree m_triangles_aabb_tree;
  AABB_segment_tree m_segments_aabb_tree;
  FT m_aabb_epsilon;

  const SizingFunction& m_sizing;

  using Incident_cells_vector = boost::container::small_vector<Cell_handle, 64>;
  std::vector<Incident_cells_vector> m_inc_cells;

  /**
  * Every finite edge of the triangulation, scanned once per smooth phase by
  * `refresh()`.
  *
  * Two of the three smooth operations -- surface and internal -- walk the
  * finite edges to accumulate their moves, and no smooth operation changes
  * topology, so the second walk re-derives a set the first one already
  * enumerated. `Tr::finite_edges()` is not a cheap range: every increment
  * circulates the cells around the candidate edge to decide whether the
  * current cell is its canonical one, so an edge of degree d is visited d
  * times to be emitted once. Walking a vector of the same edges instead reads
  * each one once, from contiguous memory.
  *
  * Only the enumeration is shared. Each operation still reads live vertex
  * positions when it computes its moves, because the operations run in
  * sequence and the earlier ones move vertices the later ones read.
  */
  std::vector<Edge> m_finite_edges;

  // The complex (1-D feature) edges, and per finite edge whether it is one.
  // The 1-D complex is constant for the whole smooth phase, and both answers
  // are wanted: the complex-edge operation walks `m_complex_edges` instead of
  // `c3t3.edges_in_complex()`, which is a filter over `finite_edges()` and so
  // pays the re-derivation described above, and the surface operation reads
  // `m_finite_edge_in_complex` instead of calling `is_in_complex()` -- a map
  // lookup -- on every finite edge to reject them. `refresh()` classifies each
  // edge once, as it scans it, for both.
  //
  // `m_complex_edges` keeps `finite_edges()` order, which is the order
  // `edges_in_complex()` emitted: the moves are floating-point sums, so a
  // different order would round differently. Both are filled only when
  // `!m_protect_boundaries`, the only case in which either consumer runs.
  std::vector<Edge> m_complex_edges;
  std::vector<bool> m_finite_edge_in_complex;

  using Vertices_surface_indices_map = std::unordered_map<Vertex_handle, std::vector<Surface_patch_index>>;
  using Vertices_normals_map =
      std::unordered_map<Vertex_handle,
                         std::unordered_map<Surface_patch_index, Vector_3, boost::hash<Surface_patch_index>>>;

  Vertices_surface_indices_map m_vertices_surface_indices;

  // The per-vertex, per-patch normals, indexed by `vertex_id()` rather than
  // stored in a `Vertices_normals_map`. A vertex carries as many normals as it
  // has incident surface patches -- one, almost always, and three at a corner
  // -- so the patch is found by a linear scan of that short vector, which is a
  // contiguous read, and never by hashing a `Surface_patch_index`.
  //
  // The nested map it replaces cost two lookups per level per accumulation and
  // allocated one inner map per surface vertex, which `clear()` then freed,
  // every time `refresh()` ran. Here the outer vector and each inner vector
  // keep their capacity between calls: after the first `refresh()` the fill
  // allocates nothing.
  using Vertex_patch_normals = std::vector<std::pair<Surface_patch_index, Vector_3>>;
  std::vector<Vertex_patch_normals> m_vertices_normals;

  // The vertices `compute_vertices_normals()` actually wrote, in the order it
  // first touched them. Only surface vertices carry a normal -- 102 k of the
  // 774 k facets of `1146193_cdt_0.5` are boundary facets -- so emptying and
  // normalizing through this list keeps both passes proportional to the
  // SURFACE, where walking `m_vertices_normals` end to end would make them
  // proportional to the whole mesh, once per `refresh()`.
  std::vector<std::pair<std::size_t, Vertex_handle>> m_vertices_with_normals;

  // Scratch for the parallel normals fan-out (`compute_vertices_normals()`).
  // All three keep their capacity between calls and are indexed so that every
  // pass costs the SURFACE, not the mesh: `m_nrm_facet_count` and
  // `m_nrm_vertex_slot` are indexed by vertex id but are only ever written for
  // a vertex that carries a normal, and `m_nrm_facet_count` is zeroed through
  // `m_vertices_with_normals`, so it is all zeros on entry to every call.
  std::vector<std::uint32_t> m_nrm_facet_count;   // per vertex id
  std::vector<std::uint32_t> m_nrm_vertex_slot;   // vertex id -> position in m_vertices_with_normals
  std::vector<std::uint32_t> m_nrm_csr_offset;    // per slot, into m_nrm_csr
  std::vector<std::uint32_t> m_nrm_csr;           // facet indices, grouped by vertex

  const CellSelector& m_cell_selector;
  const bool m_protect_boundaries;

  const bool m_smooth_constrained_edges;

  // the 2 following variables become useful and valid
  // just before flip/smooth steps, when no vertices get inserted
  // nor removed anymore
  std::unordered_map<Vertex_handle, std::size_t> m_vertex_id;
  std::vector<bool> m_free_vertices{};
  bool m_flip_smooth_steps{false};

public:
  struct Move
  {
    Vector_3 move;
    int neighbors;
    FT mass;
  };
  std::vector<Move> m_moves{};
  FT m_total_move{0};

  Vertex_smoothing_context(C3t3& c3t3,
                           const SizingFunction& sizing,
                           const CellSelector& cell_selector,
                           const bool protect_boundaries,
                           const bool smooth_constrained_edges)
      : m_sizing(sizing)
      , m_cell_selector(cell_selector)
      , m_protect_boundaries(protect_boundaries)
      , m_smooth_constrained_edges(smooth_constrained_edges)
  {
    refresh(c3t3);
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    if (m_protect_boundaries)
    {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }
    createMLSSurfaces(subdomain_FMLS,
                      subdomain_FMLS_indices,
                      vertices_normals_map(c3t3.triangulation()),
                      m_vertices_surface_indices,
                      c3t3);
#else
    build_aabb_trees(c3t3);
#endif
  }

  void refresh(C3t3& c3t3)
  {
    // The id map comes first: `compute_vertices_normals()` stores its result by
    // `vertex_id()`. The steps are independent of one another, so which one
    // runs first is free to choose.
    reset_vertex_id_map(c3t3.triangulation());
    if (!m_protect_boundaries)
    {
      collect_vertices_surface_indices(c3t3);
      compute_vertices_normals(c3t3);
    }
    reset_free_vertices(c3t3.triangulation());
    collect_incident_cells(c3t3.triangulation());
    collect_finite_edges(c3t3);
  }

  void start_flip_smooth_steps(const C3t3& c3t3)
  {
    CGAL_assertion(!m_flip_smooth_steps);
    reset_vertex_id_map(c3t3.triangulation());
    reset_free_vertices(c3t3.triangulation());

    // once this variable is set to true,
    // m_vertex_id becomes constant and
    // m_free_vertices can refer to it safely
    m_flip_smooth_steps = true;
  }

  bool in_flip_smooth_steps() const { return m_flip_smooth_steps; }

  std::size_t vertex_id(const Vertex_handle v) const
  {
    CGAL_expensive_assertion(m_vertex_id.find(v) != m_vertex_id.end());
    return m_vertex_id.at(v);
  }

  bool is_free(const Vertex_handle v) const  { return m_free_vertices[vertex_id(v)]; }
  bool is_free(const std::size_t& vid) const { return m_free_vertices[vid]; }

  // The normal of `v` on patch `si`, as `compute_vertices_normals()` left it.
  // Reading a normal a vertex does not carry is a precondition violation, as it
  // was when this was a map and the read was `.at(v).at(si)`.
  const Vector_3& vertex_normal(const std::size_t vid,
                                const Surface_patch_index& si) const
  {
    for (const auto& [patch, n] : m_vertices_normals[vid])
      if (patch == si)
        return n;
    CGAL_error_msg("no normal stored for this (vertex, surface patch)");
    return m_vertices_normals[vid].front().second;
  }
  const Vector_3& vertex_normal(const Vertex_handle v,
                                const Surface_patch_index& si) const
  {
    return vertex_normal(vertex_id(v), si);
  }

  const Incident_cells_vector& incident_cells(const Vertex_handle v) const
  {
    return m_inc_cells[vertex_id(v)];
  }
  const Incident_cells_vector& incident_cells(const std::size_t& vid) const
  {
    return m_inc_cells[vid];
  }

private:
  bool is_selected(const Cell_handle c) const
  {
    return get(m_cell_selector, c);
  }

  // this function can be used iff m_vertex_id
  // has already been initialized
  void reset_free_vertices(const Tr& tr)
  {
    // when flip-smooth steps start,
    // m_free_vertices should already be initialized
    // by the last smoothing step.
    // Then, it does not need to be recomputed
    // because no vertices are inserted nor removed anymore
    if (m_flip_smooth_steps)
    {
      CGAL_assertion(m_free_vertices.size() == tr.number_of_vertices());
      return;
    }
    m_free_vertices.clear();
    m_free_vertices.resize(tr.number_of_vertices(), false);

    for (const Cell_handle c : tr.finite_cell_handles())
    {
      if (!is_selected(c))
        continue;

      for (auto vi : tr.vertices(c))
      {
        const std::size_t idi = vertex_id(vi);
        const int dim = vi->in_dimension();

        switch (dim)
        {
        case 3:
          m_free_vertices[idi] = true;
          break;
        case 2:
          m_free_vertices[idi] = !m_protect_boundaries;
          break;
        case 1:
          m_free_vertices[idi] = !m_protect_boundaries && m_smooth_constrained_edges;
          break;
        case 0:
          m_free_vertices[idi] = false;
          break;
        default:
          CGAL_unreachable();
        }
      }
    }
  }

  // The segments of the 1D complex, for the AABB tree above.
  //
  // `c3t3.edges_in_complex()` is a filter over `finite_edges()`: walking it
  // enumerates every finite edge of the triangulation to find the complex
  // ones. The complex's own storage answers the same question in O(complex
  // edges) -- but in ITS order, and the order here is not free: it is the
  // order the primitives enter the AABB tree, which decides the tree's
  // structure and so which of two equidistant primitives a query returns.
  // Under `Parallel_tag` that is already a scheduling artefact and there is no
  // output to preserve; under `Sequential_tag` the mesh must not move, so the
  // filtered walk stays.
  void collect_complex_segments(const C3t3& c3t3)
  {
    const Tr& tr = c3t3.triangulation();
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (is_parallel)
    {
      Tetrahedral_remeshing::internal::for_each_edge_in_complex(c3t3,
        [this, &tr](const Vertex_handle& v1, const Vertex_handle& v2,
                    const auto&)
        {
          m_aabb_segments.push_back(
            tr.construct_segment(v1->point(), v2->point()));
        });
      return;
    }
#endif
    for (const Edge& e : c3t3.edges_in_complex())
      m_aabb_segments.push_back(tr.segment(e));
  }

  void build_aabb_trees(const C3t3& c3t3)
  {
    // build AABB tree of facets in complex
    for (const Facet& f : c3t3.facets_in_complex())
    {
      m_aabb_triangles.push_back(c3t3.triangulation().triangle(f));
    }
    m_triangles_aabb_tree.rebuild(m_aabb_triangles.begin(), m_aabb_triangles.end());
    m_triangles_aabb_tree.accelerate_distance_queries();

    // build AABB tree of edges in complex
    collect_complex_segments(c3t3);
    m_segments_aabb_tree.rebuild(m_aabb_segments.begin(), m_aabb_segments.end());
    m_segments_aabb_tree.accelerate_distance_queries();

    // compute epsilon for AABB tree of facets
    const CGAL::Bbox_3& bb = m_triangles_aabb_tree.bbox();
    m_aabb_epsilon = 1e-3 * (std::min)(bb.xmax() - bb.xmin(),
                            (std::min)(bb.ymax() - bb.ymin(),
                                       bb.zmax() - bb.zmin()));
  }

  // Topology is constant for the whole smooth phase, so this holds until the
  // next refresh(). The vector keeps its capacity between phases.
  void collect_finite_edges(const C3t3& c3t3)
  {
    const Tr& tr = c3t3.triangulation();

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (is_parallel)
    {
      // This walk is one of the serial sections in front of a phase that then
      // runs on every thread. `parallel_collect_from_finite_edges()` is the
      // same cells x 6-slots scan with the same smallest-incident-cell
      // ownership rule that `Finite_edges_iterator::operator++` performs, made
      // splittable -- so the edges come out in `finite_edges()` order, which
      // `m_complex_edges` and the index into `m_finite_edge_in_complex` both
      // depend on. The `is_in_complex()` lookup rides along in the same pass
      // (the 1-D complex is not written during collection, and the parallel
      // storage is a `boost::concurrent_flat_map`), so the classification is
      // parallel too rather than a second serial walk over the result.
      //
      // The fan-out below is serial but is a copy over contiguous memory, with
      // no circulation and no lookup.
      const bool classify_par = !m_protect_boundaries;
      if (!classify_par)
      {
        m_complex_edges.clear();
        m_finite_edge_in_complex.clear();
        m_finite_edges = parallel_collect_from_finite_edges<Edge>(tr,
          [](const Edge& e, std::vector<Edge>& out) { out.push_back(e); });
        return;
      }

      const std::vector<Classified_edge> classified
        = parallel_collect_from_finite_edges<Classified_edge>(tr,
            [&c3t3](const Edge& e, std::vector<Classified_edge>& out)
            { out.push_back(Classified_edge{e, c3t3.is_in_complex(e)}); });

      m_finite_edges.clear();
      m_finite_edges.reserve(classified.size());
      m_finite_edge_in_complex.clear();
      m_finite_edge_in_complex.reserve(classified.size());
      m_complex_edges.clear();
      m_complex_edges.reserve(c3t3.number_of_edges_in_complex());

      for (const Classified_edge& ce : classified)
      {
        m_finite_edges.push_back(ce.e);
        m_finite_edge_in_complex.push_back(ce.in_complex);
        if (ce.in_complex)
          m_complex_edges.push_back(ce.e);
      }
      return;
    }
#endif

    // `number_of_vertices() + number_of_cells()` is an upper bound on the
    // number of finite edges, and an O(1) one -- both are container sizes.
    // Euler on the triangulated 3-sphere the TDS holds (the infinite vertex
    // included) gives V - E + F - C = 0 with F = 2C, hence E = V + C; the
    // finite edges are that less the edges to the infinite vertex. Measured
    // over 56 smooth phases on four configs the bound held every time, with
    // 0.3% to 6.9% of slack.
    //
    // Reserving matters only for the FIRST smooth phase: `clear()` keeps the
    // capacity, so later phases reuse it and allocate only when the mesh has
    // grown. But that first phase sets the high-water mark, and growing by
    // doubling would leave the vector holding up to twice what it needs at the
    // moment it reallocates last. Peak memory is a gated metric here.
    m_finite_edges.clear();
    m_finite_edges.reserve(tr.number_of_vertices() + tr.number_of_cells());

    // `number_of_edges_in_complex()` is the exact count of the complex edges,
    // and O(1) -- it is the size of the 1-D complex's own container.
    const bool classify = !m_protect_boundaries;
    m_complex_edges.clear();
    m_finite_edge_in_complex.clear();
    if (classify)
    {
      m_complex_edges.reserve(c3t3.number_of_edges_in_complex());
      m_finite_edge_in_complex.reserve(tr.number_of_vertices() + tr.number_of_cells());
    }

    for (const Edge& e : tr.finite_edges())
    {
      m_finite_edges.push_back(e);

      if (!classify)
        continue;

      const bool in_complex = c3t3.is_in_complex(e);
      m_finite_edge_in_complex.push_back(in_complex);
      if (in_complex)
        m_complex_edges.push_back(e);
    }
  }

  void collect_incident_cells(const Tr& tr)
  {
    m_inc_cells.clear();
    const std::size_t nbv = tr.number_of_vertices();
    m_inc_cells.resize(nbv, Incident_cells_vector{});
    for (const Cell_handle c : tr.finite_cell_handles())
    {
      for (auto vi : tr.vertices(c))
      {
        const std::size_t idi = vertex_id(vi);
        if(is_free(idi))
          m_inc_cells[idi].push_back(c);
      }
    }
  }

  std::optional<Facet>
  find_adjacent_facet_on_surface(const Facet& f,
                                 const Edge& edge,
                                 const C3t3& c3t3)
  {
    CGAL_expensive_assertion(is_boundary(c3t3, f, m_cell_selector));

    typedef typename Tr::Facet_circulator Facet_circulator;

    if (c3t3.is_in_complex(edge))
      return {}; //do not "cross" complex edges
    //they are likely to be sharp and not to follow the > 0 dot product criterion

    const Surface_patch_index& patch = c3t3.surface_patch_index(f);
    const Facet& mf = c3t3.triangulation().mirror_facet(f);

    Facet_circulator fcirc = c3t3.triangulation().incident_facets(edge);
    Facet_circulator fend = fcirc;
    do
    {
      const Facet fi = *fcirc;
      if (f != fi
          && mf != fi
          && is_boundary(c3t3, fi, m_cell_selector)
          && patch == c3t3.surface_patch_index(fi))
      {
        return canonical_facet(fi); //"canonical" is important
      }
    } while (++fcirc != fend);

    return {};
  }

  template<typename Gt_>
  Vector_3 compute_normal(const Facet& f,
                          const Vector_3& reference_normal,
                          const Gt_& gt)
  {
    typename Gt_::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();
    typename Gt_::Compute_scalar_product_3
      scalar_product = gt.compute_scalar_product_3_object();

    Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, gt);
    if (scalar_product(n, reference_normal) < 0.)
      n = opp(n);

    return n;
  }

  template<typename Patch_index>
  std::string debug_to_string(const Patch_index i)
  {
    return std::to_string(i);
  }

  template<typename Patch_index>
  std::string debug_to_string(const std::pair<Patch_index, Patch_index>& pi)
  {
    std::string str = std::to_string(pi.first);
    str.append("_").append(std::to_string(pi.second));
    return str;
  }

#ifdef CGAL_LINKED_WITH_TBB
  /**
  * The parallel form of `compute_vertices_normals()` below. Both halves of the
  * routine run on all threads, and neither hash map survives.
  *
  * Pass 1 collects every boundary facet's area-weighted normal, with the ids
  * and handles of its three vertices, on all threads. The serial form puts
  * these in an `unordered_map<Facet, Vector_3>` and then walks that map --
  * 102 k random-access reads on `1146193_cdt_0.5` -- for no reason other than
  * that it is where the first loop happened to leave them. Here they stay in
  * the flat vector the collector returns.
  *
  * Passes 2 and 3 group the facets by vertex: a count per vertex, a prefix sum
  * over the vertices that have one, and a fill. Both are linear passes over
  * ~3 x 102 k entries of `std::uint32_t` and stay serial -- they are ~1% of the
  * routine, and a parallel fill would need atomics to do the same work.
  *
  * Pass 4 is the one that costs. Each vertex is summed and normalized by ONE
  * thread, so no accumulator is shared and there is nothing to reduce: a
  * `parallel_for` over vertices is the whole of it.
  *
  * WHAT MOVES. The per-vertex sums are floating point and a vertex's facets
  * are now added in facet order rather than in the order an `unordered_map`
  * happened to iterate, so the normals differ from the serial path's in the
  * last bits, and the mesh follows. That is inside what this path already
  * does: at more than one thread the remeshing is nondeterministic by
  * construction -- which operations succeed depends on the order locks are
  * taken -- so the output is not reproducible run to run whatever this routine
  * does. The `Sequential_tag` build does not execute any of this code and its
  * output is unchanged, which is the gate that matters (POLICY 5.5).
  */
  void compute_vertices_normals_parallel(const C3t3& c3t3)
  {
    const Tr& tr = c3t3.triangulation();
    const typename Tr::Geom_traits& gt = tr.geom_traits();
    typename Tr::Geom_traits::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();

    struct Facet_normal
    {
      std::array<std::size_t, 3> vids;
      std::array<Vertex_handle, 3> vhs;
      Surface_patch_index patch;
      Vector_3 n;
    };

    // ---- pass 1: the facet normals, on all threads -------------------------
    const std::vector<Facet_normal> fnormals
      = Tetrahedral_remeshing::internal::parallel_collect_from_finite_facets<Facet_normal>(
          tr,
          Tetrahedral_remeshing::internal::gather_all_cells(tr),
          [this, &c3t3, &tr, &gt, &opp](const Facet& trf, std::vector<Facet_normal>& out)
          {
            if (!is_boundary(c3t3, trf, m_cell_selector))
              return;

            const Facet f = canonical_facet(trf);
            const Cell_handle c = f.first;
            const Cell_handle neigh = f.first->neighbor(f.second);

            Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, gt);
            if (tr.is_infinite(neigh)
             || c3t3.subdomain_index(neigh) < c3t3.subdomain_index(c))
              n = opp(n);

            Facet_normal fn;
            fn.patch = c3t3.surface_patch_index(f);
            fn.n = n;
            int i = 0;
            for (const Vertex_handle vi : tr.vertices(f))
            {
              fn.vhs[i] = vi;
              fn.vids[i] = vertex_id(vi);
              ++i;
            }
            out.push_back(fn);
          });

    // ---- reset, through the list the last call left ------------------------
    for (const auto& [vid, v] : m_vertices_with_normals)
    {
      CGAL_USE(v);
      m_vertices_normals[vid].clear();
      m_nrm_facet_count[vid] = 0;
    }
    m_vertices_with_normals.clear();

    const std::size_t nbv = tr.number_of_vertices();
    m_vertices_normals.resize(nbv);
    m_nrm_facet_count.resize(nbv, 0);
    m_nrm_vertex_slot.resize(nbv);

    // ---- pass 2: count, and name the vertices that carry a normal ----------
    for (const Facet_normal& fn : fnormals)
    {
      for (int i = 0; i < 3; ++i)
      {
        const std::size_t vid = fn.vids[i];
        if (m_nrm_facet_count[vid]++ == 0)
        {
          m_nrm_vertex_slot[vid] = static_cast<std::uint32_t>(m_vertices_with_normals.size());
          m_vertices_with_normals.emplace_back(vid, fn.vhs[i]);
        }
      }
    }

    // ---- pass 3: prefix sum and fill ---------------------------------------
    const std::size_t nb_slots = m_vertices_with_normals.size();
    m_nrm_csr_offset.assign(nb_slots + 1, 0);
    for (std::size_t k = 0; k < nb_slots; ++k)
      m_nrm_csr_offset[k + 1] = m_nrm_csr_offset[k]
                              + m_nrm_facet_count[m_vertices_with_normals[k].first];

    m_nrm_csr.resize(m_nrm_csr_offset[nb_slots]);
    std::vector<std::uint32_t> cursor(m_nrm_csr_offset.begin(), m_nrm_csr_offset.end() - 1);
    for (std::size_t fi = 0; fi < fnormals.size(); ++fi)
    {
      for (int i = 0; i < 3; ++i)
      {
        const std::uint32_t k = m_nrm_vertex_slot[fnormals[fi].vids[i]];
        m_nrm_csr[cursor[k]++] = static_cast<std::uint32_t>(fi);
      }
    }

    // ---- pass 4: one vertex, one thread ------------------------------------
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb_slots),
      [&](const tbb::blocked_range<std::size_t>& range)
      {
        for (std::size_t k = range.begin(); k != range.end(); ++k)
        {
          Vertex_patch_normals& vpn = m_vertices_normals[m_vertices_with_normals[k].first];

          for (std::uint32_t j = m_nrm_csr_offset[k]; j != m_nrm_csr_offset[k + 1]; ++j)
          {
            const Facet_normal& fn = fnormals[m_nrm_csr[j]];
            auto patch_it = std::find_if(vpn.begin(), vpn.end(),
                                         [&fn](const auto& pn) { return pn.first == fn.patch; });
            if (patch_it == vpn.end())
              vpn.emplace_back(fn.patch, fn.n);
            else
              patch_it->second += fn.n;
          }

          for (auto& [surf_i, n] : vpn)
          {
            CGAL_USE(surf_i);
            CGAL::Tetrahedral_remeshing::normalize(n, gt);
          }
        }
      });
  }
#endif // CGAL_LINKED_WITH_TBB

  void compute_vertices_normals(const C3t3& c3t3)
  {
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (is_parallel)
    {
      compute_vertices_normals_parallel(c3t3);
      return;
    }
#endif

    // Emptied without giving the storage back: `resize()` only ever grows the
    // outer vector, and clearing an inner vector keeps its capacity. Only the
    // entries the last call filled need emptying, so the reset costs the
    // surface, not the mesh.
    for (const auto& [vid, v] : m_vertices_with_normals)
    {
      CGAL_USE(v);
      m_vertices_normals[vid].clear();
    }
    m_vertices_with_normals.clear();
    m_vertices_normals.resize(c3t3.triangulation().number_of_vertices());

    typename Tr::Geom_traits gt = c3t3.triangulation().geom_traits();
    typename Tr::Geom_traits::Construct_opposite_vector_3
      opp = gt.construct_opposite_vector_3_object();

    const Tr& tr = c3t3.triangulation();

    //collect all facet normals
    std::unordered_map<Facet, Vector_3, boost::hash<Facet>> fnormals;
    for (const Facet& trf : tr.finite_facets())
    {
      if(!is_boundary(c3t3, trf, m_cell_selector))
        continue;

      const Facet f = canonical_facet(trf);
      const Cell_handle c = f.first;
      const Cell_handle neigh = f.first->neighbor(f.second);

      Vector_3 n = CGAL::Tetrahedral_remeshing::normal(f, tr.geom_traits());
      if (c3t3.triangulation().is_infinite(neigh)
       || c3t3.subdomain_index(neigh) < c3t3.subdomain_index(c))
        n = opp(n);

      fnormals[f] = n; // n has length equal to the area of the facet
    }

    // accumulate the normals in normals_map
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream osf("dump_facet_normals.polylines.txt");
#endif
    for (const auto& [f, n] : fnormals)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      typename Tr::Geom_traits::Point_3 fc = CGAL::centroid(tr.triangle(f));
      osf << "2 " << fc << " " << (fc + n) << std::endl;
#endif
      const Surface_patch_index& surf_i = c3t3.surface_patch_index(f);

      for (const Vertex_handle vi : tr.vertices(f))
      {
        // Same additions, in the same order -- `fnormals` is still walked in
        // its own iteration order and a patch still gets its first normal
        // assigned and the rest added -- so the sums round exactly as before.
        const std::size_t vid = vertex_id(vi);
        Vertex_patch_normals& vpn = m_vertices_normals[vid];

        if (vpn.empty())
        {
          m_vertices_with_normals.emplace_back(vid, vi);
          vpn.emplace_back(surf_i, n);
          continue;
        }

        auto patch_it = std::find_if(vpn.begin(), vpn.end(),
                                     [&surf_i](const auto& pn) { return pn.first == surf_i; });
        if (patch_it == vpn.end())
          vpn.emplace_back(surf_i, n);
        else
          patch_it->second += n;
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    osf.close();
    std::ofstream os("dump_normals.polylines.txt");
    std::unordered_map<Surface_patch_index,
          std::vector<typename Tr::Geom_traits::Segment_3 >, boost::hash<Surface_patch_index> > ons_map;
#endif

    //normalize the computed normals
    // Only the vertices that got one, and each is normalized in place, so the
    // order this list happens to be in cannot change the result.
    for (const auto& [vid, v] : m_vertices_with_normals)
    {
      CGAL_USE(v);
      //value type is vector<pair<Surface_patch_index, Vector_3>>
      for (auto& [surf_i, n] : m_vertices_normals[vid])
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        auto p = point(v->point());
        os << "2 " << p << " " << (p + n) << std::endl;
#endif

        CGAL::Tetrahedral_remeshing::normalize(n, gt);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        const Surface_patch_index si = surf_i;
        if (ons_map.find(si) == ons_map.end())
          ons_map[si] = std::vector<typename Tr::Geom_traits::Segment_3>();
        ons_map[si].push_back(typename Tr::Geom_traits::Segment_3(p, p + n));
#endif
      }
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    os.close();
    for (auto& kv : ons_map)
    {
      std::ostringstream oss;
      oss << "dump_normals_normalized_["
        << debug_to_string(kv.first) << "].polylines.txt";
      std::ofstream ons(oss.str());
      for (auto s : kv.second)
        ons << "2 " << s.source() << " " << s.target() << std::endl;
      ons.close();
    }
#endif
  }

#ifdef CGAL_LINKED_WITH_TBB
  /**
  * The parallel form of the incident-patch cache below.
  *
  * `c3t3.facets_in_complex()` is a filter over `finite_facets()`, so reaching
  * the complex facets enumerates every facet of the mesh; that scan, and the
  * `is_in_complex()` test it applies, is what the threads take. The patch
  * index is read in the same pass, where the facet is already in hand.
  *
  * The fan-out stays serial -- it writes per-vertex vectors that adjacent
  * facets share. Its order does not matter: the only consumer outside the MLS
  * build reads `size()`, and element 0 when the size is 1, so a vertex's
  * patches are a set and nothing downstream depends on how they are arranged.
  */
  void collect_vertices_surface_indices_parallel(const C3t3& c3t3)
  {
    const Tr& tr = c3t3.triangulation();

    using Patch_facet = std::pair<Facet, Surface_patch_index>;
    const std::vector<std::vector<Patch_facet> > per_chunk
      = Tetrahedral_remeshing::internal::parallel_collect_chunks_from_finite_facets<Patch_facet>(
          tr,
          Tetrahedral_remeshing::internal::gather_all_cells(tr),
          [&c3t3](const Facet& f, std::vector<Patch_facet>& out)
          {
            if (c3t3.is_in_complex(f))
              out.emplace_back(f, c3t3.surface_patch_index(f));
          });

    for (const std::vector<Patch_facet>& chunk : per_chunk)
    {
      for (const Patch_facet& pf : chunk)
      {
        for (const Vertex_handle vi : tr.vertices(pf.first))
        {
          std::vector<Surface_patch_index>& v_surface_indices
            = m_vertices_surface_indices[vi];
          if (std::find(v_surface_indices.begin(), v_surface_indices.end(), pf.second)
              == v_surface_indices.end())
            v_surface_indices.push_back(pf.second);
        }
      }
    }
  }
#endif

  void collect_vertices_surface_indices(const C3t3& c3t3)
  {
    m_vertices_surface_indices.clear();
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (is_parallel)
    {
      collect_vertices_surface_indices_parallel(c3t3);
      return;
    }
#endif
    for (Facet fit : c3t3.facets_in_complex())
    {
      const Surface_patch_index& surface_index = c3t3.surface_patch_index(fit);

      for (const Vertex_handle vi : c3t3.triangulation().vertices(fit))
      {
        std::vector<Surface_patch_index>& v_surface_indices = m_vertices_surface_indices[vi];
        if (std::find(v_surface_indices.begin(), v_surface_indices.end(), surface_index) == v_surface_indices.end())
          v_surface_indices.push_back(surface_index);
      }
    }
  }

  // `createMLSSurfaces()` reads the normals as `.at(v).at(si)`. It runs once,
  // at construction, so it is given a map built from the dense storage rather
  // than the dense storage being shaped around it.
  Vertices_normals_map vertices_normals_map(const Tr& tr) const
  {
    CGAL_USE(tr);
    Vertices_normals_map map;
    for (const auto& [vid, v] : m_vertices_with_normals)
    {
      for (const auto& [surf_i, n] : m_vertices_normals[vid])
        map[v][surf_i] = n;
    }
    return map;
  }

  void reset_vertex_id_map(const Tr& tr)
  {
    // when flip-smooth steps start,
    // m_vertex_id should already be initialized,
    // done by the last smoothing step.
    // Then, it does not need to be recomputed
    // because no vertices are inserted nor removed anymore
    if(m_flip_smooth_steps)
      return;
    m_vertex_id.clear();
    std::size_t id = 0;
    for (const Vertex_handle v : tr.finite_vertex_handles())
    {
      m_vertex_id[v] = id++;
    }
  }
};


template <typename C3t3, typename SizingFunction, typename CellSelector>
class Vertex_smooth_operation_base
    : public Elementary_operation<C3t3,
                                  typename C3t3::Triangulation::Vertex_handle,
                                  typename C3t3::Triangulation::Finite_vertex_handles>
{
protected:
  using Base_operation = Elementary_operation<C3t3,
                                              typename C3t3::Triangulation::Vertex_handle,
                                              typename C3t3::Triangulation::Finite_vertex_handles>;

public:
  // the executors read these off the concrete operation
  using Element_type = typename Base_operation::Element_type;
  using Element_range = typename Base_operation::Element_range;
  static_assert(std::is_same_v<Element_type, typename C3t3::Triangulation::Vertex_handle>,
                "Element_type should be Vertex_handle");

protected:
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Facet Facet;

  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  typedef typename Gt::FT FT;

  using Triangle_vec = std::vector<typename Tr::Triangle>;
  using Triangle_iter = typename Triangle_vec::iterator;
  using Triangle_primitive = CGAL::AABB_triangle_primitive_3<Gt, Triangle_iter>;
  using AABB_triangle_traits = CGAL::AABB_traits_3<Gt, Triangle_primitive>;
  using AABB_triangle_tree = CGAL::AABB_tree<AABB_triangle_traits>;

  using Context = Vertex_smoothing_context<C3t3, SizingFunction, CellSelector>;

  std::shared_ptr<Context> m_context{nullptr};

public:
  Vertex_smooth_operation_base(std::shared_ptr<Context> context)
      : m_context(context) {}

protected:
  /**
  * The position this operation would move `v` to, or nothing if it declines
  * `v`. The only geometry that differs between the three smooth operations;
  * everything around it is shared and lives in this class.
  */
  virtual std::optional<Point_3> compute_target_position(const Vertex_handle v, const C3t3& c3t3) = 0;

  /**
  * Fills `Vertex_smoothing_context::m_moves` for every vertex this operation can
  * move. One serial pass, run from `get_elements()` before the candidates are
  * handed out.
  */
  virtual void compute_vertex_moves(const C3t3& c3t3) const = 0;

public:
  /**
  * Every finite vertex. All three operations decline the ones they do not
  * own, in `compute_target_position()`.
  */
  Element_range get_elements(const C3t3& c3t3) const override
  {
    compute_vertex_moves(c3t3);
    return c3t3.triangulation().finite_vertex_handles();
  }

  /**
  * Move `v` to the position this operation wants, unless that inverts a cell
  * of its star -- `check_inversion_and_move()` backtracks along the move and
  * restores `v` if it cannot find a valid fraction.
  *
  * Under `Parallel_tag` the position was already computed by `lock_zone()`,
  * which had to know it to lock the grid cell `v` lands in; recomputing it
  * here would mean a second AABB-tree projection for two of the three
  * operations. The handoff is accepted only for the vertex it was built for,
  * so the sequential executor -- which never calls `lock_zone()` -- computes
  * it here instead and is unaffected.
  */
  bool execute_operation(const Element_type& v, C3t3& c3t3) override
  {
    const Smooth_destination& handoff = destination_handoff();
    const std::optional<Point_3> target = (handoff.valid && handoff.v == v)
                                        ? handoff.position
                                        : compute_target_position(v, c3t3);
    if (target == std::nullopt)
      return false;

    return check_inversion_and_move(v, target.value(), m_context->incident_cells(v),
                                    c3t3.triangulation(), m_context->m_total_move);
  }

  // vertices are independent of one another: shuffling spreads the threads out
  static constexpr bool requires_ordered_processing = false;

  /**
  * The lock zone of all three smooth operations: `star(v)`, plus every
  * position `v` may be written to.
  *
  * Smoothing moves `v` and re-checks the orientation of every cell incident to
  * it, so the star is the whole zone and cannot be made smaller. The three
  * operations share this one definition: `Elementary_operation` names
  * `lock_zone()` in its contract but never declares it, so the executor's
  * `op.lock_zone(element, c3t3)` finds this by ordinary name lookup, and an
  * operation that needs a different zone defines its own.
  *
  * The star is NOT walked.
  *
  * `try_lock_and_get_incident_cells()` would traverse the star through the
  * neighbour pointers, marking and unmarking every cell's `tds_data()` on the
  * way, and then throw the vector away -- `execute_operation()` reads the star
  * from `Vertex_smoothing_context::m_inc_cells` instead. That cache is exactly
  * `star(v)` for the whole phase: `Adaptive_remesher::smooth()` rebuilds it in
  * `refresh()` on entry, and no smooth operation changes topology, only
  * `set_point()`. So the walk is pure loss, and the same cells can be locked
  * straight out of the cache.
  *
  * The set of grid cells taken for the star is the same, with one exception in
  * this function's favour: the cache holds the FINITE star, so the infinite
  * vertex's meaningless position is not taken. Nothing here reads it either --
  * `check_inversion_and_move()` iterates the same cache.
  *
  * A vertex that is not free has an EMPTY cache (`collect_incident_cells()`
  * fills only free ones) and all three `execute_operation()` decline it
  * immediately, so locking `v` alone is sufficient. `get_elements()` hands out
  * every finite vertex, so this is the common case.
  */
  bool lock_zone(const Vertex_handle v, const C3t3& c3t3)
  {
    Smooth_destination& handoff = destination_handoff();
    handoff.valid = false;

    const Tr& tr = c3t3.triangulation();

    if (!tr.try_lock_vertex(v))
      return false;

    for (const Cell_handle c : m_context->incident_cells(v))
    {
      if (!tr.try_lock_cell(c))
        return false;
    }

    const std::optional<Point_3> target = compute_target_position(v, c3t3);
    if (target != std::nullopt && !lock_move_destinations(v, target.value(), tr))
      return false;

    handoff.v = v;
    handoff.position = target;
    handoff.valid = true;
    return true;
  }

private:
  /**
  * Where `lock_zone()` leaves the position it computed, for
  * `execute_operation()` to consume.
  *
  * The destination has to be known INSIDE the zone, because the lock is keyed
  * on position and the grid cell `v` lands in must be held before `v` moves.
  * For the surface and complex-edge operations that position is an AABB-tree
  * projection, and computing it twice would cost more than the lock it makes
  * safe -- so it is computed once, here, and handed down. Same pattern as
  * `Located_edge` (flip, split) and `Locked_stars` (collapse).
  *
  * Thread-local, and `valid` is set LAST, so a zone abandoned part-way through
  * a failed lock attempt never matches. `execute_operation()` accepts the
  * handoff only for the vertex it was built for; anything else -- the
  * sequential executor, which never calls `lock_zone()` -- computes the
  * position itself, so the sequential path is unchanged by construction.
  */
  struct Smooth_destination
  {
    Vertex_handle v{};
    std::optional<Point_3> position{};
    bool valid{false};
  };

  static Smooth_destination& destination_handoff()
  {
    static thread_local Smooth_destination d;
    return d;
  }

  /**
  * Every position `check_inversion_and_move()` may write to `v`.
  *
  * It tries `pv + frac * move` and halves `frac` on a failed orientation or a
  * worsened angle, while `frac > 0.1`. So `frac` takes 1, 1/2, 1/4 and 1/8 and
  * then stops -- 1/16 fails the guard -- and those four points, all on one
  * segment, are the whole destination set. The restore to `pv` needs no lock
  * of its own: it is `v`'s own position, already held.
  *
  * Taking a point twice is cheap, the grid cell is already this thread's, and
  * on a grid coarser than the move they are one cell.
  */
  bool lock_move_destinations(const Vertex_handle v, const Point_3& final_pos, const Tr& tr) const
  {
    const Point_3 pv = point(v->point());
    const Vector_3 move(pv, final_pos);

    for (double frac = 1.0; frac > 0.1; frac = 0.5 * frac)
    {
      if (!tr.try_lock_point(pv + frac * move))
        return false;
    }
    return true;
  }

protected:
  Point_3 project_on_tangent_plane(const Point_3& gi, const Point_3& pi, const Vector_3& normal)
  {
    Vector_3 diff(gi, pi);
    return gi + (normal * diff) * normal;
  }

  FT density_along_segment(const Edge& e, const C3t3& c3t3, bool boundary_edge = false) const
  {
    const auto [pt, dim, index] = midpoint_with_info(e, boundary_edge, c3t3);
    const FT s = sizing_at_midpoint(e, pt, dim, index, m_context->m_sizing, c3t3, m_context->m_cell_selector);
    return 1. / s;
  }

  bool is_selected(const Cell_handle c) const { return get(m_context->m_cell_selector, c); }

#ifdef CGAL_LINKED_WITH_TBB
  /**
  * The parallel form of `accumulate_edge_moves()` below: the same pass, with
  * `keep_edge` answered on every thread and the accumulation left serial.
  *
  * The split is what makes this safe. `keep_edge` is read-only -- the surface
  * operation's is a lookup in `m_finite_edge_in_complex` and an `is_boundary()`
  * facet circulation, the internal one an `is_outside()` test -- and the
  * triangulation is not written anywhere in `get_elements()`, so the predicate
  * can run on all the edges at once. The four `moves[i] += ...` cannot: two
  * edges sharing an endpoint write the same entry, and the sum is in floating
  * point, so a reduce would change the result. Pass 2 therefore keeps the
  * original loop, in the original edge order, and reads the answers pass 1
  * stored. Same additions, same order, bit for bit.
  *
  * The filter is where the time is. On `1146193_cdt_0.5` at 4 threads the
  * surface operation scans 8.06 M edges to keep 1.41 M, and its loop costs
  * 0.88 s against the internal operation's 0.44 s over the same edges for 4.5x
  * as many kept -- the `is_boundary()` circulation, not the accumulation.
  */
  template <typename EdgeRange, typename KeepEdge, typename MovesVertex>
  void accumulate_edge_moves_parallel_filter(const EdgeRange& edges,
                                             const C3t3& c3t3,
                                             const bool boundary_edge,
                                             KeepEdge keep_edge,
                                             MovesVertex moves_vertex) const
  {
    auto& moves = m_context->m_moves;
    using Move = typename Context::Move;
    const Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
    moves.assign(c3t3.triangulation().number_of_vertices(), default_move);

    // One byte per edge rather than the kept edges themselves: the flags are
    // written by index, so no chunk has to be concatenated to put the kept
    // edges back into range order, and pass 2 walks `edges` exactly as the
    // serial loop does.
    const std::size_t nb_edges = edges.size();
    std::vector<char> kept(nb_edges, 0);

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb_edges),
      [&](const tbb::blocked_range<std::size_t>& range)
      {
        for (std::size_t ei = range.begin(); ei != range.end(); ++ei)
          kept[ei] = keep_edge(edges[ei], ei) ? 1 : 0;
      });

    std::size_t edge_index = 0;
    for (const Edge& e : edges)
    {
      if (!kept[edge_index++])
        continue;

      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t i0 = m_context->vertex_id(vh0);
      const std::size_t i1 = m_context->vertex_id(vh1);

      const bool vh0_moving = moves_vertex(vh0, i0);
      const bool vh1_moving = moves_vertex(vh1, i1);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3, boundary_edge);

      if (vh0_moving)
      {
        moves[i0].move += density * Vector_3(p0, p1);
        moves[i0].mass += density;
        ++moves[i0].neighbors;
      }
      if (vh1_moving)
      {
        moves[i1].move += density * Vector_3(p1, p0);
        moves[i1].mass += density;
        ++moves[i1].neighbors;
      }
    }
  }
#endif // CGAL_LINKED_WITH_TBB

  /**
  * The Laplacian-style move each movable endpoint of `edges` gets pulled by,
  * accumulated into `Vertex_smoothing_context::m_moves`.
  *
  * All three smooth operations do exactly this and differ only in three
  * places, which are the arguments: which edges they walk, which endpoints
  * they consider movable, and whether the edge counts as a boundary edge for
  * the sizing field. `keep_edge` is a predicate rather than a pre-filtered
  * range because the cached ranges are shared between operations that keep
  * different subsets of them. `keep_edge` gets the edge and its index in the
  * range, so that it can read a per-edge answer `refresh()` already computed.
  */
  template <typename EdgeRange, typename KeepEdge, typename MovesVertex>
  void accumulate_edge_moves(const EdgeRange& edges,
                             const C3t3& c3t3,
                             const bool boundary_edge,
                             KeepEdge keep_edge,
                             MovesVertex moves_vertex) const
  {
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (Context::is_parallel)
    {
      accumulate_edge_moves_parallel_filter(edges, c3t3, boundary_edge,
                                            keep_edge, moves_vertex);
      return;
    }
#endif

    auto& moves = m_context->m_moves;
    using Move = typename Context::Move;
    const Move default_move{CGAL::NULL_VECTOR, 0 /*neighbors*/, 0. /*mass*/};
    moves.assign(c3t3.triangulation().number_of_vertices(), default_move);

    std::size_t edge_index = 0;
    for (const Edge& e : edges)
    {
      if (!keep_edge(e, edge_index++))
        continue;

      const Vertex_handle vh0 = e.first->vertex(e.second);
      const Vertex_handle vh1 = e.first->vertex(e.third);

      const std::size_t i0 = m_context->vertex_id(vh0);
      const std::size_t i1 = m_context->vertex_id(vh1);

      const bool vh0_moving = moves_vertex(vh0, i0);
      const bool vh1_moving = moves_vertex(vh1, i1);

      if (!vh0_moving && !vh1_moving)
        continue;

      const Point_3& p0 = point(vh0->point());
      const Point_3& p1 = point(vh1->point());
      const FT density = density_along_segment(e, c3t3, boundary_edge);

      if (vh0_moving)
      {
        moves[i0].move += density * Vector_3(p0, p1);
        moves[i0].mass += density;
        ++moves[i0].neighbors;
      }
      if (vh1_moving)
      {
        moves[i1].move += density * Vector_3(p1, p0);
        moves[i1].mass += density;
        ++moves[i1].neighbors;
      }
    }
  }

  template <typename CellRange>
  Dihedral_angle_cosine max_cosine(const Tr& tr, const CellRange& cells) const
  {
    Dihedral_angle_cosine max_cos_dh = cosine_of_90_degrees();// = 0.
    for (Cell_handle c : cells)
    {
      if(!is_selected(c))
        continue;
      Dihedral_angle_cosine cos_dh = max_cos_dihedral_angle(tr, c, false);
      if (max_cos_dh < cos_dh)
        max_cos_dh = cos_dh;
    }
    return max_cos_dh;
  }

  // in flip-smooth steps, this function also checks that it improves
  // dihedral angles
  template<typename CellRange, typename Tr>
  bool check_inversion_and_move(const typename Tr::Vertex_handle v,
                                const typename Tr::Geom_traits::Point_3& final_pos,
                                const CellRange& inc_cells,
                                const Tr& tr,
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
                                FT& total_move) const
#else
                                FT&) const
#endif
  {
    const typename Tr::Point backup = v->point();//backup v's position
    const typename Tr::Geom_traits::Point_3 pv = point(backup);

    bool valid_orientation = false;
    bool angles_improved = true;
    double frac = 1.0;
    typename Tr::Geom_traits::Vector_3 move(pv, final_pos);

    const Dihedral_angle_cosine curr_max_cos = m_context->m_flip_smooth_steps
      ? max_cosine(tr, inc_cells)
      : Dihedral_angle_cosine(CGAL::ZERO, 0., 1.);//Dummy unused value

    bool valid_try = true;
    do
    {
      v->set_point(typename Tr::Point(pv + frac * move));

      valid_try = true;
      valid_orientation = true;
      angles_improved = true;

      for (const typename Tr::Cell_handle& ci : inc_cells)
      {
        if (CGAL::POSITIVE != CGAL::orientation(point(ci->vertex(0)->point()),
                                                point(ci->vertex(1)->point()),
                                                point(ci->vertex(2)->point()),
                                                point(ci->vertex(3)->point())))
        {
          frac = 0.5 * frac;
          valid_try = false;
          valid_orientation = false;
          break;
        }
        else if (m_context->m_flip_smooth_steps) //check that dihedral angles get improved
        {
          if(is_selected(ci))
          {
            Dihedral_angle_cosine max_cos_ci = max_cos_dihedral_angle(tr, ci, false);
            if (curr_max_cos < max_cos_ci)
            {
              // keep move only if new cosine is smaller than previous one
              // i.e. if angle is larger
              frac = 0.5 * frac;
              valid_try = false;
              angles_improved = false;
              break;
            }
          }
        }
      }
    }
    while(!valid_try && frac > 0.1);

    // if move failed, cancel move
    bool valid_move = valid_orientation && angles_improved;

    if(!valid_move)
      v->set_point(backup);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      total_move += CGAL::approximate_sqrt(CGAL::squared_distance(pv, point(v->point())));
#endif

    return valid_move;
  }
};

template <typename C3t3, typename SizingFunction, typename CellSelector>
class Complex_edge_vertex_smooth_operation
    : public Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>
{
public:
  using BaseClass = Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>;
  using Surface_patch_index = typename C3t3::Surface_patch_index;

  using typename BaseClass::Element_type;
  using typename BaseClass::Element_range;

protected:
  using typename BaseClass::Vertex_handle;
  using BaseClass::m_context;

  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Point_3;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;

public:
  Complex_edge_vertex_smooth_operation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

private:
  std::optional<Point_3> compute_target_position(const Vertex_handle v, const C3t3& c3t3) override
  {
    CGAL_USE(c3t3);
    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid) || !is_on_feature(v))
      return std::nullopt;

    const Point_3 current_pos = point(v->point());
    const auto& moves = m_context->m_moves;

    const std::size_t nb_neighbors = moves[vid].neighbors;
    if (nb_neighbors == 0)
      return std::nullopt;

    CGAL_assertion(moves[vid].mass > 0);
    const Vector_3 move = (nb_neighbors > 0)
                        ? moves[vid].move / moves[vid].mass
                        : CGAL::NULL_VECTOR;
    const Point_3 smoothed_position = current_pos + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
    Vector_3 sum_projections = CGAL::NULL_VECTOR;
    Point_3 tmp_pos = current_pos;

#ifndef CGAL_TET_REMESHING_EDGE_SMOOTHING_DISABLE_PROJECTION
    const std::vector<Surface_patch_index>& v_surface_indices = m_context->m_vertices_surface_indices.at(v);
    for (const Surface_patch_index& si : v_surface_indices)
    {
      Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position, current_pos,
                                                                      m_context->vertex_normal(v, si));
      sum_projections += Vector_3(tmp_pos, normal_projection);
      tmp_pos = normal_projection;
    }
#endif

    return current_pos + sum_projections;
#else
    return m_context->m_segments_aabb_tree.closest_point(smoothed_position);
#endif
  }

public:
  std::string operation_name() const override { return "Vertex Smooth (Complex Edge Vertices)"; }

private:
  void compute_vertex_moves(const C3t3& c3t3) const override
  {
    BaseClass::accumulate_edge_moves(
        m_context->m_complex_edges, c3t3, true /*boundary_edge*/,
        [](const Edge& e, const std::size_t) {
          CGAL_expensive_assertion(is_on_feature(e.first->vertex(e.second)));
          CGAL_expensive_assertion(is_on_feature(e.first->vertex(e.third)));
          CGAL_USE(e);
          return true;
        },
        [this](const Vertex_handle, const std::size_t vid) {
          return m_context->is_free(vid);
        });
  }
};

template <typename C3t3, typename SizingFunction, typename CellSelector>
class Surface_vertex_smooth_operation
    : public Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>
{
public:
  using BaseClass = Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>;

  using typename BaseClass::Element_type;
  using typename BaseClass::Element_range;

protected:
  using BaseClass::m_context;

  using typename BaseClass::AABB_triangle_tree;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Gt;
  using typename BaseClass::Point_3;
  using typename BaseClass::Surface_patch_index;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Vertex_handle;

private:
private:
  void compute_vertex_moves(const C3t3& c3t3) const override
  {
    BaseClass::accumulate_edge_moves(
        m_context->m_finite_edges, c3t3, true /*boundary_edge*/,
        [&c3t3, this](const Edge& e, const std::size_t ei) {
          CGAL_assertion(m_context->m_finite_edge_in_complex.size()
                         == m_context->m_finite_edges.size());
          return !m_context->m_finite_edge_in_complex[ei]
                 && is_boundary(c3t3, e, m_context->m_cell_selector);
        },
        [this](const Vertex_handle vh, const std::size_t vid) {
          return !is_on_feature(vh) && m_context->is_free(vid);
        });
  }

  std::optional<Point_3> project(const Surface_patch_index& si, const Point_3& gi)
  {
    CGAL_expensive_assertion(m_context->subdomain_FMLS_indices.find(si) != m_context->subdomain_FMLS_indices.end());
    CGAL_assertion(!std::isnan(gi.x()) && !std::isnan(gi.y()) && !std::isnan(gi.z()));

    Vector_3 point_vec(gi.x(), gi.y(), gi.z());
    Vector_3 res_normal = CGAL::NULL_VECTOR;
    Vector_3 result(CGAL::ORIGIN, gi);

    const typename BaseClass::Context::FMLS& fmls = m_context->subdomain_FMLS[m_context->subdomain_FMLS_indices.at(si)];

    int it_nb = 0;
    const int max_it_nb = 5;
    const double epsilon = fmls.getPNScale() / 1000.;
    const double sq_eps = CGAL::square(epsilon);

    do
    {
      point_vec = result;
      fmls.fastProjectionCPU(point_vec, result, res_normal);
      if(std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "MLS error detected si "
                  << "\t(size : " << fmls.getPNSize() << ")"
                  << "\t(point = " << point_vec << " )" << std::endl;
#endif
        return {};
      }
    } while((result - point_vec).squared_length() > sq_eps && ++it_nb < max_it_nb);

    return Point_3(result[0], result[1], result[2]);
  }

public:
  Surface_vertex_smooth_operation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

private:
  std::optional<Point_3> compute_target_position(const Vertex_handle v, const C3t3& c3t3) override
  {
    auto& tr = c3t3.triangulation();
    auto& moves = m_context->m_moves;
    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid) || v->in_dimension() != 2)
      return std::nullopt;

    const std::size_t nb_neighbors = moves[vid].neighbors;
    const Point_3 current_pos = point(v->point());

    CGAL_assertion(m_context->m_vertices_surface_indices.find(v) != m_context->m_vertices_surface_indices.end());
    const auto& incident_surface_patches = m_context->m_vertices_surface_indices.at(v);

    if (incident_surface_patches.size() > 1)
      return std::nullopt;

    const Surface_patch_index si = incident_surface_patches[0];
    CGAL_assertion(si != Surface_patch_index());
    CGAL_expensive_assertion_code(auto siv = surface_patch_index(v, c3t3));
    CGAL_expensive_assertion(si == siv);

    Point_3 new_pos;

    if (nb_neighbors > 1)
    {
      const Vector_3 move = moves[vid].move / moves[vid].mass;
      const Point_3 smoothed_position = point(v->point()) + move;

#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      Point_3 normal_projection = BaseClass::project_on_tangent_plane(smoothed_position, current_pos,
                                                                      m_context->vertex_normal(v, si));
      std::optional<Point_3> mls_projection = project(si, normal_projection);
      new_pos = (mls_projection != std::nullopt) ? *mls_projection : smoothed_position;
#else
      if(m_context->m_triangles_aabb_tree.squared_distance(smoothed_position) < m_context->m_aabb_epsilon)
      {
        new_pos = m_context->m_triangles_aabb_tree.closest_point(smoothed_position);
      }
      else
      {
        using Ray = typename Tr::Geom_traits::Ray_3;
        using Projection =
            std::optional<typename AABB_triangle_tree::template Intersection_and_primitive_id<Ray>::Type>;

        auto get_intersection_point = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          if(const Point_3* p = std::get_if<Point_3>(&intersection))
            return *p;
          return std::nullopt;
        };

        auto get_intersection_midpoint = [](const Projection& proj) -> std::optional<Point_3> {
          const auto intersection = proj.value().first;
          using Segment = typename Tr::Geom_traits::Segment_3;
          if(const Segment* s = std::get_if<Segment>(&intersection))
            return CGAL::midpoint(s->source(), s->target());
          CGAL_assertion(false);
          return std::nullopt;
        };

        const auto n = m_context->vertex_normal(v, si);
        const Ray ray = tr.geom_traits().construct_ray_3_object()(current_pos, n);
        const Projection proj = m_context->m_triangles_aabb_tree.first_intersection(ray);
        const Projection proj_opp = m_context->m_triangles_aabb_tree.first_intersection(
            tr.geom_traits().construct_opposite_ray_3_object()(ray));

        if(proj != std::nullopt && proj_opp == std::nullopt)
        {
          const auto p = get_intersection_point(proj);
          new_pos = (p != std::nullopt) ? p.value() : get_intersection_midpoint(proj).value();
        }
        else if(proj == std::nullopt && proj_opp != std::nullopt)
        {
          const auto p = get_intersection_point(proj_opp);
          new_pos = (p != std::nullopt)
                    ? p.value()
                    : get_intersection_midpoint(proj_opp).value();
        }
        else if(proj != std::nullopt && proj_opp != std::nullopt)
        {
          const auto op1 = get_intersection_point(proj);
          const auto op2 = get_intersection_point(proj_opp);

          const FT sqd1 = (op1 == std::nullopt) ? 0.
            : CGAL::squared_distance(smoothed_position, op1.value());
          const FT sqd2 = (op2 == std::nullopt) ? 0.
            : CGAL::squared_distance(smoothed_position, op2.value());

          if (sqd1 != 0. && sqd1 < sqd2)
            new_pos = op1.value();
          else if (sqd2 != 0)
            new_pos = op2.value();
          else
            new_pos = smoothed_position;
        }
        else //no valid projection
        {
          new_pos = smoothed_position;
        }
      }
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      return new_pos;
    }
    else if (nb_neighbors > 0)
    {
#ifdef CGAL_TET_REMESHING_SMOOTHING_WITH_MLS
      std::optional<Point_3> mls_proj = project(si, current_pos);
      if(mls_proj == std::nullopt)
        return std::nullopt;

      new_pos = *mls_proj;
#else // AABB_tree projection
      new_pos = m_context->m_segments_aabb_tree.closest_point(current_pos);
#endif //CGAL_TET_REMESHING_SMOOTHING_WITH_MLS

      return new_pos;
    }

    return std::nullopt;
  }

public:
  std::string operation_name() const override { return "Vertex Smooth (Surface Vertices)"; }
};

template <typename C3t3, typename SizingFunction, typename CellSelector>
class Internal_vertex_smooth_operation
    : public Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>
{
public:
  using BaseClass = Vertex_smooth_operation_base<C3t3, SizingFunction, CellSelector>;

  using typename BaseClass::Element_type;
  using typename BaseClass::Element_range;

protected:
  using BaseClass::m_context;

  using typename BaseClass::Cell_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::FT;
  using typename BaseClass::Point_3;
  using typename BaseClass::Tr;
  using typename BaseClass::Vector_3;
  using typename BaseClass::Vertex_handle;

public:
  Internal_vertex_smooth_operation(std::shared_ptr<typename BaseClass::Context> context)
      : BaseClass(context) {}

private:
  void compute_vertex_moves(const C3t3& c3t3) const override
  {
    /*for dim 3 vertices, start counting neighbors directly from 0*/
    BaseClass::accumulate_edge_moves(
        m_context->m_finite_edges, c3t3, false /*boundary_edge*/,
        [&c3t3, this](const Edge& e, const std::size_t) {
          return !is_outside(e, c3t3, m_context->m_cell_selector);
        },
        [&c3t3, this](const Vertex_handle vh, const std::size_t vid) {
          return c3t3.in_dimension(vh) == 3 && m_context->is_free(vid);
        });
  }

private:
  std::optional<Point_3> compute_target_position(const Vertex_handle v, const C3t3& c3t3) override
  {
    auto& moves = m_context->m_moves;

    const std::size_t vid = m_context->vertex_id(v);
    if (!m_context->is_free(vid))
      return std::nullopt;

    if (c3t3.in_dimension(v) == 3 && moves[vid].neighbors > 1)
    {
      const Vector_3 move = moves[vid].move / moves[vid].mass;
      return point(v->point()) + move;
    }
    return std::nullopt;
  }

public:
  std::string operation_name() const override { return "Vertex Smooth (Internal Vertices)"; }
};

}//namespace internal
}//namespace Tetrahedral_adaptive_remeshing
}//namespace CGAL

#endif //CGAL_INTERNAL_SMOOTH_VERTICES_H
