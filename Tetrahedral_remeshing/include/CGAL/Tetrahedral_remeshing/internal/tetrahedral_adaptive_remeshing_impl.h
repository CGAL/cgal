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

#ifndef TETRAHEDRAL_REMESHING_IMPL_H
#define TETRAHEDRAL_REMESHING_IMPL_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#endif

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/split_long_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/collapse_short_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/flip_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/smooth_vertices.h>
#include <CGAL/Tetrahedral_remeshing/internal/peel_slivers.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Tetrahedral_remeshing/internal/compute_c3t3_statistics.h>

#include <optional>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
class Default_remeshing_visitor
{
public:
  template<typename Tr>
  void before_split(const Tr& /* tr */, const typename Tr::Edge& /* e */) {}
  template<typename Tr>
  void after_split(const Tr& /* tr */, const typename Tr::Vertex_handle /* new_v */) {}

  template<typename CellHandleOld, typename CellHandleNew>
  void after_add_cell(CellHandleOld /* co */, CellHandleNew /* cn */) const {}

  template<typename CellHandle>
  void before_flip(const CellHandle /* c */) {}
  template<typename CellHandle>
  void after_flip(CellHandle /* c */) {}
};

template<typename Tr>
struct All_cells_selected
{
  using key_type = typename Tr::Cell_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;

  friend value_type get(const All_cells_selected&, const key_type& c)
  {
    using SI = typename Tr::Cell::Subdomain_index;
    return c->subdomain_index() != SI();
  }
  friend void put(All_cells_selected&, const key_type&, const value_type)
  {} //nothing to do : subdomain indices are updated in remeshing};
};


template<typename Triangulation
         , typename SizingFunction
         , typename EdgeIsConstrainedMap
         , typename FacetIsConstrainedMap
         , typename CellSelector
         , typename Visitor
         , typename CornerIndex = int
         , typename CurveIndex = int
         >
class Adaptive_remesher
{
  typedef Triangulation Tr;
  typedef typename Tr::Geom_traits::FT FT;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr, CornerIndex, CurveIndex> C3t3;

  typedef typename C3t3::Cell_handle         Cell_handle;
  typedef typename C3t3::Vertex_handle       Vertex_handle;
  typedef typename C3t3::Edge                Edge;
  typedef typename C3t3::Subdomain_index     Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Curve_index         Curve_index;
  typedef typename C3t3::Corner_index        Corner_index;

  typedef Tetrahedral_remeshing_smoother<C3t3> Smoother;

private:
  C3t3 m_c3t3;
  const SizingFunction& m_sizing;
  const bool m_protect_boundaries;
  CellSelector m_cell_selector;
  Visitor& m_visitor;
  Smoother m_vertex_smoother;//initialized with initial surface

  C3t3* m_c3t3_pbackup;
  Triangulation* m_tr_pbackup; //backup to re-swap triangulations when done

public:
  Adaptive_remesher(Triangulation& tr
                    , const SizingFunction& sizing
                    , const bool protect_boundaries
                    , EdgeIsConstrainedMap ecmap
                    , FacetIsConstrainedMap fcmap
                    , bool smooth_constrained_edges
                    , CellSelector cell_selector
                    , Visitor& visitor
                   )
    : m_c3t3()
    , m_sizing(sizing)
    , m_protect_boundaries(protect_boundaries)
    , m_cell_selector(cell_selector)
    , m_visitor(visitor)
    , m_c3t3_pbackup(NULL)
    , m_tr_pbackup(&tr)
  {
    m_c3t3.triangulation().swap(tr);

    init_c3t3(ecmap, fcmap);
    m_vertex_smoother.init(m_c3t3, m_cell_selector, smooth_constrained_edges);

#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "00-init");
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "00-facets_in_complex_after_init.off");
#endif
  }

  Adaptive_remesher(C3t3& c3t3
                    , const SizingFunction& sizing
                    , const bool protect_boundaries
                    , EdgeIsConstrainedMap ecmap
                    , FacetIsConstrainedMap fcmap
                    , bool smooth_constrained_edges
                    , CellSelector cell_selector
                    , Visitor& visitor
                   )
    : m_c3t3()
    , m_sizing(sizing)
    , m_protect_boundaries(protect_boundaries)
    , m_cell_selector(cell_selector)
    , m_visitor(visitor)
    , m_c3t3_pbackup(&c3t3)
    , m_tr_pbackup(NULL)
  {
    m_c3t3.swap(c3t3);

    init_c3t3(ecmap, fcmap);
    m_vertex_smoother.init(m_c3t3, m_cell_selector, smooth_constrained_edges);

#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "00-init");
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "00-facets_in_complex_after_init.off");
#endif
  }

  bool input_is_c3t3() const
  {
    return m_c3t3_pbackup != NULL;
  }

  void split()
  {
    CGAL_assertion(check_vertex_dimensions());

    const FT target_edge_length = m_sizing(CGAL::ORIGIN);
    const FT emax = FT(4)/FT(3) * target_edge_length;
    split_long_edges(m_c3t3, emax, m_protect_boundaries,
                     m_cell_selector, m_visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3,
      "1-facets_in_complex_after_split.off");
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "1-c3t3_vertices_after_split");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "1-split");
#endif
  }

  void collapse()
  {
    CGAL_assertion(check_vertex_dimensions());

    const FT target_edge_length = m_sizing(CGAL::ORIGIN);
    FT emin = FT(4)/FT(5) * target_edge_length;
    FT emax = FT(4)/FT(3) * target_edge_length;
    collapse_short_edges(m_c3t3, emin, emax, m_protect_boundaries,
                         m_cell_selector, m_visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "2-c3t3_vertices_after_collapse");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "2-collapse");
#endif
  }

  void flip()
  {
    flip_edges(m_c3t3, m_protect_boundaries,
               m_cell_selector, m_visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "3-c3t3_vertices_after_flip");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "3-flip");
#endif
  }

  void smooth()
  {
    m_vertex_smoother.smooth_vertices(m_c3t3, m_protect_boundaries, m_cell_selector);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "4-c3t3_vertices_after_smooth");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "4-smooth");
#endif
  }

  bool resolution_reached()
  {
    const FT target_edge_length = m_sizing(CGAL::ORIGIN);

    FT emax = FT(4) / FT(3) * target_edge_length;
    FT emin = FT(4) / FT(5) * target_edge_length;

    FT sqmax = emax * emax;
    FT sqmin = emin * emin;

    for (const Edge& e : tr().finite_edges())
    {
      // skip protected edges
      if (m_protect_boundaries)
      {
        if(  m_c3t3.is_in_complex(e)
             || is_boundary(m_c3t3, e, m_cell_selector))
          continue;
      }

      FT sqlen = tr().segment(e).squared_length();
      if (sqlen < sqmin || sqlen > sqmax)
        return false;
    }
    std::cout << "Resolution reached" << std::endl;
    return true;
  }

  //peel off slivers
  std::size_t postprocess(const double sliver_angle = 2.)
  {
    if (m_protect_boundaries)
      return 0;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Postprocess...";
    std::cout.flush();
#endif

    const std::size_t nb_peeled
      = CGAL::Tetrahedral_remeshing::peel_slivers(m_c3t3, sliver_angle, m_cell_selector);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL_assertion(tr().tds().is_valid(true));
    CGAL_assertion(debug::are_cell_orientations_valid(tr()));
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
    CGAL::Tetrahedral_remeshing::debug::dump_c3t3(m_c3t3, "99-postprocess");
#endif

    return nb_peeled;
  }

  void finalize()
  {
    if (m_c3t3_pbackup != NULL)
      m_c3t3_pbackup->swap(m_c3t3);
    else
      m_tr_pbackup->swap(m_c3t3.triangulation());
  }

private:
  void init_c3t3(const EdgeIsConstrainedMap& ecmap,
                 const FacetIsConstrainedMap& fcmap)
  {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    debug_c3t3();
    std::size_t nbc = 0;
    std::size_t nbf = 0;
    std::size_t nbe = 0;
    std::size_t nbv = 0;
#endif
    //update number_of_cells and number_of_facets in c3t3
    m_c3t3.rescan_after_load_of_triangulation();

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "00-c3t3_vertices_before_init_");
#endif

    const Subdomain_index default_subdomain = default_subdomain_index();

    //tag cells
    for (Cell_handle cit : tr().finite_cell_handles())
    {
      if (get(m_cell_selector, cit))
      {
        const Subdomain_index index = cit->subdomain_index();
        if (m_c3t3.is_in_complex(cit))
          m_c3t3.remove_from_complex(cit);

        const Subdomain_index new_index = (Subdomain_index() != index)
          ? index
          : default_subdomain;
        m_c3t3.add_to_complex(cit, new_index);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbc;
#endif
      }

      for (Vertex_handle vi : CGAL::Tetrahedral_remeshing::vertices(cit, tr()))
        set_dimension(vi, 3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      //else
        if (input_is_c3t3() && m_c3t3.is_in_complex(cit))
        ++nbc;
#endif
    }

    //tag facets
    typedef typename Tr::Facet Facet;
    for (const Facet& f : tr().finite_facets())
    {
      const Facet mf = tr().mirror_facet(f);
      const Subdomain_index s1 = f.first->subdomain_index();
      const Subdomain_index s2 = mf.first->subdomain_index();
      if (s1 != s2
          || get(fcmap, f)
          || get(fcmap, mf)
          || (m_c3t3_pbackup == NULL && f.first->is_facet_on_surface(f.second)))
      {
        Surface_patch_index patch = f.first->surface_patch_index(f.second);
        if(patch == Surface_patch_index())
          make_surface_patch_index(s1, s2, patch);

        if(m_c3t3.is_in_complex(f))
          m_c3t3.remove_from_complex(f);
        m_c3t3.add_to_complex(f, patch);

        for (Vertex_handle vij : CGAL::Tetrahedral_remeshing::vertices(f, tr()))
          set_dimension(vij, 2);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbf;
#endif
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3, "facets_in_complex.off");
#endif

    //tag edges
    const Curve_index default_curve_id = default_curve_index();
    for (const Edge& e : tr().finite_edges())
    {
      if (get(ecmap, CGAL::Tetrahedral_remeshing::make_vertex_pair<Tr>(e))
          || m_c3t3.is_in_complex(e)
          || nb_incident_subdomains(e, m_c3t3) > 2
          || nb_incident_surface_patches(e, m_c3t3) > 1)
      {
        const bool in_complex = m_c3t3.is_in_complex(e);
        typename C3t3::Curve_index curve_id = in_complex
          ? m_c3t3.curve_index(e)
          : default_curve_id;

        if (in_complex)
          m_c3t3.remove_from_complex(e);
        m_c3t3.add_to_complex(e, curve_id);

        for (Vertex_handle v : CGAL::Tetrahedral_remeshing::vertices(e, tr()))
          set_dimension(v, 1);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbe;
#endif
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3, "edges_in_complex.polylines.txt");
#endif

    //tag vertices
    Corner_index corner_id = 0;
    for (Vertex_handle vit : tr().finite_vertex_handles())
    {
      if ( vit->in_dimension() == 0
           || nb_incident_complex_edges(vit, m_c3t3) > 2)
      {
        if (!m_c3t3.is_in_complex(vit))
          m_c3t3.add_to_complex(vit, ++corner_id);

        set_dimension(vit, 0);
        vit->set_index(corner_id);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ++nbv;
#endif
      }
    }

    // set all indices depending on underlying dimension
    for (Vertex_handle v : tr().finite_vertex_handles())
      set_index(v, m_c3t3);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::cout << "C3t3 ready :" << std::endl;
    std::cout << "\t cells    = " << nbc << std::endl;
    std::cout << "\t facets   = " << nbf << std::endl;
    std::cout << "\t edges    = " << nbe << std::endl;
    std::cout << "\t vertices = " << nbv << std::endl;

    CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
      m_c3t3.triangulation(), "0-c3t3_vertices_after_init_");
    CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
    CGAL::Tetrahedral_remeshing::debug::count_far_points(m_c3t3);
#endif
  }

private:
  bool dimension_is_modifiable(const Vertex_handle& v, const int new_dim) const
  {
    if (!input_is_c3t3())
      return true;

    const int vdim = v->in_dimension();
    // feature edges and tip/endpoints vertices are kept
    switch (vdim)
    {
    case -1: return false;//far points are not modified
    case 3 : return true;
    case 2 : return true;//surface vertices may not be part of a triangulation surface
                        // in this case, we want to be able to set it
    case 1 : return new_dim == 0; //features can be modified to corners
    case 0 : return false;// corners remain corners
    default:
      return true;
    }
    CGAL_unreachable();
    return true;
  }

  void set_dimension(Vertex_handle v, const int new_dim)
  {
    if (dimension_is_modifiable(v, new_dim))
      v->set_dimension(new_dim);
  }

  bool check_vertex_dimensions()
  {
    for (Vertex_handle vit : tr().finite_vertex_handles())
    {
      // dimension is -1 for Mesh_3 "far points"
      // for other vertices, it is in [0; 3]
      if (vit->in_dimension() < -1 || vit->in_dimension() > 3)
        return false;
    }
    return true;
  }
  void debug_c3t3()
  {
    for (typename Tr::Facet f : tr().finite_facets())
    {
      typename Tr::Facet mf = tr().mirror_facet(f);
      CGAL_assertion(m_c3t3.is_in_complex(f) == m_c3t3.is_in_complex(mf));
    }
  }

  template<typename PatchIndex>
  void make_surface_patch_index(const Subdomain_index& s1,
                                const Subdomain_index& s2,
                                PatchIndex& patch)
  {
    patch = (s1 < s2) ? (s1 * 1000 + s2) : (s2 * 1000 + s1);
  }

  void make_surface_patch_index(const Subdomain_index& s1,
                                const Subdomain_index& s2,
                                std::pair<Subdomain_index, Subdomain_index>& patch)
  {
    patch = (s1 < s2) ? std::make_pair(s1, s2) : std::make_pair(s2, s1);
  }

  Subdomain_index max_subdomain_index() const
  {
    Subdomain_index max_index
      = (std::numeric_limits<Subdomain_index>::min)();
    for (Cell_handle cit : tr().finite_cell_handles())
    {
      const Subdomain_index cid = cit->subdomain_index();
      if (cid > max_index && cid != Subdomain_index())
        max_index = cid;
    }
    return max_index;
  }

  Subdomain_index default_subdomain_index() const
  {
    return max_subdomain_index() + 1;
  }

  Curve_index max_curve_index() const
  {
    Curve_index max_index = (std::numeric_limits<Curve_index>::min)();
    for (const Edge& e : m_c3t3.edges_in_complex())
    {
      const Curve_index cid = m_c3t3.curve_index(e);
      if (cid > max_index)
        max_index = cid;
    }
    return max_index + 1;
  }

  Curve_index default_curve_index() const
  {
    return max_curve_index() + 1;
  }

public:
  Tr& tr()
  {
    return m_c3t3.triangulation();
  }
  const Tr& tr() const
  {
    return m_c3t3.triangulation();
  }

  void remesh(const std::size_t& max_it,
              const std::size_t& nb_extra_iterations)
  {
    std::size_t it_nb = 0;
    while (it_nb < max_it)
    {
      ++it_nb;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " #" << std::endl;
#endif
      if (!resolution_reached())
      {
        split();
        collapse();
      }
      flip();
      smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " done : "
                << tr().number_of_vertices()
                << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
        tr(), m_cell_selector, ossi.str().c_str());
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      CGAL::Tetrahedral_remeshing::debug::check_surface_patch_indices(m_c3t3);
#endif
    }

    while (it_nb < max_it + nb_extra_iterations)
    {
      ++it_nb;

      flip();
      smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "# Iteration " << it_nb << " (flip and smooth only) done : "
                << tr().number_of_vertices()
                << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
      std::ostringstream ossi;
      ossi << "statistics_" << it_nb << ".txt";
      Tetrahedral_remeshing::internal::compute_statistics(
        tr(),  m_cell_selector, ossi.str().c_str());
#endif
    }

    postprocess(); //peel off boundary slivers

    finalize();
    //Warning : triangulation() is now empty
  }

};//end class Adaptive_remesher

}//end namespace internal
}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //TETRAHEDRAL_REMESHING_IMPL_H
