// Copyright (c) 2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois

#ifndef TETRAHEDRAL_REMESHING_IMPL_H
#define TETRAHEDRAL_REMESHING_IMPL_H

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#endif

#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Triangulation_utils_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/add_imaginary_layer.h>
#include <CGAL/Tetrahedral_remeshing/internal/split_long_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/collapse_short_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/flip_edges.h>
#include <CGAL/Tetrahedral_remeshing/internal/smooth_vertices.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

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
    void before_split(const Tr& tr, const typename Tr::Edge& e) {}
    template<typename Tr>
    void after_split(const Tr& tr, const typename Tr::Vertex_handle new_v) {}

    template<typename CellHandleOld, typename CellHandleNew>
    void after_add_cell(CellHandleOld co, CellHandleNew cn) const {}

    template<typename CellHandle>
    void before_flip(const CellHandle c) {}
    template<typename CellHandle>
    void after_flip(CellHandle c) {}
  };

  template<typename Tr>
  struct All_cells_selected
  {
    typedef typename Tr::Cell_handle argument_type;
    typedef typename Tr::Cell::Subdomain_index Subdomain_index;

    typedef bool                     result_type;

    result_type operator()(const argument_type c) const
    {
      return c->subdomain_index() != Subdomain_index();
    }
  };

  template<typename Primitive>
  struct No_constraint_pmap
  {
  public:
    typedef Primitive                           key_type;
    typedef bool                                value_type;
    typedef value_type&                         reference;
    typedef boost::read_write_property_map_tag  category;

    friend bool get(const No_constraint_pmap&, const key_type&) {
      return false;
    }
    friend void put(No_constraint_pmap&, const key_type&, const bool) {}
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
    typedef typename C3t3::Subdomain_index     Subdomain_index;
    typedef typename C3t3::Surface_patch_index Surface_patch_index;

  private:
    C3t3 m_c3t3;
    const SizingFunction& m_sizing;
    const bool m_protect_boundaries;
    CellSelector m_cell_selector;
    Subdomain_index m_imaginary_index;
    Visitor& m_visitor;

    Triangulation* m_tr_pbackup; //backup to re-swap triangulations when done
    C3t3* m_c3t3_pbackup;

  public:
    Adaptive_remesher(Triangulation& tr
      , const SizingFunction& sizing
      , const bool protect_boundaries
      , EdgeIsConstrainedMap ecmap
      , FacetIsConstrainedMap fcmap
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

#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_without_imaginary(m_c3t3.triangulation(),
        "00-init-no-imaginary.mesh", m_imaginary_index);
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "00-init.binary.cgal");
#endif
    }

    Adaptive_remesher(C3t3& c3t3
      , const SizingFunction& sizing
      , const bool protect_boundaries
      , EdgeIsConstrainedMap ecmap
      , FacetIsConstrainedMap fcmap
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

#ifdef CGAL_DUMP_REMESHING_STEPS
      //CGAL::Tetrahedral_remeshing::debug::dump_without_imaginary(m_c3t3.triangulation(),
      //  "00-init-no-imaginary.mesh", m_imaginary_index);
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "00-init.binary.cgal");
#endif
    }

    const Subdomain_index& imaginary_index() const
    {
      return m_imaginary_index;
    }

//    void preprocess()
//    {
//#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//      std::cout << "Preprocess...";
//      std::cout.flush();
//#endif
//
//      add_layer_of_imaginary_tets(tr(), m_imaginary_index);
//      CGAL_assertion(tr().tds().is_valid(true));
//
//#ifdef CGAL_DUMP_REMESHING_STEPS
//      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(), "0-preprocess.mesh");
//      CGAL::Tetrahedral_remeshing::debug::dump_without_imaginary(tr(),
//        "0-preprocess-no-imaginary.mesh", m_imaginary_index);
//      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "0-preprocess.binary.cgal");
//#endif
//#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//      std::cout << "done." << std::endl;
//#endif
//    }

    void split()
    {
      CGAL_assertion(check_vertex_dimensions());

      const FT target_edge_length = m_sizing(CGAL::ORIGIN);
      const FT emax = FT(4)/FT(3) * target_edge_length;
      split_long_edges(m_c3t3, emax, m_protect_boundaries,
                       m_cell_selector, m_visitor);

      CGAL_assertion(tr().tds().is_valid(true));
#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(), "1-split.mesh");
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "1-split.binary.cgal");
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

      CGAL_assertion(tr().tds().is_valid(true));
#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(),
        "2-collapse.mesh");
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "2-collapse.binary.cgal");
#endif
    }

    void flip()
    {
      flip_edges(m_c3t3, m_protect_boundaries,
                 m_cell_selector, m_visitor);

      CGAL_assertion(tr().tds().is_valid(true));
#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(), "3-flip.mesh");
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "3-flip.binary.cgal");
#endif
    }

    void smooth()
    {
      smooth_vertices_new(m_c3t3, m_imaginary_index, m_protect_boundaries,
                      m_cell_selector);

      CGAL_assertion(tr().tds().is_valid(true));
#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(),
        "4-smooth.mesh");
      CGAL::Tetrahedral_remeshing::debug::dump_without_imaginary(tr(),
        "4-smooth-no-imaginary.mesh", m_imaginary_index);
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "4-smooth.binary.cgal");
#endif
    }

    bool resolution_reached()
    {
      const FT target_edge_length = m_sizing(CGAL::ORIGIN);

      FT emax = FT(4) / FT(3) * target_edge_length;
      FT emin = FT(4) / FT(5) * target_edge_length;

      FT sqmax = emax * emax;
      FT sqmin = emin * emin;

      typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
      for (Finite_edges_iterator eit = tr().finite_edges_begin();
           eit != tr().finite_edges_end();
           ++eit)
      {
        typename Tr::Edge e = *eit;
        // skip protected edges
        if (m_protect_boundaries)
        {
          if(  m_c3t3.is_in_complex(e)
            || is_boundary(m_c3t3, e, m_cell_selector))
            continue;
        }
        // skip imaginary edges
        if (is_imaginary(e, m_c3t3, m_imaginary_index))
          continue;

        FT sqlen = tr().segment(e).squared_length();
        if (sqlen < sqmin || sqlen > sqmax)
          return false;
      }
      std::cout << "Resolution reached" << std::endl;
      return true;
    }

    void postprocess()
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "Postprocess...";
      std::cout.flush();
#endif
      //unset imaginary cells
      typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
      for (Finite_cells_iterator cit = tr().finite_cells_begin();
           cit != tr().finite_cells_end(); ++cit)
      {
        if (cit->subdomain_index() == m_imaginary_index)
        {
          m_c3t3.remove_from_complex(cit);
        }
      }

      CGAL_assertion(tr().tds().is_valid(true));
#ifdef CGAL_DUMP_REMESHING_STEPS
      CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(tr(),
        "99-postprocess.mesh");
      CGAL::Tetrahedral_remeshing::debug::dump_without_imaginary(tr(),
        "99-postprocess-no-imaginary.mesh", m_imaginary_index);
      CGAL::Tetrahedral_remeshing::debug::dump_binary(m_c3t3, "99-postprocess.binary.cgal");
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "done." << std::endl;
#endif
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
      std::size_t nbc = 0;
      std::size_t nbf = 0;
      std::size_t nbe = 0;
      std::size_t nbv = 0;
#endif

      Subdomain_index max_si = 0;

      //tag cells (no imaginary cell yet)
      typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
      for (Finite_cells_iterator cit = tr().finite_cells_begin();
           cit != tr().finite_cells_end();
           ++cit)
      {
        if (m_cell_selector(cit))
        {
          m_c3t3.add_to_complex(cit, cit->subdomain_index());
          max_si = (std::max)(max_si, cit->subdomain_index());
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
          ++nbc;
#endif
        }
        for (int i = 0; i < 4; ++i)
        {
          if (cit->vertex(i)->in_dimension() == -1)
            cit->vertex(i)->set_dimension(3);
        }
      }
      m_imaginary_index = max_si + 1;
      if(max_si == 0)
        std::cerr << "Warning : Maximal subdomain index is 0" << std::endl
                  << "          Remeshing is likely to fail." << std::endl;

      //tag facets
      typedef typename Tr::Facet                  Facet;
      typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
      for (Finite_facets_iterator fit = tr().finite_facets_begin();
           fit != tr().finite_facets_end();
           ++fit)
      {
        Facet f = *fit;
        Facet mf = tr().mirror_facet(f);
        Subdomain_index s1 = f.first->subdomain_index();
        Subdomain_index s2 = mf.first->subdomain_index();
        if ( s1 != s2
          || get(fcmap, f)
          || get(fcmap, mf)
          || (m_c3t3_pbackup == NULL && f.first->is_facet_on_surface(f.second)))
        {
          m_c3t3.add_to_complex(f, 1);

          const int i = f.second;
          for (int j = 0; j < 3; ++j)
          {
            Vertex_handle vij = f.first->vertex(Tr::vertex_triple_index(i, j));
            if (vij->in_dimension() == -1 || vij->in_dimension() > 2)
              vij->set_dimension(2);
          }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
          ++nbf;
#endif
        }
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(m_c3t3, "facets_in_complex.off");
#endif

      //tag edges
      typedef typename Tr::Edge                  Edge;
      typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
      for (Finite_edges_iterator eit = tr().finite_edges_begin();
           eit != tr().finite_edges_end();
           ++eit)
      {
        Edge e = *eit;
        if (get(ecmap, CGAL::Tetrahedral_remeshing::make_vertex_pair<Tr>(e))
          || nb_incident_subdomains(e, m_c3t3) > 2)
        {
          m_c3t3.add_to_complex(e, 1);

          Vertex_handle v = e.first->vertex(e.second);
          if(v->in_dimension() == -1 || v->in_dimension() > 1)
            v->set_dimension(1);

          v = e.first->vertex(e.third);
          if (v->in_dimension() == -1 || v->in_dimension() > 1)
            v->set_dimension(1);
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
          ++nbe;
#endif
        }
      }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      CGAL::Tetrahedral_remeshing::debug::dump_edges_in_complex(m_c3t3, "edges_in_complex.polylines.txt");
#endif

      //tag vertices
      typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
      unsigned int corner_id = 0;
      for (Finite_vertices_iterator vit = tr().finite_vertices_begin();
           vit != tr().finite_vertices_end();
           ++vit)
      {
        if ( vit->in_dimension() == 0
          || nb_incident_complex_edges(vit, m_c3t3) > 2)
        {
          m_c3t3.add_to_complex(vit, ++corner_id);

          if (vit->in_dimension() == -1 || vit->in_dimension() > 0)
            vit->set_dimension(0);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
          ++nbv;
#endif
        }
      }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      std::cout << "C3t3 ready :" << std::endl;
      std::cout << "\t cells    = " << nbc << std::endl;
      std::cout << "\t facets   = " << nbf << std::endl;
      std::cout << "\t edges    = " << nbe << std::endl;
      std::cout << "\t vertices = " << nbv << std::endl;

      CGAL::Tetrahedral_remeshing::debug::dump_vertices_by_dimension(
        m_c3t3.triangulation(), "c3t3_vertices_");
#endif
    }

  private:

    bool check_vertex_dimensions()
    {
      typename Tr::Finite_vertices_iterator vit;
      for (vit = tr().finite_vertices_begin();
           vit != tr().finite_vertices_end(); ++vit)
      {
        if (vit->in_dimension() < 0 || vit->in_dimension() > 3)
          return false;
      }
      return true;
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
//      preprocess();

      std::size_t it_nb = 0;
      while (it_nb++ < max_it)
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "# Iteration " << it_nb << " #" << std::endl;
#endif
        if (!resolution_reached())
        {
          split();
          collapse();
        }
        flip();
//        smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "# Iteration " << it_nb << " done : "
          << tr().number_of_vertices()
          << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
        std::ostringstream ossi;
        ossi << "statistics_" << it_nb << ".txt";
        Tetrahedral_remeshing::internal::compute_statistics(
          tr(), imaginary_index(), m_cell_selector, ossi.str().c_str());
#endif
      }

      while (it_nb++ < max_it + nb_extra_iterations)
      {
        //      flip();
        //      smooth();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        std::cout << "# Iteration " << it_nb << " (flip and smooth only) done : "
          << tr().number_of_vertices()
          << " vertices #" << std::endl;
#endif
#ifdef CGAL_DUMP_REMESHING_STEPS
        std::ostringstream ossi;
        ossi << "statistics_" << it_nb << ".txt";
        Tetrahedral_remeshing::internal::compute_statistics(
          tr(), imaginary_index(), m_cell_selector, ossi.str().c_str());
#endif
      }

//      postprocess(); //remove imaginary cells

      finalize();
      //triangulation() is now empty
    }

  };//end class Adaptive_remesher
}//end namespace internal
}//end namespace Tetrahedral_remeshing
}//end namespace CGAL

#endif //TETRAHEDRAL_REMESHING_IMPL_H
