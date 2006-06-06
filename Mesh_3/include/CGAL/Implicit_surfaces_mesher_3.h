// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_IMPLICIT_SURFACES_MESHER_3_H
#define CGAL_IMPLICIT_SURFACES_MESHER_3_H

#include <CGAL/Mesh_3/Refine_tets.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesher/Surface_mesher_visitor.h>
#include <CGAL/Mesh_3/Implicit_surface_mesher_visitor.h>
#include <CGAL/Mesh_3/Refine_tets_visitor.h>

namespace CGAL {

  namespace Mesh_3 
  {
    template <class C2T3>
    bool check_c2t3(const C2T3& c2t3)
    {
      typedef typename C2T3::Triangulation Tr;
      typedef typename Tr::Cell_handle Cell_handle;

      const Tr& tr = c2t3.triangulation();

      bool result = true;
      // first check that all vertices on the c2t3 have index > 0.
      for (typename Tr::Finite_facets_iterator fit = tr.finite_facets_begin();
           fit != tr.finite_facets_end(); ++fit)
        if (fit->first->is_facet_on_surface (fit->second))
        {
          const typename Tr::Cell_handle& c = fit->first;
          const int& index = fit->second;
          for(int i = 0; i < 4; ++i)
            if(i != index && c->vertex(i)->point().surface_index() == 0)
            {
              std::cerr << "error: vertex " << &*(c->vertex(i)) 
                        << " should have surface_index() > 0\n";
              result = false;
            }
        }

      typename Tr::size_type point_on_surface_count = 0;
      for(typename Tr::Finite_vertices_iterator vit = 
            tr.finite_vertices_begin();
          vit != tr.finite_vertices_end();
          ++vit)
        if( vit->point().surface_index() == 0 )
        {
          std::vector<Cell_handle> cells;
          cells.reserve(64);
        
          tr.incident_cells(vit, std::back_inserter(cells));
          for(typename std::vector<Cell_handle>::iterator
                cit = cells.begin(); 
              cit != cells.end();
              ++cit)
            if( ! (*cit)->is_in_domain() )
            {
              std::cerr << "error: aroung vertex "<< &*vit << ", cell " << &**cit
                        << " should have is_in_domain()==true... (";
              int count = 0;
              for(int i = 0; i < 4; ++i)
                if((*cit)->vertex(i)->point().surface_index()>0)
                  ++count;
              std::cerr << count << " vertices on surface)\n";
              result = false;
            }
        }
        else
          ++point_on_surface_count;
      std::cerr << point_on_surface_count << " vertices on surface\n"
                << tr.number_of_vertices() << " vertices.\n";
      return result;
    } // end check_c2t3()

  } // end namespace Mesh_3

template <
  typename C2T3,
  typename Surface,
  typename Facets_criteria,
  typename Tets_criteria,
  typename VolumeMeshTraits = 
    typename CGAL::Surface_mesh_traits_generator_3<Surface>::type
  >
class Implicit_surfaces_mesher_3
{
public:
  // ** C2T3 **
  typedef C2T3 C2t3;
  typedef typename C2t3::Triangulation Tr;

  // ** two mesher levels **/

#ifdef CGAL_SURFACE_MESHER_MANIFOLD
  typedef typename Surface_mesher::Surface_mesher_regular_edges_base<
    C2t3,
    Surface,
    VolumeMeshTraits,
    Facets_criteria> Facets_level_re_base;

  typedef typename Surface_mesher::Surface_mesher_manifold_base<
    C2t3,
    Surface,
    VolumeMeshTraits,
    Facets_criteria,
    Facets_level_re_base> Facets_level_base;
#else
  typedef typename Surface_mesher::Surface_mesher_base<
    C2t3,
    Surface,
    VolumeMeshTraits,
    Facets_criteria> Facets_level_base;
#endif

#ifdef NDEBUG
  static const Surface_mesher::Debug_flag debug_flag = Surface_mesher::NO_DEBUG;
#else
  static const Surface_mesher::Debug_flag debug_flag = Surface_mesher::DEBUG;
#endif

#ifdef CGAL_SURFACE_MESHER_VERBOSE
  static const Surface_mesher::Verbose_flag verbose_flag = 
    Surface_mesher::VERBOSE;
#else
  static const Surface_mesher::Verbose_flag verbose_flag = 
    Surface_mesher::NOT_VERBOSE;
#endif

  typedef typename 
  Surface_mesher::Surface_mesher<Facets_level_base,
    verbose_flag,
    debug_flag> Facets_level;

  typedef typename Mesh_3::Refine_tets<Tr,
                                       Tets_criteria,
                                       Surface,
                                       VolumeMeshTraits,
                                       Mesh_3::Refine_tets_with_oracle_base<
                                         Tr,
                                         Tets_criteria,
                                         Surface,
                                         VolumeMeshTraits>,
                                       Facets_level> Tets_level;

  // ** visitors **
  typedef typename Mesh_3::tets::Refine_facets_visitor<Tr,
     Tets_level, Null_mesh_visitor> Tets_facets_visitor;
  typedef typename Mesh_3::Visitor_for_surface <
    Tr,
    Null_mesh_visitor
    > Surface_facets_visitor;

  typedef Combine_mesh_visitor<Surface_facets_visitor,
    Tets_facets_visitor> Facets_visitor;

  typedef Surface_mesher::Visitor<Tr, Facets_level, 
    Facets_visitor> Tets_visitor;


private:
  // ** private data members **
  Null_mesher_level null_mesher_level;
  Null_mesh_visitor null_visitor;

  C2t3& c2t3;
  VolumeMeshTraits& oracle;
  Surface& surface;
  Facets_level facets;
  Tets_level tets;

  Surface_facets_visitor surface_facets_visitor;
  Tets_facets_visitor tets_facets_visitor;
  Facets_visitor facets_visitor;
  Tets_visitor tets_visitor;

  bool initialized;

  // ** types used in code **
  typedef typename Tr::Point Point;

public:
  Implicit_surfaces_mesher_3(C2t3& c2t3, Surface& surface,
                             Facets_criteria& facets_criteria,
                             Tets_criteria tets_crit,
                             VolumeMeshTraits volume_mesh_traits = 
                               VolumeMeshTraits())
    : c2t3(c2t3), oracle(volume_mesh_traits), surface(surface),
      facets(c2t3, surface, volume_mesh_traits, facets_criteria),
      tets(c2t3.triangulation(), tets_crit, surface, volume_mesh_traits, facets),
      surface_facets_visitor(&null_visitor),
      tets_facets_visitor(&tets, &null_visitor),
      facets_visitor(Facets_visitor(&surface_facets_visitor,
				    &tets_facets_visitor)),
      tets_visitor(&facets, &facets_visitor),
      initialized(false)
  {}

  void init()
  {
    Tr& tr = tets.triangulation_ref_impl();

    typedef typename Tr::Geom_traits Geom_traits;
    typename Geom_traits::Construct_circumcenter_3 circumcenter = 
      tr.geom_traits().construct_circumcenter_3_object();

    for(typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
	cit != tr.finite_cells_end();
	++cit)
      {
	const Point& p = cit->vertex(0)->point();
	const Point& q = cit->vertex(1)->point();
	const Point& r = cit->vertex(2)->point();
	const Point& s = cit->vertex(3)->point();

	cit->set_in_domain(oracle.is_in_volume(surface, 
                                               circumcenter(p,q,r,s)));
      }

    facets.scan_triangulation();
    tets.scan_triangulation();
    initialized = true;

    CGAL_assertion_code(
                        std::cerr << "checking c2t3 after init...\n";
                        CGAL::Mesh_3::check_c2t3(c2t3);
                        std::cerr << "done.\n";
                        )


  }

  C2t3& complex_2_in_triangulation_3()
  {
    return c2t3;
  }

  void refine_mesh()
  {
    if(!initialized)
      init();
    tets.refine(tets_visitor);
  }

  void refine_surface()
  {
    if(!initialized)
      init();
    std::cerr << "Starting refine_surface()\n";
    
    while( ! facets.is_algorithm_done() )
      facets.one_step(tets_visitor.previous_level());
  }

  void step_by_step()
  {
    if(!initialized)
      init();
    tets.try_to_insert_one_point(tets_visitor);
  }

  bool done()
  {
    return tets.no_longer_element_to_refine();
    
  }
}; // end Implicit_surfaces_mesher_3

} // end namespace CGAL

#endif // CGAL_IMPLICIT_SURFACES_MESHER_3_H
