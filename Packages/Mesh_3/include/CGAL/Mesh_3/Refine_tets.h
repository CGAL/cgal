// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_3_REFINE_TETS_H
#define CGAL_MESH_3_REFINE_TETS_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_3/Triangulation_mesher_level_traits_3.h>

#include <CGAL/Mesh_3/Refine_edges.h>
#include <CGAL/Mesh_3/Refine_facets.h>
#include <CGAL/Mesh_3/Double_map_container.h>

#include <list>

namespace CGAL {
namespace Mesh_3 {

template <class Tr,
          class Criteria,
          class Container = Double_map_container<typename Tr::Cell_handle,
                                                 typename Criteria::Quality>
>
class Refine_tets_base :
    public Container,
    public Triangulation_mesher_level_traits_3<Tr>,
    public No_test_point_conflict,
    public No_before_insertion,
    public No_after_insertion
{
protected:
  typedef typename Tr::Point Point;
  typedef typename Tr::Point Point;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef Triangulation_mesher_level_traits_3<Tr> Triangulation_traits;
  typedef typename Triangulation_traits::Zone Zone;

  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Facet_circulator Facet_circulator;

  typedef typename Tr::Facet Facet;

public:
  typedef typename Criteria::Quality Quality;

  using Triangulation_mesher_level_traits_3<Tr>::triangulation_ref_impl;

public:
  /** \name CONSTRUCTORS */

  Refine_tets_base(Tr& t, Criteria crit) 
    : Triangulation_mesher_level_traits_3<Tr>(t), criteria(crit) {}

protected:
  /* --- protected datas --- */
  //  Tr& tr; /**< The triangulation itself. */
  Criteria criteria; /**< Meshing criteria for tetrahedra. */

protected:
  /* --- protected functions --- */

  bool should_be_refined(const Cell_handle c, Quality& qual) const
  {
    return criteria.is_bad_object()(c,qual);
  }

  bool should_be_refined(const Cell_handle c) const
  {
    Quality q;
    return should_be_refined(c, q);
  }

  bool test_if_cell_is_bad(const Cell_handle c)
  {
    Quality q;
    if( c->is_in_domain() && should_be_refined(c, q) )
      {
	this->add_element(c, q);
	return true;
      }
    return false;
  }

public:
  /** \name Functions that this level must declare. */

  void scan_triangulation_impl()
  {
    for(typename Tr::Finite_cells_iterator cit = 
	  triangulation_ref_impl().finite_cells_begin();
        cit != triangulation_ref_impl().finite_cells_end();
        ++cit)
      test_if_cell_is_bad(cit);
  }

  Point refinement_point_impl(const Cell_handle& c) const
  {
    typename Geom_traits::Construct_circumcenter_3 circumcenter = 
      triangulation_ref_impl().geom_traits().construct_circumcenter_3_object();

    const Point& p = c->vertex(0)->point();
    const Point& q = c->vertex(1)->point();
    const Point& r = c->vertex(2)->point();
    const Point& s = c->vertex(3)->point();

    return circumcenter(p, q, r, s);
  }

  void before_conflicts_impl(const Cell_handle&, const Point& 
#if CGAL_MESH_3_DEBUG_BEFORE_CONFLICTS
			   p)
  {
    std::cerr << "Refine_tets: before conflicts of " << p;
#else
    ) {
#endif
  }

  void after_no_insertion_impl(const Cell_handle&, const Point&,
			       const Zone& )
  {
#if CGAL_MESH_3_DEBUG_AFTER_NO_INSERTION
    std::cerr << "  REJECTED!" << std::endl;
#endif
  }
}; // end Refine_tets_base  

template <class Tr,
          class Criteria,
          class Oracle,
          class Container = Double_map_container<typename Tr::Cell_handle,
                                                 typename Criteria::Quality>
>
class Refine_tets_with_oracle_base 
  : public Refine_tets_base<Tr,
                            Criteria,
                            Container>
{
public:
  typedef Refine_tets_base<Tr, Criteria, Container> Base;
  typedef Refine_tets_with_oracle_base<Tr, 
                                       Criteria,
                                       Oracle,
                                       Container> Self;
  
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::Point Point;
  typedef typename Base::Zone Zone;
  
  

  /** \name CONSTRUCTORS */

  Refine_tets_with_oracle_base(Tr& t, Criteria crit, Oracle& o)
    : Base(t, crit), oracle(o) {}

public:
  /* \name Overriden functions of this level */
  void do_before_insertion(const Cell_handle& c, const Point&,
                           Zone& zone)
  {
    remove_star_from_bad_cells(c, zone); // FIXME: name
  }

  void remove_star_from_bad_cells(const Cell_handle& c, Zone& zone)
  {
    for(typename Zone::Cells_iterator cit = zone.cells.begin();
	cit != zone.cells.end();
	++cit)
      if(*cit != c)
        this->remove_element(*cit);
  }

  void do_after_insertion(const Vertex_handle& v)
  {
    std::cerr << ".";
#if CGAL_MESH_3_DEBUG_AFTER_INSERTION
    std::cerr << "  INSERTED." << std::endl;
#endif
    scan_star(v);
  }

  void scan_star(const Vertex_handle& v)
  {
    // scan tets
    typedef std::list<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;

    triangulation_ref_impl().
      incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator cit = incident_cells.begin();
        cit != incident_cells.end();
        ++cit)
      {
	const Point& p = (*cit)->vertex(0)->point();
	const Point& q = (*cit)->vertex(1)->point();
	const Point& r = (*cit)->vertex(2)->point();
	const Point& s = (*cit)->vertex(3)->point();

	(*cit)->set_in_domain(oracle.surf_equation(circumcenter(p,q,r,s))<0.);

	test_if_cell_is_bad(*cit);
      }
  }

private:
  /* --- private datas --- */
  Oracle& oracle;

}; // end Refine_tets_with_oracle_base

  namespace tets {

    template <typename Tr, typename Refine_tets>
    class Refine_facets_visitor {
      Refine_tets* refine_tets;
      Null_mesh_visitor null_visitor;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Facet Facet;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      Refine_facets_visitor(Refine_tets* refine_tets_)
        : refine_tets(refine_tets_) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      void before_insertion(const Facet&,
                            const Point&,
                            Zone& zone) 
      {
        refine_tets->remove_star_from_bad_cells(Cell_handle(), zone); // FIXME: beurk
      }

      void after_insertion(const Vertex_handle& v)
      {
	refine_tets->scan_star(v);
      }

      template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Null_mesh_visitor& previous_level()
      {
        return null_visitor;
      }
      
    }; // end Refine_facets_visitor

    template <class Tr, typename Refine_tets>
    class Refine_tets_visitor {
    public:
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_iterator Cell_iterator;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Refine_facets_visitor<Tr, Refine_tets> Previous_visitor;


    private:
      Previous_visitor previous;

    public:
      Refine_tets_visitor(Refine_tets* refine_tets_)
        : previous(refine_tets_) {}

      template <typename E, typename P>
      void before_conflicts(E, P) const {}

      template <typename E, typename P, typename Z>
      void before_insertion(E, P, Z) const {}

      template <typename V>
      void after_insertion(V) const {}

     template <typename E, typename P, typename Z>
      void after_no_insertion(E, P, Z) const {}

      Previous_visitor& previous_level()
      {
        return previous;
      }
    }; // end Refine_tets_visitor

  }; // end namespace Mesh_3::tets
  
template <typename Tr,
          typename Criteria,
          typename Oracle,
          typename Base = 
            Refine_tets_with_oracle_base<Tr, Criteria, Oracle>,
          typename Facets_level = Refine_facets<Tr>
 >
class Refine_tets : 
  public Base, 
  public Mesher_level <
    Tr,
    Refine_tets<Tr, Criteria, Oracle, Base, Facets_level>,
    typename Tr::Cell_handle,
    Facets_level,
    Triangulation_mesher_level_traits_3<Tr>
  >
{
  Facets_level& f_level;
public:
  typedef Refine_tets<Tr, Criteria, Oracle, Base, Facets_level> Self;
  typedef Mesher_level <
    Tr,
    Refine_tets<Tr, Criteria, Oracle, Base, Facets_level>,
    typename Tr::Cell_handle,
    Facets_level,
    Triangulation_mesher_level_traits_3<Tr>
  > Mesher;
  Refine_tets(Tr& t, Criteria crit, Oracle& oracle, Facets_level& facets_level)
    : Base(t, crit, oracle), Mesher(facets_level), f_level(facets_level)
  {}


void do_before_insertion(const typename Base::Cell_handle& c,
			 const typename Base::Point& p,
			 typename Base::Zone& zone)
{
  f_level.do_before_insertion(typename Tr::Facet (),
			      p,
			      zone);
  Base::do_before_insertion(c, p, zone);
}

void do_after_insertion(const typename Base::Vertex_handle& v)
{
  f_level.restore_restricted_Delaunay(v);
  Base::do_after_insertion(v);
}

}; // end class Refine_tets

}; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_TETS_H
