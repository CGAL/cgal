// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H
#define CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Voronoi_traits_functors.h>
#include <list>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//===========================================================================
//===========================================================================

template<class Site_t, class DG>
struct Default_site_inserter
{
  typedef DG                                        Delaunay_graph;
  typedef Site_t                                    Site_2;
  typedef typename Delaunay_graph::Vertex_handle    result_type;
  typedef Arity_tag<2>                              Arity;

  Default_site_inserter() {}

  //  template<typename T> Default_site_inserter(T t) {}

  result_type operator()(Delaunay_graph& dg, const Site_2& t) const {
    return dg.insert(t);
  }

  // helpful for the case of sites that create multiple vertices
  template<class OutputIterator>
  OutputIterator operator()(Delaunay_graph& dg, const Site_2& t,
			    OutputIterator oit) const {
    *oit++ = operator()(dg, t);
    return oit;
  }
};

//===========================================================================

template<class VT, class SI>
class Default_caching_site_inserter
{
private:
  typedef VT  Voronoi_traits;
  typedef SI  Site_inserter;

public:
  typedef typename Voronoi_traits::Delaunay_graph   Delaunay_graph;
  typedef typename Site_inserter::Site_2            Site_2;
  typedef typename Delaunay_graph::Vertex_handle    result_type;
  typedef Arity_tag<2>                              Arity;

public:
  Default_caching_site_inserter(const Voronoi_traits* vt = NULL) : vt_(vt) {}

  result_type operator()(Delaunay_graph& dg, const Site_2& t) const
  {
    // THERE IS POTENTIAL PROBLEM IN THIS METHOD; WE DO NOT ACCOUNT
    // FOR THE VERTICES THAT BECOME HIDDEN ONCE A NEW SITE IN
    // INSERTED; THIS AFFECTS THE CORRECTNESS OF THE THE CACHED FACE
    // DEGENERACY TESTER, SINCE IN THIS CASE, I.E., WHEN HIDDEN SITES
    // ARE CREATED, WE NEED TO DELETE THEM FROM THE CACHED DEGENERACY
    // TESTER.
    typedef typename Delaunay_graph::Edge          Dual_edge;
    typedef typename Delaunay_graph::Face_handle   Dual_face_handle;
    typedef std::list<Dual_edge>                   Dual_edge_list;
    typedef std::list<Dual_face_handle>            Dual_face_handle_list;

    if ( dg.dimension() != 2 ) { return dg.insert(t); }

    Dual_edge_list        e_list;
    Dual_face_handle_list f_list;
    dg.get_conflicts_and_boundary(t, std::back_inserter(f_list),
				  std::back_inserter(e_list));

    for (typename Dual_edge_list::iterator it = e_list.begin();
	 it != e_list.end(); ++it) {
      vt_->edge_degeneracy_tester_object().erase(*it);
    }

    for (typename Dual_face_handle_list::iterator it = f_list.begin();
	 it != f_list.end(); ++it) {
      Dual_face_handle f = *it;
      for (int j = 0; j < 3; j++) {
	Dual_edge e(f, j);
	vt_->edge_degeneracy_tester_object().erase(e);
      }
    }
    return Site_inserter()(dg, t);
  }

  template<class OutputIterator>
  OutputIterator operator()(Delaunay_graph& dg, const Site_2& t,
			    OutputIterator oit) const {
    *oit++ = operator()(dg, t);
    return oit;
  }

private:
  const Voronoi_traits* vt_;
};

//===========================================================================

template<class VT>
struct Default_caching_site_inserter<VT,Null_functor>
{
  Default_caching_site_inserter() {}
  template<typename T> Default_caching_site_inserter(T t) {}
};

//===========================================================================
//===========================================================================

template<class VT, class Site_inserter = Default_site_inserter<VT,int> >
class Default_aggregate_site_inserter
{
private:
  typedef VT  Voronoi_traits;

public:
  typedef typename Voronoi_traits::Delaunay_graph   Delaunay_graph;
  typedef typename Voronoi_traits::Site_2           Site_2;
  typedef int                                       result_type;
  typedef Arity_tag<3>                              Arity;

public:
  Default_aggregate_site_inserter() : site_inserter() {}

  Default_aggregate_site_inserter(const Voronoi_traits* vt)
    : site_inserter(vt) {}

  template<class Iterator>
  int operator()(Delaunay_graph& dg,
		 Iterator first, Iterator beyond) const {
    int n = dg.number_of_vertices();
    for (Iterator it = first; it != beyond; ++it) {
      site_inserter(dg, *it);
    }
    return dg.number_of_vertices() - n;
  }

  template<class Iterator, class OutputIterator>
  OutputIterator operator()(Delaunay_graph& dg,
			    Iterator first, Iterator beyond,
			    OutputIterator oit) const {
    for (Iterator it = first; it != beyond; ++it) {
      oit = site_inserter(dg, *it, oit);
    }
    return oit;
  }  

private:
  Site_inserter site_inserter;
};

//===========================================================================
//===========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H
