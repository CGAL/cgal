// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H
#define CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H 1

#include <CGAL/license/Voronoi_diagram_2.h>


#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Adaptation_traits_functors.h>
#include <list>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//===========================================================================
//===========================================================================

template<class Site_t, class DG>
struct Default_site_inserter
{
  typedef DG                                        Delaunay_graph;
  typedef Site_t                                    Site_2;
  typedef typename Delaunay_graph::Vertex_handle    result_type;

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

template<class AT, class SI>
class Default_caching_site_inserter
{
private:
  typedef AT  Adaptation_traits;
  typedef SI  Site_inserter;

public:
  typedef typename Adaptation_traits::Delaunay_graph  Delaunay_graph;
  typedef typename Site_inserter::Site_2              Site_2;
  typedef typename Delaunay_graph::Vertex_handle      result_type;

public:
  Default_caching_site_inserter(const Adaptation_traits* at = NULL) : at_(at) {}

  result_type operator()(Delaunay_graph& dg, const Site_2& t) const
  {
    // THERE IS POTENTIAL PROBLEM IN THIS METHOD; WE DO NOT ACCOUNT
    // FOR THE VERTICES THAT BECOME HIDDEN ONCE A NEW SITE IN
    // INSERTED; THIS AFFECTS THE CORRECTNESS OF THE THE CACHED FACE
    // DEGENERACY TESTER, SINCE IN THIS CASE, I.E., WHEN HIDDEN SITES
    // ARE CREATED, WE NEED TO DELETE THEM FROM THE CACHED DEGENERACY
    // TESTER.
    typedef typename Delaunay_graph::Edge          Edge;
    typedef typename Delaunay_graph::Face_handle   Face_handle;
    typedef std::list<Edge>                        Edge_list;
    typedef std::list<Face_handle>                 Face_handle_list;

    if ( dg.dimension() != 2 ) { return dg.insert(t); }

    Edge_list        e_list;
    Face_handle_list f_list;
    dg.get_conflicts_and_boundary(t, std::back_inserter(f_list),
				  std::back_inserter(e_list));

    for (typename Edge_list::iterator it = e_list.begin();
	 it != e_list.end(); ++it) {
      at_->edge_rejector_object().erase(*it);
    }

    for (typename Face_handle_list::iterator it = f_list.begin();
	 it != f_list.end(); ++it) {
      Face_handle f = *it;
      for (int j = 0; j < 3; j++) {
	Edge e(f, j);
	at_->edge_rejector_object().erase(e);
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
  const Adaptation_traits* at_;
};

//===========================================================================

template<class AT>
class Default_caching_site_inserter<AT,Null_functor>
{
public:
  Default_caching_site_inserter() {}
  template<typename T> Default_caching_site_inserter(T /*t*/) {}
};

//===========================================================================
//===========================================================================

template<class AT, class Site_inserter = Default_site_inserter<AT,int> >
class Default_aggregate_site_inserter
{
private:
  typedef AT  Adaptation_traits;

public:
  typedef typename Adaptation_traits::Delaunay_graph   Delaunay_graph;
  typedef typename Adaptation_traits::Site_2           Site_2;
  typedef int                                          result_type;

public:
  Default_aggregate_site_inserter() : site_inserter() {}

  Default_aggregate_site_inserter(const Adaptation_traits* at)
    : site_inserter(at) {}

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

} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL

#endif // CGAL_VORONOI_DIAGRAM_2_DEFAULT_SITE_INSERTERS_H
