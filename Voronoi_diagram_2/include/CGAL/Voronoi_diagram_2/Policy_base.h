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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_POLICY_BASE_H
#define CGAL_VORONOI_DIAGRAM_2_POLICY_BASE_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Identity_rejectors.h>
#include <CGAL/Voronoi_diagram_2/Cached_degeneracy_testers.h>
#include <CGAL/Voronoi_diagram_2/Default_site_inserters.h>
#include <CGAL/Voronoi_diagram_2/Default_site_removers.h>

namespace CGAL {

namespace VoronoiDiagram_2 { namespace Internal {

//=========================================================================
//=========================================================================

template<class DG, class ER, class FR>
class Policy_base_base
{
private:
  typedef Policy_base_base<DG,ER,FR>  Self;

public:
  typedef DG                                          Delaunay_graph;
  typedef typename Delaunay_graph::Vertex_handle      Delaunay_vertex_handle;
  typedef typename Delaunay_graph::Face_handle        Delaunay_face_handle;
  typedef typename Delaunay_graph::Edge               Delaunay_edge;
  typedef typename Delaunay_graph::Edge_circulator    Delaunay_edge_circulator;

  typedef typename Delaunay_graph::All_edges_iterator
  All_Delaunay_edges_iterator;
  typedef typename Delaunay_graph::Finite_edges_iterator
  Finite_Delaunay_edges_iterator;


  typedef ER   Edge_rejector;
  typedef FR   Face_rejector;

  const Edge_rejector& edge_rejector_object() const {
    return e_rejector_;
  }

  const Face_rejector& face_rejector_object() const {
    return f_rejector_;
  }

  void clear() {
    e_rejector_.clear();
    f_rejector_.clear();
  }

  void swap(Self& other) {
    e_rejector_.swap(other.e_rejector_);
    f_rejector_.swap(other.f_rejector_);
  }

  bool is_valid() const {
    return e_rejector_.is_valid() && f_rejector_.is_valid();
  }

  bool is_valid(const Delaunay_graph& ) const {
    return e_rejector_.is_valid() && f_rejector_.is_valid();
  }

protected:
  Edge_rejector e_rejector_;
  Face_rejector f_rejector_;
};

//=========================================================================
//=========================================================================

template<class DG, class ER, class FR, class SI, class SR>
class Policy_base
  : public Policy_base_base<DG,ER,FR>
{
private:
  typedef Policy_base<DG,ER,FR,SI,SR>   Self;
  typedef Policy_base_base<DG,ER,FR>    Base;

public:
  typedef SI   Site_inserter;
  typedef SR   Site_remover;

  typedef typename Functor_exists<Site_inserter>::Value  Has_site_inserter;
  typedef typename Functor_exists<Site_remover>::Value   Has_site_remover;

  Site_inserter site_inserter_object() const {
    return Site_inserter();
  }

  Site_remover site_remover_object() const {
    return Site_remover();
  }
};

//=========================================================================
//=========================================================================

template<class DG, class ERB, class FRB, class SIB, class SRB>
class Caching_policy_base
  : public Policy_base_base<DG,
			    Cached_edge_rejector<ERB>,
			    Cached_face_rejector<FRB> >
{
protected:
  typedef ERB   Edge_rejector_base;
  typedef FRB   Face_rejector_base;
  typedef SIB   Site_inserter_base;
  typedef SRB   Site_remover_base;

  typedef Caching_policy_base<DG,ERB,FRB,SIB,SRB>  Self;

public:
  typedef Cached_edge_rejector<Edge_rejector_base>                Edge_rejector;
  typedef Cached_face_rejector<Face_rejector_base>                Face_rejector;
  typedef Default_caching_site_inserter<Self,Site_inserter_base>  Site_inserter;
  typedef Default_caching_site_remover<Self,Site_remover_base>    Site_remover;

  typedef typename Functor_exists<Site_inserter>::Value   Has_site_inserter;
  typedef typename Functor_exists<Site_remover>::Value    Has_site_remover;

  typedef DG  Delaunay_graph;

  Site_inserter site_inserter_object() const {
    return Site_inserter(this);
  }

  Site_remover site_remover_object() const {
    return Site_remover(this);
  }

  bool is_valid() const {
    return this->e_rejector_.is_valid() && this->f_rejector_.is_valid();
  }

  bool is_valid(const Delaunay_graph& dg) const {
    return this->e_rejector_.is_valid(dg) && this->f_rejector_.is_valid(dg);
  }
};


//=========================================================================
//=========================================================================


} } //namespace VoronoiDiagram_2::Internal

} //namespace CGAL


#endif // CGAL_VORONOI_DIAGRAM_2_POLICY_BASE_2_H
