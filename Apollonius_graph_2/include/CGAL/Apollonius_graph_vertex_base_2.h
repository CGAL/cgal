// Copyright (c) 2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
#define CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#include <list>

#include <CGAL/Triangulation_ds_vertex_base_2.h>
#include <CGAL/triangulation_assertions.h>

namespace CGAL {

template <class AGVB2_Iterator>
struct Apollonius_graph_vertex_base_nested_iterator_traits
{
  typedef AGVB2_Iterator                              Base_iterator;
  typedef typename
  Base_iterator::value_type::Hidden_sites_iterator    Iterator;

  Iterator begin(Base_iterator it) const
  {
    return it->hidden_sites_begin();
  }

  Iterator end(Base_iterator it) const
  {
    return it->hidden_sites_end();
  }

};



template <class Gt,
          bool StoreHidden = true,
          class Vb = Triangulation_ds_vertex_base_2<> >
class Apollonius_graph_vertex_base_2
  : public Vb
{
private:
  typedef typename Vb::Triangulation_data_structure   AGDS;
public:
  // TYPES
  //------
  typedef Gt                             Geom_traits;
  typedef Vb                             Base;
  typedef typename Gt::Point_2           Point;
  typedef typename Gt::Site_2            Site_2;
  typedef AGDS                                 Apollonius_graph_data_structure_2;
  typedef typename AGDS::Face_handle     Face_handle;
  typedef typename AGDS::Vertex_handle   Vertex_handle;

  enum {Store_hidden = StoreHidden};

  template < typename AGDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<AGDS2>::Other      Vb2;
    typedef Apollonius_graph_vertex_base_2<Gt,StoreHidden,Vb2>  Other;
  };


private:
  // local types
  typedef std::list<Site_2>         Container;

public:
  // TYPES (continued)
  //------------------
  //  typedef Container                        Hidden_sites_container;
  typedef typename Container::iterator     Hidden_sites_iterator;

public:
  // CREATION
  //---------
  Apollonius_graph_vertex_base_2() : Vb() {}
  Apollonius_graph_vertex_base_2(const Site_2& p) : Vb(), _p(p) {}
  Apollonius_graph_vertex_base_2(const Site_2& p, Face_handle f)
    : Vb(f), _p(p) {}

  ~Apollonius_graph_vertex_base_2()
  {
    clear_hidden_sites_container();
  }


  // ACCESS METHODS
  //---------------
  const Site_2& site() const { return _p; }
  Site_2& site() { return _p; }

  Face_handle face() const { return Vb::face(); }

  std::size_t number_of_hidden_sites() const {
    return hidden_site_list.size();
  }

  Hidden_sites_iterator hidden_sites_begin() {
    return hidden_site_list.begin();
  }

  Hidden_sites_iterator hidden_sites_end() {
    return hidden_site_list.end();
  }

public:
  // SETTING AND UNSETTING
  //----------------------
  void set_site(const Site_2& p) { _p = p; }


  void add_hidden_site(const Site_2& p)
  {
    if ( StoreHidden ) {
      hidden_site_list.push_back(p);
    }
  }

  void clear_hidden_sites_container()
  {
    hidden_site_list.clear();
  }

public:
  // VALIDITY CHECK
  bool is_valid(bool verbose = false, int level = 0) const {
    return Vb::is_valid(verbose, level);
  }

private:
  // class variables
  Container hidden_site_list;
  Site_2 _p;
};

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_VERTEX_BASE_2_H
