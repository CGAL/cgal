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

#ifndef CGAL_DELAUNAY_GRAPH_CONCEPT_H
#define CGAL_DELAUNAY_GRAPH_CONCEPT_H 1

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <iostream>
#include <fstream>
#include <CGAL/Voronoi_diagram_2/Dummy_iterator.h>


namespace CGAL {


template<class K>
class Delaunay_graph_concept
{
  typedef Delaunay_graph_concept<K>  Self;

 public:
  typedef unsigned int   size_type;
  typedef K              Geom_traits;

  typedef typename Geom_traits::Point_2  Site_2;
  typedef typename Geom_traits::Point_2  Point_2;

  struct Vertex_handle
    : public CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Vertex_handle>
  {
    typedef CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Vertex_handle>  Base;

    Vertex_handle() {}
    Vertex_handle(const Base&) {}

    Site_2 site() const { return Site_2(); }
  };

  struct Face_handle
    : public CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Face_handle>
  {
    typedef CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Face_handle>  Base;

    Face_handle() {}
    Face_handle(const Base&) {}

    Vertex_handle vertex(int) const { return Vertex_handle(); }

    Face_handle neighbor(int) const { return Face_handle(); }
    //      int mirror_index(int) const { return 0; }
  };

  typedef Vertex_handle Vertex;
  typedef Face_handle   Face;

  template<class V, class F>
  struct Data_structure_t {
    typedef F Face_handle;
    typedef V Vertex_handle;

    int mirror_index(Face_handle, int) const { return 0; }
    Vertex_handle mirror_vertex(Face_handle, int) const {
      return Vertex_handle();
    }

    bool is_edge(Face_handle, int) const { return false; }

    std::size_t number_of_edges() const { return 0; }
  };

  typedef Data_structure_t<Vertex_handle,Face_handle> Data_structure;
  typedef Data_structure        Triangulation_data_structure;

  //  typedef typename Data_structure::Face_handle     Face_handle;
  //  typedef typename Data_structure::Vertex_handle   Vertex_handle;

  static const Vertex_handle infinite_vertex() {
    static Vertex_handle inf_v;
    return inf_v;
  }

  static const Vertex_handle finite_vertex() {
    static Vertex_handle fin_v;
    return fin_v;
  }

  static const Face_handle infinite_face() {
    static Face_handle inf_f;
    return inf_f;
  }

  static const Data_structure& data_structure() {
    static Data_structure ds;
    return ds;
  }

  static const Data_structure& tds() {
    return data_structure();
  }

  static const Geom_traits& geom_traits() {
    static Geom_traits gt;
    return gt;
  }

  typedef std::pair<Face_handle,int>  Edge;

 private:
  class Dummy_edge_iterator : public Emptyset_iterator
  {
  public:
    typedef std::size_t   size_type;
    typedef Edge          value_type;
    typedef const Edge*   pointer;
    typedef Edge&         reference;
    typedef const Edge*   const_pointer;
    typedef const Edge&   const_reference;

    Dummy_edge_iterator() {}
    template<class T>
    Dummy_edge_iterator(const T&) {}

    template< class T >
    Dummy_edge_iterator& operator=(const T&) { return *this; }

    Dummy_edge_iterator& operator++()        { return *this; }
    Dummy_edge_iterator& operator++(int)     { return *this; }

    Dummy_edge_iterator& operator--()        { return *this; }
    Dummy_edge_iterator& operator--(int)     { return *this; }

    reference operator*()   { return e; }
    pointer   operator->()  { return &e; }

    const_reference operator*()  const  { return e; }
    const_pointer   operator->() const  { return (&e); }

    bool operator==(const Dummy_edge_iterator&) const {
      return true;
    }

    bool operator!=(const Dummy_edge_iterator&) const {
      return false;
    }

    bool operator<(const Dummy_edge_iterator& other) const {
      return this < &other;
    }

  protected:
    Edge e;
  };

  static const Dummy_edge_iterator& dummy_edge_it() {
    static Dummy_edge_iterator dummy_edge_it_static;
    return dummy_edge_it_static;
  }

 public:
  typedef Dummy_edge_iterator  All_edges_iterator;
  struct Finite_edges_iterator : public Dummy_edge_iterator {
    Finite_edges_iterator() {}
    Finite_edges_iterator(const Dummy_edge_iterator&) {}
  };

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Face>
  All_faces_iterator;

  typedef CGAL_VORONOI_DIAGRAM_2_INS::Dummy_iterator<Vertex>
  All_vertices_iterator;

  typedef All_faces_iterator      Finite_faces_iterator;
  typedef All_faces_iterator      Face_circulator;

  typedef All_vertices_iterator   Finite_vertices_iterator;
  typedef All_vertices_iterator   Vertex_circulator;

  struct Edge_circulator : public Dummy_edge_iterator
  {
    Edge_circulator() {}

    template< class T >
    Edge_circulator& operator=(const T&) { return *this; }

    Edge_circulator& operator++()        { return *this; }
    Edge_circulator& operator++(int)     { return *this; }

    Edge& operator*()         { return this->e; }
    Edge* operator->()        { return &this->e; }

    const Edge& operator*() const   { return this->e; }
    const Edge* operator->() const  { return &this->e; }

    bool operator==(const Edge_circulator&) const {
      return true;
    }

    bool operator!=(const Edge_circulator&) const {
      return false;
    }
  };

  

 private:
  static const Edge_circulator& dummy_edge_circulator() {
    static Edge_circulator dummy_edge_circulator_static;
    return dummy_edge_circulator_static;
  }

 public:
  Delaunay_graph_concept(const Geom_traits& /* gt */ = Geom_traits()) {}

  template<class Iterator>
  Delaunay_graph_concept(Iterator /* first */, Iterator /* beyond */,
			 const Geom_traits& /* gt */ = Geom_traits()) {}


  void insert(const Site_2&) {}

  template<class Iterator>
  int insert(Iterator, Iterator) { return 0; }

  size_type number_of_vertices() const { return 0; }
  size_type number_of_faces() const { return 0; }
  size_type number_of_edges() const { return 0; }

  size_type degree(Vertex_handle ) const { return 0; }

  int dimension() const { return -1; }

  template<class T> bool is_infinite(const T&) const { return false; }
  template<class T> bool is_infinite(const T&,int) const { return false; }

  All_edges_iterator all_edges_begin() const {
    return dummy_edge_it();
  }

  All_edges_iterator all_edges_end() const {
    return dummy_edge_it();
  }

  Finite_edges_iterator finite_edges_begin() const {
    return dummy_edge_it();
  }

  Finite_edges_iterator finite_edges_end() const {
    return dummy_edge_it();
  }

  All_vertices_iterator all_vertices_begin() const {
    return All_vertices_iterator::dummy_reference();
  }

  All_vertices_iterator all_vertices_end() const {
    return All_vertices_iterator::dummy_reference();
  }

  Finite_vertices_iterator finite_vertices_begin() const {
    return Finite_vertices_iterator::dummy_reference();
  }

  Finite_vertices_iterator finite_vertices_end() const {
    return Finite_vertices_iterator::dummy_reference();
  }

  All_faces_iterator all_faces_begin() const {
    return All_faces_iterator::dummy_reference();
  }

  All_faces_iterator all_faces_end() const {
    return All_faces_iterator::dummy_reference();
  }

  Finite_faces_iterator finite_faces_begin() const {
    return Finite_faces_iterator::dummy_reference();
  }

  Finite_faces_iterator finite_faces_end() const {
    return Finite_faces_iterator::dummy_reference();
  }

  Edge_circulator incident_edges(const Vertex_handle&,
				 const Face_handle& = Face_handle()) const
  {
    return dummy_edge_circulator();
  }

  Vertex_circulator incident_vertices(const Vertex_handle&,
				      const Face_handle& = Face_handle()) const
  {
    return Vertex_circulator::dummy_reference();
  }

  Face_circulator incident_faces(const Vertex_handle&,
				 const Face_handle& = Face_handle()) const
  {
    return Face_circulator::dummy_reference();
  }


  Point_2 dual(const Face_handle&) const { return Point_2(); }

  bool is_valid(bool = true,int = 0) const { return true; }

  void swap(Self& /* other */) {}
  void clear() {}
};


template<class K>
std::istream&
operator>>(std::istream& is, Delaunay_graph_concept<K>&)
{
  return is;
}

template<class K>
std::ostream&
operator<<(std::ostream& os, const Delaunay_graph_concept<K>&)
{
  return os;
}


template<class K>
std::istream&
operator>>(std::istream& is,
	   const typename Delaunay_graph_concept<K>::Site_2&)
{
  double d;
  is >> d;
  return is;
}

template<class K>
std::ostream&
operator<<(std::ostream& os,
	   const typename Delaunay_graph_concept<K>::Site_2&)
{
  return os;
}

} //namespace CGAL


#endif // CGAL_DELAUNAY_GRAPH_CONCEPT_H
