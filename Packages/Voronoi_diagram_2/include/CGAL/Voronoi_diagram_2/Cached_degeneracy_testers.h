// Copyright (c) 2005 University of Crete (Greece).
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

#ifndef CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <cstdlib>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================


template<class Edge_tester_t>
class Cached_edge_degeneracy_tester
{
 public:
  typedef Edge_tester_t    Edge_degeneracy_tester;
  typedef typename Edge_degeneracy_tester::Dual_graph Dual_graph;
  typedef typename Edge_degeneracy_tester::Edge Edge;
  typedef typename Edge_degeneracy_tester::Face_handle Face_handle;
  typedef typename Edge_degeneracy_tester::Edge_circulator Edge_circulator;

  typedef typename Edge_degeneracy_tester::All_edges_iterator
  All_edges_iterator;

  typedef typename Edge_degeneracy_tester::Finite_edges_iterator
  Finite_edges_iterator;

  typedef typename Edge_degeneracy_tester::result_type  result_type;
  typedef typename Edge_degeneracy_tester::Arity        Arity;

 private:
  class Edge_hash_function
    : public Handle_hash_function
  {
  private:
    typedef Handle_hash_function     Base;

  public:
    typedef Base::result_type        result_type;

    template<class Edge>
    result_type operator()(const Edge& e) const
    {
      return (Base::operator()(e.first)) << e.second;
    }
  };

  // true if degenerate, false otherwise
  typedef Unique_hash_map<Edge,bool,Edge_hash_function>  Edge_map;

  Edge opposite(const Dual_graph& dual, const Edge& e) const {
    int i_mirror = dual.tds().mirror_index(e.first, e.second);
    return Edge( e.first->neighbor(e.second), i_mirror );
  }

 public:
  bool operator()(const Dual_graph& dual, const Edge& e) const {
    if ( emap.is_defined(e) ) { return emap[e]; }
    bool b = e_tester(dual, e);
    emap[e] = b;
    emap[opposite(dual, e)] = b;
    return b;
  }

  bool operator()(const Dual_graph& dual, const Face_handle& f, int i) const {
    return operator()(dual, Edge(f,i));
  }

  bool operator()(const Dual_graph& dual, const Edge_circulator& ec) const {
    return operator()(dual, *ec);
  }

  bool operator()(const Dual_graph& dual,
		  const All_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Dual_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

 private:
  Edge_degeneracy_tester e_tester;
  mutable Edge_map emap;
};


//=========================================================================
//=========================================================================


template<class Face_degeneracy_t>
class Cached_face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef Face_degeneracy_t                     Face_degeneracy_tester;

  typedef typename Face_degeneracy_tester::Dual_graph     Dual_graph;
  typedef typename Face_degeneracy_tester::Vertex_handle  Vertex_handle;

  typedef typename Face_degeneracy_tester::Vertex_circulator
  Vertex_circulator;

  typedef typename Face_degeneracy_tester::Finite_vertices_iterator
  Finite_vertices_iterator;

  typedef typename Face_degeneracy_tester::All_vertices_iterator
  All_vertices_iterator;

  typedef typename Face_degeneracy_tester::result_type  result_type;
  typedef typename Face_degeneracy_tester::Arity        Arity;

 private:
  typedef Unique_hash_map<Vertex_handle,bool>   Vertex_map;

 public:
  bool operator()(const Dual_graph& dual, const Vertex_handle& v) const {
    if ( vmap.is_defined(v) ) { return vmap[v]; }
    bool b = f_tester(dual, v);
    vmap[v] = b;
    return b;
  }
 
  bool operator()(const Dual_graph& dual, const Vertex_circulator& vc) const {
    return operator()(dual, Vertex_handle(vc));
  }

  bool operator()(const Dual_graph& dual,
		  const Finite_vertices_iterator& vit) const {
    return operator()(dual, Vertex_handle(vit));
  }

#ifndef CGAL_T2_USE_ITERATOR_AS_HANDLE
  bool operator()(const Dual_graph& dual,
		  const All_vertices_iterator& vit) const {
    return operator()(dual, Vertex_handle(vit));
  }
#endif

 private:
  Face_degeneracy_tester f_tester;
  mutable Vertex_map vmap;
};


//=========================================================================
//=========================================================================


template<class Tester_handle_t, class Types_t>
struct Handle_to_tester_adaptor : public Types_t
{
 private:
  typedef Tester_handle_t                       Tester_handle;
  typedef typename Tester_handle::element_type  Tester;

 public:
  typedef typename Tester::result_type          result_type;

  Handle_to_tester_adaptor() {
    Tester rc_tester;
    h_.initialize_with(rc_tester);    
  }

  template<typename argument_type_1, typename argument_type_2>
  result_type operator()(const argument_type_1& arg1,
			 const argument_type_2& arg2) const {
    return h_.Ptr()->operator()(arg1, arg2);
  }

  template<typename argument_type_1, typename argument_type_2,
	   typename argument_type_3>
  result_type operator()(const argument_type_1& arg1,
			 const argument_type_2& arg2,
			 const argument_type_3& arg3) const {
    return h_.Ptr()->operator()(arg1, arg2, arg3);
  }

 private:
  Tester_handle h_;
};


//=========================================================================
//=========================================================================

template<class Edge_tester_t>
class Ref_counted_edge_degeneracy_tester_base
  : public Cached_edge_degeneracy_tester<Edge_tester_t>,
    public Ref_counted_virtual
{
 private:
  typedef Cached_edge_degeneracy_tester<Edge_tester_t>  Base;

 public:
  ~Ref_counted_edge_degeneracy_tester_base() {}
};

//=========================================================================

template<class Face_tester_t>
class Ref_counted_face_degeneracy_tester_base
  : public Cached_face_degeneracy_tester<Face_tester_t>,
    public Ref_counted_virtual
{
 private:
  typedef Cached_face_degeneracy_tester<Face_tester_t>  Base;

 public:
  ~Ref_counted_face_degeneracy_tester_base() {}
};


//=========================================================================
//=========================================================================



template<class T>
struct Edge_degeneracy_tester_types
{
  typedef typename T::result_type            result_type;
  typedef typename T::Arity                  Arity;

  typedef typename T::Dual_graph             Dual_graph;
  typedef typename T::Edge                   Edge;
  typedef typename T::Edge_circulator        Edge_circulator;
  typedef typename T::All_edges_iterator     All_edges_iterator;
  typedef typename T::Finite_edges_iterator  Finite_edges_iterator;
};

//=========================================================================

template<class T>
struct Face_degeneracy_tester_types
{
  typedef typename T::result_type               result_type;
  typedef typename T::Arity                     Arity;

  typedef typename T::Dual_graph                Dual_graph;
  typedef typename T::Vertex_handle             Vertex_handle;
  typedef typename T::Vertex_circulator         Vertex_circulator;
  typedef typename T::All_vertices_iterator     All_vertices_iterator;
  typedef typename T::Finite_vertices_iterator  Finite_vertices_iterator;
};

//=========================================================================

template<class Edge_tester_t>
class Ref_counted_edge_degeneracy_tester
  : public Handle_to_tester_adaptor
  <Handle_for_virtual
   <Ref_counted_edge_degeneracy_tester_base<Edge_tester_t> >,
   Edge_degeneracy_tester_types<Edge_tester_t>
   >
{
 private:
  typedef Edge_degeneracy_tester_types<Edge_tester_t>  Edge_tester_types;

  typedef
  Ref_counted_edge_degeneracy_tester_base<Edge_tester_t>
  Ref_counted_tester_base;

  typedef Handle_for_virtual<Ref_counted_tester_base>
  Ref_counted_tester_base_handle;

  typedef Handle_to_tester_adaptor<Ref_counted_tester_base_handle,
				   Edge_tester_types>
  Base;
};

//=========================================================================

template<class Face_tester_t>
class Ref_counted_face_degeneracy_tester
  : public Handle_to_tester_adaptor
  <Handle_for_virtual
   <Ref_counted_face_degeneracy_tester_base<Face_tester_t> >,
   Face_degeneracy_tester_types<Face_tester_t>
  >
{
 private:
  typedef Face_degeneracy_tester_types<Face_tester_t>  Face_tester_types;

  typedef
  Ref_counted_face_degeneracy_tester_base<Face_tester_t>
  Ref_counted_tester_base;

  typedef Handle_for_virtual<Ref_counted_tester_base>
  Ref_counted_tester_base_handle;

  typedef Handle_to_tester_adaptor<Ref_counted_tester_base_handle,
				   Face_tester_types>
  Base;
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
