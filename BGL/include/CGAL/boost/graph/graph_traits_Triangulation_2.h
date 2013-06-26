// Copyright (c) 2007  GeometryFactory (France).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Andreas Fabri, Fernando Cacciola

#ifndef CGAL_GRAPH_TRAITS_TRIANGULATION_2_H
#define CGAL_GRAPH_TRAITS_TRIANGULATION_2_H

#include <boost/config.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

// The functions and classes in this file allows the user to
// treat a CGAL Triangulation_2 object as a boost graph "as is". No
// wrapper is needed for the Triangulation_2 object.


namespace CGAL {

  namespace detail {

template < class T, class EdgeBase >
class Edge : public EdgeBase {
  typedef typename T::Face_handle Face_handle ;
  public:
  
  Edge()
  {}

  Edge(Face_handle  fh, int i)
    : EdgeBase(fh,i)
  {}
  
  Edge(const EdgeBase& e)
    : EdgeBase(e)
  {}

  Edge(const Edge& e)
    : EdgeBase(e)
  {}

  Edge&
  operator=(const Edge& e)
  {
    this->first = e.first;
    this->second = e.second;
    return *this;
}

  bool operator==(const Edge& other) const
  {
    if((this->first == other.first)&&(this->second == other.second)) return true;
    Face_handle fh = this->first->neighbor(this->second);
    if(other.first != fh) return false;
    int i = fh->index(this->first);
    return (other.second == i);
  }

  bool operator!=(Edge& other) const
  {
    return ! (*this == other);
  }
};

template <class Circ, class E>
class Out_edge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  Out_edge_circulator()
    : Circ()
  {}

  Out_edge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
  {
    E ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed.first->neighbor(ed.second), ed.first->neighbor(ed.second)->index(ed.first));
    return e;
  }
};
 
template <class Circ, class E>
class In_edge_circulator : public Circ
{
private:
  mutable E e;

public:

  typedef E value_type;
  typedef E* pointer;
  typedef E& reference;

  In_edge_circulator()
    : Circ()
  {}

  In_edge_circulator(Circ c)
    : Circ(c)
  {}

  const E& operator*() const
{
    typename Circ::value_type ed = static_cast<const Circ*>(this)->operator*();
    e = E(ed);
    return e;
  }
};
  

  //  The vertex iterator of the bgl must evaluate to a vertex handle, not to a vertex
template < class T>
class boost_all_vertices_iterator {
protected:
 typename T::All_vertices_iterator nt;
public:
  typedef typename T::All_vertices_iterator  Iterator;
  typedef boost_all_vertices_iterator<T> Self;

  typedef typename std::iterator_traits<Iterator>::iterator_category iterator_category;
  typedef typename T::Vertex_handle  value_type;
  typedef typename std::iterator_traits<Iterator>::difference_type           difference_type;
  typedef value_type      reference;
  typedef value_type      pointer;

  // CREATION
  // --------

  boost_all_vertices_iterator()
  {}

  boost_all_vertices_iterator( Iterator j) : nt(j) {}

  // OPERATIONS Forward Category
  // ---------------------------


  bool operator==( const Self& i) const { return ( nt == i.nt); }
  bool operator!=( const Self& i) const { return !(nt == i.nt );   }
  value_type  operator*() const  { return nt; }
  value_type    operator->()  { return nt; }

  Self& operator++() {
    ++nt;
    return *this;
  }

  Self  operator++(int) {
    Self tmp = *this;
    ++*this;
    return tmp;
  }
};


  } // namespace detail
} // namespace CGAL

namespace boost { 

  template <class GT, class TDS>
  struct graph_traits< CGAL::Triangulation_2<GT,TDS> > {

    struct T2_graph_traversal_category : 
      public virtual bidirectional_graph_tag,
      public virtual adjacency_graph_tag,
      public virtual edge_list_graph_tag,
      public virtual vertex_list_graph_tag { };

    typedef CGAL::Triangulation_2<GT,TDS> Triangulation;

    typedef typename CGAL::Triangulation_2<GT,TDS>::Vertex_handle vertex_descriptor;
    typedef CGAL::detail::Edge<CGAL::Triangulation_2<GT,TDS>, typename CGAL::Triangulation_2<GT,TDS>::Edge>  edge_descriptor;
    typedef typename CGAL::Triangulation_2<GT,TDS>::All_edges_iterator  edge_iterator;

    typedef CGAL::detail::boost_all_vertices_iterator<Triangulation> vertex_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::Out_edge_circulator<typename Triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > out_edge_iterator;
    typedef CGAL::Counting_iterator<CGAL::detail::In_edge_circulator<typename Triangulation::Edge_circulator, edge_descriptor>, edge_descriptor > in_edge_iterator;
    typedef CGAL::Counting_iterator<typename Triangulation::Vertex_circulator> Incident_vertices_iterator;
    typedef Incident_vertices_iterator adjacency_iterator;

    typedef undirected_tag directed_category;
    typedef disallow_parallel_edge_tag edge_parallel_category; 
    typedef T2_graph_traversal_category traversal_category;
    typedef typename Triangulation::size_type size_type;
    typedef size_type vertices_size_type;
    typedef size_type edges_size_type;
    typedef size_type degree_size_type;
  };


} // namespace boost


namespace boost {

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor
  source(typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::edge_descriptor e,
         const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.ccw(e.second));
  }

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor
  target(typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::edge_descriptor e,
         const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    return e.first->vertex(g.cw(e.second));
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_iterator,
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_iterator >  
  vertices(const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typedef typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_iterator
      Iter;
    return std::make_pair( Iter(g.all_vertices_begin()), Iter(g.all_vertices_end()) );
  }


  template <class Gt, class Tds>
  inline std::pair<
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::edge_iterator,
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::edge_iterator >  
  edges(const CGAL::Triangulation_2<Gt,Tds>& g)
  {    
    return std::make_pair(g.all_edges_begin(), g.all_edges_end());
  }

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type
  out_degree(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::out_edge_iterator,
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::out_edge_iterator >  
  out_edges(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >
      ::out_edge_iterator Iter;
    
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::in_edge_iterator,
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::in_edge_iterator >  
  in_edges(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Triangulation_2<Gt,Tds>::Edge_circulator ec(u,u->face());
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >
      ::in_edge_iterator Iter;
    return std::make_pair( Iter(ec), Iter(ec,out_deg) );
  }

  template <class Gt, class Tds>
  inline std::pair<
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::adjacency_iterator,
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::adjacency_iterator >  
  adjacent_vertices(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename CGAL::Triangulation_2<Gt,Tds>::Vertex_circulator vc = out_edge_iterator(u,u.face());
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type out_deg = out_degree(u,g);
    typedef typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >
      ::adjacency_iterator Iter;
    return std::make_pair( Iter(vc), Iter(vc,out_deg) );
  }

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertices_size_type
  num_vertices(const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    return g.number_of_vertices()+1;
  }  

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::edges_size_type
  num_edges(const CGAL::Triangulation_2<Gt,Tds>& g)
  {

    return  g.number_of_vertices() + 1 + g.number_of_faces() + degree(g.infinite_vertex(), g) - 2;
  }  

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type
  in_degree(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }

  template <class Gt, class Tds>
  typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type
  degree(
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::vertex_descriptor u, 
    const CGAL::Triangulation_2<Gt,Tds>& g)
  {
    typename graph_traits< CGAL::Triangulation_2<Gt,Tds> >::degree_size_type deg = 0;
    typename CGAL::Triangulation_2<Gt,Tds>::Edge_circulator c = g.incident_edges(u), done(c);
    if ( c != 0) {
        do {
            ++deg;
        } while (++c != done);
    }
    return deg;
  }


  // property maps
  template <class Gt, class Tds>
  class T2_vertex_id_map
    : public put_get_helper<int, T2_vertex_id_map<Gt,Tds> >
  {
  public:
    typedef readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Vertex_handle key_type;
    
    T2_vertex_id_map()
    {}
    
    long operator[](key_type vh) const {
      return vh->id(); 
    }
  };

  template <class Gt, class Tds>
  class T2_edge_id_map
    : public put_get_helper<int, T2_edge_id_map<Gt,Tds> >
  {
  public:
    typedef readable_property_map_tag category;
    typedef int value_type;
    typedef int reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Edge key_type;
    
    T2_edge_id_map()
    {}
    
    long operator[](key_type e) const {
      return (3 * e.first.id()) + e.second; 
    }
  };

  template <class Gt, class Tds>
  class T2_edge_weight_map
    : public put_get_helper<typename Gt::FT, T2_edge_weight_map<Gt, Tds> >
  {
  private:
    const CGAL::Triangulation_2<Gt,Tds>& tr;
  public:
    typedef readable_property_map_tag category;
    typedef typename Gt::FT value_type;
    typedef value_type reference;
    typedef typename CGAL::Triangulation_2<Gt,Tds>::Edge key_type;

    T2_edge_weight_map(const CGAL::Triangulation_2<Gt,Tds>& tr_) 
      : tr(tr_) 
    { }

    value_type operator[](key_type e) const {
      return tr.segment(e).squared_length();
    }
  };


  template <class Gt, class Tds>
  inline T2_vertex_id_map<Gt,Tds>
  get(vertex_index_t, const CGAL::Triangulation_2<Gt,Tds>&) {
    T2_vertex_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline T2_edge_id_map<Gt,Tds>
  get(edge_index_t, const CGAL::Triangulation_2<Gt,Tds>&) {
    T2_edge_id_map<Gt,Tds> m;
    return m;
  }

  template <class Gt, class Tds>
  inline T2_edge_weight_map<Gt,Tds>
  get(edge_weight_t, const CGAL::Triangulation_2<Gt,Tds>& g) {
    T2_edge_weight_map<Gt,Tds> m(g);
    return m;
  }

  template <class Tag>
  struct T2_property_map { };

  template <>
  struct T2_property_map<vertex_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_vertex_id_map<Gt,Tds> type;
      typedef T2_vertex_id_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct T2_property_map<edge_index_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_edge_id_map<Gt,Tds> type;
      typedef T2_edge_id_map<Gt,Tds> const_type;
    };
  };


  template <>
  struct T2_property_map<edge_weight_t> {
    template <class Gt, class Tds>
    struct bind_ {
      typedef T2_edge_weight_map<Gt,Tds> type;
      typedef T2_edge_weight_map<Gt,Tds> const_type;
    };
  };

  // g++ 'enumeral_type' in template unification not implemented workaround
  template <class Gt, class Tds, class Tag>
  struct property_map<CGAL::Triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
      T2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  // see struct property_map in Polyehdron for an explanation
  template <class Gt, class Tds, class Tag>
  struct property_map<const CGAL::Triangulation_2<Gt,Tds>, Tag> {
    typedef typename 
      T2_property_map<Tag>::template bind_<Gt,Tds> map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
  };

  template <class Gt, class Tds, class PropertyTag, class Key>
  inline
  typename boost::property_traits<
    typename boost::property_map<CGAL::Triangulation_2<Gt,Tds>,PropertyTag>::const_type>::value_type
  get(PropertyTag p, const CGAL::Triangulation_2<Gt,Tds>& g, const Key& key) {
    return get(get(p, g), key);
  }
  
  template <class Gt, class Tds, class PropertyTag, class Key,class Value>
  inline void
  put(PropertyTag p, CGAL::Triangulation_2<Gt,Tds>& g, 
      const Key& key, const Value& value)
  {
    typedef typename property_map<CGAL::Triangulation_2<Gt,Tds>, PropertyTag>::type Map;
    Map pmap = get(p, g);
    put(pmap, key, value);
  }


  // What are those needed for ???
  template <typename Gt, typename Tds>
  struct edge_property_type<CGAL::Triangulation_2<Gt,Tds> > {
    typedef void type;
  };  

  template <typename Gt, typename Tds>
  struct vertex_property_type<CGAL::Triangulation_2<Gt,Tds> > {
    typedef void type;
  };
} // namespace boost

//#include <CGAL/graph_traits_Delaunay_triangulation_2.h>

#endif // CGAL_GRAPH_TRAITS_TRIANGULATION_2_H
