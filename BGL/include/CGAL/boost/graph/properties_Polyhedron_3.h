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

#ifndef CGAL_BOOST_GRAPH_PROPERTIES_POLYHEDRON_3_H
#define CGAL_BOOST_GRAPH_PROPERTIES_POLYHEDRON_3_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Unique_hash_map.h>

#define CGAL_HDS_PARAM_ template < class Traits, class Items, class Alloc> class HDS

namespace CGAL {

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_edge_weight_map : public boost::put_get_helper<double, Polyhedron_edge_weight_map<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef boost::readable_property_map_tag                                category;
  typedef double                                                          value_type;
  typedef double                                                          reference;
  typedef typename boost::graph_traits<Polyhedron const>::edge_descriptor key_type;

  Polyhedron_edge_weight_map( Polyhedron const& ) {}

  reference operator[](key_type const& e) const
  {
    return CGAL::squared_distance(e->vertex()->point(), e->opposite()->vertex()->point());
  }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_edge_is_border_map : public boost::put_get_helper<bool, Polyhedron_edge_is_border_map<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef boost::readable_property_map_tag                                category;
  typedef bool                                                            value_type;
  typedef bool                                                            reference;
  typedef typename boost::graph_traits<Polyhedron const>::edge_descriptor key_type;

  Polyhedron_edge_is_border_map( Polyhedron const& ) {}

  reference operator[](key_type const& e) const { return e->is_border(); }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_edge_index_map_stored : public boost::put_get_helper<std::size_t, Polyhedron_edge_index_map_stored<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef boost::readable_property_map_tag                                category;
  typedef std::size_t                                                     value_type;
  typedef std::size_t                                                     reference;
  typedef typename boost::graph_traits<Polyhedron const>::edge_descriptor key_type;

  Polyhedron_edge_index_map_stored( Polyhedron const& ) {}

  reference operator[](key_type const& e) const { return e->id(); }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_edge_index_map_external : public boost::put_get_helper<std::size_t, Polyhedron_edge_index_map_external<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef boost::readable_property_map_tag                                category;
  typedef std::size_t                                                     value_type;
  typedef std::size_t                                                     reference;
  typedef typename boost::graph_traits<Polyhedron const>::edge_descriptor key_type;

  Polyhedron_edge_index_map_external( Polyhedron const& p) 
    : map( remove_const(p).halfedges_begin()
         , remove_const(p).halfedges_end()
         , 0
         , std::size_t(-1)
         , p.size_of_halfedges() 
         )
  {}

  reference operator[](key_type const& e) const { return map[e]; }
  
private:

  static Polyhedron& remove_const ( Polyhedron const& p ) { return const_cast<Polyhedron&>(p) ; }
  
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};


template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_vertex_point_map : public boost::put_get_helper< typename Gt::Point_3&, Polyhedron_vertex_point_map<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;

public:

  typedef typename Gt::Point_3 Point_3 ;
  
  typedef boost::lvalue_property_map_tag                              category;
  typedef Point_3                                                     value_type;
  typedef Point_3&                                                    reference;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor key_type;

  Polyhedron_vertex_point_map( Polyhedron& ) {}

  reference operator[](key_type const& v) const { return v->point(); }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_vertex_point_const_map : public boost::put_get_helper< typename Gt::Point_3 const&
                                                                      , Polyhedron_vertex_point_const_map<Gt, I, HDS, A> 
                                                                      >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;

public:

  typedef typename Gt::Point_3 Point_3 ;
  
  typedef boost::readable_property_map_tag                                  category;
  typedef Point_3                                                           value_type;
  typedef Point_3 const&                                                    reference;
  typedef typename boost::graph_traits<Polyhedron const>::vertex_descriptor key_type;

  Polyhedron_vertex_point_const_map( Polyhedron const& ) {}

  reference operator[](key_type const& v) const { return v->point(); }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_vertex_index_map_stored : public boost::put_get_helper<std::size_t, Polyhedron_vertex_index_map_stored<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;

public:

  typedef boost::readable_property_map_tag                                  category;
  typedef std::size_t                                                       value_type;
  typedef std::size_t                                                       reference;
  typedef typename boost::graph_traits<Polyhedron const>::vertex_descriptor key_type;

  Polyhedron_vertex_index_map_stored( Polyhedron const& ) {}

  reference operator[](key_type const& v) const { return v->id(); }
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
class Polyhedron_vertex_index_map_external : public boost::put_get_helper<std::size_t, Polyhedron_vertex_index_map_external<Gt, I, HDS, A> >
{
private:

  typedef CGAL::Polyhedron_3<Gt,I,HDS,A> Polyhedron ;
  
public:

  typedef boost::readable_property_map_tag                                  category;
  typedef std::size_t                                                       value_type;
  typedef std::size_t                                                       reference;
  typedef typename boost::graph_traits<Polyhedron const>::vertex_descriptor key_type;

  Polyhedron_vertex_index_map_external( Polyhedron const& p) 
    : map( remove_const(p).vertices_begin()
         , remove_const(p).vertices_end()
         , 0
         , std::size_t(-1)
         , p.size_of_vertices() 
         )
  {}

  reference operator[](key_type const& v) const { return map[v]; }
  
private:

  static Polyhedron& remove_const ( Polyhedron const& p ) { return const_cast<Polyhedron&>(p) ; }
  
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};


template <class Tag>
struct Polyhedron_property_map {};

template <>
struct Polyhedron_property_map<boost::edge_weight_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_edge_weight_map<Gt,I,HDS,A> type;
    typedef Polyhedron_edge_weight_map<Gt,I,HDS,A> const_type;
  };
};


template <>
struct Polyhedron_property_map<boost::edge_index_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_edge_index_map_stored<Gt,I,HDS,A> type;
    typedef Polyhedron_edge_index_map_stored<Gt,I,HDS,A> const_type;
  };
};

template <>
struct Polyhedron_property_map<edge_external_index_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_edge_index_map_external<Gt,I,HDS,A> type;
    typedef Polyhedron_edge_index_map_external<Gt,I,HDS,A> const_type;
  };
};

template <>
struct Polyhedron_property_map<edge_is_border_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_edge_is_border_map<Gt,I,HDS,A> type;
    typedef Polyhedron_edge_is_border_map<Gt,I,HDS,A> const_type;
  };
};

template <>
struct Polyhedron_property_map<vertex_point_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_vertex_point_map      <Gt,I,HDS,A> type;
    typedef Polyhedron_vertex_point_const_map<Gt,I,HDS,A> const_type;
  };
};

template <>
struct Polyhedron_property_map<boost::vertex_index_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_vertex_index_map_stored<Gt,I,HDS,A> type;
    typedef Polyhedron_vertex_index_map_stored<Gt,I,HDS,A> const_type;
  };
};

template <>
struct Polyhedron_property_map<vertex_external_index_t> 
{
  template<class Gt, class I, CGAL_HDS_PARAM_, class A>
  struct bind_ 
  {
    typedef Polyhedron_vertex_index_map_external<Gt,I,HDS,A> type;
    typedef Polyhedron_vertex_index_map_external<Gt,I,HDS,A> const_type;
  };
};

} //namespace CGAL

namespace boost
{

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_edge_weight_map<Gt,I,HDS,A> get( edge_weight_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_edge_weight_map<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_edge_is_border_map<Gt,I,HDS,A> get( CGAL::edge_is_border_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_edge_is_border_map<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_edge_index_map_stored<Gt,I,HDS,A> get( edge_index_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_edge_index_map_stored<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_edge_index_map_external<Gt,I,HDS,A> get( CGAL::edge_external_index_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_edge_index_map_external<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline
CGAL::Polyhedron_vertex_point_map<Gt,I,HDS,A> get(CGAL::vertex_point_t, CGAL::Polyhedron_3<Gt,I,HDS,A>& p) 
{
  CGAL::Polyhedron_vertex_point_map<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline
CGAL::Polyhedron_vertex_point_const_map<Gt,I,HDS,A> get(CGAL::vertex_point_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_vertex_point_const_map<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_vertex_index_map_stored<Gt,I,HDS,A> get(vertex_index_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_vertex_index_map_stored<Gt,I,HDS,A> m(p);
  return m;
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
inline 
CGAL::Polyhedron_vertex_index_map_external<Gt,I,HDS,A> get(CGAL::vertex_external_index_t, CGAL::Polyhedron_3<Gt,I,HDS,A> const& p) 
{
  CGAL::Polyhedron_vertex_index_map_external<Gt,I,HDS,A> m(p);
  return m;
}


template<class Gt, class I, CGAL_HDS_PARAM_, class A, class Tag>
struct property_map<CGAL::Polyhedron_3<Gt,I,HDS,A>, Tag> 
{
  typedef typename CGAL::Polyhedron_property_map<Tag>::
      template bind_<Gt,I,HDS,A> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

// This partial specialization shouldn't be needed but is due to a bug in Boost 1.51.
template<class Gt, class I, CGAL_HDS_PARAM_, class A, class Tag>
struct property_map<const CGAL::Polyhedron_3<Gt,I,HDS,A>, Tag> 
{
  typedef typename CGAL::Polyhedron_property_map<Tag>::
      template bind_<Gt,I,HDS,A> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

template<class Gt, class I, CGAL_HDS_PARAM_, class A, class PropertyTag, class Key>
inline
typename property_traits<typename property_map<CGAL::Polyhedron_3<Gt,I,HDS,A>,PropertyTag>::type>::reference
get(PropertyTag p, CGAL::Polyhedron_3<Gt,I,HDS,A>& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A, class PropertyTag, class Key>
inline
typename property_traits<typename property_map<CGAL::Polyhedron_3<Gt,I,HDS,A>,PropertyTag>::const_type>::reference
get(PropertyTag p, CGAL::Polyhedron_3<Gt,I,HDS,A> const& g, const Key& key) 
{
  return get(get(p, g), key);
}

template<class Gt, class I, CGAL_HDS_PARAM_, class A, class PropertyTag, class Key,class Value>
inline 
void put(PropertyTag p, CGAL::Polyhedron_3<Gt,I,HDS,A>& g, const Key& key, const Value& value)
{
  typedef typename property_map<CGAL::Polyhedron_3<Gt,I,HDS,A>, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}


// What are those needed for ???
template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct edge_property_type<CGAL::Polyhedron_3<Gt,I,HDS,A> > 
{
  typedef edge_weight_t type;
};  

template<class Gt, class I, CGAL_HDS_PARAM_, class A>
struct vertex_property_type<CGAL::Polyhedron_3<Gt,I,HDS,A> > 
{
  typedef CGAL::vertex_point_t type;
};
        
} // namespace boost

#undef CGAL_HDS_PARAM_

#endif // CGAL_BOOST_GRAPH_PROPERTIES_POLYHEDRON_3_H
