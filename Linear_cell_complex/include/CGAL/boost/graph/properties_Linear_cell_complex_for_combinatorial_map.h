// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
// All rights reserved.
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H
#define CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H

#include <CGAL/boost/graph/properties.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#define CGAL_LCC_ARGS unsigned int d_, unsigned int ambient_dim,        \
             class Traits_, \
             class Items_, \
             class Alloc_, \
             template<unsigned int,class,class,class,class> \
             class CMap, \
             class Storage_

#define CGAL_NAME_LCC_ARGS d_, ambient_dim,        \
             Traits_, \
             Items_, \
             Alloc_, \
             CMap, \
             Storage_

#define CGAL_LCC_TYPE CGAL::Linear_cell_complex_for_combinatorial_map\
           <d_, ambient_dim, Traits_, Items_, Alloc_, CMap , Storage_>

namespace CGAL {

  /*template<class LCC>
class LCC_edge_weight_map :
    public boost::put_get_helper<double, LCC_edge_weight_map<LCC> >
{
public:

  typedef boost::readable_property_map_tag                   category;
  typedef double                                             value_type;
  typedef double                                             reference;
  typedef typename boost::graph_traits<LCC>::edge_descriptor key_type;

  LCC_edge_weight_map(LCC const& lcc) : m_lcc(lcc) {}

  //  reference operator[](key_type const& e)
  //{ return CGAL::squared_distance(m_lcc.point(e), m_lcc.point(m_lcc.other_extremity(e))); }
  reference operator[](key_type const& e) const
  { return CGAL::squared_distance(m_lcc.point(e), m_lcc.point(m_lcc.other_extremity(e))); }

protected:
  LCC const& m_lcc;
};

template <class CMap>
class CMap_dart_index_map_external :
    public boost::put_get_helper<std::size_t&, CMap_dart_index_map_external<CMap> >
{
public:
  typedef boost::read_write_property_map_tag                  category;
  typedef std::size_t                                         value_type;
  typedef std::size_t&                                        reference;
  typedef typename boost::graph_traits<CMap>::edge_descriptor key_type;
  typedef CGAL::Unique_hash_map<key_type,std::size_t>         Map;

  CMap_dart_index_map_external(): map(new Map)
  {}

  CMap_dart_index_map_external(CMap const& cm): map(new Map)
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::Dart_range::iterator Iter;
    Iter b=cmap.template darts().begin();
    Iter e=cmap.template darts().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      (*map)[b] = i;
    }
    b=cmap.template darts().begin();
    for(; b != e; ++b)
    {
      value_type v1 = (*map)[b];
      value_type v2 = (*map)[cmap.template beta<2>(b)];
      CGAL_assertion(v1==v2+1 || v2==v1+1);
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<Map&>(*map)[e]; }

private:
  boost::shared_ptr<Map> map;
};

template <class CMap>
class CMap_vertex_index :
    public boost::put_get_helper<std::size_t&, CMap_vertex_index<CMap> >
{
public:
  typedef boost::read_write_property_map_tag                    category;
  typedef std::size_t                                           value_type;
  typedef std::size_t&                                          reference;
  typedef typename boost::graph_traits<CMap>::vertex_descriptor key_type;
  typedef CGAL::Unique_hash_map<key_type,std::size_t>           Map;


  CMap_vertex_index(): map(new Map)
  {}

  CMap_vertex_index(CMap const& cm): map(new Map)
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template Attribute_range<0>::type::iterator Iter;
    Iter b=cmap.template attributes<0>().begin();
    Iter e=cmap.template attributes<0>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      (*map)[b] = i;
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<Map&>(*map)[e]; }

private:
  boost::shared_ptr<Map> map;
};

template <class CMap>
class CMap_face_index :
    public boost::put_get_helper<std::size_t&, CMap_face_index<CMap> >
{
public:
  typedef boost::read_write_property_map_tag                    category;
  typedef std::size_t                                           value_type;
  typedef std::size_t&                                          reference;
  typedef typename boost::graph_traits<CMap>::face_descriptor key_type;

  CMap_face_index()
  {}

  CMap_face_index(CMap const& cm)
  {
    CMap &cmap = const_cast<CMap&>(cm);
    typedef typename CMap::template Attribute_range<2>::type::iterator Iter;
    Iter b=cmap.template attributes<2>().begin();
    Iter e=cmap.template attributes<2>().end();
    for(value_type i=0; b != e; ++b, ++i)
    {
      map[b] = i;
    }
  }

  reference operator[](key_type const& e) const
  { return const_cast<CGAL::Unique_hash_map<key_type,std::size_t>&>(map)[e]; }

private:
  CGAL::Unique_hash_map<key_type,std::size_t> map ;
};


template<class LCC>
class LCC_vertex_point_map :
    public boost::put_get_helper< typename LCC::Point&, LCC_vertex_point_map<LCC> >
{
public:

  typedef typename LCC::Point Point;

  typedef boost::lvalue_property_map_tag category;
  typedef Point value_type;
  typedef Point& reference;
  typedef typename boost::graph_traits<LCC>::vertex_descriptor key_type;

  LCC_vertex_point_map(LCC& lcc) : m_lcc(lcc) {}

  reference operator[](key_type const& v) const
  { return m_lcc.point_of_vertex_attribute(v); }

protected:
  LCC & m_lcc;
};

template<class LCC>
class LCC_vertex_point_const_map :
    public boost::put_get_helper< typename LCC::Point const&, LCC_vertex_point_const_map<LCC> >
{
public:
  typedef typename LCC::Point Point;

  typedef boost::readable_property_map_tag category;
  typedef Point value_type;
  typedef Point const& reference;
  typedef typename boost::graph_traits<LCC const>::vertex_descriptor key_type;

  LCC_vertex_point_const_map(LCC const& lcc) : m_lcc(lcc) {}

  reference operator[](key_type const& v) const
  { return m_lcc.point_of_vertex_attribute(v); }

protected:
  const LCC & m_lcc;
  }; */

template<typename LCC, typename FT>
struct Wrap_squared_lcc
  : boost::put_get_helper< double, Wrap_squared_lcc<LCC, FT> >
{
  typedef typename boost::graph_traits<LCC>::edge_descriptor Handle;
  typedef FT value_type;
  typedef FT reference;
  typedef Handle key_type;
  typedef boost::readable_property_map_tag category;

  Wrap_squared_lcc(const LCC& alcc): m_lcc(alcc)
  {}

  template<typename E>
  FT operator[](const E& e) const
  {
    return approximate_sqrt(CGAL::squared_distance
                            (m_lcc.point(e.first_halfedge()),
                             m_lcc.point(e.second_halfedge())));
  }
private:
  const LCC& m_lcc;
};


// the tag we dispatch on from property_map<G, Property>
template <class Tag>
struct LCC_property_map {};

// generalized 2-ary get functions
template<CGAL_LCC_ARGS, class PropertyTag>
typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::const_type
get(PropertyTag, CGAL_LCC_TYPE const&)
{ return typename boost::property_map<CGAL_LCC_TYPE, PropertyTag >::const_type(); }

template<CGAL_LCC_ARGS, class PropertyTag>
typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::type
get(PropertyTag, CGAL_LCC_TYPE &)
{ return typename boost::property_map<CGAL_LCC_TYPE, PropertyTag >::type(); }

// generalized 3-ary get functions
template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename boost::property_traits<typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::type>::reference
get(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key)
{
  return get(get(p, g), key);
}

template<CGAL_LCC_ARGS, class PropertyTag, class Key>
typename boost::property_traits<typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::const_type>::reference
get(PropertyTag p, CGAL_LCC_TYPE const& g, const Key& key)
{
  return get(get(p, g), key);
}

// generalized put
template<CGAL_LCC_ARGS, class PropertyTag, class Key, class Value>
void put(PropertyTag p, CGAL_LCC_TYPE& g, const Key& key, const Value& value)
{
  typedef typename boost::property_map<CGAL_LCC_TYPE, PropertyTag>::type Map;
  Map pmap = get(p, g);
  put(pmap, key, value);
}

} // namespace CGAL

// specialization needs to be repeated for halfedge, vertex, face
#define CGAL_LCC_INDEX_PM(ENTITY, TAG, ACCESSOR)               \
  namespace CGAL {                                             \
  template<> struct LCC_property_map<boost::ENTITY##TAG> {     \
  template<CGAL_LCC_ARGS>                                      \
  struct bind_ {                                               \
    typedef internal::ACCESSOR##_accessor<                     \
    CGAL_LCC_TYPE,                                             \
    typename boost::graph_traits<CGAL_LCC_TYPE>                \
                                 ::ENTITY##_descriptor > type; \
    typedef type const_type;                                   \
  };                                                           \
  };                                                           \
  }

CGAL_LCC_INDEX_PM(halfedge, _index_t, Index)
CGAL_LCC_INDEX_PM(vertex, _index_t, Index)
CGAL_LCC_INDEX_PM(face, _index_t, Index)

#undef CGAL_LCC_INDEX_PM

namespace CGAL {
// not done with macros, because LCC_edge::id does not return a
// reference
template <>
struct LCC_property_map<boost::edge_index_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Edge_index_accessor<
      typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor> type;
    typedef type const_type;
  };
};

template <>
struct LCC_property_map<boost::edge_weight_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef typename Traits_::FT FT;
    typedef typename boost::graph_traits<CGAL_LCC_TYPE>::edge_descriptor
               edge_descriptor;
    typedef Wrap_squared_lcc<CGAL_LCC_TYPE, FT> type;
    typedef type const_type;
  };
};

template <>
struct LCC_property_map<vertex_point_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Point_accessor<
      typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor,
      typename Traits_::Point_3, typename Traits_::Point_3&> type;

    typedef internal::Point_accessor<
      typename boost::graph_traits<CGAL_LCC_TYPE>::vertex_descriptor,
      typename Traits_::Point_3, const typename Traits_::Point_3&> const_type;
  };
};

//
// external indices
//
/*
  template <>
struct LCC_property_map<edge_external_index_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Polyhedron_edge_index_map_external<CGAL_LCC_TYPE> type;
    typedef type const_type;
  };
};

template <>
struct LCC_property_map<halfedge_external_index_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Polyhedron_index_map_external<
        typename boost::graph_traits<CGAL_LCC_TYPE >::halfedge_descriptor> type;
    typedef type const_type;
  };
};


template <>
struct LCC_property_map<vertex_external_index_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Polyhedron_index_map_external<
      typename boost::graph_traits<CGAL_LCC_TYPE >::vertex_descriptor> type;
    typedef type const_type;
  };
};

template <>
struct LCC_property_map<face_external_index_t>
{
  template<CGAL_LCC_ARGS>
  struct bind_
  {
    typedef internal::Polyhedron_index_map_external<
      typename boost::graph_traits<CGAL_LCC_TYPE >::face_descriptor> type;
    typedef type const_type;
  };
};
*/

} // namespace CGAL

namespace CGAL{

/*
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t >::const_type
get(boost::halfedge_external_index_t, CGAL_LCC_TYPE const&)
{
  CGAL_LCC_TYPE& ncmap=const_cast<CGAL_LCC_TYPE&>(cmap);
  return typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t>::
    const_type(halfedges(ncmap).begin(), halfedges(ncmap).end(), num_halfedges(ncmap)); 
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t >::const_type
get(boost::vertex_external_index_t, CGAL_LCC_TYPE const&)
{
 CGAL_LCC_TYPE& ncmap=const_cast<CGAL_LCC_TYPE&>(cmap);
  return typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t>::
    const_type(vertices(ncmap).begin(), vertices(ncmap).end(), num_vertices(ncmap)); 
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t >::const_type
get(boost::edge_external_index_t, CGAL_LCC_TYPE const&)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t >::const_type
get(boost::face_external_index_t, CGAL_LCC_TYPE const&)
{
  CGAL_LCC_TYPE& ncmap=const_cast<CGAL_LCC_TYPE&>(cmap);
  return typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t>::
    const_type(faces(ncmap).begin(), faces(ncmap).end(), num_faces(ncmap));
}
*/

/*
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_index_t >::const_type
get(boost::edge_index_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::edge_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_index_t >::const_type
get(boost::halfedge_index_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::halfedge_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_index_t >::const_type
get(boost::vertex_index_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_index_t >::const_type
get(boost::face_index_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::face_index_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

*/

/*template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_point_t >::const_type
get(boost::vertex_point_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE,boost::vertex_point_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}
*/

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t >::const_type
get(boost::edge_weight_t, CGAL_LCC_TYPE const& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t>::
    const_type(const_cast<CGAL_LCC_TYPE&>(cmap));
}

// the same blurb for non-const
/*
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t >::type
get(boost::halfedge_external_index_t, CGAL_LCC_TYPE&)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::halfedge_external_index_t>::
    type(halfedges(cmap).begin(), halfedges(cmap).end(), num_halfedges(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t >::type
get(boost::vertex_external_index_t, CGAL_LCC_TYPE&)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::vertex_external_index_t>::
    type(vertices(cmap).begin(), vertices(cmap).end(), num_vertices(cmap));
}

template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t >::type
get(boost::edge_external_index_t, CGAL_LCC_TYPE& cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::edge_external_index_t>::
    type(const_cast<CGAL_LCC_TYPE&>(cmap));
}


template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t >::type
get(boost::face_external_index_t, CGAL_LCC_TYPE&)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::face_external_index_t>::
    type(faces(cmap).begin(), faces(cmap).end(), num_faces(cmap));
}
*/
template<CGAL_LCC_ARGS>
typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t >::type
get(boost::edge_weight_t, CGAL_LCC_TYPE & cmap)
{
  return typename boost::property_map<CGAL_LCC_TYPE, boost::edge_weight_t>::
    type(cmap);
}

} // namespace CGAL

namespace boost {

// property_map dispatcher into LCC
template<CGAL_LCC_ARGS, class Tag>
struct property_map<CGAL_LCC_TYPE, Tag>
{
  typedef typename CGAL::LCC_property_map<Tag>::
      template bind_<CGAL_NAME_LCC_ARGS> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};

// property_map dispatcher into const LCC
template<CGAL_LCC_ARGS, class Tag>
struct property_map<const CGAL_LCC_TYPE, Tag>
{
  typedef typename CGAL::LCC_property_map<Tag>::
      template bind_<CGAL_NAME_LCC_ARGS> map_gen;
  typedef typename map_gen::type       type;
  typedef typename map_gen::const_type const_type;
};
  
  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, vertex_point_t>: CGAL::Tag_true {};
  
  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, edge_weight_t>: CGAL::Tag_true {};
  
  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, edge_index_t>: CGAL::Tag_true {};
  
  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, face_index_t>: CGAL::Tag_true {};

  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, halfedge_index_t>: CGAL::Tag_true {};

  template<CGAL_LCC_ARGS>
  struct graph_has_property<CGAL_LCC_TYPE, vertex_index_t>: CGAL::Tag_true {};  

} // namespace boost

#undef CGAL_LCC_ARGS
#undef CGAL_LCC_TYPE

#endif // CGAL_BOOST_GRAPH_PROPERTIES_LINEAR_CELL_COMPLEX_H

