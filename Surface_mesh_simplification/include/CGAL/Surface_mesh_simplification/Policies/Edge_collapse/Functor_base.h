// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/include/CGAL/Surface_mesh_simplification/Policies/Edge_length_cost.h $
// $Id: Edge_length_cost.h 32188 2006-07-03 17:19:46Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FUNCTOR_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FUNCTOR_BASE_H

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{


template<class Collapse_data_, class Base_>
class Functor_base : protected Base_
{
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Collapse_data::TSM TSM ;
  
  typedef typename grah_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename grah_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;
  
protected :

  typedef Base_ Base ;
      
  Point_3 const& get_point ( vertex_descriptor const& v, TSM const& aSurface ) const
  {
    vertex_point_t vertex_point ;
    return get(vertex_point,aSurface,v) ;
  }
  
  tuple<vertex_descriptor,vertex_descriptor> get_vertices ( edge_descriptor const& aEdge, TSM const& aSurface ) const
  {
    vertex_descriptor p = source(aEdge,aSurface);
    vertex_descriptor q = target(aEdge,aSurface);
    return make_tuple(p,q);
  }
};

template<class Collapse_data_, class Base_>
class Cost_functor_base : protected Functor_base<Collapse_data_,Base_>
{
  typedef Functor_base<Collapse_data,Base_> Base ;
  
  typedef typename Base::Base Base_base ;
  
  typedef typename Collapse_data::Has_cached_cost  Has_cached_cost ;
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Base::TSM             TSM ;
  typedef typename Base::edge_descriptor edge_descriptor ;
  typedef typename Base::FT              FT ;
  
  typedef optional<FT> result_type ;
    
public :
    
  result_type operator()( edge_descriptor const& aEdge
                        , TSM const&             aSurface
                        , Collapse_data const&   aData
                        ) const
  {
    return get_cost(aEdge,aSurface,aData,Has_cached_cost());
  }
  
private :
    
  result_type get_cost( edge_descriptor const& aEdge
                      , TSM const&             aSurface
                      , Collapse_data const&   aData
                      , Tag_true
                      ) const
  {
    return aData.cost();  
  }
  
  result_type get_cost( edge_descriptor const& aEdge
                      , TSM const&             aSurface
                      , Collapse_data const&   aData
                      , Tag_false 
                      ) const
  {
    return this->compute_cost(aEdge,aSurface);
  }
};

template<class Collapse_data_, class Base_>
class Placement_functor_base : protected Functor_base<Collapse_data_,Base_>
{
  typedef Functor_base<Collapse_data,Base_> Base ;
  
  typedef typename Base::Base Base_base ;
  
  typedef typename Collapse_data::Has_cached_placement Has_cached_placement ;
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  
  typedef typename Base::TSM             TSM ;
  typedef typename Base::edge_descriptor edge_descriptor ;
  typedef typename Base::Point_3         Point_3 ;
  
  
  typedef optional<Point_3> result_type ;
    
public :
    
  result_type operator()( edge_descriptor const& aEdge
                        , TSM const&             aSurface
                        , Collapse_data const&   aData
                        ) const
  {
    return get_placement(aEdge,aSurface,aData,Has_cached_placement());
  }
  
private :
    
  result_type get_placement( edge_descriptor const& aEdge
                           , TSM const&             aSurface
                           , Collapse_data const&   aData
                           , Tag_true
                           ) const
  {
    return aData.placement();  
  }
  
  result_type get_placement( edge_descriptor const& aEdge
                           , TSM const&             aSurface
                           , Collapse_data const&   aData
                           , Tag_false 
                           ) const
  {
    return this->compute_placement(aEdge,aSurface);
  }
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FUNCTOR_BASE_H
// EOF //
 
