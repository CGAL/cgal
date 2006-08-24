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

template<class TSM_>
class Functor_base 
{
public:
    
  typedef TSM_ TSM ;
  
  typedef typename boost::graph_traits<TSM>::vertex_descriptor vertex_descriptor ;
  typedef typename boost::graph_traits<TSM>::edge_descriptor   edge_descriptor ;
  
  typedef typename Surface_geometric_traits<TSM>::Point_3 Point_3 ;
  typedef typename Surface_geometric_traits<TSM>::FT      FT ;
  
protected :

  virtual ~Functor_base() {}
  
  Point_3 const& get_point ( vertex_descriptor const& v, TSM& aSurface ) const
  {
    vertex_point_t vertex_point ;
    return get(vertex_point,aSurface,v) ;
  }
  
  tuple<vertex_descriptor,vertex_descriptor> get_vertices ( edge_descriptor const& aEdge, TSM& aSurface ) const
  {
    vertex_descriptor p = source(aEdge,aSurface);
    vertex_descriptor q = target(aEdge,aSurface);
    return make_tuple(p,q);
  }
};

template<class Collapse_data_, class Params_>
class Cost_functor_base : public Functor_base< typename Collapse_data_::TSM >
{
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  typedef Params_        Params ;
  
  typedef typename Collapse_data::TSM             TSM ;
  typedef typename Collapse_data::Has_cached_cost Has_cached_cost ;
  
  typedef Functor_base<TSM> Base ;
  
  typedef typename Base::edge_descriptor edge_descriptor ;
  typedef typename Base::FT              FT ;
  
  typedef optional<FT> result_type ;

public :
    
  result_type operator()( edge_descriptor const& aEdge
                        , TSM&                   aSurface
                        , Collapse_data const&   aData
                        , Params const*          aParams
                        ) const
  {
    return get_cost(aEdge,aSurface,aData,aParams,Has_cached_cost());
  }
  
private :
    
  result_type get_cost( edge_descriptor const& //aEdge
                      , TSM&                   //aSurface
                      , Collapse_data const&   aData
                      , Params const*          //aParams
                      , Tag_true
                      ) const
  {
    return aData.cost();  
  }
  
  result_type get_cost( edge_descriptor const& aEdge
                      , TSM&                   aSurface
                      , Collapse_data const&   //aData
                      , Params const*          aParams
                      , Tag_false 
                      ) const
  {
    return compute_cost(aEdge,aSurface,aParams);
  }
  
protected :

  virtual result_type compute_cost( edge_descriptor const& aEdge
                                  , TSM&                   aSurface
                                  , Params const*          aParams
                                  ) const = 0 ;
                      
};

template<class Collapse_data_, class Params_>
class Placement_functor_base : public Functor_base<typename Collapse_data_::TSM>
{
  
public:
    
  typedef Collapse_data_ Collapse_data ;
  typedef Params_        Params ;
  
  typedef typename Collapse_data::TSM                  TSM ;
  typedef typename Collapse_data::Has_cached_placement Has_cached_placement ;
  
  typedef Functor_base<TSM> Base ;
  
  typedef typename Base::edge_descriptor edge_descriptor ;
  typedef typename Base::Point_3         Point_3 ;
  
  typedef optional<Point_3> result_type ;
  
public :
    
  result_type operator()( edge_descriptor const& aEdge
                        , TSM&                   aSurface
                        , Collapse_data const&   aData
                        , Params const*          aParams
                        ) const
  {
    return get_placement(aEdge,aSurface,aData,aParams,Has_cached_placement());
  }
  
private :
    
  result_type get_placement( edge_descriptor const& //aEdge
                           , TSM&                   //aSurface
                           , Collapse_data const&   aData
                           , Params const*          //aParams
                           , Tag_true
                           ) const
  {
    return aData.placement();  
  }
  
  result_type get_placement( edge_descriptor const& aEdge
                           , TSM&                   aSurface
                           , Collapse_data const&   //aData
                           , Params const*          aParams
                           , Tag_false 
                           ) const
  {
    return compute_placement(aEdge,aSurface,aParams);
  }
  
protected :

  virtual result_type compute_placement( edge_descriptor const& aEdge
                                       , TSM&                   aSurface
                                       , Params const*          aParams
                                       ) const = 0 ;
};

} } } // namespace Triangulated_surface_mesh::Simplification::Edge_collapse


CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_EDGE_COLLAPSE_FUNCTOR_BASE_H
// EOF //
 
