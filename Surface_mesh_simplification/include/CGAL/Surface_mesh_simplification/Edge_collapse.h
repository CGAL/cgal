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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H 1

#include <vector>
#include <set>

#include <boost/config.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_array.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Collapse_operator.h>
#include <CGAL/Modifiable_priority_queue.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification { namespace Edge_collapse
{

//
// Implementation of the vertex-pair collapse triangulated surface mesh simplification algorithm
//
template<class TSM_
        ,class Params_ 
        ,class SetCollapseData_
        ,class GetCost_
        ,class GetPlacement_
        ,class ShouldStop_
        ,class VisitorT_
        >
class EdgeCollapse
{
public:

  typedef TSM_             TSM ;
  typedef Params_          Params ;
  typedef SetCollapseData_ SetCollapseData ;
  typedef GetCost_         GetCost ;
  typedef GetPlacement_    GetPlacement ;
  typedef ShouldStop_      ShouldStop ;
  typedef VisitorT_        VisitorT ;
  
  typedef EdgeCollapse Self ;
  
  typedef boost::graph_traits           <TSM>       GraphTraits ;
  typedef boost::undirected_graph_traits<TSM>       UndirectedGraphTraits ;
  typedef boost::graph_traits           <TSM const> ConstGraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor  vertex_descriptor ;
  typedef typename GraphTraits::vertex_iterator    vertex_iterator ;
  typedef typename GraphTraits::edge_descriptor    edge_descriptor ;
  typedef typename GraphTraits::edge_iterator      edge_iterator ;
  typedef typename GraphTraits::out_edge_iterator  out_edge_iterator ;
  typedef typename GraphTraits::in_edge_iterator   in_edge_iterator ;
  typedef typename GraphTraits::traversal_category traversal_category ;
  typedef typename GraphTraits::edges_size_type    size_type ;
  
  typedef typename ConstGraphTraits::vertex_descriptor const_vertex_descriptor ;
  typedef typename ConstGraphTraits::edge_descriptor   const_edge_descriptor ;
  
  typedef typename UndirectedGraphTraits::edge_iterator undirected_edge_iterator ;

  typedef typename GetCost     ::result_type Optional_cost_type ;
  typedef typename GetPlacement::result_type Optional_placement_type ;
  
  typedef typename SetCollapseData::Collapse_data Collapse_data ;
    
  typedef Surface_geometric_traits<TSM> Traits ;
  
  typedef typename Traits::Point_3 Point_3 ;

  typedef typename Kernel_traits<Point_3>::Kernel Kernel ;
  
  typedef typename Kernel::Equal_3 Equal_3 ;

  class Edge_data ;
  
  typedef Edge_data* Edge_data_ptr ;
  
  struct Compare_cost
  {
    Compare_cost( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}
    
    Comparison_result operator() ( edge_descriptor a, edge_descriptor b ) const
    {
      return mAlgorithm->compare_cost(a,b);
    }
    
    Self const* mAlgorithm ;
  } ;
  
  typedef Modifiable_priority_queue<edge_descriptor,Compare_cost> PQ ;
  
  typedef typename PQ::iterator pq_handle ;
  
  // An Edge_data is associated with EVERY edge in the mesh (collapsable or not).
  // It relates the edge with the PQ handle needed to update the priority queue
  // It also relates the edge with a policy-based collapse aData.
  class Edge_data
  {
  public :
  
    Edge_data() : mPQHandle() {}
    
    Collapse_data const& data() const { return mData ; }
    Collapse_data &      data()       { return mData ; }
    
    pq_handle PQ_handle() const { return mPQHandle ;}
    
    bool is_in_PQ() const { return mPQHandle != null_PQ_handle() ; }
    
    void set_PQ_handle( pq_handle h ) { mPQHandle = h ; }
    
    void reset_PQ_handle() { mPQHandle = null_PQ_handle() ; }
    
  private:
    
    static pq_handle null_PQ_handle() { pq_handle h ; return h ; }
    
  private:  
    
    Collapse_data mData ;
    pq_handle     mPQHandle ;
    
  } ;
    
  typedef boost::scoped_array<Edge_data> Edge_data_array ;
  
public:

  EdgeCollapse( TSM&                   aSurface
              , Params const*          aParams // Can be NULL
              , SetCollapseData const& aSetCollapseData
              , GetCost const&         aGetCost
              , GetPlacement const&    aGetPlacement
              , ShouldStop const&      aShouldStop 
              , VisitorT*              aVisitor
              ) ;
  
  int run() ;
  
private:
  
  void Collect();
  void Loop();
  bool Is_collapsable( edge_descriptor const& aEdge ) ;
  vertex_descriptor Collapse( edge_descriptor const& aEdge ) ;
  void Update_neighbors( vertex_descriptor const& aKeptV ) ;
  
  Point_3 const& get_point ( const_vertex_descriptor const& v ) const
  {
    vertex_point_t vertex_point ;
    return get(vertex_point,const_cast<TSM const&>(mSurface),v) ;
  }
  
  bool is_vertex_fixed ( const_vertex_descriptor const& v ) const
  {
    cgal_tsms_is_vertex_fixed_t is_vertex_fixed ;
    return get(is_vertex_fixed,mSurface,v) ;
  }
  
  bool is_border ( const_vertex_descriptor const& v ) const 
  {
    vertex_is_border_t is_border_vertex_property ;
    return get(is_border_vertex_property,mSurface,v) ;
  }    
  
  bool is_border ( const_edge_descriptor const& aEdge ) const
  {
    edge_is_border_t is_border_edge_property ;
    return get(is_border_edge_property,mSurface,aEdge) ;
  }    
  
  bool is_undirected_edge_a_border ( const_edge_descriptor const& aEdge ) const
  {
    return is_border(aEdge) || is_border(opposite_edge(aEdge,mSurface)) ;
  }    
  
  Edge_data_ptr get_data ( edge_descriptor const& aEdge ) const
  {
    cgal_tsms_edge_cached_pointer_t edge_cached_pointer ;
    return static_cast<Edge_data_ptr>(get(edge_cached_pointer,mSurface,aEdge)) ;
  }
  
  void set_data ( edge_descriptor const& aEdge, Edge_data_ptr aData )
  {
    cgal_tsms_edge_cached_pointer_t edge_cached_pointer ;
    put(edge_cached_pointer,mSurface,aEdge,aData) ;
    put(edge_cached_pointer,mSurface,opposite_edge(aEdge,mSurface),aData) ;
  }
  
  tuple<const_vertex_descriptor,const_vertex_descriptor> get_vertices( const_edge_descriptor const& aEdge ) const
  {
    const_vertex_descriptor p,q ;
    p = source(aEdge,mSurface);
    q = target(aEdge,mSurface);
    return make_tuple(p,q);
  }
  
  tuple<vertex_descriptor,vertex_descriptor> get_vertices( edge_descriptor const& aEdge ) 
  {
    vertex_descriptor p,q ;
    p = source(aEdge,mSurface);
    q = target(aEdge,mSurface);
    return make_tuple(p,q);
  }
  
  std::string vertex_to_string( const_vertex_descriptor const& v ) const
  {
    Point_3 const& p = get_point(v);
    return boost::str( boost::format("[V%1%:%2%]") % v->ID % xyz_to_string(p) ) ;
  }
    
  std::string edge_to_string ( const_edge_descriptor const& aEdge ) const
  {
    const_vertex_descriptor p,q ; tie(p,q) = get_vertices(aEdge);
    return boost::str( boost::format("{E%1% %2%->%3%}") % aEdge->ID % vertex_to_string(p) % vertex_to_string(q) ) ;
  }
  
  Optional_cost_type get_cost ( edge_descriptor const& aEdge ) const
  {
    Edge_data_ptr lData = get_data(aEdge);
    CGAL_assertion(lData);
    return Get_cost(aEdge,mSurface,lData->data(),mParams);
  }
  
  Optional_placement_type get_placement( edge_descriptor const& aEdge ) const
  {
    Edge_data_ptr lData = get_data(aEdge);
    CGAL_assertion(lData);
    return Get_placement(aEdge,mSurface,lData->data(),mParams);
  }
  
  Comparison_result compare_cost( edge_descriptor const& aEdgeA, edge_descriptor const& aEdgeB ) const
  {
    // NOTE: A cost is an optional<> value.
    // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
    // In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
    return CGAL::compare(get_cost(aEdgeA),get_cost(aEdgeB));
  }
  
  void insert_in_PQ( edge_descriptor const& aEdge, Edge_data_ptr aData ) 
  {
    CGAL_precondition(!aData->is_in_PQ());
    pq_handle h = mPQ.push(aEdge); 
    aData->set_PQ_handle(h);
  }
  
  void update_in_PQ( Edge_data_ptr aData )
  {
    CGAL_precondition(aData->is_in_PQ());
    aData->set_PQ_handle( mPQ.update(aData->PQ_handle()) ) ; 
  }   
  
  void remove_from_PQ( Edge_data_ptr aData )
  {
    CGAL_precondition(aData->is_in_PQ());
    mPQ.erase(aData->PQ_handle());
    aData->reset_PQ_handle();
  }   
  
  edge_descriptor pop_from_PQ() 
  {
    edge_descriptor rEdge ;
    if ( !mPQ.empty() )
    {
      rEdge = mPQ.top(); 
      mPQ.pop();
      get_data(rEdge)->reset_PQ_handle();
    }
    return rEdge ;
  }
    
private:

  TSM&          mSurface ;
  Params const* mParams ; // Can be NULL

  SetCollapseData const&  Set_collapse_data;   
  GetCost         const&  Get_cost ;
  GetPlacement    const&  Get_placement ;
  ShouldStop      const&  Should_stop ;
  VisitorT*               Visitor ;
  
private:

  Collapse_triangulation_edge<TSM> Collapse_triangulation_edge ;  

  Edge_data_array mEdgeDataArray ;
  
  PQ  mPQ ;
  
  std::size_t mInitialEdgeCount ;
  std::size_t mCurrentEdgeCount ; 

  CGAL_TSMS_DEBUG_CODE ( unsigned mStep ; )
} ;

} } } // namespace Triangulated_surface_mesh::Simplification::edge_collapse

CGAL_END_NAMESPACE

#include <CGAL/Surface_mesh_simplification/Edge_collapse_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H //
// EOF //
 
