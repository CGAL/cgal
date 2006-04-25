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
// $URL: $
// $Id: $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H 1

#include <vector>

#include <boost/config.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/iterator_adaptors.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Handle_hash_function.h>

#include <CGAL/Surface_mesh_simplification/TSMS_common.h>
#include <CGAL/Surface_mesh_simplification/Collapse_operator.h>
#include <CGAL/Modifiable_priority_queue.h>

CGAL_BEGIN_NAMESPACE

namespace Triangulated_surface_mesh { namespace Simplification 
{

//
// Implementation of the vertex-pair collapse triangulated surface mesh simplification algorithm
//
template<class TSM_
        ,class SelectMap_
        ,class CostMap_
        ,class VertexPlacement_
        ,class StopPred_
        >
class VertexPairCollapse
{
public:

  typedef TSM_             TSM ;
  typedef SelectMap_       SelectMap ;
  typedef CostMap_         CostMap ;
  typedef VertexPlacement_ VertexPlacement ;
  typedef StopPred_        StopPred ;
  
  typedef VertexPairCollapse Self ;
  
  typedef boost::graph_traits<TSM> GraphTraits ;
  
  typedef typename GraphTraits::vertex_descriptor  vertex_descriptor ;
  typedef typename GraphTraits::vertex_iterator    vertex_iterator ;
  typedef typename GraphTraits::edge_descriptor    edge_descriptor ;
  typedef typename GraphTraits::edge_iterator      edge_iterator ;
  typedef typename GraphTraits::out_edge_iterator  out_edge_iterator ;
  typedef typename GraphTraits::in_edge_iterator   in_edge_iterator ;
  typedef typename GraphTraits::traversal_category traversal_category ;
  typedef typename GraphTraits::edges_size_type    size_type ;

  typedef typename boost::property_traits<CostMap>::value_type optional_cost_type ;
  
  typedef typename VertexPlacement::result_type optional_vertex_point_type ;
     
  typedef Surface_geometric_traits<TSM> Traits ;
  
  typedef typename Traits::Point_3 Point_3 ;

  class vertex_pair ;
    
  typedef boost::shared_ptr<vertex_pair> vertex_pair_ptr ;
  
  typedef Modifiable_priority_queue<vertex_pair_ptr> PQ ;
  
  typedef typename PQ::iterator pq_handle ;
  
  // The algoritm is centered around vertex-pairs, encapsulated in this type.
  // For each pair there is a cost value provided by the external CostMap.
  // Vertex-pairs are stored in a priority queue based on their cost.
  // 
  // Each edge in the TSM contributes one vertex_pair.
  // If the user chooses to collapse non-edge pairs too, some vertex-pairs won't be edges.
  //
  // The algortithm needs to update the costs of certain vertex-pairs (for those edges affected by a collapse).
  // For that reason, the vertex-pairs are kept in a relaxed_queue data structure that supports the update operation.
  // This relaxed_queue data structure needs an index property map. For that reason, an integer ID is stored in each
  // pair.
  //
  class vertex_pair
  {
  public :
  
    vertex_pair( size_type                aID
               , vertex_descriptor const& aP
               , vertex_descriptor const& aQ
               , edge_descriptor   const& aEdge
               , TSM*                     aSurface
               , CostMap const*           aCost_map 
               )
       : 
       mID(aID)
     , mP(aP)
     , mQ(aQ)
     , mEdge(aEdge)
     , mSurface(aSurface)
     , Cost_map(aCost_map)
     , mCostStored(false)
     , mIsFixed(false)
     , mMark(0)
    {}
    
    size_type id() const { return mID ; } 

    // The cost of collapsing a vertex-pair is cached in this record.
    // When calling cost() for the first time, the cached cost is taken from the external CostMap.
    // Subsequent calls to cost() returns the same cached cost, that is, CostMap[] IS NOT called each time.
    //
    // Such a cache cost is an optional<> value. This is because the CostMap can return "none" for too high or incomputable costs.
    // Therefore, even if the cost is cached, it might be "absent" (that is, == boost::none)
    //
    // OTOH, the algorithm can invalidate the cached cost of any given pair. 
    // If invalidate_cost() is called, the next call to cost() will take it again from the external CostMap.
    // NOTE: Whether the cost is cached or not is independent from whether it is absent or none.
    // The former is controlled by the algorithm by explicitely calling InvalidateCost() while the later is defined by the CostMap
    // which can retiurn boost::none for incomputable or logically invalid costs.
    //
    optional_cost_type cost() const { return UpdateCost() ; }
    
    void invalidate_cost() { mCostStored = false ; }
    
    vertex_descriptor p() const { return mP ; }
    vertex_descriptor q() const { return mQ ; }
    
    edge_descriptor edge() const { return mEdge ; }

    bool is_edge() const
    {
      edge_descriptor null ;
      return mEdge != null ;
    }

    // A pair is fixed if it cannot be collpased.
    bool  is_fixed() const { return mIsFixed ; }
    bool& is_fixed()       { return mIsFixed ; }
    
    int  mark() const { return mMark ; }
    int& mark()       { return mMark ; }

    pq_handle PQ_handle() const { return mPQHandle ;}
    
    bool is_in_PQ() const { return mPQHandle != null_PQ_handle() ; }
    
    void set_PQ_handle( pq_handle h ) { mPQHandle = h ; }
    
    void reset_PQ_handle() { mPQHandle = null_PQ_handle() ; }
    
    void update_vertex ( vertex_descriptor oldv, vertex_descriptor newv )
    {
      CGAL_assertion( mP == oldv || mQ == oldv );
      CGAL_assertion( oldv != newv );
      
      if ( mP == oldv )
           mP = newv ;
      else mQ = newv ;
      
      invalidate_cost();
    }
        
    friend bool operator< ( boost::shared_ptr<vertex_pair> const& a, boost::shared_ptr<vertex_pair> const& b ) 
    {
      // NOTE: cost() is an optional<> value.
      // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
      // In consequence, vertex-pairs with undefined costs will be promoted to the top of the priority queue and poped out first.
      return a->cost() < b->cost() ;
    }
    friend bool operator== ( boost::shared_ptr<vertex_pair> const& a, boost::shared_ptr<vertex_pair> const& b ) 
    {
      return a->cost() == b->cost() ;
    }

#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE

    bool is_p_in_surface   () const { return handle_exists(mSurface->vertices_begin(),mSurface->vertices_end(),p()) ; }
    bool is_q_in_surface   () const { return handle_exists(mSurface->vertices_begin(),mSurface->vertices_end(),q()) ; }
    bool is_edge_in_surface() const { return handle_exists(mSurface->halfedges_begin(),mSurface->halfedges_end(),edge()) ; }
         
    bool is_valid() const { return is_p_in_surface() && is_q_in_surface() && is_edge_in_surface() ; }
    
    friend std::ostream& operator<< ( std::ostream& out, vertex_pair const& vp ) 
    {
      out << "VP" << vp.mID << " {" ;
      if ( vp.is_p_in_surface() )
           out << "V" << vp.p()->ID << " (" << vp.p()->point().x() << "," << vp.p()->point().y() << "," << vp.p()->point().z() << ")" ;
      else out << "##p() has been erased## " ;
      out << "->" ;
      if ( vp.is_q_in_surface() )
           out << " V" << vp.q()->ID << "(" << vp.q()->point().x() << "," << vp.q()->point().y() << "," << vp.q()->point().z() << ")" ;
      else out << "##q() has been erased## " ;
      if ( vp.is_edge_in_surface() )
           out << " E" << vp.edge()->ID ;
      else out << "##e() has been erased## " ;
          
      if ( vp.mCostStored )
      {
        if ( vp.mCost )
             out << " [" << *vp.mCost << "]" ;
        else out << " [<none>]" ;
      }
      else out << " [<not cached>]" ;
      out << "}" ;
      return out ;
    }
#endif
    
  private:
    
    // Caches the cost if not currently cached, and returns it.
    optional_cost_type UpdateCost() const
    { 
      if ( !mCostStored )
      {
        mCost = boost::get(*Cost_map, boost::make_tuple(p(),q(),is_edge(),mSurface)) ; 
        mCostStored = true ;
      }  
      return mCost ;  
    }
    
    static pq_handle null_PQ_handle() { pq_handle h ; return h ; }
    
  private:  
    
    size_type          mID ;
    vertex_descriptor  mP ;
    vertex_descriptor  mQ ;
    edge_descriptor    mEdge ;
    TSM*               mSurface ;
    CostMap const*     Cost_map ;
    
    mutable bool               mCostStored ;
    mutable optional_cost_type mCost ;
    
    bool      mIsFixed ;
    pq_handle mPQHandle ;
    int       mMark ;
    
  } ;
  
  // Mapping from edges to vertex-pairs (as indexed into the Pairs_vector sequence)
  typedef Unique_hash_map<edge_descriptor,vertex_pair_ptr,Handle_hash_function> edges_2_pairs_map ;

  typedef std::vector<vertex_pair_ptr> vertex_pair_vector ;
  
  typedef typename vertex_pair_vector::iterator vertex_pair_vector_iterator ;
    
public:

  VertexPairCollapse( TSM&                   aSurface
                    , SelectMap       const& aSelectMap
                    , CostMap         const& aCostMap
                    , VertexPlacement const& aVertexPlacement
                    , StopPred        const& aStopPred 
                    , bool                   aIncludeNonEdgePairs 
                    ) ;
  
  int run() ;
  
private:
  
  void Collect();
  void Loop();
  bool Is_collapsable( vertex_pair_ptr const& aPair ) ;
  void Collapse( vertex_pair_ptr const& aPair ) ;
  void Update_neighbors( vertex_descriptor const& v ) ;
  
  vertex_pair_ptr get_pair ( edge_descriptor const& e ) 
  {
    CGAL_precondition( mEdgesToPairsMap.is_defined(e) ) ;
    return mEdgesToPairsMap[e] ;  
  }
  
  void set_pair ( edge_descriptor const& e, vertex_pair_ptr aPair )
  {
    edge_descriptor o = boost::opposite_edge(e,mSurface) ;
    
    mEdgesToPairsMap[e] = aPair ;  
    mEdgesToPairsMap[o] = aPair ;  

    CGAL_postcondition( get_pair(e)->id() == aPair->id() ) ;
    CGAL_postcondition( get_pair(o)->id() == aPair->id() ) ;
  }
  
  void insert_in_PQ( vertex_pair_ptr const& aPair ) 
  {
    CGAL_precondition(aPair);
    CGAL_precondition(!aPair->is_in_PQ());
    pq_handle h = mPQ.push(aPair); 
    aPair->set_PQ_handle(h);
  }
  
  void update_in_PQ( vertex_pair_ptr const& aPair )
  {
    CGAL_precondition(aPair);
    CGAL_precondition(aPair->is_in_PQ());
    aPair->set_PQ_handle( mPQ.update(aPair->PQ_handle()) ) ; 
  }   
  
  void remove_from_PQ( vertex_pair_ptr const& aPair )
  {
    CGAL_precondition(aPair);
    CGAL_precondition(aPair->is_in_PQ());
    mPQ.erase(aPair->PQ_handle());
    aPair->reset_PQ_handle();
  }   
  
  vertex_pair_ptr pop_from_PQ() 
  {
    vertex_pair_ptr rR ;
    if ( !mPQ.empty() )
    {
      rR = mPQ.top(); 
      mPQ.pop();
      rR->reset_PQ_handle();
    }
    
    return rR ;
  }

private:

  TSM& mSurface ;

  SelectMap const&       Select_map ;   
  CostMap const&         Cost_map ;
  VertexPlacement const& construct_new_vertex_point ;
  StopPred const&        stop_simplification ;

  Collapse_triangulation_edge<TSM> collapse_triangulation_edge ;  
  
  bool mIncludeNonEdgePairs;

  edges_2_pairs_map mEdgesToPairsMap ;
  PQ                mPQ ;
  
  std::size_t mInitialPairCount ;
  std::size_t mCurrentPairCount ; 

} ;

} } // namespace Triangulated_surface_mesh::Simplification

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Surface_mesh_simplification/Vertex_pair_collapse.C>
#endif

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_VERTEX_PAIR_COLLAPSE_IMPL_H //
// EOF //
 
