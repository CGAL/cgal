// Copyright (c) 2006  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Fernando Cacciola <fernando.cacciola@geometryfactory.com>
//
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_H 1

#include <CGAL/Surface_mesh_simplification/Detail/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>

namespace CGAL {

namespace Surface_mesh_simplification
{

//
// Implementation of the vertex-pair collapse triangulated surface mesh simplification algorithm
//
template<class ECM_
        ,class ShouldStop_
        ,class VertexIndexMap_
        ,class EdgeIndexMap_
        ,class EdgeIsBorderMap_
        ,class GetCost_
        ,class GetPlacement_
        ,class VisitorT_
        >
class EdgeCollapse
{
public:

  typedef ECM_              ECM ;
  typedef ShouldStop_       ShouldStop ;
  typedef VertexIndexMap_   VertexIndexMap ;
  typedef EdgeIndexMap_     EdgeIndexMap ;
  typedef EdgeIsBorderMap_  EdgeIsBorderMap ;
  typedef GetCost_          GetCost ;
  typedef GetPlacement_     GetPlacement ;
  typedef VisitorT_         VisitorT ;
  
  typedef EdgeCollapse Self ;
  
  typedef Edge_profile<ECM> Profile ;
  
  typedef boost::graph_traits  <ECM>       GraphTraits ;
  typedef boost::graph_traits  <ECM const> ConstGraphTraits ;
  typedef halfedge_graph_traits<ECM>       HalfedgeGraphTraits ; 
  
  typedef typename GraphTraits::vertex_descriptor      vertex_descriptor ;
  typedef typename GraphTraits::vertex_iterator        vertex_iterator ;
  typedef typename GraphTraits::edge_descriptor        edge_descriptor ;
  typedef typename GraphTraits::edge_iterator          edge_iterator ;
  typedef typename GraphTraits::out_edge_iterator      out_edge_iterator ;
  typedef typename GraphTraits::in_edge_iterator       in_edge_iterator ;
  typedef typename GraphTraits::traversal_category     traversal_category ;
  typedef typename GraphTraits::edges_size_type        size_type ;
  
  typedef typename ConstGraphTraits::vertex_descriptor const_vertex_descriptor ;
  typedef typename ConstGraphTraits::edge_descriptor   const_edge_descriptor ;
  typedef typename ConstGraphTraits::in_edge_iterator  const_in_edge_iterator ;
  
  typedef typename HalfedgeGraphTraits::undirected_edge_iterator undirected_edge_iterator ;
  typedef typename HalfedgeGraphTraits::Point                    Point ;

  typedef typename GetCost     ::result_type Cost_type ;
  typedef typename GetPlacement::result_type Placement_type ;
  
  typedef typename Kernel_traits<Point>::Kernel Kernel ;
  
  typedef typename Kernel::Equal_3 Equal_3 ;
  
  typedef typename Kernel::Vector_3 Vector ;
  typedef typename Kernel::FT       FT ;

  struct Compare_id
  {
    Compare_id() : mAlgorithm(0) {}
    
    Compare_id( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}
    
    bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const 
    {
      return mAlgorithm->get_directed_edge_id(a) < mAlgorithm->get_directed_edge_id(b);
    }
    
    Self const* mAlgorithm ;
  } ;
  
  struct Compare_cost
  {
    Compare_cost() : mAlgorithm(0) {}
    
    Compare_cost( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}
    
    bool operator() ( edge_descriptor const& a, edge_descriptor const& b ) const
    { 
      // NOTE: A cost is an optional<> value.
      // Absent optionals are ordered first; that is, "none < T" and "T > none" for any defined T != none.
      // In consequence, edges with undefined costs will be promoted to the top of the priority queue and poped out first.
      return mAlgorithm->get_data(a).cost() < mAlgorithm->get_data(b).cost();
    }
    
    Self const* mAlgorithm ;
  } ;
  
  struct Undirected_edge_id : boost::put_get_helper<size_type, Undirected_edge_id>
  {
    typedef boost::readable_property_map_tag category;
    typedef size_type                        value_type;
    typedef size_type                        reference;
    typedef edge_descriptor                  key_type;
    
    Undirected_edge_id() : mAlgorithm(0) {}
    
    Undirected_edge_id( Self const* aAlgorithm ) : mAlgorithm(aAlgorithm) {}
    
    size_type operator[] ( edge_descriptor const& e ) const { return mAlgorithm->get_undirected_edge_id(e); }
    
    Self const* mAlgorithm ;
  } ;
  
  
  typedef Modifiable_priority_queue<edge_descriptor,Compare_cost,Undirected_edge_id> PQ ;
  typedef typename PQ::handle pq_handle ;
  
  // An Edge_data is associated with EVERY _undirected_ edge in the mesh (collapsable or not).
  // It relates the edge with the PQ-handle needed to update the priority queue
  // It also relates the edge with a policy-based cache
  class Edge_data
  {
  public :
  
    Edge_data() : mPQHandle() {}
    
    Cost_type const& cost() const { return mCost ; }
    Cost_type      & cost()       { return mCost ; }
    
    pq_handle PQ_handle() const { return mPQHandle ;}
    
    bool is_in_PQ() const { return mPQHandle != PQ::null_handle() ; }
    
    void set_PQ_handle( pq_handle h ) { mPQHandle = h ; }
    
    void reset_PQ_handle() { mPQHandle = PQ::null_handle() ; }
    
  private:  
    
    Cost_type mCost ;
    pq_handle mPQHandle ;
  } ;
  typedef Edge_data* Edge_data_ptr ;
  typedef boost::scoped_array<Edge_data> Edge_data_array ;
  
  
public:

  EdgeCollapse( ECM&                    aSurface
              , ShouldStop       const& aShouldStop 
              , VertexIndexMap   const& aVertex_index_map 
              , EdgeIndexMap     const& aEdge_index_map 
              , EdgeIsBorderMap  const& aEdge_is_border_map 
              , GetCost          const& aGetCost
              , GetPlacement     const& aGetPlacement
              , VisitorT                aVisitor
              ) ;
  
  int run() ;
  
private:
  
  void Collect();
  void Loop();
  bool Is_collapse_topologically_valid( Profile const& aProfile ) ;
  bool Is_tetrahedron( edge_descriptor const& h1 ) ;
  bool Is_open_triangle( edge_descriptor const& h1 ) ;
  bool Is_collapse_geometrically_valid( Profile const& aProfile, Placement_type aPlacement ) ;
  void Collapse( Profile const& aProfile, Placement_type aPlacement ) ;
  void Update_neighbors( vertex_descriptor const& aKeptV ) ;
  
  Profile create_profile ( edge_descriptor const& aEdge )
  { 
    return Profile(aEdge,mSurface,Vertex_index_map,Edge_index_map,Edge_is_border_map);
  }  
  
  size_type get_directed_edge_id   ( const_edge_descriptor const& aEdge ) const { return Edge_index_map[aEdge]; }
  size_type get_undirected_edge_id ( const_edge_descriptor const& aEdge ) const { return get_directed_edge_id(aEdge) / 2 ; }

  bool is_primary_edge ( const_edge_descriptor const& aEdge ) const { return ( get_directed_edge_id(aEdge) % 2 ) == 0 ; }
  
  edge_descriptor primary_edge ( edge_descriptor const& aEdge ) 
  { 
    return is_primary_edge(aEdge) ? aEdge : opposite_edge(aEdge,mSurface) ;
  }  
    
  bool is_border ( const_edge_descriptor const& aEdge ) const { return Edge_is_border_map[aEdge] ; }    
  
  bool is_undirected_edge_a_border ( const_edge_descriptor const& aEdge ) const
  {
    return is_border(aEdge) || is_border(opposite_edge(aEdge,mSurface)) ;
  }    
  
  bool is_border ( const_vertex_descriptor const& aV ) const ;
  
  bool are_shared_triangles_valid( Point const& p0, Point const& p1, Point const& p2, Point const& p3 ) const ;
  
  edge_descriptor find_connection ( const_vertex_descriptor const& v0, const_vertex_descriptor const& v1 ) const ;
  
  vertex_descriptor find_exterior_link_triangle_3rd_vertex ( const_edge_descriptor const& e, const_vertex_descriptor const& v0, const_vertex_descriptor const& v1 ) const ;
  
  Edge_data& get_data ( edge_descriptor const& aEdge ) const 
  { 
    CGAL_assertion( is_primary_edge(aEdge) ) ;
    return mEdgeDataArray[get_undirected_edge_id(aEdge)];
  }
  
  Point const& get_point ( const_vertex_descriptor const& aV ) const
  {
    return get(vertex_point,mSurface,aV);
  }
  
  boost::tuple<const_vertex_descriptor,const_vertex_descriptor> get_vertices( const_edge_descriptor const& aEdge ) const
  {
    const_vertex_descriptor p,q ;
    p = boost::source(aEdge,mSurface);
    q = boost::target(aEdge,mSurface);
    return boost::make_tuple(p,q);
  }
  
  boost::tuple<vertex_descriptor,vertex_descriptor> get_vertices( edge_descriptor const& aEdge ) 
  {
    vertex_descriptor p,q ;
    p = boost::source(aEdge,mSurface);
    q = boost::target(aEdge,mSurface);
    return boost::make_tuple(p,q);
  }
  
  std::string vertex_to_string( const_vertex_descriptor const& v ) const
  {
    Point const& p = get_point(v);
    return boost::str( boost::format("[V%1%:%2%]") % v->id() % xyz_to_string(p) ) ;
  }
    
  std::string edge_to_string ( const_edge_descriptor const& aEdge ) const
  {
    const_vertex_descriptor p,q ; boost::tie(p,q) = get_vertices(aEdge);
    return boost::str( boost::format("{E%1% %2%->%3%}%4%") % aEdge->id() % vertex_to_string(p) % vertex_to_string(q) % ( is_border(aEdge) ? " (BORDER)" : ( is_border(aEdge->opposite()) ? " (~BORDER)": "" ) ) ) ;
  }
  
  Cost_type get_cost ( Profile const& aProfile ) const
  {
    return Get_cost(aProfile, get_placement(aProfile) );
  }
  
  Placement_type get_placement( Profile const& aProfile ) const
  {
    return Get_placement(aProfile);
  }
  
  void insert_in_PQ( edge_descriptor const& aEdge, Edge_data& aData ) 
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge)) ;
    CGAL_SURF_SIMPL_TEST_assertion(!aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(aEdge) ) ;

    aData.set_PQ_handle(mPQ->push(aEdge));

    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge) ) ;
  }
  
  void update_in_PQ( edge_descriptor const& aEdge, Edge_data& aData )
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge)) ;
    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge) ) ;

    aData.set_PQ_handle(mPQ->update(aEdge,aData.PQ_handle())) ; 

    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge) ) ;
  }   
  
  void remove_from_PQ( edge_descriptor const& aEdge, Edge_data& aData )
  {
    CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(aEdge)) ;
    CGAL_SURF_SIMPL_TEST_assertion(aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(mPQ->contains(aEdge) ) ;

    aData.set_PQ_handle(mPQ->erase(aEdge,aData.PQ_handle()));

    CGAL_SURF_SIMPL_TEST_assertion(!aData.is_in_PQ());
    CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(aEdge) ) ;
  }   
  
  optional<edge_descriptor> pop_from_PQ() 
  {
    optional<edge_descriptor> rEdge = mPQ->extract_top();
    if ( rEdge )
    {
      CGAL_SURF_SIMPL_TEST_assertion(is_primary_edge(*rEdge) ) ;
      CGAL_SURF_SIMPL_TEST_assertion(get_data(*rEdge).is_in_PQ()) ;

      get_data(*rEdge).reset_PQ_handle();

      CGAL_SURF_SIMPL_TEST_assertion(!get_data(*rEdge).is_in_PQ() ) ;
      CGAL_SURF_SIMPL_TEST_assertion(!mPQ->contains(*rEdge)) ;
    }  
    return rEdge ;  
  }
   
private:

  ECM&                   mSurface ;
  
  ShouldStop      const& Should_stop ;
  VertexIndexMap  const& Vertex_index_map ;
  EdgeIndexMap    const& Edge_index_map ;
  EdgeIsBorderMap const& Edge_is_border_map ;
  GetCost         const& Get_cost ;
  GetPlacement    const& Get_placement ;
  VisitorT               Visitor ; 
  
  
private:

  Edge_data_array mEdgeDataArray ;
  
  boost::scoped_ptr<PQ> mPQ ;
    
  std::size_t mInitialEdgeCount ;
  std::size_t mCurrentEdgeCount ; 

  FT          mcMaxDihedralAngleCos2 ;
  
  CGAL_ECMS_DEBUG_CODE ( unsigned mStep ; )
} ;

} // namespace Surface_mesh_simplification

} //namespace CGAL

#include <CGAL/Surface_mesh_simplification/Detail/Edge_collapse_impl.h>

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_I_EDGE_COLLAPSE_H //
// EOF //
 
