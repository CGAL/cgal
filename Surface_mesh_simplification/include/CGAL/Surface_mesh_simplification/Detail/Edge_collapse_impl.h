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
#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H

namespace CGAL {

namespace Surface_mesh_simplification 
{

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::EdgeCollapse( ECM&                    aSurface
                                                    , ShouldStop       const& aShould_stop
                                                    , VertexIndexMap   const& aVertex_index_map 
                                                    , EdgeIndexMap     const& aEdge_index_map 
                                                    , EdgeIsBorderMap  const& aEdge_is_border_map 
                                                    , GetCost          const& aGet_cost
                                                    , GetPlacement     const& aGet_placement
                                                    , VisitorT                aVisitor
                                                    )
  : 
   mSurface           (aSurface)
  ,Should_stop        (aShould_stop) 
  ,Vertex_index_map   (aVertex_index_map)
  ,Edge_index_map     (aEdge_index_map)
  ,Edge_is_border_map (aEdge_is_border_map)
  ,Get_cost           (aGet_cost)
  ,Get_placement      (aGet_placement)
  ,Visitor            (aVisitor)
{
  const FT cMaxDihedralAngleCos = std::cos( 1.0 * CGAL_PI / 180.0 ) ;
  
  mcMaxDihedralAngleCos2 = cMaxDihedralAngleCos * cMaxDihedralAngleCos ;
  
  CGAL_ECMS_TRACE(0,"EdgeCollapse of ECM with " << (num_edges(aSurface)/2) << " edges" ); 
  
  CGAL_ECMS_DEBUG_CODE ( mStep = 0 ; )
  
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
    vertex_iterator vb, ve ;     
    for ( boost::tie(vb,ve) = boost::vertices(mSurface) ; vb != ve ; ++ vb )  
      CGAL_ECMS_TRACE(1, vertex_to_string(*vb) ) ;
      
    undirected_edge_iterator eb, ee ;
    for ( boost::tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
      CGAL_ECMS_TRACE(1, edge_to_string(*eb) ) ;
#endif
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
int EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::run()
{
  CGAL_SURF_SIMPL_TEST_assertion( mSurface.is_valid() && mSurface.is_pure_triangle() ) ;

  Visitor.OnStarted(mSurface);
   
  // First collect all candidate edges in a PQ
  Collect(); 
  
  // Then proceed to collapse each edge in turn
  Loop();

  CGAL_ECMS_TRACE(0,"Finished: " << (mInitialEdgeCount - mCurrentEdgeCount) << " edges removed." ) ;

  int r = (int)(mInitialEdgeCount - mCurrentEdgeCount) ;
    
  Visitor.OnFinished(mSurface);
    
  return r ;
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
void EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Collect()
{
  CGAL_ECMS_TRACE(0,"Collecting edges...");

  //
  // Loop over all the _undirected_ edges in the surface putting them in the PQ
  //
  
  Equal_3 equal_points = Kernel().equal_3_object();
    
  size_type lSize = num_edges(mSurface) / 2 ;
  
  CGAL_SURF_SIMPL_TEST_assertion( ( lSize * 2 ) == mSurface.size_of_halfedges() ) ;

  mInitialEdgeCount = mCurrentEdgeCount = lSize;
  
  mEdgeDataArray.reset( new Edge_data[lSize] ) ;
  
  mPQ.reset( new PQ (lSize, Compare_cost(this), Undirected_edge_id(this) ) ) ;
  
  std::size_t id = 0 ;
  
  CGAL_SURF_SIMPL_TEST_assertion_code ( size_type lInserted    = 0 ) ;
  CGAL_SURF_SIMPL_TEST_assertion_code ( size_type lNotInserted = 0 ) ;

  std::vector<Profile> zero_length_edges;

  undirected_edge_iterator eb, ee ;
  for ( boost::tie(eb,ee) = undirected_edges(mSurface); eb!=ee; ++eb )
  {
    edge_descriptor lEdge = *eb ;
  
    CGAL_assertion( get_directed_edge_id(lEdge) == id ) ;
    CGAL_assertion( get_directed_edge_id(opposite_edge(lEdge,mSurface)) == id+1 ) ;

    Profile const& lProfile = create_profile(lEdge);
          
    if ( !equal_points(lProfile.p0(),lProfile.p1()) )
    {
      Edge_data& lData = get_data(lEdge);

      lData.cost() = get_cost(lProfile) ;
      insert_in_PQ(lEdge,lData);
      
      Visitor.OnCollected(lProfile,lData.cost());

      CGAL_SURF_SIMPL_TEST_assertion_code ( ++ lInserted ) ;
    }
    else
    {
      zero_length_edges.push_back(lProfile);
      CGAL_SURF_SIMPL_TEST_assertion_code ( ++ lNotInserted ) ;
    }

    
    CGAL_ECMS_TRACE(2,edge_to_string(lEdge));
    
    id += 2 ;
  }
 
  CGAL_SURF_SIMPL_TEST_assertion ( lInserted + lNotInserted == mInitialEdgeCount ) ;

  for (typename std::vector<Profile>::iterator it=zero_length_edges.begin(),it_end=zero_length_edges.end();it!=it_end;++it)
  {
    Placement_type lPlacement = get_placement(*it);
    Collapse(*it,lPlacement);
  }
  
  CGAL_ECMS_TRACE(0,"Initial edge count: " << mInitialEdgeCount ) ;
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
void EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Loop()
{
  CGAL_ECMS_TRACE(0,"Collapsing edges...") ;

  CGAL_SURF_SIMPL_TEST_assertion_code ( size_type lLoop_watchdog = 0 ) ;
  CGAL_SURF_SIMPL_TEST_assertion_code ( size_type lNonCollapsableCount = 0 ) ;

  //
  // Pops and processes each edge from the PQ
  //
  optional<edge_descriptor> lEdge ;
  while ( (lEdge = pop_from_PQ()) )
  {
    CGAL_SURF_SIMPL_TEST_assertion ( lLoop_watchdog ++ < mInitialEdgeCount ) ;

    CGAL_ECMS_TRACE(1,"Poped " << edge_to_string(*lEdge) ) ;
    
    Profile const& lProfile = create_profile(*lEdge);
      
    Cost_type lCost = get_data(*lEdge).cost();
    
    Visitor.OnSelected(lProfile,lCost,mInitialEdgeCount,mCurrentEdgeCount);
      
    if ( lCost ) 
    {
      if ( Should_stop(*lCost,lProfile,mInitialEdgeCount,mCurrentEdgeCount) )
      {
        Visitor.OnStopConditionReached(lProfile);
          
        CGAL_ECMS_TRACE(0,"Stop condition reached with InitialEdgeCount=" << mInitialEdgeCount
                       << " CurrentEdgeCount=" << mCurrentEdgeCount
                       << " Current Edge: " << edge_to_string(*lEdge)
                       );
        break ;
      }
        
      if ( Is_collapse_topologically_valid(lProfile) )
      {
        // The external function Get_new_vertex_point() is allowed to return an absent point if there is no way to place the vertex
        // satisfying its constrians. In that case the remaining vertex is simply left unmoved.
        Placement_type lPlacement = get_placement(lProfile);
        
        if ( Is_collapse_geometrically_valid(lProfile,lPlacement) )
          Collapse(lProfile,lPlacement);
      }
      else
      {
        CGAL_SURF_SIMPL_TEST_assertion_code ( lNonCollapsableCount++ ) ;

        Visitor.OnNonCollapsable(lProfile);
          
        CGAL_ECMS_TRACE(1,edge_to_string(*lEdge) << " NOT Collapsable"  );
      }  
    }
    else
    {
      CGAL_ECMS_TRACE(1,edge_to_string(*lEdge) << " uncomputable cost."  );
    }
  }
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::is_border( const_vertex_descriptor const& aV ) const
{
  bool rR = false ;
  
  const_in_edge_iterator eb, ee ; 
  for ( boost::tie(eb,ee) = in_edges(aV,mSurface) ; eb != ee ; ++ eb )
  {
    const_edge_descriptor lEdge = *eb ;
    if ( is_undirected_edge_a_border(lEdge) )
    {
      rR = true ;
      break ;
    }
  }  
    
  return rR ;  
}

// Some edges are NOT collapsable: doing so would break the topological consistency of the mesh.
// This function returns true if a edge 'p->q' can be collapsed.
//
// An edge p->q can be collapsed iff it satisfies the "link condition"
// (as described in the "Mesh Optimization" article of Hoppe et al (1993))
//
// The link conidition is as follows: for every vertex 'k' adjacent to both 'p and 'q',
// "p,k,q" is a facet of the mesh.
//
template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Is_collapse_topologically_valid( Profile const& aProfile )
{
  bool rR = true ;

  CGAL_ECMS_TRACE(3,"Testing topological collapsabilty of p_q=V" << aProfile.v0()->id() << "(%" << aProfile.v0()->vertex_degree() << ")"
                 << "->V" << aProfile.v1()->id() << "(%" << aProfile.v1()->vertex_degree() << ")" 
                 );
                 
  CGAL_ECMS_TRACE(4, "is p_q border:" << aProfile.is_v0_v1_a_border() );
  CGAL_ECMS_TRACE(4, "is q_q border:" << aProfile.is_v1_v0_a_border() );

  out_edge_iterator eb1, ee1 ; 
  out_edge_iterator eb2, ee2 ; 

  CGAL_ECMS_TRACE(4,"  t=V" 
                   << ( aProfile.left_face_exists() ? aProfile.vL()->id() : -1 )
                   << "(%" 
                   << ( aProfile.left_face_exists() ? aProfile.vL()->vertex_degree() : 0 ) 
                   << ")" 
                 );
  CGAL_ECMS_TRACE(4,"  b=V" 
                   << ( aProfile.right_face_exists() ? aProfile.vR()->id() : -1 )
                   << "(%" 
                   << ( aProfile.right_face_exists() ? aProfile.vR()->vertex_degree() :0 ) 
                   << ")" 
                   );

  // The following loop checks the link condition for v0_v1.
  // Specifically, that for every vertex 'k' adjacent to both 'p and 'q', 'pkq' is a face of the mesh.
  // 
  for ( boost::tie(eb1,ee1) = out_edges(aProfile.v0(),mSurface) ; rR && eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor v0_k = *eb1 ;
    
    if ( v0_k != aProfile.v0_v1() )
    {
      vertex_descriptor k = target(v0_k,mSurface);
      
      for ( boost::tie(eb2,ee2) = out_edges(k,mSurface) ; rR && eb2 != ee2 ; ++ eb2 )
      {
        edge_descriptor k_v1 = *eb2 ;

        if ( target(k_v1,mSurface) == aProfile.v1() )
        {
          // At this point we know p-q-k are connected and we need to determine if this triangle is a face of the mesh.
          //
          // Since the mesh is known to be triangular there are at most two faces sharing the edge p-q.
          //
          // If p->q is NOT a border edge, the top face is p->q->t where t is target(next(p->q))
          // If q->p is NOT a border edge, the bottom face is q->p->b where b is target(next(q->p))
          //
          // If k is either t or b then p-q-k *might* be a face of the mesh. It won't be if k==t but p->q is border
          // or k==b but q->b is a border (because in that case even though there exists triangles p->q->t (or q->p->b)
          // they are holes, not faces)
          // 
     
          bool lIsFace =   ( aProfile.vL() == k && aProfile.left_face_exists () )
                        || ( aProfile.vR() == k && aProfile.right_face_exists() ) ;
                        
          CGAL_SURF_SIMPL_TEST_assertion_code
          ( 
            if ( lIsFace )
            {
              // Is k_v1 the halfedge bounding the face 'k-v1-v0'?
              if ( !k_v1->is_border() && k_v1->next()->vertex() == aProfile.v0() )
              {
                CGAL_SURF_SIMPL_TEST_assertion( !k_v1->is_border() ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  k_v1                ->vertex() == aProfile.v1() ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  k_v1->next()        ->vertex() == aProfile.v0() ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  k_v1->next()->next()->vertex() == k ) ;
              }
              else // or is it the opposite?
              {
                edge_descriptor v1_k = k_v1->opposite();
                CGAL_SURF_SIMPL_TEST_assertion( !v1_k->is_border() ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  v1_k                ->vertex() == k ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  v1_k->next()        ->vertex() == aProfile.v0() ) ;
                CGAL_SURF_SIMPL_TEST_assertion(  v1_k->next()->next()->vertex() == aProfile.v1() ) ;
              }
            }
          );

          if ( !lIsFace )
          {
            CGAL_ECMS_TRACE(3,"  k=V" << k->id() << " IS NOT in a face with p-q. NON-COLLAPSABLE edge." ) ;
            rR = false ;
            break ;
          }  
          else 
          {
            CGAL_ECMS_TRACE(4,"  k=V" << k->id() << " is in a face with p-q") ;
          }
        }
      }  
    }
  }   
     
  if ( rR )
  {
    if ( aProfile.is_v0_v1_a_border() )
    {
      if ( Is_open_triangle(aProfile.v0_v1()) )
      {
        rR = false ;
        CGAL_ECMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge." ) ;
      }
    }
    else if ( aProfile.is_v1_v0_a_border() )
    {
      if ( Is_open_triangle(aProfile.v1_v0()) )
      {
        rR = false ;
        CGAL_ECMS_TRACE(3,"  p-q belongs to an open triangle. NON-COLLAPSABLE edge." ) ;
      }
    }
    else
    {
      if ( is_border(aProfile.v0()) && is_border(aProfile.v1()) )
      {
        rR = false ;
        CGAL_ECMS_TRACE(3,"  both p and q are boundary vertices but p-q is not. NON-COLLAPSABLE edge." ) ;
      }  
      else
      {
        bool lTetra = Is_tetrahedron(aProfile.v0_v1());

        CGAL_SURF_SIMPL_TEST_assertion( lTetra == mSurface.is_tetrahedron(aProfile.v0_v1()) ) ;

        if ( lTetra )
        {
          rR = false ;
          CGAL_ECMS_TRACE(3,"  p-q belongs to a tetrahedron. NON-COLLAPSABLE edge." ) ;
        }
      }
    }
  }
  
  return rR ;
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Is_tetrahedron( edge_descriptor const& h1 )
{
  //
  // Code copied from Polyhedron_3::is_tetrahedron()
  //
  edge_descriptor h2 = next_edge(h1,mSurface);
  edge_descriptor h3 = next_edge(h2,mSurface);
  
  edge_descriptor h1o = opposite_edge(h1,mSurface) ;
  edge_descriptor h2o = opposite_edge(h2,mSurface) ;
  edge_descriptor h3o = opposite_edge(h3,mSurface) ;
  
  edge_descriptor h4 = next_edge(h1o,mSurface);
  edge_descriptor h5 = next_edge(h2o,mSurface);
  edge_descriptor h6 = next_edge(h3o,mSurface);
  
  edge_descriptor h4o = opposite_edge(h4,mSurface) ;
  edge_descriptor h5o = opposite_edge(h5,mSurface) ;
  edge_descriptor h6o = opposite_edge(h6,mSurface) ;
  
  // check halfedge combinatorics.
  // at least three edges at vertices 1, 2, 3.
  if ( h4 == h3o ) return false;
  if ( h5 == h1o ) return false;
  if ( h6 == h2o ) return false;
  
  // exact three edges at vertices 1, 2, 3.
  if ( next_edge(h4o,mSurface) != h3o ) return false;
  if ( next_edge(h5o,mSurface) != h1o ) return false;
  if ( next_edge(h6o,mSurface) != h2o ) return false;
  
  // three edges at v4.
  if ( opposite_edge(next_edge(h4,mSurface),mSurface) != h5) return false;
  if ( opposite_edge(next_edge(h5,mSurface),mSurface) != h6) return false;
  if ( opposite_edge(next_edge(h6,mSurface),mSurface) != h4) return false;
  
  // All facets are triangles.
  if ( next_edge(next_edge(next_edge(h1,mSurface),mSurface),mSurface) != h1) return false;
  if ( next_edge(next_edge(next_edge(h4,mSurface),mSurface),mSurface) != h4) return false;
  if ( next_edge(next_edge(next_edge(h5,mSurface),mSurface),mSurface) != h5) return false;
  if ( next_edge(next_edge(next_edge(h6,mSurface),mSurface),mSurface) != h6) return false;
  
  // all edges are non-border edges.
  if ( is_border(h1)) return false;  // implies h2 and h3
  if ( is_border(h4)) return false;
  if ( is_border(h5)) return false;
  if ( is_border(h6)) return false;

  return true;
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Is_open_triangle( edge_descriptor const& h1 )
{
  bool rR = false ;
  
  edge_descriptor h2 = next_edge(h1,mSurface);
  edge_descriptor h3 = next_edge(h2,mSurface);
  
  // First check if it is a triangle 
  if ( next_edge(h3,mSurface) == h1 )
  {
    // Now check if it is open
    CGAL_ECMS_TRACE(4,"  p-q is a border edge... checking E" << h2->id() << " and E" << h3->id() ) ;
    
    rR = is_border(h2) && is_border(h3);  
  
    CGAL_SURF_SIMPL_TEST_assertion( rR == (  h1->is_border()
                                          && h1->next()->is_border()
                                          && h1->next()->next()->is_border()
                                          ) 
                                  ) ;
  }
  
  return rR ;
}

// Given triangles 'p0,p1,p2' and 'p0,p2,p3', both shared along edge 'v0-v2',
// determine if they are geometrically valid: that is, the ratio of their
// respective areas is no greater than a max value and the internal
// dihedral angle formed by their supporting planes is no greater than
// a given threshold
template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::are_shared_triangles_valid( Point const& p0, Point const& p1, Point const& p2, Point const& p3 ) const 
{
  bool rR = false ;
  
  Vector e01 = p1 - p0 ;
  Vector e02 = p2 - p0 ;
  Vector e03 = p3 - p0 ;
  
  Vector n012 = cross_product(e01,e02);
  Vector n023 = cross_product(e02,e03);
  
  FT l012 = n012 * n012 ;
  FT l023 = n023 * n023 ;
  
  FT larger  = (std::max)(l012,l023);
  FT smaller = (std::min)(l012,l023);
  
  const double cMaxAreaRatio = 1e8 ;
  
  CGAL_ECMS_TRACE(4,"    Testing validity of shared triangles:"
                  << "\n      p0=" << xyz_to_string(p0) << "\n      p1=" << xyz_to_string(p1) << "\n      p2=" << xyz_to_string(p2) << "\n      p3=" << xyz_to_string(p3)
                  << "\n      e01=" << xyz_to_string(e01) << "\n      e02=" << xyz_to_string(e02) << "\n      e03=" << xyz_to_string(e03)
                  << "\n      n012=" << xyz_to_string(n012) << "\n      n023=" << xyz_to_string(n023)
                  << "\n      l012=" << n_to_string(l012) << "\n      l023=" << n_to_string(l023)
                 );
  
  if ( larger < cMaxAreaRatio * smaller )
  {
    FT l0123 = n012 * n023 ;
    CGAL_ECMS_TRACE(4,"\n      l0123=" << n_to_string(l0123) );
    
    if ( CGAL_NTS is_positive(l0123) )
    {
      rR = true ;
    }
    else
    {
      CGAL_ECMS_TRACE(4,"\n      lhs: " << n_to_string(( l0123 * l0123 ) / ( l012 * l023 )) << " <= rhs: " << mcMaxDihedralAngleCos2 ) ;
      
      if ( ( l0123 * l0123 ) <= mcMaxDihedralAngleCos2 * ( l012 * l023 ) )
      {
        rR = true ;
      }
    }
  }
  
  return rR ;
}

// Returns the directed halfedge connecting v0 to v1, if exists.
template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
typename EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::edge_descriptor
EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::find_connection ( const_vertex_descriptor const& v0, const_vertex_descriptor const& v1 ) const
{
  out_edge_iterator eb, ee ; 
  for ( boost::tie(eb,ee) = out_edges(v0,mSurface) ; eb != ee ; ++ eb )
  {
    edge_descriptor out = *eb ;
    if ( target(out,mSurface) == v1 )
      return out ;
  }
      
  return edge_descriptor() ;  
}

// Given the edge 'e' around the link for the collapsinge edge "v0-v1", finds the vertex that makes a triangle adjacent to 'e' but exterior to the link (i.e not containing v0 nor v1)
// If 'e' is a mnull handle OR 'e' is a border edge, there is no such triangle and a null handle is returned.
template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
typename EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::vertex_descriptor
EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::find_exterior_link_triangle_3rd_vertex ( const_edge_descriptor const& e, const_vertex_descriptor const& v0, const_vertex_descriptor const& v1 ) const
{
  vertex_descriptor r ;
  
  if ( handle_assigned(e) )
  {
    vertex_descriptor ra = target(next_edge(e, mSurface), mSurface);
    vertex_descriptor rb = source(prev_edge(e, mSurface), mSurface);
    
    if ( ra == rb && ra != v0 && ra != v1 )
    {
      r = ra ;
    }
    else
    {
      ra = target(next_edge(opposite_edge(e,mSurface), mSurface), mSurface);
      rb = source(prev_edge(opposite_edge(e,mSurface), mSurface), mSurface);
      
      if ( ra == rb && ra != v0 && ra != v1 )
      {
        r = ra ;
      }
    }
  }
  
  return r ;
}


// A collase is geometrically valid if, in the resulting local mesh no two adjacent triangles form an internal dihedral angle
// greater than a fixed threshold (i.e. triangles do not "fold" into each other)
//
template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
bool EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Is_collapse_geometrically_valid( Profile const& aProfile, Placement_type k0)
{
  bool rR = true ;
  
  CGAL_ECMS_TRACE(3,"Testing geometrical collapsabilty of v0-v1=E" << aProfile.v0_v1()->id() );
  if ( k0 )
  {
    //
    // Use the current link to extract all local triangles incident to 'vx' in the collapsed mesh (which at this point doesn't exist yet)
    //
    typedef typename Profile::vertex_descriptor_vector::const_iterator const_link_iterator ;
    const_link_iterator linkb = aProfile.link().begin();
    const_link_iterator linke = aProfile.link().end  ();
    const_link_iterator linkl = cpp11::prev(linke) ;
    
    for ( const_link_iterator l = linkb ; l != linke && rR ; ++ l )
    {
      const_link_iterator pv = ( l == linkb ? linkl : cpp11::prev (l) );
      const_link_iterator nx = ( l == linkl ? linkb : cpp11::next  (l) ) ;
      
      // k0,k1 and k3 are three consecutive vertices along the link.
      vertex_descriptor k1 = *pv ;
      vertex_descriptor k2 = * l ;
      vertex_descriptor k3 = *nx ;
      
      CGAL_ECMS_TRACE(4,"  Screening link vertices k1=V" << k1->id() << " k2=V" << k2->id() << " k3=V" << k3->id() ) ;
      
      edge_descriptor e12 = find_connection(k1,k2);
      edge_descriptor e23 = k3 != k1 ? find_connection(k2,k3) : edge_descriptor() ;
      
      // If 'k1-k2-k3' are connected there will be two adjacent triangles 'k0,k1,k2' and 'k0,k2,k3' after the collapse.
      if ( handle_assigned(e12) && handle_assigned(e23) )
      {
        CGAL_ECMS_TRACE(4,"    Link triangles shared" ) ;
        
        if ( !are_shared_triangles_valid( *k0, get_point(k1), get_point(k2), get_point(k3) ) )
        {
          CGAL_ECMS_TRACE(3,"    Triangles VX-V" << k1->id() << "-V" << k2->id() << " and VX-V" << k3->id() << " are not geometrically valid. Collapse rejected");
          rR = false ;
        }
      }
      
      if ( rR )
      {
        // Also check the triangles 'k0,k1,k2' and it's adjacent along e12: 'k4,k2,k1', if exist
        vertex_descriptor k4 = find_exterior_link_triangle_3rd_vertex(e12,aProfile.v0(),aProfile.v1());
          
        // There is indeed a triangle shared along e12
        if ( handle_assigned(k4) )
        {
          CGAL_ECMS_TRACE(4,"    Found exterior link triangle shared along E" << e12->id() << " with third vertex: V" << k4->id() ) ;
          
          if ( !are_shared_triangles_valid( get_point(k1), get_point(k4), get_point(k2), *k0 ) )
          {
            CGAL_ECMS_TRACE(3,"    Triangles V" << k1->id() << "-V" << k4->id() << " and V" << k2->id() << "-VX are not geometrically valid. Collapse rejected");
            rR = false ;
          }
        }
      }
      
      if ( rR )
      {
        // And finally, check the triangles 'k0,k2,k3' and it's adjacent e23: 'k5,k3,k2' if exist
        vertex_descriptor k5 = find_exterior_link_triangle_3rd_vertex(e23,aProfile.v0(),aProfile.v1());
        
        // There is indeed a triangle shared along e12
        if ( handle_assigned(k5) )
        {
          CGAL_ECMS_TRACE(4,"    Found exterior link triangle shared along E" << e23->id() << " with third vertex: V" << k5->id() ) ;
          
          if ( !are_shared_triangles_valid( get_point(k2), get_point(k5), get_point(k3), *k0 ) )
          {
            CGAL_ECMS_TRACE(3,"    Triangles V" << k2->id() << "-V" << k5->id() << " and V" << k3->id() << "-VX are not geometrically valid. Collapse rejected");
            rR = false ;
          }
        }
      }
    } 
    
  }
  
  return rR ;
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
void EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Collapse( Profile const& aProfile, Placement_type aPlacement )
{

  CGAL_ECMS_TRACE(1,"S" << mStep << ". Collapsig " << edge_to_string(aProfile.v0_v1()) ) ;
  
  vertex_descriptor rResult ;
    
  Visitor.OnCollapsing(aProfile,aPlacement);

  -- mCurrentEdgeCount ;
      
  CGAL_SURF_SIMPL_TEST_assertion_code( size_type lResultingVertexCount = mSurface.size_of_vertices() ;
                                       size_type lResultingEdgeCount   = mSurface.size_of_halfedges() / 2 ;
                                     ) ; 

  // If the top/bottom facets exists, they are removed and the edges v0vt and Q-B along with them.
  // In that case their corresponding pairs must be pop off the queue
  
  if ( aProfile.left_face_exists() )
  {
    edge_descriptor lV0VL = primary_edge(aProfile.vL_v0());
    
    CGAL_ECMS_TRACE(3,"V0VL E" << lV0VL->id()
                   << "(V" <<  lV0VL->vertex()->id() << "->V" << lV0VL->opposite()->vertex()->id() << ")"
                   ) ;
                   
    Edge_data& lData = get_data(lV0VL) ;
    if ( lData.is_in_PQ() )
    {
      CGAL_ECMS_TRACE(2,"Removing E" << lV0VL->id() << " from PQ" ) ;
      remove_from_PQ(lV0VL,lData) ;
    }

    -- mCurrentEdgeCount ;
    CGAL_SURF_SIMPL_TEST_assertion_code( -- lResultingEdgeCount ) ; 
  }
  
  if ( aProfile.right_face_exists() )
  {
    edge_descriptor lVRV1 = primary_edge(aProfile.vR_v1());
    
    CGAL_ECMS_TRACE(3,"V1VRE" << lVRV1->id()
                   << "(V" <<  lVRV1->vertex()->id() << "->V" << lVRV1->opposite()->vertex()->id() << ")"
                   ) ;
                   
    Edge_data& lData = get_data(lVRV1) ;
    if ( lData.is_in_PQ() )
    {
      CGAL_ECMS_TRACE(2,"Removing E" << lVRV1->id() << " from PQ") ;
      remove_from_PQ(lVRV1,lData) ;
    }
    -- mCurrentEdgeCount ;
    CGAL_SURF_SIMPL_TEST_assertion_code( -- lResultingEdgeCount ) ; 
  }

  CGAL_ECMS_TRACE(1,"Removing:\n  v0v1: E" << aProfile.v0_v1()->id() << "(V" << aProfile.v0()->id() << "->V" << aProfile.v1()->id() << ")" );
  
    
  // Perform the actuall collapse.
  // This is an external function.
  // It's REQUIRED to remove ONLY 1 vertex (P or Q) and edges PQ,PT and QB
  // (PT and QB are removed if they are not null).
  // All other edges must be kept.
  // All directed edges incident to vertex removed are relink to the vertex kept.
  rResult = halfedge_collapse(aProfile.v0_v1(),mSurface);

  CGAL_SURF_SIMPL_TEST_assertion_code( -- lResultingEdgeCount ) ; 

  CGAL_SURF_SIMPL_TEST_assertion_code( -- lResultingVertexCount ) ; 

  CGAL_SURF_SIMPL_TEST_assertion( lResultingEdgeCount * 2 == mSurface.size_of_halfedges() ) ;

  CGAL_SURF_SIMPL_TEST_assertion( lResultingVertexCount == mSurface.size_of_vertices() ) ;

  CGAL_SURF_SIMPL_TEST_assertion( mSurface.is_valid() && mSurface.is_pure_triangle() ) ;
  
  CGAL_ECMS_TRACE(1,"V" << rResult->id() << " kept." ) ;
                 
#ifdef CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE
  out_edge_iterator eb1, ee1 ;     
  for ( boost::tie(eb1,ee1) = out_edges(rResult,mSurface) ; eb1 != ee1 ; ++ eb1 )  
    CGAL_ECMS_TRACE(2, edge_to_string(*eb1) ) ;
#endif

  if ( aPlacement )  
  {
    CGAL_ECMS_TRACE(1,"New vertex point: " << xyz_to_string(*aPlacement) ) ;
    put(vertex_point,mSurface,rResult,*aPlacement) ;
  }  

  Visitor.OnCollapsed(aProfile,rResult);

  Update_neighbors(rResult) ;
  
  CGAL_ECMS_DEBUG_CODE ( ++mStep ; )
}

template<class M,class SP, class VIM,class EIM,class EBM, class CF,class PF,class V>
void EdgeCollapse<M,SP,VIM,EIM,EBM,CF,PF,V>::Update_neighbors( vertex_descriptor const& aKeptV )
{
  CGAL_ECMS_TRACE(3,"Updating cost of neighboring edges..." ) ;

  //
  // (A) Collect all edges to update its cost: all those around each vertex adjacent to the vertex kept
  //  
  
  typedef std::set<edge_descriptor,Compare_id> edges ;
  
  edges lToUpdate(Compare_id(this)) ;
  
  // (A.1) Loop around all vertices adjacent to the vertex kept
  in_edge_iterator eb1, ee1 ; 
  for ( boost::tie(eb1,ee1) = in_edges(aKeptV,mSurface) ; eb1 != ee1 ; ++ eb1 )
  {
    edge_descriptor lEdge1 = *eb1 ;
    
    vertex_descriptor lAdj_k = source(lEdge1,mSurface);
    
    // (A.2) Loop around all edges incident on each adjacent vertex
    in_edge_iterator eb2, ee2 ; 
    for ( boost::tie(eb2,ee2) = in_edges(lAdj_k,mSurface) ; eb2 != ee2 ; ++ eb2 )
    {
      edge_descriptor lEdge2 = primary_edge(*eb2) ;
      
      Edge_data& lData2 = get_data(lEdge2);
      CGAL_ECMS_TRACE(4,"Inedge around V" << lAdj_k->id() << edge_to_string(lEdge2) ) ;
    
      // Only those edges still in the PQ _and_ not already collected are updated.
      if ( lData2.is_in_PQ() && lToUpdate.find(lEdge2) == lToUpdate.end() )
        lToUpdate.insert(lEdge2) ;
    } 
  }  
  
  //
  // (B) Proceed to update the costs.
  //
  
  for ( typename edges::iterator it = lToUpdate.begin(), eit = lToUpdate.end() ; it != eit ; ++ it )
  {
    edge_descriptor lEdge = *it;
    
    Edge_data& lData = get_data(lEdge);
    
    Profile const& lProfile = create_profile(lEdge);
    
    lData.cost() = get_cost(lProfile) ;
    
    CGAL_ECMS_TRACE(3, edge_to_string(lEdge) << " updated in the PQ") ;
    
    update_in_PQ(lEdge,lData);
  }
    
}


} // namespace Surface_mesh_simplification

} //namespace CGAL

#endif // CGAL_SURFACE_MESH_SIMPLIFICATION_DETAIL_EDGE_COLLAPSE_IMPL_H //
// EOF //
 
