#ifndef CGAL_INTERNAL_MEAN_CURVATURE_SKELETON_EDGE_COLLAPSE_H
#define CGAL_INTERNAL_MEAN_CURVATURE_SKELETON_EDGE_COLLAPSE_H


/**
\todo This function is coming from the branch BGL-redesign-GF and should be removed
      once that branch is integrated in master.

Let `v0` be the source and `v1` be the target vertices, and let `e` and `e'` be the halfedges of edge `v0v1`. 

For `e`, let `en` and `ep` be the next and previous 
halfedges, that is `en = next(e, g)`, `ep = prev(e, g)`, and let 
`eno` and `epo` be their opposite halfedges, that is 
`eno = opposite(en, g)` and `epo = opposite(ep, g)`.
Analoguously, for `e'` define  `en'`, `ep'`, `eno'`, and  `epo'`.

Then, after the collapse of edge `v0v1` the following holds for `e` (and analoguously for `e'`)

<UL> 
<LI>The edge `v0v1` is no longer in `g`. 
<LI>Either `v0`, or `v1` is no longer in `g` while the other remains. 
    Let `vgone` be the removed vertex and `vkept` be the remaining vertex. 
<LI>If `e` was a border halfedge, that is `get(bm, e, g) == true`, then `next(ep,g) == en`, and `prev(en,g) == ep`. 
<LI>If `e` was not a border halfedge, that is `get(bm, e, g) == false`, then `ep` and `epo` are no longer in `g` while `en` and `eno` are kept in `g`. 
<LI>For all halfedges `hv` in `halfedges_around_vertex(vgone, g)`, `target(*hv, g) == vkept` and `source(opposite(*hv, g), g) == vkept`. 
<LI>No other incidence information has changed in `g`. 
</UL> 

\returns vertex `vkept` (which can be either `v0` or `v1`). 

\pre For both halfedges of `v0v1` if it is not a border halfedge, the incident face must be triangular. 
\pre The edge `v0v1` must satisfy the *link condition* [\cite degn-tpec-98], which guarantees that the surface is also 2-manifold after the edge collapse. 

*/
template<typename Graph, typename BorderMap>
typename boost::graph_traits<Graph>::vertex_descriptor
collapse_edge(Graph& g,
              const BorderMap& bm,
              typename boost::graph_traits<Graph>::edge_descriptor v0v1)
{
  typedef boost::graph_traits< Graph > Traits;
  typedef typename Traits::vertex_descriptor          vertex_descriptor;
  typedef typename Traits::edge_descriptor            edge_descriptor;
  typedef typename Traits::halfedge_descriptor            halfedge_descriptor;

  halfedge_descriptor pq = halfedge(v0v1,g);
  halfedge_descriptor qp = opposite(pq, g);
  halfedge_descriptor pt = opposite(prev(pq, g), g);
  halfedge_descriptor qb = opposite(prev(qp, g), g);
  
  bool lTopFaceExists         = !bm[pq];
  bool lBottomFaceExists      = !bm[qp];
  bool lTopLeftFaceExists     = lTopFaceExists    && !bm[pt];
  bool lBottomRightFaceExists = lBottomFaceExists && !bm[qb];

  CGAL_precondition( !lTopFaceExists    || (lTopFaceExists    && ( degree(target(pt, g), g) > 2 ) ) ) ;
  CGAL_precondition( !lBottomFaceExists || (lBottomFaceExists && ( degree(target(qb, g), g) > 2 ) ) ) ;

  vertex_descriptor q = target(pq, g);
  vertex_descriptor p = source(pq, g);

  bool lP_Erased = false, lQ_Erased = false ;

  if ( lTopFaceExists )
  { 
    CGAL_precondition( !bm[opposite(pt, g)] ) ; // p-q-t is a face of the mesh
    if ( lTopLeftFaceExists )
    {
      //CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V" 
      //                << p.idx() << "->V" << target(pt, g).idx() 
      //                << ") by joining top-left face" ) ;

      join_face(g, pt);
    }
    else
    {
      //CGAL_ECMS_TRACE(3, "Removing p-t E" << pt.idx() << " (V" << p.idx() 
      //                << "->V" << target(pt, g).idx() << ") by erasing top face" ) ;

      remove_face(g, opposite(pt, g), bm);

      if ( !lBottomFaceExists )
      {
        //CGAL_ECMS_TRACE(3, "Bottom face doesn't exist so vertex P already removed" ) ;

        lP_Erased = true ;
      }  
    } 
  }

  if ( lBottomFaceExists )
  {   
    CGAL_precondition( !bm[opposite(qb, g)] ) ; // p-q-b is a face of the mesh
    if ( lBottomRightFaceExists )
    {
      //CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V" 
      //                << q.idx() << "->V" << target(qb, g).idx() 
      //                << ") by joining bottom-right face" ) ;

      join_face(g, qb);
    }
    else
    {
      //CGAL_ECMS_TRACE(3, "Removing q-b E" << qb.idx() << " (V" 
      //                << q.idx() << "->V" << target(qb, g).idx() 
      //                << ") by erasing bottom face" ) ;

      remove_face(g, opposite(qb, g), bm);

      if ( !lTopFaceExists )
      {
        //CGAL_ECMS_TRACE(3, "Top face doesn't exist so vertex Q already removed" ) ;
        lQ_Erased = true ;
      }  
    }
  }

  CGAL_assertion( !lP_Erased || !lQ_Erased ) ;

  if ( !lP_Erased && !lQ_Erased )
  {
    //CGAL_ECMS_TRACE(3, "Removing vertex P by joining pQ" ) ;

    join_vertex(g, pq);
    lP_Erased = true ;
  }    
  
  CGAL_assertion(g.is_valid());

  return lP_Erased ? q : p ;
}

#endif //CGAL_INTERNAL_MEAN_CURVATURE_SKELETON_EDGE_COLLAPSE_H
