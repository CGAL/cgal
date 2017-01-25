// Copyright (c) 2011 GeometryFactory (France).
// All rights reserved.
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
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_COPLANAR_H
#define CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_COPLANAR_H

#include <CGAL/license/Polygon_mesh_processing.h>



namespace CGAL{
  namespace internal_IOP{

//enum Intersection_type {FACET,EDGE,VERTEX,EMPTY,COPLNR};

template<class Polyhedron,class Is_const>
struct Intersection_point_coplanar{
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle Halfedge_handle;
  
  Intersection_type type;
  Halfedge_handle info;
  Halfedge_handle segment_vertex;
  
  Intersection_point_coplanar(){};
  
  Intersection_point_coplanar(Intersection_type type_,
                              Halfedge_handle info_,
                              Halfedge_handle segment_vertex_
  )     :type(type_),info(info_),segment_vertex(segment_vertex_){}
};
  
//test q vs segment bc
template<class Inter_pt_coplanar,class Halfedge_handle>
void point_vs_segment(Inter_pt_coplanar& pt,
                      Halfedge_handle bc, Halfedge_handle ab,
                      const Orientation& pqc, const Orientation& pqb)
{
  if (pqb==COLLINEAR){
    pt.type=VERTEX;
    pt.info=ab;
  }
  else{
    if (pqc==COLLINEAR){
      pt.type=VERTEX;
      pt.info=bc;
    }
    else{
      pt.type=EDGE;
      pt.info=bc;
    }
  }
}


/*
//      c
//     / \
// ___/___\____   with q-------p
//   /     \
//  a_______b
//supporting line of segment pq intersects ac and bc
*/

template<class Inter_pt_coplanar,class Point_3,class Halfedge_handle>
std::pair<Inter_pt_coplanar,Inter_pt_coplanar>
decision_tree(const Point_3* a,const Point_3* b,const Point_3* c,
              const Point_3* p,const Point_3* q,
              const Orientation& pqa,const Orientation& pqb,const Orientation& pqc,
              Halfedge_handle pq,
              Halfedge_handle ca,Halfedge_handle ab,Halfedge_handle bc)
{
  CGAL_precondition(pqa!=NEGATIVE);
  CGAL_precondition(pqb!=NEGATIVE);
  CGAL_precondition(pqc!=POSITIVE);
  
  Inter_pt_coplanar pt1(EMPTY,NULL,NULL),pt2(EMPTY,NULL,NULL);
  
  //the segment supporting line intersects ac and bc
  const Orientation bcq = coplanar_orientation(*b,*c,*q);
  switch(bcq){
    case NEGATIVE: //segment does not intersect the triangle
      return std::make_pair(pt1,pt2);
    case COLLINEAR: //q \in [bc]
      point_vs_segment(pt1,bc,ab,pqc,pqb);
      pt1.segment_vertex=pq;
      return std::make_pair(pt1,pt2);
    default:
      break;
  }

  //===> q is to the left of bc
  
  const Orientation cap = coplanar_orientation(*c,*a,*p);
  switch(cap){
    case NEGATIVE: //segment does not intersect the triangle
      return false; 
    case COLLINEAR: //p \in [ca]
      point_vs_segment(pt1,ca,bc,pqa,pqc);
      pt1.segment_vertex=pq->opposite();
      return std::make_pair(pt1,pt2);
    default:
      break;
  }
  //===> p is to the right of ac
  
   
  const Orientation cqa = coplanar_orientation(*c,*q,*a);
  const Orientation cbp = coplanar_orientation(*c,*b,*p);
  
  if (pqc==COLLINEAR && cqa==POSITIVE && cbp==POSITIVE){
    //special case when c is inside the segment pq
    pt1.type=VERTEX;
    pt1.info=bc;
    return std::make_pair(pt1,pt2);
  }
  
  //where is located q?
  switch (cqa){
    case POSITIVE: //q is outside the triangle
      point_vs_segment(pt1,ca,bc,pqa,pqc);
    break;
    case NEGATIVE: //q is inside the triangle
      pt1.type=FACET;
      pt1.info=pq;      
    break;
    case COLLINEAR: //q \in [ca]
      point_vs_segment(pt1,ca,bc,pqa,pqc);
      pt1.segment_vertex=pq;
    break;
  }
  
  //where is located p?
  switch (cbp){
    case POSITIVE: //p is outside the triangle
      point_vs_segment(pt2,bc,ab,pqc,pqb);
    break;
    case NEGATIVE: //p is inside the triangle
      pt2.type=FACET;
      pt2.info=pq->opposite();  
    break;
    case COLLINEAR: //p \in [bc]
      point_vs_segment(pt2,bc,ab,pqc,pqb);
      pt2.segment_vertex=pq->opposite();
    break;
  }
  return std::make_pair(pt1,pt2);
}


/*
//         c
//        / \
//       /   \       
//      /     \
// ____a_______b______  with q-------p
//supporting lines of segments pq and ab are the same.
*/

template<class Inter_pt_coplanar,class Point_3,class Halfedge_handle>
std::pair<Inter_pt_coplanar,Inter_pt_coplanar>
collinear_decision_tree(const Point_3* a,const Point_3* b,const Point_3* c,
                        const Point_3* p,const Point_3* q,
                        const Orientation& pqa,const Orientation& pqb,const Orientation& pqc,
                        Halfedge_handle pq,
                        Halfedge_handle ca,Halfedge_handle ab,Halfedge_handle bc)
{
  CGAL_precondition(pqa==COLLINEAR);
  CGAL_precondition(pqb==COLLINEAR);
  CGAL_precondition(pqc==NEGATIVE);
  
  Inter_pt_coplanar pt1(EMPTY,NULL,NULL),pt2(EMPTY,NULL,NULL);
  
  //the segment supporting line intersects ac and bc
  const Orientation bcq = coplanar_orientation(*b,*c,*q);
  switch(bcq){
    case NEGATIVE: //segment does not intersect the triangle
      return std::make_pair(pt1,pt2);
    case COLLINEAR: // q = b
      pt1.type=VERTEX;
      pt1.info=ab;
      pt1.segment_vertex=pq;
      return std::make_pair(pt1,pt2);
    default:
      break;
  }

  //===> q is to the left of b
  
  const Orientation cap = coplanar_orientation(*c,*a,*p);
  switch(cap){
    case NEGATIVE: //segment does not intersect the triangle
      return false; 
    case COLLINEAR: //p = a
      pt1.type=VERTEX;
      pt1.info=ca;
      pt1.segment_vertex=pq->opposite();
      return std::make_pair(pt1,pt2);
    default:
      break;
  }
  
  //===> p is to the right of a
  
   
  const Orientation cqa = coplanar_orientation(*c,*q,*a);
  const Orientation cbp = coplanar_orientation(*c,*b,*p);
 
  //where is located q?
  switch (cqa){
    case POSITIVE: //q is to the left of a
      pt1.type=VERTEX;
      pt1.info=ca;
    break;
    case NEGATIVE: //q is to the right of a
      pt1.type=EDGE;
      pt1.info=ab;      
      pt1.segment_vertex=pq;
    break;
    case COLLINEAR: //q = a
      pt1.type=VERTEX;
      pt1.info=ca;
      pt1.segment_vertex=pq;
    break;
  }
  
  //where is located p?
  switch (cbp){
    case POSITIVE: //p is to the right of b
    {
      pt2.type=VERTEX;
      pt2.info=ab;
    }
    break;
    case NEGATIVE: //p is to the left of b
    {
      pt2.type=EDGE;
      pt2.info=ab;
      pt2.segment_vertex=pq->opposite();
    }      
    break;
    case COLLINEAR: //p = b
      pt2.type=VERTEX;
      pt2.info=ab;
      pt2.segment_vertex=pq->opposite();
    break;
  }
  
  return std::make_pair(pt1,pt2);
}

//std::pair<Intersection_point_coplanar<Polyhedron,Is_const>,Intersection_point_coplanar<Polyhedron,Is_const> >
template<class Polyhedron,class Kernel,class Is_const>
std::pair<Intersection_point_coplanar<Polyhedron,Is_const>,Intersection_point_coplanar<Polyhedron,Is_const> >
do_intersect_coplanar(typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle pq,
                      typename Polyhedron_types<Polyhedron,Is_const>::Face_handle fh)
{
  typedef Intersection_point_coplanar<Polyhedron,Is_const> Inter_pt_coplanar;
  typedef std::pair<Inter_pt_coplanar,Inter_pt_coplanar> Return_type;
  
  
  typedef typename Polyhedron_types<Polyhedron,Is_const>::Halfedge_handle Halfedge_handle; 
  
  const typename Kernel::Point_3 & A = fh->halfedge()->vertex()->point();
  const typename Kernel::Point_3 & B = fh->halfedge()->next()->vertex()->point();
  const typename Kernel::Point_3 & C = fh->halfedge()->next()->next()->vertex()->point();
  const typename Kernel::Point_3 & P = pq->vertex()->point();
  const typename Kernel::Point_3 & Q = pq->opposite()->vertex()->point();

  const typename Kernel::Point_3 *  a = &A;
  const typename Kernel::Point_3 *  b = &B;
  const typename Kernel::Point_3 *  c = &C;
  
  const typename Kernel::Point_3*   p = &P;
  const typename Kernel::Point_3*   q = &Q;

  Halfedge_handle ca=fh->halfedge();
  Halfedge_handle ab=fh->halfedge()->next();
  Halfedge_handle bc=fh->halfedge()->next()->next();
  
  
  // Determine the orientation of the triangle in the common plane

  if (coplanar_orientation(A,B,C) != POSITIVE){
    // The triangle is not counterclockwise oriented swap two vertices.
    b = &C;
    c = &B;
    std::swap(bc,ab);
  }

  // Test whether the segment's supporting line intersects the
  // triangle in the common plane
  
  Orientation pqa = coplanar_orientation(*p,*q,*a);
  
  //ensure pqa >= 0
  if (pqa == NEGATIVE){
    std::swap(p,q);
    pqa=POSITIVE;
    pq=pq->opposite();
  }  
  
  const Orientation pqb = coplanar_orientation(*p,*q,*b);
  const Orientation pqc = coplanar_orientation(*p,*q,*c);

  
  
  //Handle first case pq collinear with a triangle edge
  if (pqa == COLLINEAR){
    if (pqb == COLLINEAR){
      //ab and pq are on the same line
      if (pqc == NEGATIVE)
        return collinear_decision_tree<Inter_pt_coplanar>
              (a,b,c,p,q,pqa,pqb,pqc,pq,ca,ab,bc);
      else
        return collinear_decision_tree<Inter_pt_coplanar>
              (a,b,c,q,p,pqa,pqb,NEGATIVE,pq->opposite(),ca,ab,bc);        
    }
    if (pqc == COLLINEAR){
      //ac and pq are on the same line
      if (pqb == NEGATIVE)
        return collinear_decision_tree<Inter_pt_coplanar>
              (c,a,b,p,q,pqc,pqa,pqb,pq,bc,ca,ab);
      else
        return collinear_decision_tree<Inter_pt_coplanar>
              (c,a,b,q,p,pqc,pqa,NEGATIVE,pq->opposite(),bc,ca,ab);        
    }
  }
  else
    if(pqb ==COLLINEAR && pqc == COLLINEAR){
      //bc and pq are on the same line
      return collinear_decision_tree<Inter_pt_coplanar>
            (b,c,a,q,p,pqb,pqc,NEGATIVE,pq->opposite(),ab,bc,ca);
    }
  
  
  CGAL_assertion(pqa!=NEGATIVE);

  switch ( pqa ) {
  case POSITIVE:
    switch ( pqb ) {
      case POSITIVE:
        if (pqc == POSITIVE)
          return std::make_pair(Inter_pt_coplanar(EMPTY,NULL,NULL),Inter_pt_coplanar(EMPTY,NULL,NULL));
        return decision_tree<Inter_pt_coplanar>(a,b,c,p,q,pqa,pqb,pqc,pq,ca,ab,bc);
        
      case COLLINEAR:
      case NEGATIVE:
        if (pqc == POSITIVE) // b is isolated on the negative side
          return decision_tree<Inter_pt_coplanar>(c,a,b,p,q,pqc,pqa,pqb,pq,bc,ca,ab);
        return decision_tree<Inter_pt_coplanar>
          (b,c,a,q,p,POSITIVE,opposite(pqc),NEGATIVE,pq->opposite(),ab,bc,ca);
      default:// should not happen.
        CGAL_assertion(false);
        return std::make_pair(Inter_pt_coplanar(EMPTY,NULL,NULL),Inter_pt_coplanar(EMPTY,NULL,NULL));       
    }
  case COLLINEAR:
    switch ( pqb ) {
      case POSITIVE:
        if (pqc == POSITIVE) // a is isolated on the negative side
          return decision_tree<Inter_pt_coplanar>(b,c,a,p,q,pqb,pqc,pqa,pq,ab,bc,ca);
        return decision_tree<Inter_pt_coplanar>(a,b,c,p,q,pqa,pqb,pqc,pq,ca,ab,bc);  
        
      case NEGATIVE:
        if (pqc == NEGATIVE) // a is isolated on the positive side
          return decision_tree<Inter_pt_coplanar>
            (b,c,a,q,p,POSITIVE,POSITIVE,pqa,pq->opposite(),ab,bc,ca);        
        return decision_tree<Inter_pt_coplanar>(c,a,b,p,q,pqc,pqa,pqb,pq,bc,ca,ab);
      
      default:// should not happen.
        CGAL_assertion(false);
        return std::make_pair(Inter_pt_coplanar(EMPTY,NULL,NULL),Inter_pt_coplanar(EMPTY,NULL,NULL));
    }
  default:// should not happen.
    CGAL_assertion(false);
    return std::make_pair(Inter_pt_coplanar(EMPTY,NULL,NULL),Inter_pt_coplanar(EMPTY,NULL,NULL));
    
  }
}

}}//namespace CGAL::internal_IOP

#endif //CGAL_INTERNAL_INTERSECTION_TRIANGLE_SEGMENT_3_COPLANAR_H
