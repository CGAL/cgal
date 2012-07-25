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

#ifndef CGAL_INTERNAL_INTERSECTION_COPLANAR_TRIANGLES_3_H
#define CGAL_INTERNAL_INTERSECTION_COPLANAR_TRIANGLES_3_H

#include <CGAL/internal/corefinement/intersection_triangle_segment_3.h> //for Intersection_type
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <bitset>

//TODO rename this file when doing proper integration
namespace CGAL{
namespace internal_IOP{


//intersection point of two coplanar triangles that keeps track of 
//the location of that point onto the triangles.
template <class IK,class Halfedge_handle>
struct Intersection_point_with_info
{
  typedef CGAL::Exact_predicates_exact_constructions_kernel     Exact_kernel;
  typedef IK                                                    Input_kernel;
  typedef CGAL::Cartesian_converter<Input_kernel,Exact_kernel>  Converter;
  
  Intersection_type type_1,type_2; //intersection type for 1st and 2nd facets
  Halfedge_handle info_1,info_2; //halfedge providing primitive indicated by type_1 and type_2
  typename Exact_kernel::Point_3 point; //the geometric embedding of the intersection
  
  //constructor from a vertex of first triangle initialized inside the second triangle
  Intersection_point_with_info(Halfedge_handle info1,Halfedge_handle info2):
    type_1(VERTEX),type_2(FACET),info_1(info1),info_2(info2)
  {
    Converter converter;
    point=converter(info_1->vertex()->point());
  }

  //constructor for intersection of edges. prev and curr are two points on an edge of the first facet (preserving the 
  //orientation of the facet). This edge is intersected by info2 from the second facet.
  //
  //The rational is the following: we first check whether curr and prev are on the same edge. I so we create
  //an intersection point between two edges. Otherwise, the point is a vertex of the second facet included into
  //the first facet.
  //
  //(V,F) : point initialy constructed
  //(V,E) : (V,F) updated by get_orientation_and_update_info_2 (i.e lies on one edge)
  //(V,V) : (V,E) updated by get_orientation_and_update_info_2 (i.e lies on two edges)
  //(E,E) : created in the following function when prev and curr lie on the same edge
  //(E,V) : (E,E) updated by get_orientation_and_update_info_2 (always done as lies on two edges)
  //(E,F) : impossible
  //(F,V) : detected when curr and prev and not on the same edge
  //(F,E) : impossible
  //(F,F) : impossible
  //
  Intersection_point_with_info(Intersection_point_with_info prev,Intersection_point_with_info curr,
                               Halfedge_handle info1,Halfedge_handle info2):
    type_2(EDGE),info_2(info2)
  {
    #ifdef CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
    std::cout << "prev: "; prev.print_debug();
    std::cout << "curr: "; curr.print_debug(); std::cout << std::endl;
    #endif //CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
    Converter converter;
    if (prev.type_1==VERTEX && prev.info_1->next() == curr.info_1){
      CGAL_assertion(curr.type_1!=FACET);
      type_1=EDGE;
      info_1=curr.info_1;
    }
    else{
      if(curr.type_1==VERTEX && prev.info_1 == curr.info_1){
        CGAL_assertion(prev.type_1!=FACET);
        type_1=EDGE;
        info_1=curr.info_1;        
      }
      else{
        if (curr.type_1==EDGE && prev.type_1==EDGE &&  curr.info_1==prev.info_1){
          type_1=EDGE;
          info_1=curr.info_1;
        }
        else{
          //curr and prev are not on the same edge of the first facet.
          //The intersection point to be computed is a VERTEX of the second facet
          type_1=FACET;
          info_1=info1;
          type_2=VERTEX;
          
          //this is used to select the correct endpoint of the edge of the second facet
          typename Exact_kernel::Collinear_3 is_collinear = Exact_kernel().collinear_3_object();
          if ( !is_collinear(prev.point,curr.point,converter(info_2->vertex()->point()) ) ){
            info_2=info_2->next()->next();
            CGAL_assertion( is_collinear(prev.point,curr.point,converter(info_2->vertex()->point()) ) );
          }
          point = converter( info_2->vertex()->point() );
          return;
        }
      }
    }
    
    //handle degenerate case when two edges overlap
    //at least one of the two vertex has already been found as a vertex of a facet. Here we set it for the second point
    if(prev.type_2!=FACET && curr.type_2!=FACET && (prev.type_1==VERTEX || prev.type_2==VERTEX) && (curr.type_1==VERTEX || curr.type_2==VERTEX)){
        typename Exact_kernel::Collinear_3 is_collinear = Exact_kernel().collinear_3_object();
        if ( is_collinear(prev.point,curr.point,converter(info_2->opposite()->vertex()->point()) ) ){
          info_2=info_2->next()->next();
          type_2=VERTEX;
          point = converter( info_2->vertex()->point() );
          return;          
        }
        if ( is_collinear(prev.point,curr.point,converter(info_2->vertex()->point()) ) ){
          type_2=VERTEX;
          point = converter( info_2->vertex()->point() );
          return;
        }
    }

    //handle regular intersection of two edges
    typename Exact_kernel::Construct_line_3 line_3=Exact_kernel().construct_line_3_object();
    typename Exact_kernel::Line_3 l1=
      line_3(converter(info_2->vertex()->point()),converter(info_2->opposite()->vertex()->point()));
    typename Exact_kernel::Line_3 l2=line_3(prev.point,curr.point);
    CGAL::Object res=Exact_kernel().intersect_3_object()(l1,l2);
    const typename Exact_kernel::Point_3* ptptr=CGAL::object_cast<typename Exact_kernel::Point_3>(&res);
    CGAL_assertion(ptptr!=NULL);
    point=*ptptr;
  }
  
  void print_debug() const{
    std::cout << " (";
    if (type_1==VERTEX) std::cout << "V";
    if (type_1==EDGE)   std::cout << "E";
    if (type_1==FACET) std::cout << "F";
    if (type_1==EMPTY) std::cout << "?";
    std::cout << "-" << &(*info_1);
    std::cout << ";";
    if (type_2==VERTEX) std::cout << "V";
    if (type_2==EDGE)   std::cout << "E";
    if (type_2==FACET)   std::cout << "F";
    if (type_2==EMPTY) std::cout << "?";
    std::cout << ")" << "[" << CGAL::to_double(point.x()) << "," << CGAL::to_double(point.y()) << "," << CGAL::to_double(point.z()) << "]";
  }
  
  int debug_unique_type_int() const{
    int res=0;
    switch (type_1){
      case VERTEX: res+=1; break;
      case EDGE: res+=3; break;
      case FACET: res+=7; break;
      default: break;
    } 
    switch (type_2){
      case VERTEX: res+=1; break;
      case EDGE: res+=3; break;
      case FACET: res+=7; break;
      default: break;
    }
    return res;
  }
  
  bool is_valid(Intersection_type type,Halfedge_handle info){
    bool valid=true;
    Converter converter;
    switch (type){
      case VERTEX: valid&= converter(info->vertex()->point())==point; break;
      case EDGE:{
        typename Exact_kernel::Segment_3 seg=
          Exact_kernel().construct_segment_3_object()( converter(info->vertex()->point()),
                                                       converter(info->opposite()->vertex()->point()) );
        valid&= Exact_kernel().has_on_3_object()(seg,point);
      }
      break;
      case FACET:{
        typename Exact_kernel::Coplanar_orientation_3 orient=Exact_kernel().coplanar_orientation_3_object();
        typename Exact_kernel::Point_3 p=converter(info->vertex()->point());
        typename Exact_kernel::Point_3 q=converter(info->next()->vertex()->point());
        typename Exact_kernel::Point_3 r=converter(info->opposite()->vertex()->point());
        valid &= orient(p,q,r,point)==POSITIVE;
        valid &= orient(q,r,p,point)==POSITIVE;
        valid &= orient(r,p,q,point)==POSITIVE;
      }
      break;
      default: valid=false;
    } 
    return valid;
  }
  
  bool is_valid(){
    return is_valid(type_1,info_1) && is_valid(type_2,info_2);
  }
  
};
  
template<class Halfedge_handle,class Inter_pt>  
CGAL::Orientation get_orientation_and_update_info_2(Halfedge_handle h,Inter_pt& p)
{
  typename Inter_pt::Exact_kernel::Coplanar_orientation_3 orient=
    typename Inter_pt::Exact_kernel().coplanar_orientation_3_object();
  typename Inter_pt::Converter converter;
  
  CGAL::Orientation res = orient(converter(h->opposite()->vertex()->point()),
                                 converter(h->vertex()->point()),
                                 converter(h->next()->vertex()->point()),
                                 p.point);
  
  if ( (p.type_1==VERTEX || p.type_1==EDGE) && res==COLLINEAR){
    if (p.type_2==FACET){ //detect a case (VERTEX,EDGE)
      p.type_2=EDGE;
      p.info_2=h;
    }
    else{
      //detect a case (VERTEX,VERTEX) or (EDGE,VERTEX)
      CGAL_assertion(p.type_2==EDGE);
      p.type_2=VERTEX;
      if (p.info_2->next()!=h){
        CGAL_assertion(h->next()==p.info_2);
        p.info_2=h;
      }
    }
  }
  
  return res;
}
  
template <class Facet_handle,class Inter_pt>
void intersection_coplanar_facets_cutoff(Facet_handle f,std::list<Inter_pt>& inter_pts,Facet_handle other)
{
  #ifdef CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
  std::cout << "cutoff: " << f->opposite()->vertex()->point() << " " << f->vertex()->point() << std::endl;
  #endif //CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
  if ( inter_pts.empty() ) return;
  typedef typename std::list<Inter_pt>::iterator Iterator;
 
  std::map<Inter_pt*,Orientation> orientations;
  for (Iterator it=inter_pts.begin();it!=inter_pts.end();++it){
    #ifdef CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
    it->print_debug();
    #endif //CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
    orientations[ &(*it) ]=get_orientation_and_update_info_2(f,*it);
  }
  #ifdef CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
  std::cout << std::endl;
  #endif //CGAL_DEBUG_COPLANAR_TRIANGLE_INTERSECTION
  
  int pt_added=0;
  
  Inter_pt* prev = &(*boost::prior(inter_pts.end()));
  bool inter_pts_size_g_2 = inter_pts.size() > 2;
  Iterator stop = inter_pts_size_g_2 ? inter_pts.end() : boost::prior(inter_pts.end());
  for (Iterator it=inter_pts.begin();it!=stop;++it)
  {
    Inter_pt* curr=&(*it);
    if (!inter_pts_size_g_2) std::swap(prev,curr);
    Orientation or_prev=orientations[prev],or_curr=orientations[curr];
    if ( (or_prev==POSITIVE && or_curr==NEGATIVE) || (or_prev==NEGATIVE && or_curr==POSITIVE) )
    {
      Iterator it_curr = inter_pts_size_g_2 ? it:boost::next(it);
      prev=&(* inter_pts.insert( it_curr,Inter_pt(*prev,*curr,other,f) ) );
      orientations[prev]=COLLINEAR;
      ++pt_added;
    }
    prev=&(*it);    
  }
  
  CGAL_kernel_assertion(pt_added<3);
  Iterator it=inter_pts.begin();
  std::size_t nb_interpt=inter_pts.size();
  //this boolean allows to reverse order of intersection points in case there were 3 remaining intersection points
  //and the point in the middle was removed. In that case the order must be reversed to preserve the orientations
  //of the last edge:
  // A---X---B  --> AB  to be consistent with the other cases this should be BA!
  // X---B---A  --> BA
  // B---A---X  --> BA
  //
  bool should_revert_list=false;
  
  while(it!=inter_pts.end())
  {
    if (orientations[&(*it)]==NEGATIVE){
      inter_pts.erase(it++);
      if (--nb_interpt == 2 && it!=inter_pts.end() && boost::next(it)==inter_pts.end()) should_revert_list=true; 
    }
    else
      ++it;
  }
  if (should_revert_list && nb_interpt==2) inter_pts.reverse();
}
  
template <class Input_kernel,class Halfedge_handle>
void
intersection_coplanar_facets(
  Halfedge_handle f1,
  Halfedge_handle f2,
  std::list<Intersection_point_with_info<Input_kernel,Halfedge_handle> >& output )
{
  typedef Intersection_point_with_info<Input_kernel,Halfedge_handle> Inter_pt;
  output.push_back( Inter_pt(f1,f2) );
  output.push_back( Inter_pt(f1->next(),f2) );
  output.push_back( Inter_pt(f1->next()->next(),f2) );

  //intersect f2 with the three half planes which intersection defines f1
  intersection_coplanar_facets_cutoff(f2,output,f1);
  intersection_coplanar_facets_cutoff(f2->next(),output,f1);
  intersection_coplanar_facets_cutoff(f2->next()->next(),output,f1);
}



} } //namespace CGAL::internal_IOP


#endif //CGAL_INTERNAL_INTERSECTION_COPLANAR_TRIANGLES_3_H

