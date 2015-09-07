// Copyright (c) 2013 Technical University Braunschweig (Germany).
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
// Author(s):  Michael Hemmer <michael.hemmer@cgal.org>
//             

#ifndef CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2_H
#define CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2_H

#include <CGAL/Arrangement_2.h>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/assertions.h>

namespace CGAL {

template<class Arrangement_2_ , class RegularizationCategory = CGAL::Tag_true >
class Triangular_expansion_visibility_2 {
  typedef typename Arrangement_2_::Geometry_traits_2    Geometry_traits_2;
  typedef typename Geometry_traits_2::Kernel            K;

  typedef Triangular_expansion_visibility_2<
    Arrangement_2_, RegularizationCategory>             Self;

public:
  typedef Arrangement_2_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2              Traits_2;
  typedef typename Arrangement_2::Halfedge              Halfedge;
  typedef typename Arrangement_2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement_2::Halfedge_handle       Halfedge_handle;
  typedef typename Arrangement_2::Edge_const_iterator   Edge_const_iterator;
  typedef typename Arrangement_2::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;
  typedef typename Arrangement_2::Ccb_halfedge_circulator
    Ccb_halfedge_circulator;
  typedef typename Arrangement_2::Face_const_handle     Face_const_handle;
  typedef typename Arrangement_2::Face_handle           Face_handle;
  typedef typename Arrangement_2::Vertex_const_handle   Vertex_const_handle;
  typedef typename Arrangement_2::Vertex_handle         Vertex_handle;

  typedef typename K::Point_2                           Point_2;
  typedef typename Geometry_traits_2::Ray_2             Ray_2;
  typedef typename Geometry_traits_2::Segment_2         Segment_2;
  typedef typename Geometry_traits_2::Line_2            Line_2;
  typedef typename Geometry_traits_2::Vector_2          Vector_2;
  typedef typename Geometry_traits_2::Direction_2       Direction_2;
  typedef typename Geometry_traits_2::FT                Number_type;
  typedef typename Geometry_traits_2::Object_2          Object_2;

  typedef RegularizationCategory                       Regularization_category;
  
  typedef CGAL::Tag_true                      Supports_general_polygon_category;
  typedef CGAL::Tag_true                      Supports_simple_polygon_category;

private:
  typedef CGAL::Triangulation_vertex_base_2<K>                     Vb;
  typedef CGAL::Constrained_triangulation_face_base_2<K>           Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb>              TDS;
  typedef CGAL::No_intersection_tag                                Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;

  typedef std::pair<Point_2,Point_2>                               Constraint;

  // Functor to create edge constraints for the CDT out of Halfedges
  struct Make_constraint
  {
      typedef Constraint                                            result_type;

      Constraint operator()(const Halfedge& edge) const {
          return std::make_pair(edge.source()->point(),
                                edge.target()->point());
      }
  };

  // Observer to track any changes of the attached arrangement.
  class Observer : public Arr_observer<Arrangement_2>
  {
      
      typedef Arr_observer<Arrangement_2>                           Base;
      typedef Observer                                              Self;


  public:
      bool has_changed;

      Observer() : Base(), has_changed(false)
      {}

      Observer(const Arrangement_2& arr)
          : Base(const_cast<Arrangement_2&>(arr)), has_changed(false)
      {}

      // Arr_observer interface

      void after_attach() { has_changed = false; }


      void after_global_change() { has_changed = true; }
      void after_create_vertex(Vertex_handle) { has_changed = true; }
      void after_create_boundary_vertex(Vertex_handle) { has_changed = true; }
      void after_create_edge(Halfedge_handle) { has_changed = true; }
      void after_modify_vertex(Vertex_handle) { has_changed = true; }
      void after_modify_edge(Halfedge_handle) { has_changed = true; }
      void after_split_edge(Halfedge_handle, Halfedge_handle) {
          has_changed = true; }
      void after_split_fictitious_edge(Halfedge_handle, Halfedge_handle) {
          has_changed = true; }
      void after_split_face(Face_handle, Face_handle, bool) {
          has_changed = true; }
      void after_split_outer_ccb(Face_handle, Ccb_halfedge_circulator,
                                 Ccb_halfedge_circulator) {
          has_changed = true; }
      void after_split_inner_ccb(Face_handle, Ccb_halfedge_circulator,
                                 Ccb_halfedge_circulator) {
          has_changed = true; }
      void after_add_outer_ccb(Ccb_halfedge_circulator) { has_changed = true; }
      void after_add_inner_ccb(Ccb_halfedge_circulator) { has_changed = true; }
      void after_add_isolated_vertex(Vertex_handle) { has_changed = true; }
      void after_merge_edge(Halfedge_handle) { has_changed = true; }
      void after_merge_fictitious_edge(Halfedge_handle) { has_changed = true; }
      void after_merge_face(Face_handle) { has_changed = true; }
      void after_merge_outer_ccb(Face_handle, Ccb_halfedge_circulator) {
          has_changed = true; }
      void after_merge_inner_ccb(Face_handle, Ccb_halfedge_circulator) {
          has_changed = true; }
      void after_move_outer_ccb(Ccb_halfedge_circulator) { has_changed = true; }
      void after_move_inner_ccb(Ccb_halfedge_circulator) { has_changed = true; }
      void after_move_isolated_vertex(Vertex_handle) { has_changed = true; }
      void after_remove_vertex() { has_changed = true; }
      void after_remove_edge() { has_changed = true; }
      void after_remove_outer_ccb(Face_handle) { has_changed = true; }
      void after_remove_inner_ccb(Face_handle) { has_changed = true; }
  };
  

private:
  const Arrangement_2* p_arr;

  // May change during visibility computation
  mutable Observer observer;
  mutable boost::shared_ptr<CDT> p_cdt;
  mutable std::vector<Segment_2> needles;

  // Copy constructor and assignment not supported
  Triangular_expansion_visibility_2(const Self&);
  Self& operator= (const Self& );


public: 
  Triangular_expansion_visibility_2() : p_arr(NULL){}

  /*! Constructor given an arrangement. */
  Triangular_expansion_visibility_2 (const Arrangement_2& arr)
    : p_arr(&arr), observer(arr)
  {
    init_cdt(); 
  }

  const std::string name() const { return std::string("T_visibility_2"); }


  bool is_attached() const {
    //std::cout << "is_attached" << std::endl;
    return (p_arr != NULL);
  }

  void attach(const Arrangement_2& arr) {
    if(p_arr != &arr){
      p_arr = &arr;
      observer.detach();
      observer.attach(const_cast<Arrangement_2&>(arr));
      init_cdt(); 
    }   
    //std::cout << "attach done" << std::endl;
  }

  void detach() {
    //std::cout << "detach" << std::endl;
    observer.detach();
    p_arr = NULL; 
    p_cdt.reset();
  }

  const Arrangement_2& arrangement_2() const {
    return *p_arr;
  }


  template <typename VARR>
  typename VARR::Face_handle
  compute_visibility(const Point_2& q,
                     const Face_const_handle face,
                     VARR& out_arr )
  const {
    //std::cout << "query in face interior" << std::endl;

    if(observer.has_changed) {
        init_cdt();
    }

    out_arr.clear();
    needles.clear();
    CGAL_assertion(!face->is_unbounded());


    std::vector<Point_2> raw_output;
    typename CDT::Face_handle fh = p_cdt->locate(q);

    raw_output.push_back(fh->vertex(1)->point());
    if(!p_cdt->is_constrained(get_edge(fh,0))){
      //std::cout<< "edge 0 is not constrained" << std::endl;
      expand_edge(
          q,
          fh->vertex(2)->point(),
          fh->vertex(1)->point(),
          fh,0,std::back_inserter(raw_output));
    }

    raw_output.push_back(fh->vertex(2)->point());
    if(!p_cdt->is_constrained(get_edge(fh,1))){
      //std::cout << "edge 1 is not constrained" << std::endl;
      expand_edge(
          q,
          fh->vertex(0)->point(),
          fh->vertex(2)->point(),
          fh,1,std::back_inserter(raw_output));
    }

    raw_output.push_back(fh->vertex(0)->point());
    if(!p_cdt->is_constrained(get_edge(fh,2))){
      //std::cout << "edge 2 is not constrained" << std::endl;
      expand_edge(
          q,
          fh->vertex(1)->point(),
          fh->vertex(0)->point(),
          fh,2,std::back_inserter(raw_output));
    }


    return output(raw_output,out_arr);
  }

  template <typename VARR>
  typename VARR::Face_handle
  compute_visibility(const Point_2& q,
                     const Halfedge_const_handle he,
                     VARR& out_arr)
  const {
    //std::cout << "visibility_region he" << std::endl;

    if(observer.has_changed) {
        init_cdt();
    }

    CGAL_assertion(!he->face()->is_unbounded());
    out_arr.clear();
    needles.clear();

    std::vector<Point_2> raw_output;
    typename CDT::Locate_type location;
    int index;
    typename CDT::Face_handle fh = p_cdt->locate(q,location,index);
    CGAL_assertion(location == CDT::EDGE || location == CDT::VERTEX);
    //the following code tries to figure out which triangle one should start in.


    if(location == CDT::EDGE){
      //std::cout << "query on edge" << std::endl;
      // this is the easy part, there are only two possible faces
      // index indicates the edge = vertex on the other side of the edge
      // the next vertex in cw order should be the target of given edge
      if(fh->vertex(p_cdt->cw(index))->point() != he->target()->point()){
        //std::cout << "need to swap face" << std::endl;
        // take face on the other side if this is not the case
        typename CDT::Face_handle nfh = fh->neighbor(index);
        index = nfh->index(fh);
        fh = nfh;
      }
      CGAL_assertion(fh->vertex(p_cdt->cw(index))->point() == he->target()->point());
      CGAL_assertion(!p_cdt->is_infinite(fh->vertex(index)));


      // output the edge the query lies on
      raw_output.push_back(he->source()->point());
      raw_output.push_back(he->target()->point());

      if(!p_cdt->is_constrained(get_edge(fh,p_cdt->ccw(index)))){
        expand_edge(
            q,
            fh->vertex(index)->point(), //left
            he->target()->point()        , //right
            fh,
            p_cdt->ccw(index),
            std::back_inserter(raw_output));
      }
      raw_output.push_back(fh->vertex(index)->point());

      if(!p_cdt->is_constrained(get_edge(fh,p_cdt->cw(index)))){
        expand_edge(
            q,
            he->source()->point()        , //left
            fh->vertex(index)->point(), //right
            fh,
            p_cdt->cw(index),
            std::back_inserter(raw_output));
      }
    }

    if(location == CDT::VERTEX){
      //std::cout << "query on vertex" << std::endl;

      //bool query_point_on_vertex_is_not_working_yet = false;
      //CGAL_assertion(query_point_on_vertex_is_not_working_yet);

      CGAL_assertion(q  ==  he->target()->point());
      CGAL_assertion(fh->vertex(index)->point() ==  he->target()->point());

      // push points that are seen anyway
      // raw_output.push_back(he->source()->point()); inserted last
      raw_output.push_back(he->target()->point());
      raw_output.push_back(he->next()->target()->point());

      // now start in the triangle that contains he->next()
      while(
            p_cdt->is_infinite(fh->vertex(p_cdt->ccw(index))) ||
            he->next()->target()->point() !=
            fh->vertex(p_cdt->ccw(index))->point()
            )
      {
        typename CDT::Face_handle nfh = fh->neighbor(p_cdt->ccw(index));
        int nindex = nfh->index(fh);
        index = p_cdt->ccw(nindex);
        fh = nfh;
        CGAL_assertion(he->target()->point() == fh->vertex(index)->point());
      }


      CGAL_assertion(he->next()->source()->point() == fh->vertex(index)->point());
      CGAL_assertion(he->next()->target()->point() ==
             fh->vertex(p_cdt->ccw(index))->point());
      CGAL_assertion(!p_cdt->is_infinite(fh));
      CGAL_assertion(p_cdt->is_constrained(get_edge(fh,p_cdt->cw(index))));

      while(he->source()->point() != fh->vertex(p_cdt->ccw(index))->point()){

        if(!p_cdt->is_constrained(get_edge(fh,index))){
          expand_edge(
              q,
              fh->vertex(p_cdt-> cw(index))->point(), //left
              fh->vertex(p_cdt->ccw(index))->point(), //right
              fh,
              index,
              std::back_inserter(raw_output));
        }
        // push left end point of edge into output
        raw_output.push_back(fh->vertex(p_cdt-> cw(index))->point());

        // take the next triangle around q in ccw order
        typename CDT::Face_handle nfh = fh->neighbor(p_cdt->ccw(index));
        int nindex = nfh->index(fh);
        index = p_cdt->ccw(nindex);
        fh = nfh;
        CGAL_assertion(fh->vertex(index)->point() ==  he->target()->point());
      }
    }
    return output(raw_output,out_arr);
  }



private:

  typename CDT::Edge get_edge(typename CDT::Face_handle fh, int i) const {
    return std::make_pair(fh,i);
  }

  Point_2 ray_seg_intersection(
      const Point_2& q, const Point_2& b, // the ray 
      const Point_2& s, const Point_2& t  // the segment
    ) const {

    Ray_2 ray(q,b);
    Segment_2 seg(s,t);
    CGAL_assertion(typename K::Do_intersect_2()(ray,seg));
    CGAL::Object obj = typename K::Intersect_2()(ray,seg); 
    Point_2 result =  object_cast<Point_2>(obj);
    return result; 
  }

  void collect_needle(
      const Point_2& q, 
      const typename CDT::Vertex_handle vh, 
      const typename CDT::Face_handle fh, 
      int index)
  const {

    // the expanded edge should not be constrained 
    CGAL_assertion(!p_cdt->is_constrained(get_edge(fh,index)));
    CGAL_assertion(!p_cdt->is_infinite(fh));
    // go into the new face  
    const typename CDT::Face_handle nfh(fh->neighbor(index)); 
    CGAL_assertion(!p_cdt->is_infinite(nfh));

    // get indices of neighbors 
    int nindex = nfh->index(fh); // index of new vertex and old face 
    int rindex = p_cdt->ccw(nindex); // index of face behind right edge 
    int lindex = p_cdt-> cw(nindex); // index of face behind left edge 
    
    // get vertices seen from entering edge 
    const typename CDT::Vertex_handle nvh(nfh->vertex(nindex));
    const typename CDT::Vertex_handle rvh(nfh->vertex(p_cdt->cw (nindex)));
    const typename CDT::Vertex_handle lvh(nfh->vertex(p_cdt->ccw(nindex)));
    CGAL_assertion(!p_cdt->is_infinite(nvh));
    CGAL_assertion(!p_cdt->is_infinite(lvh));
    CGAL_assertion(!p_cdt->is_infinite(rvh));
    
    // get edges seen from entering edge 
    typename CDT::Edge re = get_edge(nfh,p_cdt->ccw(nindex));
    typename CDT::Edge le = get_edge(nfh,p_cdt-> cw(nindex));
     
    // do orientation computation once for new vertex 
    typename K::Orientation_2 orientation = 
      p_cdt->geom_traits().orientation_2_object();
    CGAL::Orientation orient = orientation(q,vh->point(),nvh->point());
    
    
    //std::cout << "\n collect_needle" <<std::endl;
    //std::cout << "q             "<< q << std::endl ;
    //std::cout << "vh->point()  "<<  vh->point() << std::endl;  
    //std::cout << "lvh->point()  "<< lvh->point() << std::endl ;
    //std::cout << "nvh->point()  "<< nvh->point() << std::endl ;
    //std::cout << "rvh->point()  "<< rvh->point() << std::endl<< std::endl;


    switch ( orient ) {
    case CGAL::COUNTERCLOCKWISE:
      // looking on to the right edge 
      if(p_cdt->is_constrained(re)) {
        if(vh != rvh) {
          Point_2 p = ray_seg_intersection(q, vh->point(),
                                           nvh->point(), rvh->point());
          //std::cout << vh->point() <<" -1- "<< p <<std::endl; 
          needles.push_back(Segment_2(vh->point(),p));
        }
      } else {
        collect_needle(q,vh,nfh,rindex);
      }
      break;
    case CGAL::CLOCKWISE:
      // looking on to the left edge 
      if(p_cdt->is_constrained(le)){
        if(vh != lvh){
          Point_2 p = ray_seg_intersection(q, vh->point(),
                                           nvh->point(), lvh->point());
          //std::cout << vh->point() <<" -2- "<< p <<std::endl; 
          needles.push_back(Segment_2(vh->point(),p));
        }
      } else {
        collect_needle(q,vh,nfh,lindex);
      }      
      break;
    default:
      CGAL_assertion(orient == CGAL::COLLINEAR);
      // looking on nvh, so it must be reported 
      // if it wasn't already (triangles rotate around vh)    
      if(vh != nvh){
        //std::cout << vh->point() <<" -3- "<< nvh->point() <<std::endl; 
        needles.push_back(Segment_2(vh->point(),nvh->point())); 
      }
      // but we may also contiue looking along the vertex 
      if(!p_cdt->is_constrained(re)) {
        collect_needle(q,nvh,nfh,rindex);
      }
      if(!p_cdt->is_constrained(le)) {
        collect_needle(q,nvh,nfh,lindex);
      }
      break;
    }
  }

  template<class OIT> 
  OIT expand_edge(
      const Point_2& q, 
      const Point_2& left, 
      const Point_2& right, 
      typename CDT::Face_handle fh, 
      int index, 
      OIT oit)
  const {

    // the expanded edge should not be constrained 
    CGAL_assertion(!p_cdt->is_constrained(get_edge(fh,index)));
    CGAL_assertion(!p_cdt->is_infinite(fh));
    
    // go into the new face  
    const typename CDT::Face_handle nfh(fh->neighbor(index)); 
    CGAL_assertion(!p_cdt->is_infinite(nfh));

    // get indices of neighbors 
    int nindex = nfh->index(fh); // index of new vertex and old face 
    int rindex = p_cdt->ccw(nindex); // index of face behind right edge 
    int lindex = p_cdt-> cw(nindex); // index of face behind left edge 
    
    // get vertices seen from entering edge 
    const typename CDT::Vertex_handle nvh(nfh->vertex(nindex));
    const typename CDT::Vertex_handle rvh(nfh->vertex(p_cdt->cw (nindex)));
    const typename CDT::Vertex_handle lvh(nfh->vertex(p_cdt->ccw(nindex)));
    CGAL_assertion(!p_cdt->is_infinite(nvh));
    CGAL_assertion(!p_cdt->is_infinite(lvh));
    CGAL_assertion(!p_cdt->is_infinite(rvh));
    
    // get edges seen from entering edge 
    typename CDT::Edge re = get_edge(nfh,p_cdt->ccw(nindex));
    typename CDT::Edge le = get_edge(nfh,p_cdt-> cw(nindex));
     
    // do orientation computation once for new vertex 
    typename K::Orientation_2 orientation = 
      p_cdt->geom_traits().orientation_2_object();
    CGAL::Orientation ro = orientation(q,right,nvh->point());
    CGAL::Orientation lo = orientation(q,left ,nvh->point());
    
    CGAL_assertion(typename K::Orientation_2()(q,left ,lvh->point())
           != CGAL::CLOCKWISE);
    CGAL_assertion(typename K::Orientation_2()(q,right,rvh->point())
           != CGAL::COUNTERCLOCKWISE);

    //std::cout << (ro == CGAL::COUNTERCLOCKWISE) << " " <<
    //(lo == CGAL::CLOCKWISE) << std::endl;
    
    //right edge is seen if new vertex is counter clockwise of right boarder 
    if(ro == CGAL::COUNTERCLOCKWISE){
      if(p_cdt->is_constrained(re)){
        // the edge is constrained
        // report intersection with right boarder ray 
        // if it is not already the right vertex (already reported)
        if(right != rvh->point()){
          *oit++ = ray_seg_intersection(q,right,nvh->point(),rvh->point());
        }
        
        // then report intersection with left boarder if it exists
        if(lo == CGAL::COUNTERCLOCKWISE){
          *oit++ = ray_seg_intersection(q,left,nvh->point(),rvh->point());
        }
      }else{
        // the edge is not a constrained 
        if(lo == CGAL::COUNTERCLOCKWISE){
          // no split needed and return 
          //std::cout<< "h1"<< std::endl;
          oit = expand_edge(q,left,right,nfh,rindex,oit);
          //std::cout<< "h1 done"<< std::endl;
          return oit;          
        }else{
          // spliting at new vertex 
          //std::cout<< "h2"<< std::endl;
          *oit++ = expand_edge(q,nvh->point(),right,nfh,rindex,oit);
          //std::cout<< "h2 done"<< std::endl;
        }        
      }
    }
    

    //std::cout << "q             "<< q << std::endl ;
    //std::cout << "lvh->point()  "<< lvh->point() << std::endl;  
    //std::cout << "left          "<< left << std::endl  ;
    //std::cout << "nvh->point()  "<< nvh->point() << std::endl ;
    //std::cout << "right         "<< right << std::endl ;
    //std::cout << "rvh->point()  "<< rvh->point() << std::endl<< std::endl;
    
    
    // determin whether new vertex needs to be reported 
    if(ro != CGAL::CLOCKWISE && lo != CGAL::COUNTERCLOCKWISE){
      *oit++ = nvh->point(); 
    }
    if(!Regularization_category::value){
      CGAL_assertion(!(ro == CGAL::COLLINEAR && lo == CGAL::COLLINEAR));
      // we have to check whether a needle starts here. 
      if(p_cdt->is_constrained(le) && !p_cdt->is_constrained(re)
              && ro == CGAL::COLLINEAR)
        collect_needle(q,nvh,nfh,rindex);

      if(p_cdt->is_constrained(re) && !p_cdt->is_constrained(le)
              && lo == CGAL::COLLINEAR)
        collect_needle(q,nvh,nfh,lindex);    
    }

    //left edge is seen if new vertex is clockwise of left boarder 
    if(lo == CGAL::CLOCKWISE){
      if(p_cdt->is_constrained(le)){
        // the edge is constrained
        // report interesection with right boarder if exists 
        if(ro == CGAL::CLOCKWISE){
          *oit++ = ray_seg_intersection(q,right,nvh->point(),lvh->point());
        }
        // then report intersection with left boarder ray 
        // if it is not already the left vertex (already reported)
        if(left != lvh->point()){
          *oit++ = ray_seg_intersection(q,left,nvh->point(),lvh->point());          
        }
        return oit; 
      }else{
        // the edge is not a constrained 
        if(ro == CGAL::CLOCKWISE){
          // no split needed and return
          //std::cout<< "h3"<< std::endl;
          oit = expand_edge(q,left,right,nfh,lindex,oit);
          //std::cout<< "h3 done"<< std::endl;
          return oit; 
        }else{
          // spliting at new vertex 
          //std::cout<< "h4"<< std::endl;
          oit = expand_edge(q,left,nvh->point(),nfh,lindex,oit);
          //std::cout<< "h4 done"<< std::endl;
          return oit;          
        }        
      }
    }
    
    return oit;
  }


  template <typename VARR> 
  typename VARR::Face_handle 
  output(std::vector<Point_2>& raw_output, VARR& out_arr) const {

    if(!needles.empty()){
      std::vector<Segment_2> segments(needles.begin(),needles.end()); 
      for(unsigned int i = 0; i < raw_output.size(); i++){
//      //std::cout <<  raw_output[i] << " -- " 
//                <<  raw_output[(i+1)%raw_output.size()] << std::endl; 
        segments.push_back(Segment_2(raw_output[i],
                                     raw_output[(i+1) % raw_output.size()]));
      }

      CGAL::insert_non_intersecting_curves(out_arr,
                                           segments.begin(),
                                           segments.end());

    } else {
      typename VARR::Vertex_handle v_last, v_first;
      v_last = v_first = 
        out_arr.insert_in_face_interior(raw_output[0],out_arr.unbounded_face());
      
      for(unsigned int i = 0; i < raw_output.size()-1; i++){
//      std::cout <<  raw_output[i] << " -- "
//                <<  raw_output[(i+1)%raw_output.size()] << std::endl;
        if(raw_output[i] < raw_output[(i+1)]){
          v_last = out_arr.insert_from_left_vertex (
                             Segment_2(raw_output[i], raw_output[i+1]), v_last
                           )->target();
        } else {
          v_last = out_arr.insert_from_right_vertex(
                             Segment_2(raw_output[i], raw_output[i+1]), v_last
                           )->target();
        }        
      }
      out_arr.insert_at_vertices(
                Segment_2(raw_output.front(), raw_output.back()),
                v_last, v_first
              );
    }

    CGAL_assertion(out_arr.number_of_faces() == 2);

    if(out_arr.faces_begin()->is_unbounded())
      return ++out_arr.faces_begin();
    else
      return out_arr.faces_begin();
  }

  void init_cdt() const {
    //std::cout<< "==============" <<std::endl;
    //std::cout<< "Input Polygon:" <<std::endl;

    typedef typename boost::transform_iterator<Make_constraint,
                                               Edge_const_iterator>        Iter;

    Iter begin = boost::make_transform_iterator(p_arr->edges_begin(),
                                                Make_constraint());

    Iter end = boost::make_transform_iterator(p_arr->edges_end(),
                                              Make_constraint());

    //std::cout << "init_cdt new CDT" << std::endl;
    p_cdt = boost::shared_ptr<CDT>(new CDT(begin, end));
    observer.has_changed = false;
    //std::cout << "init_cdt done" << std::endl;
    //std::cout << std::endl;
  }
};

} // namespace CGAL

#endif // CGAL_TRIANGULAR_EXPANSION_VISIBILITY_2_H
