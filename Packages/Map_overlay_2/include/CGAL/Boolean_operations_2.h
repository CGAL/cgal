// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.5-I-11 $
// release_date  : $CGAL_Date: 2002/08/04 $
//
// file          : include/CGAL/Boolean_operations_2.h
// package       : Map_overlay (1.12)
// maintainer    : Efi Fogel <efif@math.tau.ac.il>
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Eti Ezra          <estere@post.tau.ac.il>
//
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef BOOLEAN_OPERATIONS_2_H
#define BOOLEAN_OPERATIONS_2_H

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <fstream.h>
#include <vector>
#include <list>

#ifndef CGAL_MAP_OVERLAY_BASE_H
#include <CGAL/Map_overlay_base.h>
#endif

#ifndef CGAL_MAP_OVERLAY_H
#include <CGAL/Map_overlay.h>
#endif

//#ifndef CGAL_ARR_2_BOP_DCEL_H
//#include <CGAL/Bop_default_dcel.h>
//#endif

CGAL_BEGIN_NAMESPACE

template <class Map_overlay_>
class Boolean_operations_2
{
  typedef Boolean_operations_2<Map_overlay_>               Self;
public:
  typedef Map_overlay_                                     Map_overlay;
 
  // typedef Halfedges_output_container_            Halfedges_output_container;
  //typedef Vertices_output_container_              Vertices_output_container;

  typedef  typename Map_overlay::Subdivision              Subdivision;
  typedef  typename Map_overlay::Change_notification      Change_notification;
  typedef  typename Map_overlay::Map_overlay_base         Map_overlay_base;
  typedef  typename Map_overlay::Point_location_base      Point_location_base; 
  
  typedef  typename Subdivision::Vertex_iterator         Vertex_iterator;
  typedef  typename Subdivision::Vertex_const_iterator   Vertex_const_iterator;
  typedef  typename Subdivision::Halfedge_iterator       Halfedge_iterator;
  typedef  typename Subdivision::Halfedge_const_iterator  
                                                     Halfedge_const_iterator;
  typedef  typename Subdivision::Face_iterator           Face_iterator;
  typedef  typename Subdivision::Face_const_iterator     Face_const_iterator;
  typedef  typename Subdivision::Vertex_handle           Vertex_handle;
  typedef  typename Subdivision::Vertex_const_handle     Vertex_const_handle;
  typedef  typename Subdivision::Halfedge_handle         Halfedge_handle;
  typedef  typename Subdivision::Halfedge_const_handle   Halfedge_const_handle;
  typedef  typename Subdivision::Face_handle             Face_handle;
  typedef  typename Subdivision::Face_const_handle       Face_const_handle;
  
  typedef  std::list<Face_const_handle>                  Faces_container;
  typedef  std::list<Halfedge_const_handle>              Halfedges_container;
  typedef  std::list<Vertex_const_handle>                Vertices_container;

  Boolean_operations_2() {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2) : 
    //arr1_( *(new Map_overlay(arr1)) ), 
    //arr2_( *(new Map_overlay(arr2)) ),
    map_overlay_( (new Map_overlay(arr1, arr2)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf),
    //Map_overlay(arr2, &ovl_notf), &ovl_notf) )
  {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Map_overlay_base *ovl_ptr) : 
    //arr1_( *(new Map_overlay(arr1, ovl_ptr)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_ptr)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, ovl_ptr)) ) {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Change_notification* ovl_notf) : 
    //arr1_( *(new Map_overlay(arr1, ovl_notf)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_notf)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, ovl_notf)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf),
    //Map_overlay(arr2, &ovl_notf), &ovl_notf) )
  {}

  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Change_notification* ovl_notf_ptr, 
                       Map_overlay_base *ovl_ptr) : 
    //arr1_( *(new Map_overlay(arr1, ovl_notf, ovl_ptr)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_notf_ptr, ovl_ptr)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, ovl_notf_ptr, ovl_ptr)) )
  {}
  
  Boolean_operations_2(const Map_overlay& map_ovl) :  
    //arr1_(*(map_ovl.get_first_subdivision())), 
    //arr2_(*(map_ovl.get_second_subdivision())), 
    map_overlay_( new Map_overlay(map_ovl))
  {}

  // More constructors getting a Point location object as an argument.
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Point_location_base* pl_ptr) : 
    //arr1_( *(new Map_overlay(arr1)) ), 
    //arr2_( *(new Map_overlay(arr2)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, pl_ptr)) )
  {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Point_location_base* pl_ptr,
                       Map_overlay_base *ovl_ptr) : 
    //arr1_( *(new Map_overlay(arr1, ovl_ptr)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_ptr)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, pl_ptr, ovl_ptr)) ) {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Point_location_base* pl_ptr,
                       Change_notification* ovl_notf) : 
    //arr1_( *(new Map_overlay(arr1, ovl_notf)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_notf)) ),
    map_overlay_( (new Map_overlay(arr1, arr2, pl_ptr, ovl_notf)) )
  {}
  
  Boolean_operations_2(const Subdivision& arr1, 
                       const Subdivision& arr2, 
                       Point_location_base* pl_ptr,
                       Change_notification* ovl_notf_ptr, 
                       Map_overlay_base *ovl_ptr) : 
    //arr1_( *(new Map_overlay(arr1, ovl_notf, ovl_ptr)) ), 
    //arr2_( *(new Map_overlay(arr2, ovl_notf_ptr, ovl_ptr)) ),
    map_overlay_( (new Map_overlay(arr1,arr2,pl_ptr,ovl_notf_ptr,ovl_ptr)) )
  {}
  
  // Copy contructor
  Boolean_operations_2(const Self& bop) : 
    map_overlay_( new Map_overlay(*(bop.map_overlay_)))
  {} 
  
  virtual ~Boolean_operations_2() 
  {
    delete map_overlay_;
  }
  

  // -------------------- Assignement operator --------------------------
  const Self& operator=(const Self& bop)
  {
    map_overlay_ = new Map_overlay(*(bop.map_overlay_)); 

    return *this;
  }
  
  //template <class Faces_output_container, class
  //Halfedges_output_container, class Vertices_output_container >
  void intersection (Faces_container& list_of_faces,
                     Halfedges_container& list_of_halfedges, 
                     Vertices_container& list_of_vertices) const 
  {
    const Subdivision&              arr = map_overlay_->subdivision();
    //Change_notification  tmp_notf;
    
    // a vertex is in the intersection if it gots two vertices above
    // it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) && 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the intersection if it gots two halfedges
    // above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) && 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //(tmp_notf.get_second_face_above(face_iter)->bop()))
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void intersection (Faces_container& list_of_faces) const 
  {
    const Subdivision&              arr = map_overlay_->subdivision();
    
    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      
      if (face_is_below_first_map(face_iter) && 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void intersection (Halfedges_container& list_of_halfedges) const
  {
    const Subdivision&              arr = map_overlay_->subdivision();
    
    // an halfedge is in the intersection if it gots two halfedges
    // above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) && 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }

  void intersection (Vertices_container& list_of_vertices) const
  {
    const Subdivision&              arr = map_overlay_->subdivision();
    
    // a vertex is in the intersection if it gots two vertices above
    // it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) && 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
  }

  void intersection (Faces_container& list_of_faces,
                     Halfedges_container& list_of_halfedges) const
  {
    const Subdivision&         arr = map_overlay_->subdivision();
    
    // an halfedge is in the intersection if it gots two halfedges
    // above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) && 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //(tmp_notf.get_second_face_above(face_iter)->bop()))
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  void intersection (Halfedges_container& list_of_halfedges, 
                     Vertices_container& list_of_vertices)  const
  {
    const Subdivision&      arr = map_overlay_->subdivision();
    
    // a vertex is in the intersection if it gots two vertices above
    // it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) && 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the intersection if it gots two halfedges
    // above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) && 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }
  
  void intersection (Faces_container& list_of_faces,
                     Vertices_container& list_of_vertices) const
  {
    const Subdivision&         arr = map_overlay_->subdivision();
    
    // a vertex is in the intersection if it gots two vertices above
    // it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) && 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    } 

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //(tmp_notf.get_second_face_above(face_iter)->bop()))
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  //template <class Faces_output_container, class
  //Halfedges_output_container, class Vertices_output_container>
  void Union (Faces_container& list_of_faces, 
              Halfedges_container& list_of_halfedges, 
              Vertices_container& list_of_vertices)  const 
  {
    const Subdivision& arr = map_overlay_->subdivision();

    // a vertex is on the union if it gots at least one vertex above
    // it, or one halfedg or one face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) || 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the union if it gots at least one halfedge
    // above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) || 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() ||
      //tmp_notf.get_second_face_above(face_iter)->bop() )
      //list_of_faces.push_back(face_iter);
      
      if (face_is_below_first_map(face_iter) || 
          face_is_below_second_map(face_iter)){
        list_of_faces.push_back(face_iter);
      }
    }
  }
  
  void Union (Faces_container& list_of_faces) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() ||
      //tmp_notf.get_second_face_above(face_iter)->bop() )
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  //template <class Faces_output_container, class
  //Halfedges_output_container, class Vertices_output_container>
  void Union (Halfedges_container& list_of_halfedges) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    // an halfedge is in the union if it gots at least one halfedge
    // above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) || 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }

  void Union (Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    // a vertex is on the union if it gots at least one vertex above
    // it, or one halfedg or one face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) || 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
  }

  void Union (Faces_container& list_of_faces, 
              Halfedges_container& list_of_halfedges) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    // an halfedge is in the union if it gots at least one halfedge
    // above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) || 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() ||
      //tmp_notf.get_second_face_above(face_iter)->bop() )
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  void Union (Halfedges_container& list_of_halfedges, 
              Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    // a vertex is on the union if it gots at least one vertex above
    // it, or one halfedg or one face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) || 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the union if it gots at least one halfedge
    // above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (halfedge_is_below_first_map(halfedge_iter) || 
          halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }
  
  void Union (Faces_container& list_of_faces, 
              Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    // a vertex is on the union if it gots at least one vertex above
    // it, or one halfedg or one face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (vertex_is_below_first_map(vertices_iter) || 
          vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() ||
      //tmp_notf.get_second_face_above(face_iter)->bop() )
      //list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || 
          face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void symmetric_difference (Faces_container& list_of_faces, 
                             Halfedges_container& list_of_halfedges, 
                             Vertices_container& list_of_vertices) const 
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if ((vertex_is_below_first_map(vertices_iter) && 
           !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && 
           vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if ((halfedge_is_below_first_map(halfedge_iter) && 
           !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && 
           halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() &&
      //!(tmp_notf.get_second_face_above(face_iter)->bop()) ) ||
      //(!(tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //tmp_notf.get_second_face_above(face_iter)->bop() ))
      //list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && 
           !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && 
           face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }
  
  void symmetric_difference (Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if ((vertex_is_below_first_map(vertices_iter) && 
           !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && 
           vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
  }
  
  void symmetric_difference (Halfedges_container& list_of_halfedges) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if ((halfedge_is_below_first_map(halfedge_iter) && 
           !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && 
           halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }
  }

  void symmetric_difference (Faces_container& list_of_faces) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() &&
      //!(tmp_notf.get_second_face_above(face_iter)->bop()) ) ||
      //(!(tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //tmp_notf.get_second_face_above(face_iter)->bop() ))
      //list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && 
           !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && 
           face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }

  void symmetric_difference (Halfedges_container& list_of_halfedges, 
                             Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if ((vertex_is_below_first_map(vertices_iter) && 
           !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && 
           vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if ((halfedge_is_below_first_map(halfedge_iter) && 
           !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && 
           halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }
  }

   void symmetric_difference (Faces_container& list_of_faces, 
                              Halfedges_container& list_of_halfedges) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if ((halfedge_is_below_first_map(halfedge_iter) && 
           !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && 
           halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() &&
      //!(tmp_notf.get_second_face_above(face_iter)->bop()) ) ||
      //(!(tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //tmp_notf.get_second_face_above(face_iter)->bop() ))
      //list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && 
           !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && 
           face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }

  void symmetric_difference (Faces_container& list_of_faces, 
                             Vertices_container& list_of_vertices) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if ((vertex_is_below_first_map(vertices_iter) && 
           !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && 
           vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() &&
      //!(tmp_notf.get_second_face_above(face_iter)->bop()) ) ||
      //(!(tmp_notf.get_first_face_above(face_iter)->bop()) &&
      //tmp_notf.get_second_face_above(face_iter)->bop() ))
      //list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && 
           !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && 
           face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }
  
  //template <class Faces_output_container, class
  //Halfedges_output_container, class Vertices_output_container>
  void difference (Faces_container& list_of_faces, 
                   Halfedges_container& list_of_halfedges, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && 
            !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && 
            vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && 
            !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && 
            halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() &&
        //!(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && 
            !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) &&
        //tmp_notf.get_second_face_above(face_iter)->bop())
        //list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && 
            face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

  void difference (Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
  

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && 
            !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && 
            vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
  }

   void difference (Halfedges_container& list_of_halfedges, 
                    bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && 
            !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && 
            halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }
  }
  
  void difference (Faces_container& list_of_faces, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() &&
        //!(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && 
            !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) &&
        //tmp_notf.get_second_face_above(face_iter)->bop())
        //list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && 
            face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

  void difference (Halfedges_container& list_of_halfedges, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && 
            !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && 
            vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && 
            !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && 
            halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }
  }

   void difference (Faces_container& list_of_faces, 
                   Halfedges_container& list_of_halfedges, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();
    
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); 
         halfedge_iter !=  arr.halfedges_end(); ++halfedge_iter){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && 
            !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && 
            halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() &&
        //!(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && 
            !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) &&
        //tmp_notf.get_second_face_above(face_iter)->bop())
        //list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && 
            face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

   void difference (Faces_container& list_of_faces, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    const Subdivision& arr = map_overlay_->subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); 
         vertices_iter !=  arr.vertices_end(); ++vertices_iter){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && 
            !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && 
            vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Face_const_iterator face_iter = arr.faces_begin(); 
         face_iter !=  arr.faces_end(); ++face_iter){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() &&
        //!(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && 
            !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) &&
        //tmp_notf.get_second_face_above(face_iter)->bop())
        //list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && 
            face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }
  
  const Map_overlay&  map_overlay() const 
  {
    return *map_overlay_;
  }
  
  
private:
  bool  vertex_is_below_first_map(Vertex_const_handle v) const
  {
    const Change_notification* notf = map_overlay_->change_notification();
      
    return ((notf->get_first_vertex_above(v) != v && 
             notf->get_first_vertex_above(v)->bop())
            ||
            (notf->get_first_halfedge_above(v) != v->incident_halfedges() && 
             notf->get_first_halfedge_above(v)->bop())
            || 
            (notf->get_first_face_above(v) != v->incident_halfedges()->face() 
             && 
             notf->get_first_face_above(v)->bop()) );
  }
  
  bool  vertex_is_below_second_map(Vertex_const_handle v) const
  {
    const Change_notification* notf = map_overlay_->change_notification();
    
    return ((notf->get_second_vertex_above(v) != v &&  
             notf->get_second_vertex_above(v)->bop())
            ||
            (notf->get_second_halfedge_above(v) != v->incident_halfedges() && 
             notf->get_second_halfedge_above(v)->bop())
            || 
            (notf->get_second_face_above(v) != v->incident_halfedges()->face() 
             && 
             notf->get_second_face_above(v)->bop()) );
  }
  
  bool  halfedge_is_below_first_map(Halfedge_const_handle h) const
  {
    const Change_notification* notf = map_overlay_->change_notification();
    
    return ((notf->get_first_halfedge_above(h) != h && 
             notf->get_first_halfedge_above(h)->bop())
            || 
            (notf->get_first_face_above(h) != h->face() &&  
             notf->get_first_face_above(h)->bop()) );
  }
  
  bool  halfedge_is_below_second_map(Halfedge_const_handle h) const
  {
    const Change_notification* notf = map_overlay_->change_notification();
    
    return ((notf->get_second_halfedge_above(h) != h && 
             notf->get_second_halfedge_above(h)->bop())
            || 
            (notf->get_second_face_above(h) != h->face() &&  
             notf->get_second_face_above(h)->bop()) );
  }

  bool  face_is_below_first_map(Face_const_handle f) const
  {
    //Change_notification  notf;
    const Change_notification* notf = map_overlay_->change_notification();
    
    return (notf->get_first_face_above(f) != f && 
            notf->get_first_face_above(f)->bop());
  }
  
  bool  face_is_below_second_map(Face_const_handle f) const
  {
    //Change_notification  notf;
    const Change_notification* notf = map_overlay_->change_notification();
    
    return (notf->get_second_face_above(f) != f && 
            notf->get_second_face_above(f)->bop());
  }
  
  //const Map_overlay &arr1_, &arr2_;
  Map_overlay  *map_overlay_;
};

CGAL_END_NAMESPACE

#endif








