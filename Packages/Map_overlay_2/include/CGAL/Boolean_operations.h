#ifndef BOOLEAN_OPERATIONS_H
#define BOOLEAN_OPERATIONS_H

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

#ifndef CGAL_ARR_2_BOP_DCEL_H
#include <CGAL/Bop_default_dcel.h>
#endif

//#ifndef CGAL_MAP_OVERLAY_MISC_H
//#include <CGAL/Map_overlay_default_notifier.h>
//#endif

//#ifndef CGAL_ARRANGEMENT_2_H
//#include <CGAL/Arrangement_2.h>
//#endif

CGAL_BEGIN_NAMESPACE

template <class Map_overlay_>
class Boolean_operations
{
public:
  typedef Map_overlay_                                             Map_overlay;
 
  // typedef Halfedges_output_container_                              Halfedges_output_container;
  //typedef Vertices_output_container_                               Vertices_output_container;

  typedef typename  Map_overlay::Arrangement                             Arrangement;
  typedef typename  Map_overlay::Map_overlay_change_notification         Map_overlay_change_notification;
  typedef Map_overlay_base<Arrangement, Map_overlay_change_notification>  Map_ovl_base;

  typedef  typename Arrangement::Vertex_iterator                    Vertex_iterator;
  typedef  typename Arrangement::Vertex_const_iterator              Vertex_const_iterator;
  typedef  typename Arrangement::Halfedge_iterator                  Halfedge_iterator;
  typedef  typename Arrangement::Halfedge_const_iterator            Halfedge_const_iterator;
  typedef  typename Arrangement::Face_iterator                      Face_iterator;
  typedef  typename Arrangement::Face_const_iterator                Face_const_iterator;
  typedef  typename Arrangement::Vertex_handle                      Vertex_handle;
  typedef  typename Arrangement::Vertex_const_handle                Vertex_const_handle;
  typedef  typename Arrangement::Halfedge_handle                    Halfedge_handle;
  typedef  typename Arrangement::Halfedge_const_handle              Halfedge_const_handle;
  typedef  typename Arrangement::Face_handle                        Face_handle;
  typedef  typename Arrangement::Face_const_handle                  Face_const_handle;

  typedef  std::list<Face_const_handle>                            Faces_container;
  typedef  std::list<Halfedge_const_handle>                        Halfedges_container;
  typedef  std::list<Vertex_const_handle>                          Vertices_container;

  Boolean_operations() {}

  Boolean_operations(const Arrangement& arr1, 
                     const Arrangement& arr2) : 
    arr1_( *(new Map_overlay(arr1)) ), 
    arr2_( *(new Map_overlay(arr2)) ),
    map_overlay( *(new Map_overlay(arr1_, arr2_)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf), Map_overlay(arr2, &ovl_notf), &ovl_notf) ) 
  {}

  Boolean_operations(const Arrangement& arr1, 
                     const Arrangement& arr2, 
                     Map_ovl_base *ovl_ptr) : 
    arr1_( *(new Map_overlay(arr1, ovl_ptr)) ), 
    arr2_( *(new Map_overlay(arr2, ovl_ptr)) ),
    map_overlay( *(new Map_overlay(arr1_, arr2_, ovl_ptr)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf), Map_overlay(arr2, &ovl_notf), &ovl_notf) ) 
  {}

  Boolean_operations(const Arrangement& arr1, 
                     const Arrangement& arr2, 
                     Map_overlay_change_notification* ovl_notf_ptr, Map_ovl_base *ovl_ptr) : 
    arr1_( *(new Map_overlay(arr1, ovl_notf, ovl_ptr)) ), 
    arr2_( *(new Map_overlay(arr2, ovl_notf_ptr, ovl_ptr)) ),
    map_overlay(*(new Map_overlay(arr1_, arr2_, ovl_notf_ptr, ovl_ptr)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf), Map_overlay(arr2, &ovl_notf), &ovl_notf) ) 
  {}

  Boolean_operations(const Arrangement& arr1, const Arrangement& arr2, Map_overlay_change_notification* ovl_notf) : 
    arr1_( *(new Map_overlay(arr1, ovl_notf)) ), arr2_( *(new Map_overlay(arr2, ovl_notf)) ),
    map_overlay( *(new Map_overlay(arr1_, arr2_, ovl_notf)) )
    //map_overlay( Map_overlay(Map_overlay(arr1, &ovl_notf), Map_overlay(arr2, &ovl_notf), &ovl_notf) ) 
  {}
  
  Boolean_operations(const Map_overlay& map_ovl) :  
    arr1_(*(map_ovl.get_first_subdivision())), arr2_(*(map_ovl.get_second_subdivision())), 
    map_overlay(map_ovl)
  {}

  
  //template <class Faces_output_container, class Halfedges_output_container, class Vertices_output_container >
  void intersection (Faces_container& list_of_faces,
                     Halfedges_container& list_of_halfedges, 
                     Vertices_container& list_of_vertices) const 
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // a vertex is in the intersection if it gots two vertices above it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the intersection if it gots two halfedges above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) && (tmp_notf.get_second_face_above(face_iter)->bop())) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void intersection (Faces_container& list_of_faces) const 
  {
    const Arrangement&              arr = map_overlay.get_arr();
    //Map_overlay_change_notification  tmp_notf;
    
    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) && (tmp_notf.get_second_face_above(face_iter)->bop())) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void intersection (Halfedges_container& list_of_halfedges) const
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // an halfedge is in the intersection if it gots two halfedges above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }

  void intersection (Vertices_container& list_of_vertices) const
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // a vertex is in the intersection if it gots two vertices above it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
  }

  void intersection (Faces_container& list_of_faces,
                     Halfedges_container& list_of_halfedges) const
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // an halfedge is in the intersection if it gots two halfedges above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) && (tmp_notf.get_second_face_above(face_iter)->bop())) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  void intersection (Halfedges_container& list_of_halfedges, 
                     Vertices_container& list_of_vertices)  const
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // a vertex is in the intersection if it gots two vertices above it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the intersection if it gots two halfedges above it or two faces.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }
  
  void intersection (Faces_container& list_of_faces,
                     Vertices_container& list_of_vertices) const
  {
    const Arrangement&              arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // a vertex is in the intersection if it gots two vertices above it or two halfedges or two face.
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    } 

    // a face is in the intersection if it gots two faces above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ((tmp_notf.get_first_face_above(face_iter)->bop()) && (tmp_notf.get_second_face_above(face_iter)->bop())) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  //template <class Faces_output_container, class Halfedges_output_container, class Vertices_output_container>
  void Union (Faces_container& list_of_faces, 
              Halfedges_container& list_of_halfedges, 
              Vertices_container& list_of_vertices)  const 
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;

    // a vertex is on the union if it gots at least one vertex above it, or one halfedg or one face. 
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) || vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the union if it gots at least one halfedge above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) || halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() || tmp_notf.get_second_face_above(face_iter)->bop() ) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void Union (Faces_container& list_of_faces) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_change_notification  tmp_notf;
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() || tmp_notf.get_second_face_above(face_iter)->bop() ) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  //template <class Faces_output_container, class Halfedges_output_container, class Vertices_output_container>
  void Union (Halfedges_container& list_of_halfedges) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    // an halfedge is in the union if it gots at least one halfedge above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) || halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }

  void Union (Vertices_container& list_of_vertices) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    // a vertex is on the union if it gots at least one vertex above it, or one halfedg or one face. 
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) || vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
  }

  void Union (Faces_container& list_of_faces, 
              Halfedges_container& list_of_halfedges) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;
    
    // an halfedge is in the union if it gots at least one halfedge above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) || halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() || tmp_notf.get_second_face_above(face_iter)->bop() ) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }

  void Union (Halfedges_container& list_of_halfedges, 
              Vertices_container& list_of_vertices) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    // a vertex is on the union if it gots at least one vertex above it, or one halfedg or one face. 
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) || vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // an halfedge is in the union if it gots at least one halfedge above it or one face.
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (halfedge_is_below_first_map(halfedge_iter) || halfedge_is_below_second_map(halfedge_iter))
        list_of_halfedges.push_back(halfedge_iter);
    }
  }
  
  void Union (Faces_container& list_of_faces, 
              Vertices_container& list_of_vertices) const
  {
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    // a vertex is on the union if it gots at least one vertex above it, or one halfedg or one face. 
    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (vertex_is_below_first_map(vertices_iter) || vertex_is_below_second_map(vertices_iter))
        list_of_vertices.push_back(vertices_iter);
    }
    
    // a face is in the union if it gots at least one face above it.
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( tmp_notf.get_first_face_above(face_iter)->bop() || tmp_notf.get_second_face_above(face_iter)->bop() ) 
      //  list_of_faces.push_back(face_iter);
      if (face_is_below_first_map(face_iter) || face_is_below_second_map(face_iter))
        list_of_faces.push_back(face_iter);
    }
  }
  
  void symmetric_difference (Faces_container& list_of_faces, 
                             Halfedges_container& list_of_halfedges, 
                             Vertices_container& list_of_vertices) const 
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if ((vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if ((halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) ) 
      //     || (!(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop() )) 
      // list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }
  
  void symmetric_difference (Vertices_container& list_of_vertices) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if ((vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
  }
  
  void symmetric_difference (Halfedges_container& list_of_halfedges) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if ((halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }
  }

  void symmetric_difference (Faces_container& list_of_faces) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) ) 
      //     || (!(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop() )) 
      // list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }

  void symmetric_difference (Halfedges_container& list_of_halfedges, 
                             Vertices_container& list_of_vertices) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if ((vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if ((halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }
  }

   void symmetric_difference (Faces_container& list_of_faces, 
                              Halfedges_container& list_of_halfedges) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if ((halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter)) || 
          (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter)) )
        list_of_halfedges.push_back(halfedge_iter);
      
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) ) 
      //     || (!(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop() )) 
      // list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }

  void symmetric_difference (Faces_container& list_of_faces, 
                             Vertices_container& list_of_vertices) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if ((vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter)) || 
          (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter)) )
        list_of_vertices.push_back(vertices_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      //if ( (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) ) 
      //     || (!(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop() )) 
      // list_of_faces.push_back(face_iter);

      if ((face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter)) || 
          (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter)) )
          list_of_faces.push_back(face_iter);
    }
  }
  
  //template <class Faces_output_container, class Halfedges_output_container, class Vertices_output_container>
  void difference (Faces_container& list_of_faces, 
                   Halfedges_container& list_of_halfedges, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //  list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop()) 
        //  list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

  void difference (Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
  }

   void difference (Halfedges_container& list_of_halfedges, 
                    bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }
  }

  void difference (Faces_container& list_of_faces, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //  list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop()) 
        //  list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

  void difference (Halfedges_container& list_of_halfedges, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }
  }

   void difference (Faces_container& list_of_faces, 
                   Halfedges_container& list_of_halfedges, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    
    
    for (Halfedge_const_iterator halfedge_iter = arr.halfedges_begin(); halfedge_iter !=  arr.halfedges_end(); halfedge_iter++){
      if (first){
        if (halfedge_is_below_first_map(halfedge_iter) && !halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
      }
      else
        if (!halfedge_is_below_first_map(halfedge_iter) && halfedge_is_below_second_map(halfedge_iter))
          list_of_halfedges.push_back(halfedge_iter);
    }

    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //  list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop()) 
        //  list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }

   void difference (Faces_container& list_of_faces, 
                   Vertices_container& list_of_vertices, 
                   bool first = true) const
  {
    //Map_overlay tmp_ovl;
    const Arrangement& arr = map_overlay.get_subdivision();
    //Map_overlay_changeNotification  tmp_notf;

    for (Vertex_const_iterator vertices_iter = arr.vertices_begin(); vertices_iter !=  arr.vertices_end(); vertices_iter++){
      if (first){
        if (vertex_is_below_first_map(vertices_iter) && !vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
      }
      else  // second minus first. 
        if (!vertex_is_below_first_map(vertices_iter) && vertex_is_below_second_map(vertices_iter))
          list_of_vertices.push_back(vertices_iter);
    }
    
    for (Face_const_iterator face_iter = arr.faces_begin(); face_iter !=  arr.faces_end(); face_iter++){
      if (first){
        //if (tmp_notf.get_first_face_above(face_iter)->bop() && !(tmp_notf.get_second_face_above(face_iter)->bop()) )
        //  list_of_faces.push_back(face_iter);
        
        if (face_is_below_first_map(face_iter) && !face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
      }
      else
        //if ( !(tmp_notf.get_first_face_above(face_iter)->bop()) && tmp_notf.get_second_face_above(face_iter)->bop()) 
        //  list_of_faces.push_back(face_iter);
        if (!face_is_below_first_map(face_iter) && face_is_below_second_map(face_iter))
          list_of_faces.push_back(face_iter);
    }
  }
  
  const Map_overlay&  get_map_overlay() const {return map_overlay;}
  
  ~Boolean_operations () {}
  
private:
  bool  vertex_is_below_first_map(Vertex_const_handle v) const
  {
    Map_overlay_change_notification  tmp_notf;
    return ((tmp_notf.get_first_vertex_above(v) != v &&  tmp_notf.get_first_vertex_above(v)->bop())
            ||
            (tmp_notf.get_first_halfedge_above(v) != v->incident_halfedges() && tmp_notf.get_first_halfedge_above(v)->bop())
            || 
            (tmp_notf.get_first_face_above(v) != v->incident_halfedges()->face() &&  tmp_notf.get_first_face_above(v)->bop()) );
  }
  
  bool  vertex_is_below_second_map(Vertex_const_handle v) const
  {
    Map_overlay_change_notification  tmp_notf;
    return ((tmp_notf.get_second_vertex_above(v) != v &&  tmp_notf.get_second_vertex_above(v)->bop())
            ||
            (tmp_notf.get_second_halfedge_above(v) != v->incident_halfedges() && tmp_notf.get_second_halfedge_above(v)->bop())
            || 
            (tmp_notf.get_second_face_above(v) != v->incident_halfedges()->face() && tmp_notf.get_second_face_above(v)->bop()) );
  }
  
  bool  halfedge_is_below_first_map(Halfedge_const_handle h) const
  {
    Map_overlay_change_notification  tmp_notf;
    return ((tmp_notf.get_first_halfedge_above(h) != h && tmp_notf.get_first_halfedge_above(h)->bop())
            || 
            (tmp_notf.get_first_face_above(h) != h->face() &&  tmp_notf.get_first_face_above(h)->bop()) );
  }
  
  bool  halfedge_is_below_second_map(Halfedge_const_handle h) const
  {
    Map_overlay_change_notification  tmp_notf;
    return ((tmp_notf.get_second_halfedge_above(h) != h && tmp_notf.get_second_halfedge_above(h)->bop())
            || 
            (tmp_notf.get_second_face_above(h) != h->face() &&  tmp_notf.get_second_face_above(h)->bop()) );
  }

  bool  face_is_below_first_map(Face_const_handle f) const
  {
    Map_overlay_change_notification  tmp_notf;
    return (tmp_notf.get_first_face_above(f) != f && tmp_notf.get_first_face_above(f)->bop());
  }
  
  bool  face_is_below_second_map(Face_const_handle f) const
  {
    Map_overlay_change_notification  tmp_notf;
    return (tmp_notf.get_second_face_above(f) != f && tmp_notf.get_second_face_above(f)->bop());
  }
  
  const Map_overlay &arr1_, &arr2_;
  const Map_overlay &map_overlay;
  //std::list<Face_const_iterator> faces;
};

CGAL_END_NAMESPACE

#endif




