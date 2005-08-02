// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Kapelushnik Lior <liorkape@post.tau.ac.il>

/*! \file
 * spherical arrangements of none intersecting arcs of great circles on a sphere
 */

#ifndef CGAL_SPHERE_POINT_LOCATION_H_
#define CGAL_SPHERE_POINT_LOCATION_H_

#include "CGAL/Spherical_map.h"
#include <CGAL/Object.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>

CGAL_BEGIN_NAMESPACE

/*
 generic class for point location, can be templated with the required
 cgm arrangement point location

 _Spherical_map - the spherical map to point locate over
 _PointLocation - the type of point location in a cgm face
*/
template<class _Spherical_map,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  template <class _T>
#endif
  class _PointLocation >
class Sphere_base_point_location {
public:  
  // public definitions
  typedef _Spherical_map                   Spherical_map;
  typedef typename Spherical_map::CGM           CGM;
  typedef typename Spherical_map::Vertex_const_handle   Vertrx_const_handle;
  typedef typename Spherical_map::Vertex_handle       Vertrx_handle;
  typedef typename Spherical_map::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Spherical_map::Halfedge_handle      Halfedge_handle;
  typedef typename Spherical_map::Face_const_handle     Face_const_handle;  
  typedef typename Spherical_map::Face_handle       Face_handle;  
  typedef typename Spherical_map::Direction_3         Direction_3;
  typedef typename Spherical_map::Point_3           Point_3;
  typedef typename CGM::Arrangement            Int_arr;
  // the templated point location
  typedef _PointLocation<Int_arr>              Point_location;
  
  /*
   constructor, based uppon a spherical map

   sphere - the spherical map to point locate
  */
  Sphere_base_point_location(Spherical_map *sphere):
    m_sphere(sphere) {}

  /*
   empty constructor
  */
  Sphere_base_point_location(): m_sphere(0) {}

  /*
   set the spherical object to point locate over

   sphere - the spherical map to point locate
  */
  void init(Spherical_map * sphere) {
    m_sphere = sphere;
  }

  /*
   locate the object at a sphere direction

   d - direction on the sphere to locate it's object
   
   return value - an object stating the element at direction d,
     the element may be a sphere Vertex_handle, Halfedge_handle or Face_handle
  */
  CGAL::Object locate(const Direction_3 &d) {
    m_sphere->update(); // make sure the search sphere has updated data
    Projected_normal pn; // projected normal with direction d
    pn.compute_projection(d.vector());

    // find a cubical arrangement id which holds a direction point on the cube
    // with the direction d
    unsigned int faces_mask = pn.get_faces_mask();
    unsigned int face_id;
    for (face_id = 0; face_id < CGM::NUM_FACES; ++face_id) {
      if (faces_mask & CGM::get_mask(face_id)) {
        break;
      }
    }    
    
    // find the 2d point on the cubical face
    const Point_3 &point3d = pn.get_projected_normal();
    Point_2 point = (m_sphere->get_cgm())->construct_point(point3d, face_id);

    // locate on the cubical map
    Point_location pl((m_sphere->get_cgm())->get_arrangement(face_id));
    CGAL::Object obj = pl.locate(point);

    // handle possible cubical elements cases
    if (const Int_vertex_const_handle *const_ver_ptr =
      CGAL::object_cast<Int_vertex_const_handle>(&obj)) {
      // in case the direction points to a cubical vertex

      Int_vertex_const_handle ver = *const_ver_ptr;
      if (ver->getReal()) {
        // a real vertex, return it
        while (!ver->is_rep()) {
          // get the representative vertex in degenerate cases
          void *tmp = ver->get_adjacent_vertex();
          ver = *((Int_vertex_const_handle *)(&tmp));
        }  
      } else {
        // not a real vertex, either a corner or part of an arc
        if (ver->degree() == 2) {
          // a corner has degree 2, all other unreal vertices are
          // located on a cubic edge and have another incoming halfedge
          // that caused their creation so they have degree > 2      
          Int_halfedge_around_vertex_const_circulator icirc;
          icirc = ver->incident_halfedges();
          int ind;
          // check if a corner that is part of a real halfedge
          for (ind = 0; ind < 2; ++ind) {
            if (icirc->get_is_real()) {
              break;        
            }
            ++icirc;        
          }
        
          if (!icirc->get_is_real()) {
            // corner, not part of a sphere edge

            // find an edge near the corner vertex
            Int_halfedge_const_handle ihedge;
            ihedge = ver->incident_halfedges();
            if ((ihedge->face())->is_unbounded()) {
              ihedge = ihedge->twin();
            }  
            // return the face of adjacent halfedge
            return CGAL::make_object
            (*((Face_handle *)((ihedge->face())->getSphereFace())));    
          } else {
            // corner, part of a sphere edge
            if ((icirc->face())->is_unbounded()) {
              icirc = icirc->twin();
            }  
            while (!icirc->target()->getReal()) {
              // find the last halfedge of the spherical halfedge that
              // holds the cubical halfedge
              icirc = icirc->next();
              while (!icirc->get_is_real()) {
                icirc = m_sphere->get_cgm()->get_adjacent_halfedge_handle(icirc);
                icirc = icirc->next();
              }
            }
            return CGAL::
            make_object(*((Halfedge_handle *)(icirc->getSphereEdge())));
          }  
        } else {
          // unreal vertex, not a corner, part of a sphere halfedge
          Int_halfedge_around_vertex_const_circulator icirc;
          icirc = ver->incident_halfedges();
          // locate a real cgm halfedge around the vertex (part of sphere halfedge)
          while (!icirc->get_is_real()) {
            ++icirc;
          }
          if ((icirc->face())->is_unbounded()) {
            icirc = icirc->twin();
          }  
          // find the last halfedge of the spherical halfedge that  
          // holds the cubical halfedge
          while (!icirc->target()->getReal()) {
            icirc = icirc->next();
            while (!icirc->get_is_real()) {
              icirc = m_sphere->get_cgm()->get_adjacent_halfedge_handle(icirc);
              icirc = icirc->next();
            }
          }
          return CGAL::
          make_object(*((Halfedge_handle *)(icirc->getSphereEdge())));
        }      
      }
      // a real vertex, return the spherical vertex
      return CGAL::
      make_object(*((Vertrx_handle *)(ver->getSphereVertex())));
    }

    
    if (const Int_halfedge_const_handle *const_hed_ptr =
      CGAL::object_cast<Int_halfedge_const_handle>(&obj)) {
      
      // in case the direction points to a cubical halfedge
      Int_halfedge_const_handle  ihedge = *const_hed_ptr;

      if (ihedge->get_is_real()) {
        // part of a real sphere halfedge
        if ((ihedge->face())->is_unbounded()) {
          ihedge = ihedge->twin();
        }
        // find the last cgm halfedge of the spherical halfedge
        while (!ihedge->target()->getReal()) {
          ihedge = ihedge->next();
          while (!ihedge->get_is_real()) {
            ihedge = m_sphere->get_cgm()->get_adjacent_halfedge_handle(ihedge);
            ihedge = ihedge->next();
          }
        }
        return CGAL::
        make_object(*((Halfedge_handle *)(ihedge->getSphereEdge())));          
      } else {            
        // part of cube unreal halfedge, return a face near halfedge  
        if ((ihedge->face())->is_unbounded()) {
          ihedge = ihedge->twin();
        }
        return CGAL::
        make_object
        (*((Face_handle *)((ihedge->face())->getSphereFace())));
      }
    }
    
    if (const Int_face_const_handle *const_face_ptr =
      CGAL::object_cast<Int_face_const_handle>(&obj)) {
      // found a cubical face, return it's spherical face handle

      Int_face_const_handle  iface = *const_face_ptr;
      return CGAL::make_object
      (*((Face_handle *)(iface->getSphereFace())));  
    }

    // getting here means that in the cubical gaussian map the located element
    // was not a face, halfedge or vertex, should not get here
    std::cerr << "internal error locating direction on cgm" << std::endl;
    CGAL_assertion(false);
    return CGAL::make_object(0);
  }
 
private:
  // private definitions
  typedef typename CGM::Projected_normal         Projected_normal;
  typedef typename CGM::Kernel::Point_2         Point_2;
  typedef typename CGM::Arr_vertex_const_handle    Int_vertex_const_handle;  
  typedef typename CGM::Arr_halfedge_const_handle     Int_halfedge_const_handle;
  typedef typename CGM::Arr_face_const_handle      Int_face_const_handle;  
  //  typedef typename CGM::Arr_ccb_halfedge_circulator   Int_halfedge_circulator;
  typedef typename CGM::Arr_halfedge_around_vertex_const_circulator Int_halfedge_around_vertex_const_circulator;

  Spherical_map *m_sphere; // the spherical map to point locate over
};

/*
 a spherical point location class that does the cgm point location using
 the naive point location method
 inherits from the general point location templated with naive point location

 _Spherical_map - the spherical map to point locate over
*/
template <class _Spherical_map>
class Sphere_naive_point_location : 
public Sphere_base_point_location<_Spherical_map, CGAL::Arr_naive_point_location>
{ 
public:
  
  /*
   constructor, based uppon a spherical map

   sphere - the spherical map to point locate
  */
  Sphere_naive_point_location(_Spherical_map *sphere):
    Sphere_base_point_location<_Spherical_map, CGAL::Arr_naive_point_location>::
  Sphere_base_point_location(sphere) {}
  
  /*
   empty constructor
  */
  Sphere_naive_point_location() {}
};

/*
 a spherical point location class that does the cgm point location using
 the walk along line point location method
 inherits from the general point location templated with walk along line point location

 _Spherical_map - the spherical map to point locate over
*/
template <class _Spherical_map>
class Sphere_walk_along_line_point_location : 
public Sphere_base_point_location<_Spherical_map, CGAL::Arr_walk_along_line_point_location>
{ 
public:
  /*
   constructor, based uppon a spherical map

   sphere - the spherical map to point locate
  */
  Sphere_walk_along_line_point_location(_Spherical_map *sphere):
    Sphere_base_point_location<_Spherical_map, CGAL::Arr_walk_along_line_point_location>::
  Sphere_base_point_location(sphere) {}
  
  /*
   empty constructor
  */
  Sphere_walk_along_line_point_location() {}
};

CGAL_END_NAMESPACE

#endif
