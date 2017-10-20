// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Naama mayer       <naamamay@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_H
#define CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Spherical_gaussian_map is a data dtructure that represents a Gaussinal map
 * embedded on the sphere.
 *
 * This file consists of the definition of the main type, namely
 * Arr_spherical_gaussian_map_2 and a service tye,
 * namely Arr_sgm_initializer, that initializes an object of the main type.
 */

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>

#if defined(CGAL_USE_LEDA)
#include <LEDA/numbers/rational.h>
#endif

namespace CGAL {

// #define CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG 1

/*! Sgm_normalizer normalizes the two coordinates of a point in tow-
 * dimensional cartesian coordinate-system. It is parameterized by the
 * field number-type of the point coordinates and by the type of the point
 * which coordinates are to be normalized. Normalization is performed through
 * a call an operator that accepts the point as a parameter.
 * The normalization of field number-types that initiate normalization
 * automatically (e.g., Gmpq) do not need to be initiated. Therefore, the
 * default instance does nothing. Normalizers for field number-types that
 * do not initiate normalization automatically, must be specialed.
 */
template <class FT, class Point_2>
class Sgm_normalizer {
public:
  /*! Normalize the coordinates of the given point, but in fact does
   * nothing.
   * \param p the point which coordinates are to be normalized
   */ 
  void operator()(Point_2 &) {}
};

#if defined(CGAL_USE_LEDA)
/*! Sgm_normalizer normalizes the two coordinates of a point in tow-
 * dimensional cartesian coordinate-system. It is parameterized by the
 * type of the point which coordinates are to be normalized, and
 * specialized for the leda::rational field number-type.
 */
template <class Point_2>
class Sgm_normalizer<leda::rational, Point_2> {
public:
  /*! Normalize the coordinates of the given point
   * \param p the point which coordinates are to be normalized
   */
  void operator()(Point_2 & p)
  {
    leda::rational x = p.x();
    x.normalize();
    leda::rational y = p.y();
    y.normalize();
    leda::rational z = p.z();
    z.normalize();
    p = Point_2(x, y, z);
  }
};
#endif

/*! Arr_sgm_initializer is an algorothmic framework that initializes a
 * Arr_spherical_gaussian_map_3 structure. It is parameterized by the SGM to
 * be initialized and by a visitor class.
 */
template <typename Sgm, typename T_Traits = typename Sgm::Traits>
class Arr_sgm_initializer {
public:
  typedef T_Traits                                        Traits;
  typedef typename Traits::Vector_3                       Vector_3;

  typedef typename Sgm::Geometry_traits_2                 Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Geometry_traits_2::Curve_2             Curve_2;

  typedef typename Sgm::Vertex_handle                     Vertex_handle;
  typedef typename Sgm::Halfedge_handle                   Halfedge_handle;
  typedef typename Sgm::Face_handle                       Face_handle;

  /*! Constructor */
  Arr_sgm_initializer(Sgm & sgm) : m_sgm(sgm) { }

  /*! Destructor */
  virtual ~Arr_sgm_initializer() {}
  
  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the SGM. Each normal defines an end point of the greate arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   */
  Halfedge_handle insert_non_intersecting(const Vector_3 & normal1,
                                          const Vector_3 & normal2)
  {
    X_monotone_curve_2 xc(normal1.direction(), normal2.direction());
    return insert_non_intersecting_curve(m_sgm, xc);
  }

  /*! Make x-monotone
   */
  template<typename OutputIterator>
  OutputIterator make_x_monotone(const Vector_3 & normal1,
                                 const Vector_3 & normal2,
                                 OutputIterator oi)
  {
    Curve_2 cv(normal1.direction(), normal2.direction());
    const Geometry_traits_2 * traits = this->m_sgm.geometry_traits();
    oi = traits->make_x_monotone_2_object()(cv, oi);
    return oi;
  }
  
  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the SGM. Each normal defines an end point of the greate arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   * \pre the SGM is empty.
   */
  template<typename OutputIterator>
  OutputIterator insert(const Vector_3 & normal1, const Vector_3 & normal2,
                        OutputIterator oi)
  {
    std::list<CGAL::Object> x_objects;
    make_x_monotone(normal1, normal2, std::back_inserter(x_objects));

    typename std::list<CGAL::Object>::iterator it = x_objects.begin();
    const X_monotone_curve_2 * xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "1.a. insert_in_face_interior(" << *xc << ")" << std::endl;
#endif
    Halfedge_handle he =
      m_sgm.insert_in_face_interior(*xc, m_sgm.faces_begin());
    if (!xc->is_directed_right()) he = he->twin();
    *oi++ = he;

    ++it;
    if (it == x_objects.end()) return oi;

    xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "1.b. insert_from_vertex(" << *xc << ")" << std::endl;
#endif
    *oi++ = (xc->is_directed_right()) ?
      m_sgm.insert_from_left_vertex(*xc, he->target()) :     
      m_sgm.insert_from_right_vertex(*xc, he->target());
    return oi;
  }

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the SGM. Each normal defines an end point of the greate arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   * \param vertex the handle of the vertex that is the source of the arc
   */
  template<typename OutputIterator>
  OutputIterator insert(const Vector_3 & normal1, Vertex_handle vertex1,
                        const Vector_3 & normal2,
                        OutputIterator oi)
  {
    std::list<CGAL::Object> x_objects;
    make_x_monotone(normal1, normal2, std::back_inserter(x_objects));

    typename std::list<CGAL::Object>::iterator it = x_objects.begin();
    const X_monotone_curve_2 * xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "2.a. insert_from_vertex(" << *xc << ", "
              << vertex1->point() << ")" << std::endl;
#endif

    Halfedge_handle he = (xc->is_directed_right()) ?
      m_sgm.insert_from_left_vertex(*xc, vertex1) :
      m_sgm.insert_from_right_vertex(*xc, vertex1);
    *oi++ = he;

    ++it;
    if (it == x_objects.end()) return oi;

    xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "2.b. insert_from_vertex(" << *xc << ")" << std::endl;
#endif
    *oi++ = (xc->is_directed_right()) ?
      m_sgm.insert_from_left_vertex(*xc, he->target()) :
      m_sgm.insert_from_right_vertex(*xc, he->target());
    return oi;
  }

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the SGM. Each normal defines an end point of the greate arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \param vertex the handle of the vertex that is the source of the arc
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   */
  template<typename OutputIterator>
  OutputIterator insert(const Vector_3 & normal1, 
                        const Vector_3 & normal2, Vertex_handle vertex2,
                        OutputIterator oi)
  {
    std::list<CGAL::Object> x_objects;
    make_x_monotone(normal1, normal2, std::back_inserter(x_objects));

    typename std::list<CGAL::Object>::iterator it = x_objects.begin();
    if (x_objects.size() == 1) {
      const X_monotone_curve_2 * xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
      std::cout << "3. insert_from_vertex(" << *xc << ")" << std::endl;
#endif
      Halfedge_handle he = (xc->is_directed_right()) ?
        m_sgm.insert_from_right_vertex(*xc, vertex2) :
        m_sgm.insert_from_left_vertex(*xc, vertex2);
      *oi++ = he->twin();
      return oi;
    }

    const X_monotone_curve_2 * xc1 = object_cast<X_monotone_curve_2>(&(*it++));
    const X_monotone_curve_2 * xc2 = object_cast<X_monotone_curve_2>(&(*it));

#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "3.a. insert_from_vertex(" << *xc2 << ")" << std::endl;
#endif
    Halfedge_handle he2 = (xc2->is_directed_right()) ?
      m_sgm.insert_from_right_vertex(*xc2, vertex2) :
      m_sgm.insert_from_left_vertex(*xc2, vertex2);
    he2 = he2->twin();
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "3.b. insert_from_vertex(" << *xc1 << ")" << std::endl;
#endif
    Halfedge_handle he1 = (xc1->is_directed_right()) ?
      m_sgm.insert_from_right_vertex(*xc1, he2->source()) :
      m_sgm.insert_from_left_vertex(*xc1, he2->source());
    he1 = he1->twin();
    *oi++ = he1;
    *oi++ = he2;
    return oi;
  }

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the SGM. Each normal defines an end point of the greate arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \param vertex1 the handle of the vertex that is the source of the arc
   * \param vertex2 the handle of the vertex that is the target of the arc
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   */
  template<typename OutputIterator>
  OutputIterator insert(const Vector_3 & normal1, Vertex_handle vertex1,
                        const Vector_3 & normal2, Vertex_handle vertex2,
                        OutputIterator oi)
  {
    std::list<CGAL::Object> x_objects;
    make_x_monotone(normal1, normal2, std::back_inserter(x_objects));
    typename std::list<CGAL::Object>::iterator it = x_objects.begin();
    if (x_objects.size() == 1) {
      const X_monotone_curve_2 * xc = object_cast<X_monotone_curve_2>(&(*it));
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
      std::cout << "4. insert_at_vertices(" << *xc << ")" << std::endl;
#endif
      *oi++ = m_sgm.insert_at_vertices(*xc, vertex1, vertex2);
      return oi;
    }

    const X_monotone_curve_2 * xc1 = object_cast<X_monotone_curve_2>(&(*it++));
    const X_monotone_curve_2 * xc2 = object_cast<X_monotone_curve_2>(&(*it));

#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "4.a. insert_from_vertex(" << *xc1
              << "," << vertex1->point() << ")" << std::endl;
#endif
    Halfedge_handle he = (xc1->is_directed_right()) ?
      m_sgm.insert_from_left_vertex(*xc1, vertex1) :
      m_sgm.insert_from_right_vertex(*xc1, vertex1);
    *oi++ = he;
#if CGAL_ARR_SPHERICAL_GAUSSIAN_MAP_3_DEBUG==1
    std::cout << "4.b. insert_at_vertices(" << *xc2 << ")" << std::endl;
#endif
    *oi++ = m_sgm.insert_at_vertices(*xc2, he->target(), vertex2);
    return oi;
  }

protected:
  /*! The SGM to initialize */
  Sgm & m_sgm;
};

/* Spherical_gaussian_map is a data dtructure that represents a Gaussinal map
 * embedded on the sphere.
 */
template <class T_Traits,
          template <class T>
          class T_Dcel = Arr_default_dcel>
class Arr_spherical_gaussian_map_3 :
  public Arrangement_on_surface_2<T_Traits,
    Arr_spherical_topology_traits_2<T_Traits,
      T_Dcel<T_Traits>
    >
  >
{
private:
  typedef Arr_spherical_gaussian_map_3<T_Traits, T_Dcel>    Self;
  
public:
  typedef T_Traits                                          Traits;
  typedef Traits                                            Geometry_traits_2;
  
  typedef Arrangement_on_surface_2<Traits, 
	Arr_spherical_topology_traits_2<Traits, T_Dcel<Traits> > >
															Base;

  /*! Parameter-less Constructor */
  Arr_spherical_gaussian_map_3() { }

  /*! Copy Constructor */
  Arr_spherical_gaussian_map_3
  (const Arr_spherical_gaussian_map_3 &)
  {
    // Not implemented yet!
    CGAL_error();
  }

  /*! Destructor */
  virtual ~Arr_spherical_gaussian_map_3() { this->clear(); }

#if 0
  /*! Clear the internal representation and auxiliary data structures */
  void clear()
  {
    CGAL_error_msg( "Not implemented yet!");
  }
  
  /*! returns true if the representation is empty */
  bool is_empty() const
  {
    CGAL_error_msg( "Not implemented yet!");
    return true;
  }

  /*! Return the degree of a vertex
   * \param vh the vertex handle
   * \return the degree
   */
  unsigned int degree(Vertex_const_handle vh) const
  {
    CGAL_error_msg( "Not implemented yet!");
    return vh->degree();
  }
#endif
  
  /*! Print statistics */
  void print_stat()
  {
    std::cout << "No. vertices: " << this->number_of_vertices()
              << ",  no. halfedges: " << this->number_of_halfedges()
              << ",  no. faces: " << this->number_of_faces()
              << std::endl;
  }

  /*! Allow the initializer to update the SGM data members */
  friend class CGAL::Arr_sgm_initializer<Self>;
};

} //namespace CGAL

#endif
