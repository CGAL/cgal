// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_FACE_EXTENDED_DCEL_H
#define CGAL_ARR_FACE_EXTENDED_DCEL_H

/*! \file
 * This class extends the face topological-feature of the Dcel. While it is
 * possible to store extra (non-geometric) data with the curves or points by
 * extending their types respectively, it is also possible to extend the
 * vertex, halfedge, or face types of the Dcel through inheritance. Many
 * times it is desired to associate extra data just with the arrangement
 * faces. For example, when an arrangement represents the subdivision of a
 * country into regions associated with their population density. In this
 * case, there is no alternative other than to extend the Dcel face. As this
 * technique is somewhat cumbersome and difficult for inexperienced users, we
 * provide the class-template Arr_face_extended_dcel<FaceData, Dcel>, which
 * extends each face in the Dcel (defaulted to Arr_default_dcel) class with a
 * FaceData object.
 */

#include <CGAL/basic.h>
#include <CGAL/Arr_dcel_base.h>

CGAL_BEGIN_NAMESPACE

/*! \class An extended Dcel face */
template <class T_Face, class T_Face_data> class Arr_extended_face :
  public T_Face
{
public:
  typedef T_Face_data                                   Face_data;
  typedef T_Face                                        Face;
  
  /*! Constructor */
  Arr_extended_face() {}

  /*! Destructor */
  virtual ~Arr_extended_face() {}

  /*! Set the data */
  void set_data(const Face_data & data) { m_data = data; }

  /*! Obtain a const version of the data */
  const Face_data & get_data() const { return m_data; }

  /*! Obtain a non const version of the data */
  Face_data & get_data() { return m_data; }
    
private:
  /*! The data the face is extended with */
  Face_data m_data;
};

/*! \class A Dcel with extended face */
template <class T_Face_data,
          class T_Traits,
          class T_Vertex_base = Arr_vertex_base<typename T_Traits::Point_2>,
          class T_Halfedge_base =
            Arr_halfedge_base<typename T_Traits::X_monotone_curve_2>,
          class T_Face_base = Arr_face_base>
class Arr_face_extended_dcel :
  public Arr_dcel_base<T_Vertex_base,
                       T_Halfedge_base,
                       Arr_extended_face<T_Face_base, T_Face_data> >
{
public:
  /*! Constructor */
  Arr_face_extended_dcel() {}
  
  /*! Destructor */
  virtual ~Arr_face_extended_dcel() {}
};

CGAL_END_NAMESPACE

#endif
