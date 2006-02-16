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
// $URL$
// $Id$
// 
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef GPS_DCEL_H
#define GPS_DCEL_H

#include <CGAL/Arr_default_dcel.h>

CGAL_BEGIN_NAMESPACE


class Gps_face : public Arr_face_base
{
  protected:
  mutable char m_info;
  
  enum
  {
    CONTAINED = 1,
    VISITED   = 2
  };


public:

  
  //Constructor
  Gps_face() : Arr_face_base(),
               m_info(0)
  {}

   /*! Assign from another face. */
  virtual void assign (const Arr_face_base& f)
  {
    Arr_face_base::assign (f);

    const Gps_face&  ex_f = static_cast<const Gps_face&>(f);
    m_info = ex_f.m_info;
  }
  
  bool contained() const
  {
    return (m_info & CONTAINED) != 0;
  }

  void set_contained(bool b)
  {
    if(b)
      m_info |= CONTAINED;
    else
      m_info &= ~CONTAINED;
  }

  bool visited() const
  {
    return (m_info & VISITED) != 0;
  }

  void set_visited(bool b) const
  {
    if(b)
      m_info |= VISITED;
    else
      m_info &= ~VISITED;
  }

};


template <class Traits_>
class Gps_dcel :
  public Arr_dcel_base<Arr_vertex_base<typename Traits_::Point_2>,
                       Arr_halfedge_base<typename Traits_::X_monotone_curve_2>,
                       Gps_face>
{
public:

  /*! Default constructor. */
  Gps_dcel()
  {}
};



CGAL_END_NAMESPACE

#endif
