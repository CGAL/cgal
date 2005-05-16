// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Tali Zvi <talizvi@post.tau.ac.il>
//                 Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_ARR_INSERT_INFO_H
#define CGAL_ARR_INSERT_INFO_H

CGAL_BEGIN_NAMESPACE

template<class _Halfedge_handle>
class Arr_insert_info
{
public:

  typedef _Halfedge_handle Halfedge_handle;

  /*! Constructor */
  Arr_insert_info() : m_halfedge(Halfedge_handle(NULL)),
                       m_right_curves_counter(0)

  {
  }

  void set_halfedge_handle(Halfedge_handle h) {
    m_halfedge = h;
  }

  Halfedge_handle get_halfedge_handle() const {
    return m_halfedge;
  }

  unsigned int dec_right_curves_counter()
  {
     return --m_right_curves_counter;
  }

  void inc_right_curves_counter()
  {
    ++m_right_curves_counter;
  }

  unsigned int get_right_curves_counter() const
  {
    return m_right_curves_counter;
  }


  void Print()
  {
    if ( m_halfedge == Halfedge_handle(NULL))
      std::cout << "halfedge: NULL\n";
    else 
      std::cout << "halfedge: " << m_halfedge->source()->point() 
		<< " " << m_halfedge->target()->point() << "\n";
  }
 
private:
  Halfedge_handle m_halfedge;
  unsigned int m_right_curves_counter;
  
};

CGAL_END_NAMESPACE

#endif // CGAL_ARR_INSERT_INFO_H

