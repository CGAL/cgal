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
#ifndef CGAL_PMWX_INSERT_INFO_H
#define CGAL_PMWX_INSERT_INFO_H

CGAL_BEGIN_NAMESPACE

template<class _Halfedge_handle>
class Pmwx_insert_info
{
public:

  typedef _Halfedge_handle Halfedge_handle;

  /*! Constructor */
  Pmwx_insert_info() : 
    m_halfedge(Halfedge_handle(NULL))
  {
  }

  void set_halfedge_handle(Halfedge_handle h) {
    m_halfedge = h;
  }

  Halfedge_handle get_halfedge_handle() const {
    return m_halfedge;
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
  
};

CGAL_END_NAMESPACE

#endif // CGAL_PMWX_INSERT_INFO_H

