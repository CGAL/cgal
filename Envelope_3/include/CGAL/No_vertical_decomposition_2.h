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
// $Id$ $Date$
// 
//
// Author(s)     : Michal Meyerovitch     <gorgymic@post.tau.ac.il>

#ifndef NO_VERTICAL_DECOMPOSITION_2_H
#define NO_VERTICAL_DECOMPOSITION_2_H

CGAL_BEGIN_NAMESPACE

// No_vertical_decomposition is a class with the interface needed for
// vertical decomposition, which doesn't do a decomposition at all
template <class ArrangementWithOverlayData>
class No_vertical_decomposition_2
{
public:
  typedef ArrangementWithOverlayData                           Pmwx;
  typedef typename Pmwx::Halfedge_iterator                     Halfedge_iterator;
  typedef typename Pmwx::Halfedge_handle                       Halfedge_handle;
  typedef typename Pmwx::Face_handle                           Face_handle;
  typedef typename Pmwx::Face_iterator                         Face_iterator;
  typedef typename Pmwx::Ccb_halfedge_circulator               Ccb_halfedge_circulator;
  typedef typename Pmwx::Hole_iterator                        Hole_iterator;
  typedef typename Pmwx::Face_const_iterator                   Face_const_iterator;
  typedef typename Pmwx::Vertex_iterator                       Vertex_iterator;
  typedef typename Pmwx::Vertex_handle                         Vertex_handle;
    
  typedef typename Pmwx::Traits_2                              Traits;
  typedef typename Traits::X_monotone_curve_2                  X_monotone_curve_2;
  typedef typename Traits::Curve_2                             Curve_2;
  typedef typename Traits::Point_2                             Point_2;

  typedef typename Pmwx::Dcel                                  Pm_dcel;
  typedef typename Pm_dcel::Face_data                          Pm_face_data;
  
  void operator()(Pmwx& pm)
  {
  }
};




CGAL_END_NAMESPACE

#endif  // OVERLAY_2_H
