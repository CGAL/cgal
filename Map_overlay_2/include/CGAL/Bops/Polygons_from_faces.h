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
// Author(s)     : Eti Ezra          <estere@post.tau.ac.il>
#ifndef CGAL_POLYGONS_FROM_FACES_H
#define CGAL_POLYGONS_FROM_FACES_H   


CGAL_BEGIN_NAMESPACE     

template <class Planar_map_, class Polygon_>
class Polygons_from_faces
{
public:
  typedef Planar_map_  Planar_map;
  typedef Polygon_     Polygon;
  //typedef typename Polygon::Point_2  Point;
  //typedef typename Planar_map::Face_const_iterator Face_const_iterator;
  
  template <class InputIterator, class OutputIterator>
  void operator()(InputIterator faces_begin,
                  InputIterator faces_end,
                  OutputIterator polygons)
  {
    //std::list<Face_const_iterator>::iterator  lit;
    InputIterator f_iter;
    
    Polygon                        poly;
    
    for (f_iter = faces_begin; f_iter != faces_end; ++f_iter) 
      {
        
        if (f_iter->is_unbounded() || f_iter->in_hole())
          continue;
        
        poly.erase(poly.vertices_begin(), poly.vertices_end());
        typename Planar_map::Ccb_halfedge_circulator 
          cc=f_iter->outer_ccb();
        do 
          {
            poly.push_back(cc->source()->point());
          } while (++cc != f_iter->outer_ccb());
      
        *polygons++ = poly;
      }
  }
  
};

CGAL_END_NAMESPACE     

#endif
