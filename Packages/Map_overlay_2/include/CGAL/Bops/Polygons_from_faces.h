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
// file          : include/CGAL/Polygons_from_faces.h
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
