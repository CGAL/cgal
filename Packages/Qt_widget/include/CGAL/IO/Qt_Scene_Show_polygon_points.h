// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// file          : include/CGAL/IO/Qt_Scene_Show_polygon.h
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_SCENE_SHOW_POLYGON_POINTS_H
#define CGAL_QT_SCENE_SHOW_POLYGON_POINTS_H

#include <CGAL/IO/Qt_Scene.h>
#include <qobject.h>


namespace CGAL {

template <class T>
class Qt_scene_show_polygon_points : public Qt_scene
{
    //Q_OBJECT
  typedef typename T::Point_2	Point_2;
public:
  

  Qt_scene_show_polygon_points(T &p) : polygon(p){};
  void draw_scene(Qt_widget &widget)
  {
    typename T::const_iterator vert_it;

    
    for (vert_it = polygon.vertices_begin(); vert_it != polygon.vertices_end(); 
        vert_it++)
    {
      
      widget << CGAL::GREEN << CGAL::PointSize (5) << CGAL::PointStyle (CGAL::DISC);
      widget << Point_2((*vert_it).x(), (*vert_it).y());
    }
  };
	
private:
  T &polygon;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_SCENE_SHOW_POLYGON_H
