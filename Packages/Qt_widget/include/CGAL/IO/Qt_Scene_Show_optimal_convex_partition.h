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
// file          : include/CGAL/IO/Qt_Scene_Show_green_approximation.h
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_SCENE_SHOW_OPTIMAL_CONVEX_H
#define CGAL_QT_SCENE_SHOW_OPTIMAL_CONVEX_H

#include <CGAL/IO/Qt_Scene.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>


namespace CGAL {

template <class T>
class Qt_scene_show_optimal_convex : public Qt_scene
{
    //Q_OBJECT
public:
  typedef typename T::FT	      FT;
  typedef CGAL::Cartesian<FT>	      K;
  typedef CGAL::Partition_traits_2<K> Traits;


  Qt_scene_show_optimal_convex(T &p) : polygon(p)
  {};
  void draw_scene(Qt_widget &widget)
  {
    optimal_convex.clear();
    Traits  partition_traits;
    
    CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                       polygon.vertices_end(),
                           std::back_inserter(optimal_convex),
                                            partition_traits);    
    
    std::list<T>::const_iterator p_it;
    for(p_it = optimal_convex.begin(); p_it != optimal_convex.end(); p_it++)
    {
      widget << CGAL::YELLOW; 
      widget << *p_it;
    }
    
  };
private:
  T		&polygon;
  std::list<T>	optimal_convex;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_WINDOW_GET_SEGMENT_H
