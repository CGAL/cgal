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
// file          : include/CGAL/IO/Qt_ScenesWidget.h
// package       : Qt_widget
// author(s)     : Laurent Rineau
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifndef CGAL_QT_SCENESWIDGET_H
#define CGAL_QT_SCENESWIDGET_H

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_scene.h>
#include <map>

namespace CGAL {

class Qt_scenes_widget : public Qt_widget {
Q_OBJECT
public:
  // constructor
  Qt_scenes_widget(QWidget *parent = 0, const char *name = 0) 
    : Qt_widget(parent,name), scenes_to_display() {};

  // destructor
  ~Qt_scenes_widget() {};

  using Qt_widget::show; // not to mask void Qt_widget::show()
  void show(Qt_scene* s);

  
	using Qt_widget::operator<<; // not to mask Qt_widget's << operators
	// the following operator add a scene and show it
  inline
  Qt_scenes_widget& operator<<(Qt_scene* s)
  {
    add_scene(s);
    show(s);
  };
  

signals:
  void redrawed(); // some programers may want to paint something on scenes
	

public slots:

  // redraw shown scenes
  // ***** Should be call when:
  //    - an editable scene is changed (should be call by tools)
  //    - ranges are changed
  void redraw();
  
  // add a scene in the list of displayable scenes
  void add_scene(Qt_scene* s);

  // remove a scene from the list of displayable scenes
  void remove_scene(Qt_scene* s);
  
private:
  typedef std::map<Qt_scene*,bool> Map_scene_bool;
  Map_scene_bool scenes_to_display;
};

} // end namespace CGAL

#endif // CGAL_QT_SCENESWIDGET_H
