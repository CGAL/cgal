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
// file          : include/CGAL/IO/Qt_ScenesWidget.C
// package       : QT_window
// author(s)     : Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_ScenesWidget.h>
#include <vector>
#include <map>

namespace CGAL {
	void Qt_scenes_widget::show(Qt_scene* s)
  {
    if(scenes_to_display.find(s)!=scenes_to_display.end())
    {
			scenes_to_display[s]=true;
			redraw();
		}
  };


  // redraw shown scenes
  // ***** Should be call when:
  //    - an editable scene is changed (should be call by tools)
  //    - ranges are changed
  void Qt_scenes_widget::redraw()
  {
    clear();
    typedef Map_scene_bool::const_iterator CI;
    for(CI it=scenes_to_display.begin(); it!=scenes_to_display.end();it++)
      if(it->second)
				(it->first)->draw_scene(this);
    emit(redrawed());
  }

  // add a scene in the list of displayable scenes
  void Qt_scenes_widget::add_scene(Qt_scene* s)
  {
    if(scenes_to_display.find(s)==scenes_to_display.end())
      scenes_to_display[s]=false;
    connect(s,SIGNAL(dying(Qt_scene*)),this,SIGNAL(remove_scene(Qt_scene*)));
  }

  // remove a scene from the list of displayable scenes
  void Qt_scenes_widget::remove_scene(Qt_scene* s)
  {
    if(scenes_to_display.find(s)!=scenes_to_display.end())
    {
			bool to_be_redrawn=!scenes_to_display[s];
			scenes_to_display.erase(s);
			if(to_be_redrawn)
			redraw();
    }
  }


} // end namespace CGAL
#include "Qt_ScenesWidget.moc"

#endif // CGAL_USE_QT
