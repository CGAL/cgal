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
// file          : include/CGAL/IO/Qt_widget_layer.C
// package       : Qt_widget
// author(s)     : Laurent Rineau & Radu Ursu
// release       : 
// release_date  : 
//
// coordinator   : Laurent Rineau <rineau@clipper.ens.fr>
//
// ============================================================================

#ifdef CGAL_USE_QT

#include <CGAL/IO/Qt_widget_layer.h>

namespace CGAL {
	void	Qt_widget_layer::attach(Qt_widget *w) {
				widget=w;
        if(activate())
				  emit(activated(this));
	}
  void Qt_widget_layer::stateChanged(int i){
    if(i==2)
      activate();
    else if(i == 0)
      deactivate();
  }
  bool Qt_widget_layer::activate(){
    if(active)
      return false;
    else {
      active = true;
      activating();
      emit(activated(this));
      return true;
    }
  }
  bool Qt_widget_layer::deactivate(){
    if(!active)
      return false;
    else {
      active = false;
      deactivating();
      emit(deactivated(this));
      return true;
    }
      
  }
} // namespace CGAL

#include "Qt_widget_layer.moc"

#endif // CGAL_USE_QT
