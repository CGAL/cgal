// Copyright (c) 2002-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Laurent Rineau and Radu Ursu

#include <CGAL/basic.h>


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

  void Qt_widget_layer::toggle(bool b) {
    if(b)
      activate();
    else
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

