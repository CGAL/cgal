// Copyright (c) 1997-2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Laurent Rineau & Radu Ursu

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
