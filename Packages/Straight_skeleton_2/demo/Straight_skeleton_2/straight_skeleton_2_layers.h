// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
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
// Author(s)     : Radu Ursu

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/Cartesian.h>


template <class T>
class Qt_layer_show_skeleton : public CGAL::Qt_widget_layer
{
public:

  Qt_layer_show_skeleton(T &h) : hds(h)
  {};
  void draw()
  {  
    *widget << CGAL::GREEN; 
    for ( Halfedge_iterator i = hds.halfedges_begin(); i != hds.halfedges_end(); ++i ){
      if(i->is_bisector()){
	*widget << i->segment();
      }
    } 
  }
private:
  T &hds;

};//end class 


template <class T>
class Qt_layer_show_polygon : public CGAL::Qt_widget_layer
{
public:
  
  Qt_layer_show_polygon(T &h) : hds(h){};
  void draw()
  {  
    *widget << CGAL::BLUE; 
    for ( Halfedge_iterator i = hds.halfedges_begin(); i != hds.halfedges_end(); ++i ){
      if(! i->is_bisector()){
	*widget << i->segment();
      }
    } 
  }
	
private:
  T &hds;
};//end class 


