// Copyright (c) 2001  Utrecht University (The Netherlands),
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
// Author(s)     : 

#include <CGAL/IO/PS_edge_3.h>

CGAL_BEGIN_NAMESPACE
  
//Constructor
PS_edge_3::PS_edge_3(Point3 &PA,Point3 &PB,Color color,VISIBILITY v)
  : _first_point(PA),_second_point(PB),_edge_color(color),_visibility(v) {}

// Copy constructor
PS_edge_3::PS_edge_3(const PS_edge_3& ar) {
_first_point=ar.first_point();
_second_point=ar.second_point();
_edge_color=ar.color();
_visibility=ar.visibility();
}

//Setting operator
PS_edge_3 PS_edge_3::operator=(const PS_edge_3& ar){
_first_point=ar.first_point();
_second_point=ar.second_point();
_edge_color=ar.color();
_visibility=ar.visibility();
return *this;
}

ostream& operator<<(ostream& os, const PS_edge_3& ar) {
  os << "(" << ar.first_point() << "),(" << ar.second_point() << ")| Color : "
     << ar.color(); 
  if (ar.visibility() == VISIBLE) cout << "| VISIBLE ";
  else cout << "| INVISIBLE ";

  cout << endl;
return os;
}

//Return the minimum or the maximum of the two points of the edge
//according the x,y,z 
coord_type PS_edge_3::xmin() {
  coord_type xa = _first_point.x();
  coord_type xb = _second_point.x();
  if (xa>xb) return xb;
  else return xa;
}
  
coord_type PS_edge_3::ymin() {
  coord_type ya = _first_point.y();
  coord_type yb = _second_point.y();
  if (ya>yb) return yb;
  else return ya;
}
  
coord_type PS_edge_3::zmin() {
  coord_type za = _first_point.z();
  coord_type zb = _second_point.z();
  if (za>zb) return zb;
  else return za;
}

coord_type PS_edge_3::xmax() {
  coord_type xa = _first_point.x();
  coord_type xb = _second_point.x();
  if (xa<xb) return xb;
  else return xa;
}
  
coord_type PS_edge_3::ymax() {
  coord_type ya = _first_point.y();
  coord_type yb = _second_point.y();
  if (ya<yb) return yb;
  else return ya;
}
  
coord_type PS_edge_3::zmax() {
  coord_type za = _first_point.z();
  coord_type zb = _second_point.z();
  if (za<zb) return zb;
  else return za;
}

//Function to transform the two points of the edge according the
//Transformation t
void PS_edge_3::transformation(Transformation &t) {
  _first_point=t(_first_point);
  _second_point=t(_second_point);
}

CGAL_END_NAMESPACE
