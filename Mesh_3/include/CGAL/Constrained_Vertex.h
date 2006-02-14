// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Laurent Rineau

#ifndef SOMMET
#define SOMMET

#include <string>
#include <CGAL/basic.h>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <utility>
#include <map> 
#include <set>
#include <list>

#include "Constrained_Element.h"
#include "Constrained_Facet.h"
#include "Constrained_Edge.h"


CGAL_BEGIN_NAMESPACE

template <class Gt>
class Constrained_Vertex 
  : public virtual Constrained_Element <Gt>

{
public :
  typedef Gt                                      Geom_traits;
  typedef typename Gt::FT                         Nb;
  typedef typename Gt::Point_3                    Point; 
  typedef typename Gt::Vector_3                   Vector;
  typedef typename Gt::Segment_3                  Segment;
  typedef typename Gt::Iso_cuboid_3               Iso_Cuboid;
  typedef typename Gt::Triangle_3                 Triangle;
  typedef typename Gt::Plane_3                    Plane;
  typedef typename Gt::Tetrahedron_3              Tetrahedron;
  typedef typename Gt::Line_3                     Line;
  typedef typename Gt::Object_3                   Object;
  typedef typename Gt::Ray_3                      Ray;

  typedef CGAL :: Bbox_3                          Bbox_3;
  typedef Constrained_Element <Gt>								C_elt;
  typedef Constrained_Facet <Gt>							    C_facet;
  typedef Constrained_Edge<Gt>										C_edge;
  typedef Constrained_Vertex<Gt>									C_vertex;

 
 
private :
  Point p;
  int num;
  Geom_traits gt;

 
public :
  //constructeur
  Constrained_Vertex(void) : C_elt("vertex"){}

  Constrained_Vertex(Point pt) : C_elt("vertex"), p(pt) { }
  
  //Declaration Des Fonctions 
  
  Bbox_3 getBBox();
  Point point(); 
 
 
  
  bool does_intersect( Iso_Cuboid c);
  bool does_intersect( Tetrahedron tet);
  bool does_intersect( Segment seg );
  bool does_intersect( C_elt* ce);
  bool does_intersect( Ray ray );
  bool does_intersect_vertex( Constrained_Vertex* s);
  bool does_intersect_edge( C_edge*  a);
  bool does_intersect_facet( C_facet* f); 

  bool operator== (Constrained_Vertex s);
  bool operator!= (Constrained_Vertex s);
  Nb   sqared_distance_to_sheet(Point q);

  Object intersection (Segment seg);
  Object intersection (Ray ray);

};

//-------------------------------------------------------------------
template <class Gt>
bool 
Constrained_Vertex<Gt> ::operator== (Constrained_Vertex s){ 
  
  typename Geom_traits::Equal_3 equal_3 =  gt.equal_3_object();
  
  return ((equal_3(p,s.point()) ) && (p,s.getNum()));
}

//-------------------------------------------------------------------
template <class Gt>
bool 
Constrained_Vertex<Gt> ::operator!= (Constrained_Vertex s){ 
  
  typename Geom_traits::Equal_3 equal_3 =  gt.equal_3_object();
 
  return !((equal_3(p,s.point()) ) && (p,s.getNum()));
}

//-------------------------------------------------------------------
// fonction qui retourne la bounding box de la contrainte 
template <class Gt>
typename Constrained_Vertex<Gt> :: Bbox_3 
Constrained_Vertex<Gt> :: getBBox() { // ????????
  return p.bbox(); //commande CGAL
  //return construct_bbox_3(this.p);
}


//-------------------------------------------------------------------
// fonction qui retourne la bounding boxe Point
template <class Gt>
typename Constrained_Vertex<Gt> :: Point
Constrained_Vertex<Gt> :: point() { // ????????
  return p;
}


//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte un
//isocuboide donne C
template <class Gt>
bool
Constrained_Vertex<Gt> :: does_intersect( Iso_Cuboid c)  
{ 
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
	=gt.has_on_bounded_side_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();

  //std :: cerr<< "Vertex does_intersect "
  //	     << ( (has_on_bounded_side(c,p))  ||  (has_on_boundary(c,p)))  <<endl;
  return ( (has_on_bounded_side(c,p))  ||  (has_on_boundary(c,p)));
} 


//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//cellule donnee c
template <class Gt>
bool
Constrained_Vertex<Gt> ::does_intersect( Tetrahedron tet){ 
 
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
	=gt.has_on_bounded_side_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();
 
  
  
  return ( (has_on_bounded_side(tet,p))  ||  (has_on_boundary(tet,p)));
} 
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//une ligne donnee
template <class Gt>
bool
Constrained_Vertex<Gt> ::does_intersect(Segment seg){ 
 
  typename Geom_traits:: Has_on_3  has_on
	=gt.has_on_3_object();

  return ( has_on(seg,p));
} 
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//une ligne donnee
template <class Gt>
bool
Constrained_Vertex<Gt> ::does_intersect(Ray ray){ 
 
  typename Geom_traits:: Has_on_3  has_on
	=gt.has_on_3_object();

  return ( has_on(ray,p));
} 

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Vertex<Gt> :: does_intersect_vertex(Constrained_Vertex* s)  
{
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();
  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();

  //std::cerr<<"le point : "<<p<<" et l'autre "<<s<<endl;
  return ( equal_3(p,s->point()) ) ;
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Vertex<Gt> :: does_intersect_edge(C_edge* a) 
{
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();
  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();

  Point 
    pp = construct_vertex(a->s,0),
    q = construct_vertex(a->s,1);
   
  //std::cerr<<"le point : "<<p<<" et le segment "<<a->s<<endl;
  return ((equal_3(p,pp) ) || (equal_3(p,q)));
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Vertex<Gt> :: does_intersect_facet( C_facet* facet ) 
{ 
  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();
  
  Triangle t = facet->getTriangle();
  Point r;

  for(int i=0;i<3;i++) {
    r = construct_vertex(t,i);
    if  (equal_3(r,p))
      return true;
  }
  // ne pas oublier de l'effacetr
  //std::cerr<<endl;
  return false;
  
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Vertex<Gt> :: does_intersect( C_elt* ce) 
{
  
  char* type = ce->getType();
  if ((strcmp(type,"vertex") == 0))
    return does_intersect_vertex(dynamic_cast<C_vertex*>(ce) );
  else  if ((strcmp(type,"edge") == 0))
    return does_intersect_edge(dynamic_cast<C_edge*>(ce) );
  else 
    return does_intersect_facet(dynamic_cast<C_facet*>(ce) );
}

//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Vertex<Gt> ::Nb 
Constrained_Vertex<Gt> :: sqared_distance_to_sheet(Point q)
{
  typename Geom_traits:: Compute_squared_distance_3
    compute_squared_distance = gt.compute_squared_distance_3_object();
  
  //std::cerr<<"sommet : "<<p<<endl;

  return compute_squared_distance(p,q);

}

template <class Gt>
typename Constrained_Vertex<Gt> :: Object
Constrained_Vertex<Gt> :: intersection (Segment seg)
{
  return Object();
}


template <class Gt>
typename Constrained_Vertex<Gt> :: Object
Constrained_Vertex<Gt> :: intersection (Ray ray)
{
  return Object();
}

CGAL_END_NAMESPACE
#endif
