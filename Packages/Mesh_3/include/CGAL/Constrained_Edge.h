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

#ifndef ARETE
#define ARETE

#include <string>
#include <CGAL/basic.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <utility> 
#include <map> 
#include <set>
#include <list>
#include <vector>
#include <fstream>
#include <CGAL/Timer.h>


#include "Constrained_Element.h"
#include "Constrained_Vertex.h"
#include "Constrained_Facet.h"

CGAL_BEGIN_NAMESPACE
template <class Gt> class Constrained_Vertex;
template <class Gt> class Constrained_Facet;

template <class Gt>
class Constrained_Edge : public virtual  Constrained_Element <Gt>
{

public :
  typedef Gt                                       Geom_traits;
  typedef typename Gt::FT                          Nb;

  typedef typename Gt::Point_3                     Point; 
  typedef typename Gt::Vector_3                    Vector;
  typedef typename Gt::Segment_3                   Segment;
  typedef typename Gt::Iso_cuboid_3                Iso_Cuboid;
  typedef typename Gt::Triangle_3                  Triangle;
  typedef typename Gt::Plane_3                     Plane;
  typedef typename Gt::Tetrahedron_3               Tetrahedron;
  typedef typename Gt::Line_3                      Line;
  typedef typename Gt::Object_3                    Object;
  typedef typename Gt::Ray_3                       Ray;
  typedef CGAL :: Bbox_3                           Bbox_3; 

  typedef Constrained_Element <Gt>                 C_elt;
  typedef Constrained_Vertex<Gt>                   C_vertex;
  typedef Constrained_Facet<Gt>                    C_facet;
 
  //typedef Triangulation_2_traits_3<Geom_traits>     Dt;
public :
  Segment s;
  Geom_traits gt;

public :
  //constructeur
  Constrained_Edge()            : C_elt("edge"){  }
  
  Constrained_Edge(Segment seg) : C_elt("edge"), s(seg) { }

 
  Segment getSegment(){return s;}

  //Declaration Des Fonctions 
  Bbox_3           getBBox() ;
  bool             does_intersect( Iso_Cuboid c) ;
  bool             does_intersect( Tetrahedron tet) ;
  bool             does_intersect( Segment seg );
  bool             does_intersect( C_elt* ce);
  bool             does_intersect( Ray ray );
  bool             does_intersect_vertex (C_vertex* som) ;
  bool             does_intersect_edge   (Constrained_Edge* a)  ;
  bool             does_intersect_facet   (C_facet* f)  ;
  bool             do_intersect2( Iso_Cuboid c, Segment seg);
  bool             do_intersect2( Tetrahedron c, Segment seg);
  
  //  Constrained_Edge create (Paire paire );
  //Constrained_Edge operator= (Constrained_Edge ce);

  Nb               sqared_distance_to_sheet(Point q);

  Object           intersection (Segment seg);
  Object           intersection (Ray ray);
  
}; 

//-------------------------------------------------------------------
//template <class Gt>
//typename Constrained_Edge<Gt> :: Constrained_Edge
//Constrained_Edge<Gt> :: operator= (Constrained_Edge ce)
//{
//  s = ce.s;
//  return Constrained_Edge(ce.s);
//}

//-------------------------------------------------------------------
// fonction qui retourne la bounding box de la contrainte 
template <class Gt>
typename Constrained_Edge<Gt> :: Bbox_3 
Constrained_Edge<Gt> :: getBBox()  { // ????????
  return s.bbox();
}


//-------------------------------------------------------------------
// fonction qui retourne true si le segment s intersecte
// l'isocuboide c
template <class Gt>
bool
Constrained_Edge<Gt> :: do_intersect2(Iso_Cuboid c, Segment seg)  
{
  // on commence par tester si une des extremites du segment est
  // dans c

  typename Geom_traits::Construct_vertex_3 
      construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
    =gt.has_on_bounded_side_3_object();
  typename Geom_traits::Construct_triangle_3
    construct_triangle =  gt.construct_triangle_3_object();

  typename Geom_traits:: Bounded_side_3  bounded_side
    =gt.bounded_side_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();
  
  Point 
    p = construct_vertex(seg,0),
    q = construct_vertex(seg,1);


  //std :: cout<< "do_intersect2 "<< p <<"   "<<q <<endl;

  /**if (bounded_side(c,p)==ON_UNBOUNDED_SIDE) 
     std :: cout<<"****ON_UNBOUNDED_SIDE "<<endl;
  */
  
  if (has_on_bounded_side(c,p) ||  has_on_boundary(c,q) ) 
    {
    return true;
    //std :: cout<< " le premier point est ds la boite "<<endl;
    }
  if (has_on_bounded_side(c,q) ||  has_on_boundary(c,q) )
    {
    return true;
    //std :: cout<< " le second point est ds la boite "<<endl;
    }

  //std :: cout<< "pas de points ds la boite "<<endl;
  //auccune des extremites n'est dans c.
  // on construit les triangles qui forment les facets de c
  Point sc[8];
  for (int i=0; i<8;i++)
    sc[i] = construct_vertex(c,i);
  
  Triangle tf;
  
  //on va regarder pour toutes les facets du cube.
  //on coupe une facet en deux triangle
    tf = construct_triangle(sc[1],sc[3],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[1],sc[3],sc[2]);
  if (do_intersect(tf,seg)) return true;
  
  tf = construct_triangle(sc[1],sc[6],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[0],sc[6],sc[5]);
  if (do_intersect(tf,seg)) return true;
  
  tf = construct_triangle(sc[1],sc[7],sc[6]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[1],sc[7],sc[2]);
  if (do_intersect(tf,seg)) return true;
  
  tf = construct_triangle(sc[2],sc[4],sc[3]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[2],sc[4],sc[7]);
  if (do_intersect(tf,seg)) return true;     
  
  tf = construct_triangle(sc[3],sc[5],sc[0]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[3],sc[5],sc[4]);
  if (do_intersect(tf,seg)) return true;
  
  tf = construct_triangle(sc[5],sc[4],sc[6]);
  if (do_intersect(tf,seg)) return true;
  tf = construct_triangle(sc[6],sc[4],sc[7]);
  if (do_intersect(tf,seg)) return true;
  
  return false;
}

//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte un
//isocuboide donne C
template <class Gt>
  bool
  Constrained_Edge<Gt> :: does_intersect( Iso_Cuboid c)  
{
  //std :: cout<< "Edge does_intersect( Iso_Cuboid c)"<<endl;
  return do_intersect2(c,s);
}

//-------------------------------------------------------------------
// fonction qui retourne true si le segment s intersecte
// l'isocuboide c
template <class Gt>
bool
Constrained_Edge<Gt> :: do_intersect2(Tetrahedron tet, Segment seg)  
{
  // on commence par tester si une des extremites du segment est
  // dans c

  typename Geom_traits::Construct_vertex_3 
      construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Has_on_bounded_side_3  has_on_bounded_side
    =gt.has_on_bounded_side_3_object();
  typename Geom_traits::Construct_triangle_3
    construct_triangle =  gt.construct_triangle_3_object();

  typename Geom_traits:: Bounded_side_3  bounded_side
    =gt.bounded_side_3_object();
  typename Geom_traits:: Has_on_boundary_3  has_on_boundary
    =gt.has_on_boundary_3_object();
  
  Point 
    p = construct_vertex(seg,0),
    q = construct_vertex(seg,1);
  
  if (has_on_bounded_side(tet,p) ||  has_on_boundary(tet,q) ) 
    {
      return true;
    }
  if (has_on_bounded_side(tet,q) ||  has_on_boundary(tet,q) )
    {
      return true;
    }

  //std :: cout<< "pas de points ds la boite "<<endl;
  //auccune des extremites n'est dans tet
  // on construit les triangles qui forment les facets de tet
  Point stet[4];
  for (int i=0; i<4;i++)
    stet[i] = construct_vertex(tet,i);
  
  Triangle tf;
  
  // facet 3
  tf = construct_triangle(stet[0],stet[1],stet[2]);
  if (do_intersect(tf,seg)) return true;
 
  // facet 2
  tf = construct_triangle(stet[0],stet[1],stet[3]);
  if (do_intersect(tf,seg)) return true;
  
  // facet 1
  tf = construct_triangle(stet[0],stet[2],stet[3]);
  if (do_intersect(tf,seg)) return true;
  
  // facet 0
  tf = construct_triangle(stet[1],stet[2],stet[3]);
  if (do_intersect(tf,seg)) return true;
  
  return false;
}
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte un
//tetraedre donne
template <class Gt>
  bool
  Constrained_Edge<Gt> :: does_intersect( Tetrahedron tet)  
{ 
  return do_intersect2(tet,s);
}
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//une ligne donnee
template <class Gt>
  bool
  Constrained_Edge<Gt> :: does_intersect(Segment seg)  
{
  return false; // c'est le cas presque tjs... mais y reflechir !
}
//-------------------------------------------------------------------
//fonction qui retourne true si la contrainte intersecte une
//une ligne donnee
template <class Gt>
  bool
  Constrained_Edge<Gt> :: does_intersect(Ray ray)  
{
  return false; // c'est le cas presque tjs... mais y reflechir !
}
//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Edge<Gt> :: does_intersect_vertex(C_vertex* som) 
{
  typename Geom_traits::Construct_vertex_3 
   construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();

  Point 
    p = construct_vertex(s,0),
    q = construct_vertex(s,1);

  //std::cerr<<"le segment : "<<s<<" et le point "<<som->point()<<endl;
  
  return ((equal_3(som->point(),p) ) || (equal_3(som->point(),q)));
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Edge<Gt> :: does_intersect_edge( Constrained_Edge* a) 
{
  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();

  Point 	
    p1 = construct_vertex(s,0),
    q1 = construct_vertex(s,1),
    p2 = construct_vertex(a->s,0),
    q2 = construct_vertex(a->s,1);
  
  //std::cerr<<"le segment : "<<s<<" et l'autre "<< a->getSegment()
  //   <<endl;
  //std::cerr<<p1<<endl<<p2<<endl<<q1<<endl<<q2<<endl;

  return (equal_3(p1,p2) || equal_3(p1,q2) || 
	  equal_3(q1,p2) || equal_3(q1,q2));
}


//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Edge<Gt> :: does_intersect_facet (C_facet* facet) 
{

  typename Geom_traits::Construct_vertex_3 
    construct_vertex =  gt.construct_vertex_3_object();
  typename Geom_traits:: Equal_3  equal_3 = gt.equal_3_object();


  Point r,
    p = construct_vertex(s,0),
    q = construct_vertex(s,1);

  Triangle t = facet->getTriangle();
  for (int i = 0; i<3;i++) {
    r =  construct_vertex(t,i);
    if ((equal_3(r,p)) ||(equal_3(r,q)) )
      return true;
  }

  return false;
}

//-------------------------------------------------------------------
template <class Gt>
bool
Constrained_Edge<Gt> :: does_intersect( C_elt* ce) 
{
	std::cout<<"entree ds does_intersect( C_elt* ce)"<<std::endl;
  char* type = ce->getType();
  if ((strcmp(type,"vertex") == 0))
    return does_intersect_vertex(dynamic_cast<C_vertex*>(ce) );
  else  if ((strcmp(type,"edge") == 0))
    return does_intersect_edge(dynamic_cast<Constrained_Edge*>(ce) );
  else 
    return does_intersect_facet(dynamic_cast<C_facet*>(ce) ); 
  
}

//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Edge<Gt> :: Nb 
Constrained_Edge<Gt> :: sqared_distance_to_sheet(Point q)
{
  // std::cerr << "Arete : "<<s<<endl;
  //std::cerr<<"point : "<<q<<endl;

  return squared_distance(q,s);

}
//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Edge<Gt> :: Object
Constrained_Edge<Gt> :: intersection (Segment seg)
{
  return Object();
}
//-------------------------------------------------------------------
template <class Gt>
typename Constrained_Edge<Gt> :: Object
Constrained_Edge<Gt> :: intersection (Ray ray )
{
  return Object();
}

CGAL_END_NAMESPACE
#endif
