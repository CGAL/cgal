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
// Author(s)     : Radu Ursu

#ifndef CGAL_apply_to_range_h
#define CGAL_apply_to_range_h

#include <CGAL/Triangulation_2.h>
#include <CGAL/Unique_hash_map.h>
#include <stack>

namespace CGAL{


template <class Tr, class Fct, class R>
void apply_to_range(const Tr &t, 
             const Point_2<R> &p1, const Point_2<R> &p2,
             Fct& fct)
{
  if (t.dimension()<2) return;
  typedef typename Tr::Point         POINT;
  typedef typename Tr::Face_handle   hFACE;
  typedef typename Tr::Vertex_handle hVERTEX;
  typedef typename Tr::Line_face_circulator      LFC;
  typedef typename Tr::Finite_vertices_iterator  FVI;
  typedef typename Tr::Finite_faces_iterator     FFI;
  typedef typename Kernel_traits<POINT>::Kernel K; 
  typedef typename K::FT FT;

  LFC l1, l2, l3, l4;         //the faces that intersect the pixmap RECTANGLE
  hFACE hface1, hface2, 
        hface3, hface4;       //the faces where we start to search
  FT    xr_left,   yr_top, 
        xr_right,  yr_bottom;//the coordinates of the screen boundaries  
  CGAL::Unique_hash_map<hFACE, bool> visited(false);//used for DFS  
  std::stack<hFACE>               face_stack; //used for DFS
  
  xr_left = p1.x(); xr_right = p2.x();
  yr_top = p1.y(); yr_bottom = p2.y();

  hface1 = t.locate(POINT(xr_left, yr_top));
  hface2 = t.locate(POINT(xr_right, yr_top));
  hface3 = t.locate(POINT(xr_right, yr_bottom));
  hface4 = t.locate(POINT(xr_left, yr_bottom));
  
  l1 = t.line_walk(POINT(xr_left, yr_top), POINT(xr_right, yr_top), hface1);
  l2 = t.line_walk(POINT(xr_right, yr_top), POINT(xr_right, yr_bottom), hface2);
  l3 = t.line_walk(POINT(xr_right, yr_bottom), POINT(xr_left, yr_bottom), hface3);
  l4 = t.line_walk(POINT(xr_left, yr_bottom), POINT(xr_left, yr_top), hface4);

  //test if everything is inside or outside
  if( (l1 == (Nullptr_t) NULL) && (l2 == (Nullptr_t) NULL) &&
      (l3 == (Nullptr_t) NULL) && (l4 == (Nullptr_t) NULL)) 
  {
    FVI v = t.finite_vertices_begin();
    if((*v).point().x() < xr_left || (*v).point().x() > xr_right || 
       (*v).point().y() < yr_bottom || (*v).point().y() > yr_top) //true if everything is outside
      return;
    else{ //everything is inside
      FFI it = t.finite_faces_begin();
      while(it != t.finite_faces_end())
      {
        fct(it);
        it++;
      }
    }
    return;
  }

  //if we are here, then a part of the triangulation is inside, the other is outside

  //put all the faces on the boundaries in the stack and the map
  if(l1 != (Nullptr_t) NULL) //found at least one face that intersect the TOP segment
  {
    while (t.is_infinite(l1)) l1++; //we should start with a finite face
    do{                             //put all of them in the stack;
      face_stack.push(l1);
      visited[l1] = true;
      l1++;
    }while(!t.is_infinite(l1) && 
	   t.triangle(l1).has_on_unbounded_side(POINT(xr_right, yr_top)));
  }
  if(l2 != (Nullptr_t) NULL) //found at least one face that intersect the RIGHT segment
  {
    while (t.is_infinite(l2)) l2++; //we should start with a finite face
    do{                             //put all of them in the stack;
      if(!visited[l2]){
        face_stack.push(l2);
        visited[l2] = true;       
      }
      l2++;
    }while(!t.is_infinite(l2) && 
	   t.triangle(l2).has_on_unbounded_side(POINT(xr_right, yr_bottom)));
  }
  if(l3 != (Nullptr_t) NULL) //found at least one face that intersect the BOTTOM segment
  {
    while (t.is_infinite(l3)) l3++; //we should start with a finite face
    do{                             //put all of them in the stack;
      if(!visited[l3]){
        face_stack.push(l3);
        visited[l3] = true;        
      }
      l3++;
    }while(!t.is_infinite(l3) && 
	   t.triangle(l3).has_on_unbounded_side(POINT(xr_left, yr_bottom)));
  }
  if(l4 != (Nullptr_t) NULL) //found at least one face that intersect the LEFT segment
  {
    while (t.is_infinite(l4)) l4++; //we should start with a finite face
    do{                             //put all of them in the stack;
      if(!visited[l4]){
        face_stack.push(l4);
        visited[l4] = true;
      }
      l4++;
    }while(!t.is_infinite(l4) && 
	   t.triangle(l4).has_on_unbounded_side(POINT(xr_left, yr_top)));
  }
  
  //HERE we begin to walk through the faces DFS
  hFACE hf;
  typename CGAL::Unique_hash_map<hFACE,bool>::Data& 
     data_ref_start(visited[hf]);
  data_ref_start = true;
  while(!face_stack.empty()){
    hf = face_stack.top();
    face_stack.pop();         //done with this face
    for (int i=0; i<3; i++){  //visit all the neighbors
      if(!visited[(*hf).neighbor(i)] )
        if(!t.is_infinite((*hf).neighbor(i))){                //true if it is not an infinite face
          hVERTEX hv = (*(*hf).neighbor(i)).vertex((*(*hf).neighbor(i)).index(hf));
          if(!((*hv).point().x() < xr_left || (*hv).point().x() > xr_right ||
               (*hv).point().y() < yr_bottom || (*hv).point().y() > yr_top)) //true if the vertex is outside
          face_stack.push((*hf).neighbor(i));
          typename CGAL::Unique_hash_map<hFACE,bool>::Data& 
              data_ref(visited[(*hf).neighbor(i)]);
          data_ref = true;
        }
    }
    fct(hf);
  }
}

}//end namespace

#endif
