// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : Triangulation/include/CGAL/Triangulation_ds_face_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================

#ifndef CGAL_TRIANGULATION_DS_FACE_2_H
#define CGAL_TRIANGULATION_DS_FACE_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>

CGAL_BEGIN_NAMESPACE 

template < class Fb >
class  Triangulation_ds_face_2
  : public Fb  
{
public:
  typedef typename Fb::Triangulation_data_structure Tds;
  typedef typename Tds::Vertex             Vertex;
  typedef typename Tds::Face               Face;
  typedef typename Tds::Vertex_handle      Vertex_handle;
  typedef typename Tds::Face_handle        Face_handle;

public :
  // creators
  Triangulation_ds_face_2()
    : Fb()
  {}
    
  Triangulation_ds_face_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2)
    :  Fb(v0,v1,v2)
  {}
    
  Triangulation_ds_face_2(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2,
			  Face_handle n0, Face_handle n1, Face_handle n2)
    :  Fb(v0,v1,v2,n0,n1,n2)
  {}

  Triangulation_ds_face_2( const Face& f)
    : Fb(f)
    {}

  static int ccw(int i) {return Triangulation_cw_ccw_2::ccw(i);}
  static int  cw(int i) {return Triangulation_cw_ccw_2::cw(i);}

  Vertex_handle mirror_vertex(int i) const;
  int mirror_index(int i) const;

  Face_handle handle() const {return const_cast<Face*>(this);}
  Face_handle handle() {return const_cast<Face*>(this);}
  bool is_valid(bool verbose = false, int level = 0) const;

};

template < class Fb >
inline
typename Triangulation_ds_face_2<Fb>::Vertex_handle
Triangulation_ds_face_2<Fb>::
mirror_vertex(int i) const
{
  CGAL_triangulation_precondition ( neighbor(i) != NULL &&  dimension() >= 1);
  //return neighbor(i)->vertex(neighbor(i)->index(this->handle()));
  return neighbor(i)->vertex(mirror_index(i));
}

template < class Fb >
inline int
Triangulation_ds_face_2<Fb>::
mirror_index(int i) const
{
  // return the index of opposite vertex in neighbor(i);
  CGAL_triangulation_precondition (neighbor(i) != NULL &&  dimension() >= 1);
  if (dimension() == 1) return neighbor(i)->index(this->handle());
  return ccw( neighbor(i)->index(vertex(ccw(i))));
}


template < class Fb >
bool
Triangulation_ds_face_2<Fb>::  
is_valid(bool verbose, int level) const
{
  bool result = Fb::is_valid(verbose, level);
  for(int i = 0; i <= dimension(); i++) {
    Face_handle n = neighbor(i);
    // the strange formulation in case dimension()==2 
    // is used to handle the cases of TDS allowing
    // two faces with two common edges
    int in;
    if (dimension() == 0) in = n->index(this->handle());
    else in = n->index(mirror_vertex(i));
    result = result && ( this->handle() == n->neighbor(in) );
    switch(dimension()) {
    case 0 : 
      break;
    case 1 :
      result = result &&  in == 1-i;
      result = result && ( this->vertex(1-i) == n->vertex(1-in));
      break;
    case 2 :
      result = result && ( this->vertex(cw(i)) == n->vertex(ccw(in)))
	              && ( this->vertex(ccw(i)) == n->vertex(cw(in)));
      break;
    }
  }
  return result;
}


CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_FACE_2_H
