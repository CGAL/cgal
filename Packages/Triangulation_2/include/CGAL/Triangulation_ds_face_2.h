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
// release       : $CGAL_Revision: CGAL-2.0-I-12 $
// release_date  : $CGAL_Date: 1999/04/28 $
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

#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE 

template <class Vb, class Fb >
class  Triangulation_ds_vertex_2 ;


template < class Vb, class Fb >
class  Triangulation_ds_face_2
  : public Fb
{
public:
  typedef Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef Triangulation_ds_face_2<Vb,Fb> Face;

  // creators
  Triangulation_ds_face_2()
    : Fb()
  {}
    
  Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2)
    :  Fb(v0,v1,v2)
  {}
    
  Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2,
			  Face* n0, Face* n1, Face* n2)
    :  Fb(v0,v1,v2,n0,n1,n2)
  {}

  Triangulation_ds_face_2( const Face * f);
 
  //setting
  void set_vertex(int i, Vertex* v) { Fb::set_vertex(i,v);}
  void set_neighbor(int i, Face* n) { Fb::set_neighbor(i,n);}
  void set_vertices() { Fb::set_vertices();}
  void set_neighbors() { Fb::set_neighbors();}
  void set_vertices(Vertex* v0, Vertex* v1, Vertex* v2);
  void set_neighbors(Face* n0, Face* n1, Face* n2);
  void reorient();
 
  //Vertex Access Member Functions
  Vertex* vertex(int i) const;
  Vertex* opposite_vertex(int i) const;
  bool has_vertex(const Vertex* v) const;
  bool has_vertex(const Vertex* v, int& i) const;
  int index(const Vertex* v) const;

  // Neighbors Access Functions
  Face* neighbor(int i) const;
  bool has_neighbor(const Face* n) const;
  bool has_neighbor(const Face* n, int& i) const;
  int index(const Face* n) const;
  int opposite_index(int i) const;

  //Miscelleanous
  int dimension() const;
  bool is_valid(bool verbose = false, int level = 0) const;
};


template < class Vb, class Fb >
Triangulation_ds_face_2<Vb,Fb>::
Triangulation_ds_face_2(const Face * f)
  :  Fb()
{
  set_vertices(f->vertex(0), f->vertex(1), f->vertex(2));
  set_neighbors(f->neighbor(0), f->neighbor(1), f->neighbor(2));
}

template < class Vb, class Fb >
inline void 
Triangulation_ds_face_2<Vb,Fb>::
set_vertices(Vertex* v0, Vertex* v1, Vertex* v2)
{
  Fb::set_vertices(v0,v1,v2);
}

template < class Vb, class Fb >
inline void 
Triangulation_ds_face_2<Vb,Fb>::
set_neighbors(Face* n0, Face* n1, Face* n2)
{
  Fb::set_neighbors(n0,n1,n2);
}

template < class Vb, class Fb >
void 
Triangulation_ds_face_2<Vb,Fb>::
reorient()
{
  Vertex* vtemp = vertex(0);
  Face* ftemp   = neighbor(0);
  set_vertex(0, vertex(1)) ; set_vertex(1,vtemp);
  set_neighbor(0, neighbor(1));set_neighbor(1,ftemp);
}

template < class Vb, class Fb >
inline
Triangulation_ds_vertex_2<Vb,Fb> *
Triangulation_ds_face_2<Vb,Fb>::
vertex(int i) const
{
  return( (Vertex*) (Fb::vertex(i)));
} 

template < class Vb, class Fb >
inline
Triangulation_ds_vertex_2<Vb,Fb> *
Triangulation_ds_face_2<Vb,Fb>::
opposite_vertex(int i) const
{
  return neighbor(i)->vertex(neighbor(i)->index(this));
}

template < class Vb, class Fb >
inline int
Triangulation_ds_face_2<Vb,Fb>::
opposite_index(int i) const
{
  return neighbor(i)->index(this);
}

template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::
has_vertex(const Vertex* v) const
{
  return (Fb::has_vertex(v));
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::    
has_vertex(const Vertex* v, int& i) const
{
  return (Fb::has_vertex(v,i));
}
    
template < class Vb, class Fb >
inline  int 
Triangulation_ds_face_2<Vb,Fb>::  
index(const Vertex* v) const
{
  return(Fb::vertex_index(v));
}

// Neighbors Access Functions
template < class Vb, class Fb >
inline   
Triangulation_ds_face_2<Vb,Fb>* 
Triangulation_ds_face_2<Vb,Fb>::  
neighbor(int i) const
{
  return ((Face*) Fb::neighbor(i));
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::  
has_neighbor(const Face* n) const
{
  return (Fb::has_neighbor(n));
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::      
has_neighbor(const Face* n, int& i) const
{
  return (Fb::has_neighbor(n,i));
}
    
template < class Vb, class Fb >
inline  int 
Triangulation_ds_face_2<Vb,Fb>::    
index(const Face* n) const
{
  return(Fb::face_index(n));
}
    
//Miscelleanoustemplate < class Vb, class Fb >
template < class Vb, class Fb >
inline  int 
Triangulation_ds_face_2<Vb,Fb>::   
dimension() const
{
  if (vertex(2) != NULL) {return 2;}
  else return( vertex(1) != NULL ? 1 : 0);
}

template < class Vb, class Fb >
bool
Triangulation_ds_face_2<Vb,Fb>::  
is_valid(bool verbose = false, int level = 0) const
{
  bool result = Fb::is_valid(verbose, level);
  for(int i = 0; i <= dimension(); i++) {
    Face* n = neighbor(i);
    int in = n->index(this);
    result = result && ( this == n->neighbor(in) );
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
