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

template <class Vb, class Fb>
class  Triangulation_ds_vertex_2 ;


template < class Vb, class Fb>
class  Triangulation_ds_face_2
  : public Fb
{
public:
  typedef Vb Vertex_base;
  typedef Fb Face_base;
  typedef Triangulation_ds_vertex_2<Vertex_base,Face_base> Vertex;
  typedef Triangulation_ds_face_2<Vertex_base,Face_base> Face;

  // creators
  Triangulation_ds_face_2()
    : Face_base()
  {}
    
  Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2)
    :  Face_base(v0,v1,v2)
  {}
    
  Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2,
			  Face* n0, Face* n1, Face* n2)
    :  Face_base(v0,v1,v2,n0,n1,n2)
  {}

  Triangulation_ds_face_2( const Face & f)
    : Face_base(f)
    {}

  //setting
  void set_vertex(int i, Vertex* v) { Face_base::set_vertex(i,v);}
  void set_neighbor(int i, Face* n) { Face_base::set_neighbor(i,n);}
  void set_vertices() { Face_base::set_vertices();}
  void set_neighbors() { Face_base::set_neighbors();}
  void set_vertices(Vertex* v0, Vertex* v1, Vertex* v2);
  void set_neighbors(Face* n0, Face* n1, Face* n2);
  //void reorient();  inherited from Face_base
 
  //Vertex Access Member Functions
  Vertex* vertex(int i) const;
  Vertex* mirror_vertex(int i) const;
  bool has_vertex(const Vertex* v) const;
  bool has_vertex(const Vertex* v, int& i) const;
  int index(const Vertex* v) const;

  // Neighbors Access Functions
  Face* neighbor(int i) const;
  bool has_neighbor(const Face* n) const;
  bool has_neighbor(const Face* n, int& i) const;
  int index(const Face* n) const;
  int mirror_index(int i) const;

  //Miscelleanous
  bool is_valid(bool verbose = false, int level = 0) const;
};



template < class Vb, class Fb >
inline void 
Triangulation_ds_face_2<Vb,Fb>::
set_vertices(Vertex* v0, Vertex* v1, Vertex* v2)
{
  Face_base::set_vertices(v0,v1,v2);
}

template < class Vb, class Fb >
inline void 
Triangulation_ds_face_2<Vb,Fb>::
set_neighbors(Face* n0, Face* n1, Face* n2)
{
  Face_base::set_neighbors(n0,n1,n2);
}

template < class Vb, class Fb >
inline
Triangulation_ds_vertex_2<Vb,Fb> *
Triangulation_ds_face_2<Vb,Fb>::
vertex(int i) const
{
  return( static_cast<Vertex*>(Face_base::vertex(i)) );
} 

template < class Vb, class Fb >
inline
Triangulation_ds_vertex_2<Vb,Fb> *
Triangulation_ds_face_2<Vb,Fb>::
mirror_vertex(int i) const
{
  CGAL_triangulation_precondition ( neighbor(i) != NULL);
  return neighbor(i)->vertex(neighbor(i)->index(this));
}

template < class Vb, class Fb >
inline int
Triangulation_ds_face_2<Vb,Fb>::
mirror_index(int i) const
{
  CGAL_triangulation_precondition (neighbor(i) != NULL);
  return neighbor(i)->index(this);
}

template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::
has_vertex(const Vertex* v) const
{
  return (Face_base::has_vertex(v));
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::    
has_vertex(const Vertex* v, int& i) const
{
  return (Face_base::has_vertex(v,i));
}
    
template < class Vb, class Fb >
inline  int 
Triangulation_ds_face_2<Vb,Fb>::  
index(const Vertex* v) const
{
  return(Face_base::vertex_index(v));
}

// Neighbors Access Functions
template < class Vb, class Fb >
inline   
Triangulation_ds_face_2<Vb,Fb>* 
Triangulation_ds_face_2<Vb,Fb>::  
neighbor(int i) const
{
  return (static_cast<Face*>(Face_base::neighbor(i)) );
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::  
has_neighbor(const Face* n) const
{
  return (Face_base::has_neighbor(n));
}
    
template < class Vb, class Fb >
inline  bool 
Triangulation_ds_face_2<Vb,Fb>::      
has_neighbor(const Face* n, int& i) const
{
  return (Face_base::has_neighbor(n,i));
}
    
template < class Vb, class Fb >
inline  int 
Triangulation_ds_face_2<Vb,Fb>::    
index(const Face* n) const
{
  return(Face_base::face_index(n));
}
    
//Miscelleanous
template < class Vb, class Fb >
bool
Triangulation_ds_face_2<Vb,Fb>::  
is_valid(bool verbose, int level) const
{
  bool result = Face_base::is_valid(verbose, level);
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
