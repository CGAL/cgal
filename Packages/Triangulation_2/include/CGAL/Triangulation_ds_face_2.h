// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : Triangulation/include/CGAL/Triangulation_ds_face_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_FACE_2_H
#define CGAL_TRIANGULATION_DS_FACE_2_H

#include <CGAL/Triangulation_short_names_2.h>

template <class Vb, class Fb >
class  CGAL_Triangulation_ds_vertex_2 ;



template < class Vb, class Fb >
class  CGAL_Triangulation_ds_face_2
  : public Fb
{
public:
  //typedef typename Fb::Triangle Triangle;
  typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;

  // creators
  CGAL_Triangulation_ds_face_2()
    : Fb()
  {}
    
  CGAL_Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2)
    :  Fb(v0,v1,v2)
  {}
    
  CGAL_Triangulation_ds_face_2(Vertex* v0, Vertex* v1, Vertex* v2,
				Face* n0, Face* n1, Face* n2)
    :  Fb(v0,v1,v2,n0,n1,n2)
  {}

   CGAL_Triangulation_ds_face_2( const Face * f)
    :  Fb()
  {
    set_vertices(f->vertex(0), f->vertex(1), f->vertex(2));
    set_neighbors(f->neighbor(0), f->neighbor(1), f->neighbor(2));
  }

  //setting
  inline 
  void set_vertex(int i, Vertex* v)
  {
    Fb::set_vertex(i,v);
  }
    
    
  inline 
   void set_neighbor(int i, Face* n)
  {
    Fb::set_neighbor(i,n);
  }

  inline
  void set_vertices() 
  {
    Fb::set_vertices();
  }
      
  inline 
  void set_vertices(Vertex* v0,
		    Vertex* v1,
		    Vertex* v2)
  {
    Fb::set_vertices(v0,v1,v2);
   }
    
  inline
  void set_neighbors() 
  {
    Fb::set_neighbors();
  }
     
  inline
  void set_neighbors(Face* n0,
		     Face* n1,
		     Face* n2)
  {
    Fb::set_neighbors(n0,n1,n2);
  }

  //Vertex Access Member Functions
  Vertex* vertex(int i) const
  {
    return( (Vertex*) (Fb::vertex(i)));
  } 

 inline 
 bool has_vertex(const Vertex* v) const
  {
    return (Fb::has_vertex(v));
  }
    
    
  inline 
  bool has_vertex(const Vertex* v, int& i) const
  {
    return (Fb::has_vertex(v,i));
  }
    
  inline 
  int index(const Vertex* v) const
  {
    return(Fb::vertex_index(v));
  }

  // Neighbors Access Functions
  inline 
  Face* neighbor(int i) const
  {
    return ((Face*) Fb::neighbor(i));
  }
    
  inline 
  bool has_neighbor(const Face* n) const
  {
    return (Fb::has_neighbor(n));
  }
    
    
  inline 
  bool has_neighbor(const Face* n, int& i) const
  {
    return (Fb::has_neighbor(n,i));
  }
    
    
  inline 
  int index(const Face* n) const
  {
    return(Fb::face_index(n));
  }
    
  //Miscelleanous
  int dimension()
  {
    if (vertex(2) != NULL) {return 2;}
    else return( vertex(1) != NULL ? 1 : 0);
  }



  //Additionnal Operations

  //the following function has been moved to the tds class
//   void insert_in_face(Vertex*& v)

//   void insert_in_edge(const Vertex* v, int i)

//   bool insert_outside(const Vertex* v, int i)
//   bool remove(Vertex* v)

//   void flip(int i)


   bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Fb::is_valid();
    for(int i = 0; i < 3; i++) {
      Face* n = neighbor(i);
            
      // The following seems natural, but it may fail if the faces
      // this and n are neighbors on two edges (1-dim triangulation,
      // with infinite faces
      // int ni = n->index(this);

      //  int ni = cw(n->index(vertex(cw(i))));
      // CGAL_triangulation_assertion( this == n->neighbor(ni) );
      // result = result && (vertex(cw(i)) == n->vertex(ccw(ni)));
      // result = result && (vertex(ccw(i)) == n->vertex(cw(ni)));

      int in;
      if (! n->has_vertex(vertex(cw(i)),in )) return false;
      in = cw(in); 
      result = result && ( this == n->neighbor(in) );
      result = result && (vertex(ccw(i)) == n->vertex(cw(in)));

    }
    return result;
  }
   

};

#endif CGAL_TRIANGULATION_DS_FACE_2_H
