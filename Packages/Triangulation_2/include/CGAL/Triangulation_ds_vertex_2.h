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
// file          : Triangulation/include/CGAL/Triangulation_ds_vertex_2.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_VERTEX_2_H
#define CGAL_TRIANGULATION_DS_VERTEX_2_H

#include <utility>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

CGAL_BEGIN_NAMESPACE

template <class Vb, class Fb >
class  Triangulation_ds_vertex_2 
  : public Vb,
    public Triangulation_cw_ccw_2
{
public:
  typedef typename Vb::Point Point;
  typedef Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef Triangulation_ds_face_2<Vb,Fb> Face;
  typedef std::pair< Face*,int> Edge;
  typedef Triangulation_ds_face_circulator_2<Vertex,Face> Face_circulator;
  typedef Triangulation_ds_vertex_circulator_2<Vertex,Face> 
                                                           Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_2<Vertex,Face> Edge_circulator;

  Triangulation_ds_vertex_2()
    : Vb()
  {}
    
  Triangulation_ds_vertex_2(const Point & p)
    :  Vb(p)
  {}
    
  Triangulation_ds_vertex_2(const Point & p, Face * f)
    :  Vb(p, f )
  {}

  // set_point()
  // point()
  //inherited from Vb

  inline 
  void set_face(Face* f)
  {
    Vb::set_face(f);
  }

  inline Face* face() const
  {
    return ( (Face *) (Vb::face()) );
  }
    
  int degree() const
  {
    Face* f = face();
    
    if (f == NULL) {
      return 0;
    }
    int i = f->index(this);
    
    Face* ptr1 = f->neighbor(ccw(i));
    Face* ptr2 = f;
    f = f->neighbor(cw(i));
    
    int count = 2;
    while(ptr1 != f){
      count++;
      i = ptr1->index(ptr2);
      ptr2 = ptr1;
      ptr1 = ptr1->neighbor(cw(i));
    }
    return count;
  }
  
  inline Vertex_circulator incident_vertices() 
  {
    return Vertex_circulator(this, face());
  }
    
  inline Face_circulator incident_faces() 
  {
    return Face_circulator(this, face());
  }
    
  inline Face_circulator incident_faces(const Face* f) 
  {
    return Face_circulator(this, f);
  }
    
  inline Edge_circulator incident_edges() 
  {
    return Edge_circulator(this, face());
  }
       
  inline Edge_circulator incident_edges(const Face* f) 
  {
    return Edge_circulator(this, f);
  }
    
   bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Vb::is_valid();
    CGAL_triangulation_assertion(result);
    CGAL_triangulation_assertion(face() != NULL);
    result = result && face()->has_vertex(this);
    return result;
  }
};

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_VERTEX_2_H
