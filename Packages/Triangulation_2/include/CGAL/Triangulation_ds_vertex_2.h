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
// author(s)     : Mariette Yvinec
//
// coordinator   : Mariette Yvinec  < Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_DS_VERTEX_2_H
#define CGAL_TRIANGULATION_DS_VERTEX_2_H

#include <CGAL/basic.h>
#include <CGAL/Triangulation_short_names_2.h>

CGAL_BEGIN_NAMESPACE

template <class Vb>
class  Triangulation_ds_vertex_2 
  : public Vb 
{
  typedef typename Vb::Triangulation_data_structure Tds;
public:
  typedef typename Tds::Vertex             Vertex;
  typedef typename Tds::Vertex_handle      Vertex_handle;
  typedef typename Tds::Face_handle        Face_handle;
  typedef typename Tds::Vertex_circulator  Vertex_circulator;
  typedef typename Tds::Face_circulator    Face_circulator;
  typedef typename Tds::Edge_circulator    Edge_circulator;
 
  //CREATORS
  Triangulation_ds_vertex_2() : Vb() {}

  //ACCESS
  Vertex_handle handle() {return const_cast<Vertex*>(this); }
  int degree(); //should be const

  // the following should be const
  // when Face_circulator, Vertex_circulator and Edge_circulator
  // are created from 
  // Face_const_handle and Face_const_vertex
  Vertex_circulator incident_vertices()     
    {return Vertex_circulator(this);}
 
  Vertex_circulator incident_vertices( Face_handle f)  
    {return Vertex_circulator(this,f);}
  
  Face_circulator incident_faces()  
    { return Face_circulator(this) ;}
  
  Face_circulator incident_faces( Face_handle f)    
    { return Face_circulator(this, f);}
  
  Edge_circulator incident_edges()   
    { return Edge_circulator(this);}
  
  Edge_circulator incident_edges( Face_handle f)  
    { return Edge_circulator(this, f);}
  
  bool is_valid(bool verbose = false, int level = 0);

};

template <class Tds>
int
Triangulation_ds_vertex_2 <Tds> ::
degree() //const
{
  int count = 0;
  Vertex_circulator vc = incident_vertices(), done(vc);
  if ( ! vc.is_empty()) {
    do { 
      count += 1;
    } while (++vc != done);
  }
  return count;
}


    
template <class Vb>
bool 
Triangulation_ds_vertex_2<Vb> ::  
is_valid(bool verbose, int level) 
{
  bool result = Vb::is_valid(verbose, level);
  CGAL_triangulation_assertion(result);
  if (face() != NULL) { // face==NULL if dim <0
    result = result && face()->has_vertex(handle());
  }
  return result;
}

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_VERTEX_2_H
