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
  int degree(); //should be const

  //Deprecated access to circulators - for bacward compatibility
  // the following should be const
  // when Face_circulator, Vertex_circulator and Edge_circulator
  // are created from 
  // Face_const_handle and Face_const_vertex
  Vertex_circulator incident_vertices()     
    {return Vertex_circulator(handle());}
 
  Vertex_circulator incident_vertices( Face_handle f)  
    {return Vertex_circulator(handle(),f);}
  
  Face_circulator incident_faces()  
    { return Face_circulator(handle()) ;}
  
  Face_circulator incident_faces( Face_handle f)    
    { return Face_circulator(handle(), f);}
  
  Edge_circulator incident_edges()   
    { return Edge_circulator(handle());}
  
  Edge_circulator incident_edges( Face_handle f)  
    { return Edge_circulator(handle(), f);}
  
  bool is_valid(bool verbose = false, int level = 0);

private:
  // used to implement deprected access to circulators
  Vertex_handle handle();
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

template <class Tds>
typename Triangulation_ds_vertex_2<Tds>::Vertex_handle
Triangulation_ds_vertex_2 <Tds> ::
handle()
{
  Face_handle fh = this->face();
  for(int i = 0 ; i < 3 ; ++i){
    if ( &*fh->vertex(i) == this) return fh->vertex(i);
  }
  return Vertex_handle();				    
}
    
template <class Vb>
bool 
Triangulation_ds_vertex_2<Vb> ::  
is_valid(bool verbose, int level) 
{
  bool result = Vb::is_valid(verbose, level);
  CGAL_triangulation_assertion(result);
  if (face() != NULL) { // face==NULL if dim <0
    result = result && ( &*face()->vertex(0) == this ||
			 &*face()->vertex(1) == this ||
			 &*face()->vertex(2) == this );
  }
  return result;
}

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_DS_VERTEX_2_H
