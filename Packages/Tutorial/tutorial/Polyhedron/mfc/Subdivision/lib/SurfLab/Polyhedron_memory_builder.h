// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_memory_builder.h
// Description   : Build the polyhedron from the input stream 
//                 (.10/.mesh format)
// Creation_date : 31 Jan 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

/** @file Polyhedron_memory_builder.h
*/

#ifndef _POLYHEDRON_MEMORY_BUILDER_H_01312002
#define _POLYHEDRON_MEMORY_BUILDER_H_01312002

#include "config.h"

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <strstream>

SURFLAB_BEGIN_NAMESPACE

/**
 * Polyhedron_memory_builder is a modifier (check CGAL ref.) creating 
 * a polyhedron from a input stream that has .10/.mesh format.
*/
template <class _P> // P should be a class type of polyhedron  
class Polyhedron_memory_builder : 
  public CGAL::Modifier_base<typename _P::HalfedgeDS> {
protected:
  typedef _P                                           Polyhedron;
  typedef typename Polyhedron::HalfedgeDS              HDS;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Polyhedron::FT                           FT;

  typedef CGAL::Polyhedron_incremental_builder_3<HDS>  PIB;

public:
  ///
  inline Polyhedron_memory_builder() : 
    num_point(0), point_buffer(NULL), num_facet(0), facet_buffer(NULL) {}
  ///
  inline Polyhedron_memory_builder(int nump, FT* pb, int numf, int** fb) : 
    num_point(nump), point_buffer(pb), num_facet(numf), facet_buffer(fb) {}
  ///
  inline void setInputMemory(int nump, FT* pb, int numf, int** fb) { 
    num_point = nump;
    point_buffer(pb);
    num_facet(numf); 
    facet_buffer(fb);
  }
  ///
  inline void operator()(HDS& hds) {
    if (point_buffer != NULL && facet_buffer != NULL) {
      PIB pib(hds, true);
      buildPolyhedron(pib);
    }
  }

private:
  inline void buildPolyhedron(PIB& B); // from input stream

protected:
  int num_point;
  FT* point_buffer;
  int num_facet;
  int** facet_buffer;
};

template <class HDS>
void Polyhedron_memory_builder<HDS>::buildPolyhedron(PIB& pb) {
  pb.begin_surface(num_point, num_facet); {
    for (int i = 0; i < num_point; ++i) 
      pb.add_vertex(Point(point_buffer[i*3+0], 
			  point_buffer[i*3+1], 
			  point_buffer[i*3+2]));	
    for (int i = 0; i < num_facet; ++i) {
      pb.begin_facet(); {
	for (int n = 0; n < facet_buffer[i][0]; ++n)
	  pb.add_vertex_to_facet(facet_buffer[i][n+1]);
      }
      pb.end_facet();
    }
  }
  pb.end_surface();
}

SURFLAB_END_NAMESPACE

#endif //_POLYHEDRON_MEMORY_BUILDER_H_01312002
