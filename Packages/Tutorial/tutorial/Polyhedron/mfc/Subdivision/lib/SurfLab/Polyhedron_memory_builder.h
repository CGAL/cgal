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

#include <CGAL/Cartesian.h>
//#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <sstream>

/**
 * Polyhedron_memory_builder is a modifier (check CGAL ref.) creating 
 * a polyhedron from a input stream that has .10/.mesh format.
*/
template <class _Poly> // P should be a class type of polyhedron  
class Polyhedron_memory_builder : 
  public CGAL::Modifier_base<typename _Poly::HalfedgeDS> {
protected:
  typedef _Poly                                        Polyhedron;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::HalfedgeDS              HDS;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Kernel::FT                          FT;

  typedef CGAL::Polyhedron_incremental_builder_3<HDS>  PIB;

public:
  ///
  Polyhedron_memory_builder() : 
    num_point(0), point_buffer(NULL), num_facet(0), facet_buffer(NULL) {}
  ///
  Polyhedron_memory_builder(int nump, FT* pb, int numf, int** fb) :
    num_point(nump), point_buffer(pb), num_facet(numf), facet_buffer(fb), 
    c_pointer(true) {}
  ///
  Polyhedron_memory_builder(int nump, Point* pb, int numf, int** fb) : 
    num_point(nump), point_buffer(pb), num_facet(numf), facet_buffer(fb),
    c_pointer(false) {}

  ///
  inline void setInputMemory(int nump, FT* pb, int numf, int** fb) { 
    num_point = nump;    point_buffer = pb;
    num_facet = numf;    facet_buffer = fb;
    c_pointer = true;
  }
  ///
  inline void setInputMemory(int nump, Point* pb, int numf, int** fb) { 
    num_point = nump;    point_buffer = pb;
    num_facet = numf;    facet_buffer = fb;
    c_pointer = false;
  }

  ///
  inline void operator()(HDS& hds) {
    if (point_buffer != NULL && facet_buffer != NULL) {
      PIB pib(hds, true);
      c_pointer ? buildPolyhedron_c(pib) : buildPolyhedron_pt(pib);
    }
  }

private:
  inline void buildPolyhedron_c(PIB& B);
  inline void buildPolyhedron_pt(PIB& B);

protected:
  int num_point;
  void* point_buffer;
  int num_facet;
  int** facet_buffer;

  bool c_pointer;
};

template <class HDS>
void Polyhedron_memory_builder<HDS>::buildPolyhedron_c(PIB& pb) {
  FT* p = (FT*) point_buffer;
  pb.begin_surface(num_point, num_facet); {
    for (int i = 0; i < num_point; ++i) 
      pb.add_vertex(Point(p[i*3+0], p[i*3+1], p[i*3+2]));	
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

template <class HDS>
void Polyhedron_memory_builder<HDS>::buildPolyhedron_pt(PIB& pb) {
  Point* p = (Point*) point_buffer;
  pb.begin_surface(num_point, num_facet); {
    for (int i = 0; i < num_point; ++i) pb.add_vertex(p[i]);	
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

#endif //_POLYHEDRON_MEMORY_BUILDER_H_01312002
