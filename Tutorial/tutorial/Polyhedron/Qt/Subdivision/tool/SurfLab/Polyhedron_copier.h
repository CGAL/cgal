// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_copier.h
// Description   : 
// Creation_date : 21 Feb 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

/** @file Polyhedron_copier.h
*/

#ifndef _POLYHEDRON_COPIER_H_02212002
#define _POLYHEDRON_COPIER_H_02212002

#include <SurfLab/config.h>

//#include <CGAL/Cartesian.h>
//#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

SURFLAB_BEGIN_NAMESPACE

/**
 * Polyhedron_copier is a modifier (check CGAL ref.) creating 
 * a polyhedron from another polyhedron (different type).
 * _P1/_P2 is the target/source type of the polyhedron  
*/
template <class _P1, class _P2> 
class Polyhedron_copier:public CGAL::Modifier_base<typename _P1::HalfedgeDS> {
protected:
  typedef _P1                                          TargetPolyhedron;
  typedef _P2                                          SourcePolyhedron;

  typedef typename TargetPolyhedron::HalfedgeDS        HDS;
  typedef typename TargetPolyhedron::Point_3           Point;
  typedef typename Point::FT                     FT;

  typedef CGAL::Polyhedron_incremental_builder_3<HDS>  PIB;

public:
  ///
  inline Polyhedron_copier() : source(NULL) {}
  ///
  inline Polyhedron_copier(SourcePolyhedron& spoly) : source(&spoly) {}
  ///
  inline void setSourcePolyhedron(SourcePolyhedron& spoly) { source = &spoly; }
  ///
  inline void operator()(HDS& hds) {
    if (source != NULL) {
      hds.clear();
      hds.reserve(source->size_of_vertices(), 
		  source->size_of_halfedges(), 
		  source->size_of_facets());
      PIB pib(hds, true);
      buildPolyhedron(pib);
    }
  }

private:
  inline void buildPolyhedron(PIB& B); // from input stream

protected:
  SourcePolyhedron* source; 
};

template <class _P1, class _P2>
void Polyhedron_copier<_P1, _P2>::buildPolyhedron(PIB& pb) {
  int num_point = source->size_of_vertices();
  int num_facet = source->size_of_facets();

  typename SourcePolyhedron::Vertex_iterator vitr = source->vertices_begin();
  typename SourcePolyhedron::Vertex_iterator v0 = source->vertices_begin();
  typename SourcePolyhedron::Facet_iterator fitr = source->facets_begin();

  pb.begin_surface(num_point, num_facet); {
    for (int i = 0; i < num_point; ++i, ++vitr) 
      pb.add_vertex(vitr->point());	
    for (int i = 0; i < num_facet; ++i, ++fitr) {
      pb.begin_facet(); {
	typename SourcePolyhedron::Halfedge_around_facet_circulator 
	  cir = fitr->facet_begin();
	int valance = CGAL::circulator_size(cir);
	for (int n = 0; n < valance; ++n, ++cir)
	  pb.add_vertex_to_facet(std::distance(v0, cir->vertex()));
      }
      pb.end_facet();
    }
  }
  pb.end_surface();
}

SURFLAB_END_NAMESPACE

#endif //_POLYHEDRON_COPIER_H_02212002
