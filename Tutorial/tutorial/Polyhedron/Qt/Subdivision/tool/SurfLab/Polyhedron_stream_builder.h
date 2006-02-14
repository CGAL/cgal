// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_stream_builder.h
// Description   : Build the polyhedron from the input stream 
//                 (.10/.mesh format)
// Creation_date : 28 Jan 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id$

/** @file Polyhedron_stream_builder.h
*/

#ifndef _POLYHEDRON_STREAM_BUILDER_H_01282002
#define _POLYHEDRON_STREAM_BUILDER_H_01282002

#include <SurfLab/config.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>

#include <strstream>

SURFLAB_BEGIN_NAMESPACE

// ======================================================================
/**
 * Polyhedron_stream_builder is a modifier (check CGAL ref.) creating 
 * a polyhedron from a input stream that has .10/.mesh format.
*/
template <class _Poly> // P should be a class type of polyhedron  
class Polyhedron_stream_builder : 
  public CGAL::Modifier_base<typename _Poly::HalfedgeDS> {
protected:
  typedef _Poly                                        Polyhedron;
  typedef typename Polyhedron::HalfedgeDS              HDS;
  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Point::FT                           FT;

  typedef CGAL::Polyhedron_incremental_builder_3<HDS>  PIB;

public:
  enum FileFormat { FF_TEN = 1, FF_OFF};

public:
  ///
  Polyhedron_stream_builder(FileFormat f = FF_TEN): 
    mpInputFile(NULL), mFormat(f) {}
  ///
  Polyhedron_stream_builder(std::istream& in, FileFormat f = FF_TEN): 
    mpInputFile(&in), mFormat(f) {}
  ///
  void setInStream(std::istream& in, FileFormat f = FF_TEN) {
    mpInputFile = &in;
    mFormat = f;
  }
  ///
  void operator()(HDS& hds) {
    if (mpInputFile != NULL) {
      PIB pib(hds, true);
      switch (mFormat) {
      case FF_TEN:
	buildPolyhedron_TEN(pib); break;
      case FF_OFF:
	buildPolyhedron_OFF(pib); break;
      }
    }
  }

private:
  ///
  void buildPolyhedron_TEN(PIB& B); // from input stream of .10 format
  ///
  void buildPolyhedron_OFF(PIB& B); // from input stream of .OFF format

protected:
  ///
  std::istream* mpInputFile;

  ///
  FileFormat mFormat;
};


#define inf (*mpInputFile)

// ----------------------------------------------------------------------
template <class HDS>
void Polyhedron_stream_builder<HDS>::buildPolyhedron_TEN(PIB& pb) {
  int numV(0), idxV(0), numF(0), degreeF(0);
  FT* xyz = NULL;

  inf >> numV;
  xyz = new FT[(numV+1)*3];
  for (int i = 0; i < numV; i++) {
    inf >> idxV;
    inf >> xyz[idxV*3+0] >> xyz[idxV*3+1] >> xyz[idxV*3+2]; 
  }

  inf >> numF;
  pb.begin_surface( numV, numF); {
    for (int i = 1; i <= numV; i++) 
      pb.add_vertex(Point(xyz[i*3+0], xyz[i*3+1], xyz[i*3+2]));
    for (int i = 0; i < numF; i++) {
      inf >> degreeF;
      pb.begin_facet(); {
	for (int n = 0; n < degreeF; n++) {
	  inf >> idxV; 
	  // .10 file starts the vertex index as 1
	  // but CGAL starts it as 0 !
	  pb.add_vertex_to_facet(idxV-1);
	}		
      }
      pb.end_facet();
    }
  }
  pb.end_surface();

  pb.remove_unconnected_vertices();

  delete[] xyz; 
}

// ----------------------------------------------------------------------
template <class HDS>
void Polyhedron_stream_builder<HDS>::buildPolyhedron_OFF(PIB& pb) {
  int numV(0), numF(0), numE(0), idxV, degreeF;

  char str[256];
  inf >> str; // should be "OFF"
  
  inf >> numV >> numF >> numE;

  pb.begin_surface( numV, numF); {
    FT x, y, z;
    for (int i = 0; i < numV; i++) {
      inf >> x >> y >> z; 
      pb.add_vertex(Point(x, y, z));
    }
    for (int i = 0; i < numF; i++) {
      inf >> degreeF;
      pb.begin_facet(); {
	for (int n = 0; n < degreeF; n++) {
	  inf >> idxV; 
	  pb.add_vertex_to_facet(idxV);
	}
	inf.getline(str, 256);
      }
      pb.end_facet();
    }
  }
  pb.end_surface();

  pb.remove_unconnected_vertices();
}

#undef inf

SURFLAB_END_NAMESPACE

#endif //_POLYHEDRON_STREAM_BUILDER_H_01282002
