// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h
// package       : Nef_S2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef NEF_POLYHEDRON_S2_OGLUT_STREAM_H
#define NEF_POLYHEDRON_S2_OGLUT_STREAM_H

#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/Nef_S2/SM_visualizor.h>
#include <string>

CGAL_BEGIN_NAMESPACE


struct OGLUT_stream { // dummy class
  OGLUT_stream() {}
  void display() { CGAL::OGL::start_viewer(); }
};

static OGLUT_stream ogl;

template <typename T>
CGAL::OGLUT_stream& operator<<(CGAL::OGLUT_stream& ogls, 
			       const Nef_polyhedron_S2<T>& P)
{
  typedef Nef_polyhedron_S2<T> Polyhedron;
  typedef typename Polyhedron::Sphere_map Sphere_map;
  typedef typename Polyhedron::Sphere_kernel Sphere_kernel;
  typedef CGAL::SM_visualizor<Sphere_map,Sphere_kernel> Visualizor;
  CGAL::OGL::add_sphere();
  Visualizor V(P.sphere_map(),CGAL::OGL::spheres_.back()); V.draw_map();
  // CGAL::OGL::spheres_.back().print();
  return ogls;
}

static CGAL::OGLUT_stream& operator<<(CGAL::OGLUT_stream& ogls, 
			              const char* s)
{
  CGAL::OGL::titles_.push_back(std::string(s));
  return ogls;
}

CGAL_END_NAMESPACE

#endif // NEF_POLYHEDRON_S2_OGLUT_STREAM_H


