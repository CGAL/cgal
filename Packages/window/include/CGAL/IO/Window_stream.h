// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// source        : leda_window.fw
// file          : include/CGAL/IO/Window_stream.h
// package       : window (2.8.0)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.8.0
// revision_date : 17 May 2001 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
// use the port of the LEDA window (CGAL::window) or the "normal" LEDA window ... 
 
#if defined(__USE_CGAL_WINDOW__)
#include <CGAL/IO/cgal_window.h>
#else
#include <CGAL/IO/leda_window.h>
#endif

#ifndef IO_TRIANGULATION_WINDOW_STREAM_H
#include <CGAL/IO/triangulation_Window_stream.h>
#endif  // IO_TRIANGULATION_WINDOW_STREAM_H
#ifndef IO_OPTIMISATION_WINDOW_STREAM_H
#include <CGAL/IO/optimisation_Window_stream.h>
#endif // IO_OPTIMISATION_WINDOW_STREAM_H
#ifndef IO_POLYGON_WINDOW_STREAM_H
#include <CGAL/IO/polygon_Window_stream.h>
#endif // IO_POLYGON_WINDOW_STREAM_H

