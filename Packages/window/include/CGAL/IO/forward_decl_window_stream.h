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
// file          : include/CGAL/IO/forward_decl_window_stream.h
// package       : window (2.8.0)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.8.0
// revision_date : 17 May 2001 
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 
#ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
#define CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H

#if defined(__USE_CGAL_WINDOW__)
#include <CGAL/LEDA/window.h>

class CGAL::window;
#else
class leda_window;
#endif

CGAL_BEGIN_NAMESPACE

#if defined(__USE_CGAL_WINDOW__)
typedef   CGAL::window    Window_stream;
#else
typedef   leda_window    Window_stream;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
