// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// file          : include/CGAL/IO/Null_output_stream.h
// package       : IO
// source        : $URL$
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================


#ifndef CGAL_NULL_OUTPUT_STREAM_H
#define CGAL_NULL_OUTPUT_STREAM_H


namespace CGAL {


struct Null_output_stream {};

#if defined(__INTEL_COMPILER)
template<class T>
Null_output_stream&
operator<<(Null_output_stream& nos, const T&);
#endif

template<class T>
Null_output_stream&
operator<<(Null_output_stream& nos, const T&)
{
  return nos;
}


} //namespace CGAL

#endif // CGAL_NULL_OUTPUT_STREAM_H
