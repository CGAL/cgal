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
// file          : misc.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_MISC_H
#define CGAL_MISC_H

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
// A helper class:
// ---------------------
template <class Target, class Source>
struct converter
{
    static inline Target do_it(const Source& s)
    { return static_cast<Target>(s); }
};

template <class Target, class Source>
inline
Target
convert_to (const Source& s)
{ return converter<Target, Source>::do_it(s); }

/*
template <class Target, class Source>
inline
Target
convert_to( const Source& s)
{ return Target(s); }
*/
#endif // CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION

template <class Target, class Source>
inline
Target
convert_from_to( const Target& t, const Source& s)
{ return Target(s); }

CGAL_END_NAMESPACE

#endif // CGAL_MISC_H
