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
// release_date  : 2000, December 13
// 
// source        : Double.fw
// file          : src/Double.C
// package       : Number_types (4.2)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 4.2
// revision_date : 13 Dec 2000 
// author(s)     : Stefan Schirra <stschirr@mpi-sb.mpg.de>
//                 Geert-Jan Giezeman
//                 Sylvain Pion
//
// coordinator   : MPI Informatik, Saarbruecken
//
// ======================================================================


#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE


#ifdef CGAL_OLD_FINITE_VALID

#if !defined(__sgi) && !defined(__sun) && !defined(__hpux) && !defined(__linux)

bool is_valid(double d)
{
    return (d == d);             /* !!! */
}

bool is_finite(double d)
{
    return (d == d) && (CGAL_is_valid(d-d));
}

#else  // custom definitions.

#ifdef __sgi

// implementation for SGI IRIX 5.3.
#include <fp_class.h>

bool is_finite(double d)
{
    switch (fp_class_d(d)) {
    case FP_POS_NORM:
    case FP_NEG_NORM:
    case FP_POS_ZERO:
    case FP_NEG_ZERO:
    case FP_POS_DENORM:
    case FP_NEG_DENORM:
        return true;
    case FP_SNAN:
    case FP_QNAN:
    case FP_POS_INF:
    case FP_NEG_INF:
        return false;
    }
    return false; // NOT REACHED
}

bool is_valid(double d)
{
    switch (fp_class_d(d)) {
    case FP_POS_NORM:
    case FP_NEG_NORM:
    case FP_POS_ZERO:
    case FP_NEG_ZERO:
    case FP_POS_INF:
    case FP_NEG_INF:
    case FP_POS_DENORM:
    case FP_NEG_DENORM:
        return true;
    case FP_SNAN:
    case FP_QNAN:
        return false;
    }
    return false; // NOT REACHED
}

#endif // __sgi

#ifdef __hpux

// implementation for HP

bool is_valid(double d)
{
    return isnan(d) == 0;
}

bool is_finite(double d)
{
    switch (fpclassify(d)) {
    case FP_PLUS_NORM:
    case FP_MINUS_NORM:
    case FP_PLUS_ZERO:
    case FP_MINUS_ZERO:
    case FP_PLUS_DENORM:
    case FP_MINUS_DENORM:
        return true;
    case FP_PLUS_INF:
    case FP_MINUS_INF:
    case FP_SNAN:
    case FP_QNAN:
        return false;
    }
    return false; // NOT REACHED
}

#endif // __hpux

#ifdef __sun

// implementation for SUN

#ifdef __SVR4
#include <ieeefp.h>
#endif // __SVR4

#ifdef __svr4__
#include <ieeefp.h>
#endif //  __svr4__

#include <CGAL/config.h>
bool is_finite(double d)
{
  if(finite(d)){
    return true;
  }
  return false; // NOT REACHED
}

bool is_valid(double d)
{
  return isnan(d) == 0;
}

#endif // __sun

#ifdef __linux

// Implementation for Linux

bool is_finite(double d)
{
  return finite(d) != 0;
}

bool is_valid(double d)
{
  return isnan(d) == 0;
}

#endif // __linux
#endif // of custom definitions
#endif // CGAL_OLD_FINITE_VALID

CGAL_END_NAMESPACE


