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
// release_date  : 2000, December 10
// 
// source        : aff_transformation_tags.fw
// file          : aff_transformation_tags.h
// package       : Kernel_basic (3.17)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 3.17
// revision_date : 10 Dec 2000 
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_AFF_TRANSFORMATION_TAGS_H
#define CGAL_AFF_TRANSFORMATION_TAGS_H

#ifndef CGAL_CONFIG_H
#include <CGAL/config.h>
#endif // CGAL_CONFIG_H

CGAL_BEGIN_NAMESPACE

class Translation {};
class Rotation {};
class Scaling {};
class Reflection {};
class Identity_transformation {};

extern  Translation              TRANSLATION;
extern  Rotation                 ROTATION;
extern  Scaling                  SCALING;
extern  Reflection               REFLECTION;
extern  Identity_transformation  IDENTITY;

CGAL_END_NAMESPACE

#endif // CGAL_AFF_TRANSFORMATION_TAGS_H
