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
// file          : aff_transformation_tags.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_AFF_TRANSFORMATION_TAGS_H
#define CGAL_AFF_TRANSFORMATION_TAGS_H

#include <CGAL/config.h>

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
