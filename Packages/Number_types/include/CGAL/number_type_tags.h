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
// release       : $CGAL_Revision: CGAL-2.4-I-65 $
// release_date  : $CGAL_Date: 2002/03/19 $
// 
// file          : include/CGAL/number_type_tags.h
// package       : Number_types (4.46)
// maintainer    : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_NUMBER_TYPE_TAGS_H
#define CGAL_NUMBER_TYPE_TAGS_H

#ifndef CGAL_ENUM_H
#include <CGAL/enum.h>
#endif // CGAL_ENUM_H



CGAL_BEGIN_NAMESPACE

struct No_number_tag{};
struct Number_tag{};
struct Quotient_tag{};
struct Ring_number_type_tag {};
struct Euclidean_ring_number_type_tag {};
struct Field_number_type_tag {};

CGAL_END_NAMESPACE


#endif // CGAL_NUMBER_TYPE_TAGS_H
