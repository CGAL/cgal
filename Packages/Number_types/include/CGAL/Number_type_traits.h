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
// file          : include/CGAL/Number_type_traits.h
// package       : Number_types (4.46)
// maintainer    : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert, Michael Hoffmann
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_NUMBER_TYPE_TRAITS_H
#define CGAL_NUMBER_TYPE_TRAITS_H

CGAL_BEGIN_NAMESPACE

template < class NT >
struct Number_type_traits {
  typedef typename NT::Has_gcd_tag       Has_gcd_tag;
  typedef typename NT::Has_division_tag  Has_division_tag;
  typedef typename NT::Has_sqrt_tag      Has_sqrt_tag;
};

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_TRAITS_H
