// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Kernel_traits_common.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

// This file is intentionally not protected for re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

#define CGAL_Kernel_pred(Y,Z) typedef CGALi::Y<R> Y; Y Z() const {return Y();}
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

#include <CGAL/Kernel/interface_macros.h>
