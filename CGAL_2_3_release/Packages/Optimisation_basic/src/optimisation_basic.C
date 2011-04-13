// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : src/optimisation_basic.C
// package       : $CGAL_Package: Optimisation_basic $
// chapter       : Geometric Optimisation
//
// source        : web/Optimisation_basic.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: basic things for optimisation algorithms
// ============================================================================

#include <CGAL/Optimisation/basic.h>

CGAL_BEGIN_NAMESPACE

// Function implementations
// ========================

// is_valid failure function
// -------------------------
bool
_optimisation_is_valid_fail( CGAL::Verbose_ostream& verr,
                             const char*            message)
{
    verr << "FAILED." << std::endl;
    verr << "  --> " << message << std::endl;
    verr << "  object is NOT valid!" << std::endl;
    return( false);
}

CGAL_END_NAMESPACE

// ===== EOF ==================================================================
