// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// package       : $CGAL_Package: MinSphere $
// file          : min_sphere_traits_3_test.C
// source        : web/Optimisation/Min_sphere_d.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
//                 Bernd Gärtner
//
// coordinator   : ETH Zurich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: dD Smallest Enclosing Sphere
// ============================================================================
// hack to overcome internal compiler error of egcs
#define FT_(x) FT(x,1)
#define RT_(x) RT(x,1)

// hack to overcome external linkage conflict warning in MIPS
#define __SGI_STL_INTERNAL_RELOPS

#include <CGAL/basic.h>

// only do something if LEDA available
#ifdef CGAL_USE_LEDA

#include <cassert>
#include<CGAL/Random.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Optimisation_d_traits_3.h>
#include<CGAL/Min_sphere_d.h>
#include<CGAL/leda_rational.h>
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_rational);
using namespace CGAL;

typedef leda_rational                               NT;
typedef NT                                          FT;
typedef NT                                          RT;
typedef Cartesian<FT>                               C;
typedef Homogeneous<RT>                     H;
typedef Optimisation_d_traits_3<C>          Cartesian_traits;
typedef Optimisation_d_traits_3<H>          Homogeneous_traits;
typedef Min_sphere_d<Cartesian_traits>      Min_sphereC;
typedef Min_sphere_d<Homogeneous_traits>    Min_sphereH;
typedef C::Point_3                          PointC;
typedef H::Point_3                          PointH;


// checked traits
Cartesian_traits                                    tC;
Homogeneous_traits                                  tH;

const int n = 10;                           // number of points
const int D = 3;
const int r = 128;                          // coordinate range
bool verbose = true;

// NT == Quotient<NT> ?
bool equals (NT x, Quotient<NT> y)
{
    return (Quotient<NT>(x) == y);
}

// PointC == Point H ?
bool equals (PointC p, PointH q)
{
    int dim =  p.dimension();
    if (q.dimension() != dim) return false;
    for (int j=0; j<dim; ++j)
       if (!equals(p.cartesian(j),q.cartesian(j)))
            return false;
    return true;
}

int main ()
{
    Point_3<C>      PC[n];          // cartesian point set
    Point_3<H>      PH[n];          // homogeneous point set
    int             coord[D];       // arrays for random coordinates
    FT              coordC[D];      // ... cartesian version
    RT              coordH[D+1];    // ... homogeneous version
    Random          my_random;              // random number generator


    for (int i=0; i<n; ++i) {
        int j;
        // random coordinates
        for (j=0; j<D; ++j)
            coord[j] = my_random (r);

        // cartesian point
        for (j=0; j<D; ++j)
            coordC[j] = FT (coord[j]);
        PC[i] = Point_3<C>(coordC[0], coordC[1], coordC[2]);

        // homogeneous point
        for (j=0; j<D; ++j)
            coordH[j] = RT(2*coord[j]);
        coordH[D] = RT(2);
        PH[i] = Point_3<H>(coordH[0], coordH[1], coordH[2], coordH[3]);
    }

    Min_sphereC     msC (PC, PC+n, tC); assert (msC.is_valid(verbose));
    Min_sphereH     msH (PH, PH+n, tH); assert (msH.is_valid(verbose));

    PointC                  centerC (msC.center());
    PointH                  centerH (msH.center());

    FT                      radiusC (msC.squared_radius());
    Quotient<RT>    radiusH (msH.squared_radius());

    assert (equals (centerC, centerH));
    assert (equals (radiusC, radiusH));

    return 0;
}

#else

int main()
{
  return 0;
}

#endif



// ===== EOF ==================================================================

