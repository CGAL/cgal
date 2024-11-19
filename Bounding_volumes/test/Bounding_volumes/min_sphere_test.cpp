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
// chapter       : $CGAL_Chapter: Geometric Optimization $
// package       : $CGAL_Package: MinSphere $
// file          : min_sphere_test.C
// source        : web/Optimisation/Min_sphere_d.aw
// revision      : $Id$
// revision_date : $Date$
// author(s)     : Sven Schönherr <sven@inf.fu-berlin.de>
//                 Bernd Gärtner
//
// coordinator   : ETH Zurich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: dD Smallest Enclosing Sphere
// ============================================================================

// hack to overcome external linkage conflict warning in MIPS
#define __SGI_STL_INTERNAL_RELOPS


#include <CGAL/Exact_rational.h>
#include <CGAL/Random.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Min_sphere_annulus_d_traits_d.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Quotient.h>

#include <cassert>
#include <sstream>

using namespace CGAL;

// leda_rational, or Gmpq, or Quotient<MP_float>
typedef CGAL::Exact_rational                NT;
typedef NT                                  FT;
typedef NT                                  RT;
typedef Cartesian_d<FT>                     C;
typedef Homogeneous_d<RT>                   H;
typedef Min_sphere_annulus_d_traits_d<C>    Cartesian_traits;
typedef Min_sphere_annulus_d_traits_d<H>    Homogeneous_traits;
typedef Min_sphere_d<Cartesian_traits>      Min_sphereC;
typedef Min_sphere_d<Homogeneous_traits>    Min_sphereH;
typedef Point_d<C>                          PointC;
typedef Point_d<H>                          PointH;


// checked traits
Cartesian_traits                                    tC;
Homogeneous_traits                                  tH;

const int n = 10;                           // number of points
const int D = 2;                            // dimension of points
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
    Point_d<C>  PC[n];          // cartesian point set
    Point_d<H>  PH[n];          // homogeneous point set
    int         coord[D];       // arrays for random coordinates
    FT          coordC[D];      // ... cartesian version
    RT          coordH[D+1];    // ... homogeneous version
    Random      my_random;      // random number generator


    for (int i=0; i<n; ++i) {
        int j;
        // random coordinates
        for (j=0; j<D; ++j)
            coord[j] = my_random (r);

        // cartesian point
        for (j=0; j<D; ++j)
            coordC[j] = FT (coord[j]);
        PC[i] = Point_d<C>(D, coordC, coordC+D);

        // homogeneous point
        for (j=0; j<D; ++j)
            coordH[j] = RT(2*coord[j]);
        coordH[D] = RT(2);
        PH[i] = Point_d<H>(D, coordH, coordH+D+1);
    }


    // test constructors

    // default
    Min_sphereC         msC_empty; assert(msC_empty.is_valid(verbose));
    Min_sphereH         msH_empty; assert(msH_empty.is_valid(verbose));

    // from range
    Min_sphereC         msC (PC, PC+n, tC); assert(msC.is_valid(verbose));
    Min_sphereH         msH (PH, PH+n, tH); assert(msH.is_valid(verbose));

    // copy
    Min_sphereC         msC1 (msC); assert(msC1.is_valid(verbose));
    Min_sphereH         msH1 (msH); assert(msH1.is_valid(verbose));

    PointC              centerC (msC.center()), centerC1 (msC1.center());
    PointH              centerH (msH.center()), centerH1 (msH1.center());

    FT                  radiusC (msC.squared_radius()),
                        radiusC1 (msC1.squared_radius());
    Quotient<RT>        radiusH (msH.squared_radius()),
                        radiusH1 (msH1.squared_radius());

    assert(equals (centerC, centerH));
    assert(centerC == centerC1); assert(centerH == centerH1);
    assert(equals (radiusC, radiusH));
    assert(radiusC == radiusC1); assert(radiusH == radiusH1);

    // assignment
    msC1 = msC; msH1 = msH;
    assert(centerC == centerC1); assert(centerH == centerH1);
    assert(radiusC == radiusC1); assert(radiusH == radiusH1);


    // test set method
    msC.set (PC, PC+n); assert(msC.is_valid(verbose));
    msH.set (PH, PH+n); assert(msH.is_valid(verbose));
    assert(centerC == msC.center());
    assert(centerH == msH.center());
    assert(radiusC == msC.squared_radius());
    assert(radiusH == msH.squared_radius());

    // test clear and insert methods
    msC.clear(); assert(msC.is_valid(verbose));
    msH.clear(); assert(msH.is_valid(verbose));
    msC.insert (PC, PC+n); assert(msC.is_valid(verbose));
    msH.insert (PH, PH+n); assert(msH.is_valid(verbose));
    assert(centerC == msC.center());
    assert(centerH == msH.center());
    assert(radiusC == msC.squared_radius());
    assert(radiusH == msH.squared_radius());

    // combined set and insert
    msC.set (PC, PC+n/2); msC.insert (PC+n/2, PC+n);
    assert(msC.is_valid(verbose));
    msH.set (PH, PH+n/2); msH.insert (PH+n/2, PH+n);
    assert(msH.is_valid(verbose));
    assert(centerC == msC.center());
    assert(centerH == msH.center());
    assert(radiusC == msC.squared_radius());
    assert(radiusH == msH.squared_radius());


    // test access functions

    // number_of_points
    assert(msC.number_of_points() == n);
    assert(msH.number_of_points() == n);

    // number_of_support_points
    assert(msC.number_of_support_points() == msH.number_of_support_points());

    // points_begin, points_end
    std::ptrdiff_t m;
    m = std::distance (msC.points_begin(), msC.points_end());
    assert(m == n);
    m = std::distance (msH.points_begin(), msH.points_end());
    assert(m == n);

    // support_points_begin, support_points_end
    m =
    std::distance (msC.support_points_begin(), msC.support_points_end());
    assert(m == msC.number_of_support_points());
    m =
    std::distance (msH.support_points_begin(), msH.support_points_end());
    assert(m == msH.number_of_support_points());

    // ambient dim
    assert(msC.ambient_dimension() == D);
    assert(msH.ambient_dimension() == D);

    // center and squared radius already tested


    // test predicates

    // bounded_side
    assert(msC.bounded_side (centerC) == ON_BOUNDED_SIDE);
    assert(msH.bounded_side (centerH) == ON_BOUNDED_SIDE);

    // has_on_bounded_side
    assert(msC.has_on_bounded_side (centerC));
    assert(msH.has_on_bounded_side (centerH));

    // has_on_boundary already tested in is_valid method
    // has_on_unbounded_side already tested in is_valid method

    // is_empty
    assert(!msC.is_empty());
    assert(!msH.is_empty());

    // is_degenerate
    assert(!msC.is_degenerate());
    assert(!msH.is_degenerate());


    std::ostringstream ost;           // output string
    IO::set_ascii_mode (ost);
    ost << msC << msH << std::endl;      // write spheres

    std::istringstream ist (ost.str().c_str());  // input string
    IO::set_ascii_mode (ist);
    ist >> msC >> msH;              // read spheres

    assert(centerC == msC.center());
    assert(centerH == msH.center());
    assert(radiusC == msC.squared_radius());
    assert(radiusH == msH.squared_radius());

    return 0;
}


// ===== EOF ==================================================================

