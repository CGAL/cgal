// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation are part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation are provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - Please check the CGAL web site http://www.cgal.org/index2.html for 
//   availability.
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-2.4
// release_date  : 2002, May 16
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres_pivot.C
// package       : Min_sphere_of_spheres_d (1.00)
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Kaspar Fischer
//
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================

#ifndef CGAL_MIN_SPHERE_OF_SPHERES_PIVOT_C
#define CGAL_MIN_SPHERE_OF_SPHERES_PIVOT_C

#include <vector>

#include <CGAL/basic.h>
#include <CGAL/Min_sphere_of_spheres_d.h>

namespace CGAL {

  // another namespace to "hide" implementation details:
  namespace Min_sphere_of_spheres_impl {

    // The following routines are used to decide, at the end of
    // routine pivot(), whether we have made any progress.  (Progress
    // is guaranteed in case of exact computation.)

    template<typename FT>
    inline bool is_better(const FT& old,const FT& now,
			  const Tag_false is_exact) {
      return now > old;
    }
    
    template<typename FT>
    inline bool is_better(const FT&,const FT&,
			  const Tag_true is_exact) {
      return true;
    }

    // When the client uses the farthest-first heuristic in
    // combination with exact arithmetic, floating-point arithmetic is
    // used to quickly find the ball "farthest away."  Due to rounding
    // errors however, this ball need not really be the ball farthest
    // away.  Consequently, we will have to run the LP-algorithm at
    // the end of the farthest-first heuristic in order to find the
    // exact miniball. --- It is for this distinction that we need the
    // following two routines:

    inline bool is_heuristic(const Farthest_first_heuristic) {
      return true;
    }

    inline bool is_heuristic(const LP_algorithm) {
      return false;
    }

  } // namespace Min_sphere_of_spheres_impl
    
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::pivot(const int D) {
    // remember old radius:
    const Result old = ss.radius();
    
    // cerr << "[Entering pivot() with radius " << old << "]" << endl;

    // reset basis to {D}:
    ss.clear();
    ss.push(*l[D]);
    
    // try all subsets:
    std::vector<bool> T(ss.dim,false);
    int pos = e;
    bool up = true;
    while (pos >= 0) {
      if (pos == e) {
	bool isEnclosingSupporting = ss.isValid();
	
	if (isEnclosingSupporting)
	  for(int i=0; i<e; ++i)
	    if (!T[i] && !ss.contains(t.center_coordinates_begin(*l[i]),
				      t.radius(*l[i]),
				      Min_sphere_of_spheres_impl::Tol,
				      Is_exact())) {
	      isEnclosingSupporting = false;
	      break;
	    }
	
	if (isEnclosingSupporting) {
	  // rearrange pointers:
	  int next = 0;
	  for(int i=0; i<e; ++i)
	    if (T[i])
	      std::swap(l[next++],l[i]);
	  std::swap(l[next++],l[D]);
	  e = next;
	  return Min_sphere_of_spheres_impl::
	    is_better(old,ss.radius(),Is_exact());
	}
	
	--pos;
	up = false;
      } else if (!T[pos]) {
	if (up)
	  ++pos;
	else
	  if (ss.push(*l[pos])) {
	    T[pos] = true;
	    ++pos;
	    up = true;
	  } else
	    --pos;
      } else {
	ss.pop();
	T[pos] = false;
	--pos;
	up = false;
      }
    }
    
    // Here, no basis has been found. --- We can only get here under
    // approximate computation:
    CGAL_assertion(Min_sphere_of_spheres_impl::is_approximate(Is_exact()));

    // Output a warning if the user is running the LP-algorithm
    // instead of the FarthestFirst heuristic (because the former
    // doesn't handle degeneracies as well as the latter does):
    CGAL_warning_msg(Min_sphere_of_spheres_impl::is_heuristic(Algorithm()),
    "LP_algorithm can't cope with degeneracies, use Farthest_first_heuristic");
    
    // revert basis:
    ss.clear();
    for (int i=0; i<e; ++i)
      ss.push(*l[i]);
    ss.isValid();
    
    // signal that we failed:
    return false;
  }

} // namespace CGAL

#endif // CGAL_MIN_SPHERE_OF_SPHERES_PIVOT_C
