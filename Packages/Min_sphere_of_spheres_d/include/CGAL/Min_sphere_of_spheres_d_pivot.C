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
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// file          : include/CGAL/Min_sphere_of_spheres_d_pivot.C
// package       : Min_sphere_of_spheres_d (1.10)
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Kaspar Fischer
// maintainer    : Kaspar Fischer <fischerk@inf.ethz.ch>
// coordinator   : ETH Zurich (Kaspar Fischer)
//
// implementation: dD Smallest Enclosing Sphere of Spheres
// email         : contact@cgal.org
// www           : http://www.cgal.org
//
// ======================================================================


#ifndef CGAL_MINIBALL_PIVOT_C
#define CGAL_MINIBALL_PIVOT_C

#include <CGAL/Min_sphere_of_spheres_d.h>

namespace CGAL_MINIBALL_NAMESPACE {

  template<typename FT>
  bool is_better(const FT& old,const FT& now,const Tag_false is_exact) {
    return now > old;
  }
  
  template<typename FT>
  bool is_better(const FT&,const FT&,const Tag_true is_exact) {
    return true;
  }
  
  template<class Traits>
  bool Min_sphere_of_spheres_d<Traits>::pivot(const int D) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // remember old radius:
    const Result old = ss.radius();
  
    // reset basis to {D}:
    ss.reset();
    ss.push(*l[D]);
  
    // try all subsets:
    std::bitset<Traits::D+1> T;
    int pos = e;
    bool up = true;
    while (pos >= 0) {
      if (pos == e) {
        bool isEnclosingSupporting = ss.is_spanning();
        if (isEnclosingSupporting)
          for(int i=0; i<e; ++i)
            if (!T.test(i) && !ss.contains(t.center_cartesian_begin(*l[i]),
                                           t.radius(*l[i]),
                                           Tol,Is_exact())) {
              isEnclosingSupporting = false;
              break;
            }
        
        if (isEnclosingSupporting) {
          // rearrange pointers:
          int next = 0;
          for(int i=0; i<e; ++i)
            if (T.test(i))
              std::swap(l[next++],l[i]);
          std::swap(l[next++],l[D]);
          e = next;
          return is_better(old,ss.radius(),Is_exact());
        }
        
        --pos;
        up = false;
      } else if (!T.test(pos)) {
        if (up)
          ++pos;
        else
          if (ss.push(*l[pos])) {
            T.set(pos,true);
            ++pos;
            up = true;
          } else
            --pos;
      } else {
        ss.pop();
        T.set(pos,false);
        --pos;
        up = false;
      }
    }
  
    // Here, no basis has been found (this only happens because
    // of rounding errors):
    #ifdef CGAL_MINIBALL_WARNINGS
    std::cerr << '!';
    #endif
    
    // revert basis:
    ss.reset();
    for (int i=0; i<e; ++i)
      ss.push(*l[i]);
    ss.is_spanning();
    
    // signal that we failed:
    return false;
  }
  
} // namespace CGAL_MINIBALL_NAMESPACE

#endif // CGAL_MINIBALL_PIVOT_CC
