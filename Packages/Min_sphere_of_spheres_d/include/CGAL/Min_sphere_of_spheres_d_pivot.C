// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Kaspar Fischer


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
  bool Min_sphere_of_spheres_d<Traits>::pivot(const int d) {
    using namespace Min_sphere_of_spheres_d_impl;
  
    // remember old radius:
    const Result old = ss.radius();
  
    // reset basis to {d}:
    ss.reset();
    ss.push(*l[d]);
  
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
          // rearrange balls:
          int next = 0;
          for(int i=0; i<e; ++i)
            if (T.test(i))
              std::swap(l[next++],l[i]);
          std::swap(l[next++],l[d]);
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
