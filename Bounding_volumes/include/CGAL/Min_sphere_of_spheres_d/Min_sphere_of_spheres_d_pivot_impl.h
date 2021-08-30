// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Kaspar Fischer


#ifndef CGAL_MINIBALL_PIVOT_C
#define CGAL_MINIBALL_PIVOT_C

#include <CGAL/license/Bounding_volumes.h>


#include <CGAL/Min_sphere_of_spheres_d.h>

namespace CGAL_MINIBALL_NAMESPACE {

  namespace Min_sphere_of_spheres_d_impl {

    template<typename FT>
    bool is_better(const FT& old,const FT& now,const Tag_false /* is_exact */) {
      return now > old;
    }

    template<typename FT>
    bool is_better(const FT&,const FT&,const Tag_true /* is_exact */) {
      return true;
    }

    template<class Traits>
    bool Support_set<Traits>::pivot(std::vector<const typename
                                    Traits::Sphere *>& l,
                                    int& e,
                                    const int d) {

      // remember old radius:
      const Result old = radius();

      // reset basis to {d}:
      reset();
      push(*l[d]);

      // try all subsets:
      std::bitset<Traits::D+1> T;
      int pos = e;
      bool up = true;
      while (pos >= 0) {
        if (pos == e) {
          bool isEnclosingSupporting = is_spanning();
          if (isEnclosingSupporting)
            for(int i=0; i<e; ++i)
              if (!T.test(i) && !contains(t.center_cartesian_begin(*l[i]),
                                          t.radius(*l[i]),
                                          Tol<FT>::result(),Is_exact())) {
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
            return is_better(old,radius(),Is_exact());
          }

          --pos;
          up = false;
        } else if (!T.test(pos)) {
          if (up)
            ++pos;
          else
            if (push(*l[pos])) {
              T.set(pos,true);
              ++pos;
              up = true;
            } else
              --pos;
        } else {
          pop();
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
      reset();
      for (int i=0; i<e; ++i)
        push(*l[i]);
      is_spanning();

      // signal that we failed:
      return false;
    }

  } // namespace Min_sphere_of_spheres_d_impl
} // namespace CGAL_MINIBALL_NAMESPACE

#endif // CGAL_MINIBALL_PIVOT_CC
