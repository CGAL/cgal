// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_VERTEX_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_VERTEX_CONFLICT_2_H

#include <CGAL/license/Apollonius_graph_2.h>



// FIXME: We include the old traits class file for now to get the functors.
#include <CGAL/Uncertain.h>
#include <CGAL/Apollonius_graph_traits_2.h>

namespace CGAL {

//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Uncertain_vertex_conflict_new_2
{
public:
    typedef typename K::Site_2                Site_2;
    typedef typename K::RT                    RT;
    typedef Uncertain<Sign>                   result_type;

private:

    inline
    Uncertain<bool> is_less (const Site_2 &p0, const Site_2 &p1) const
    {
      static const Uncertain<bool> uncertain_bool =
        Uncertain<bool>::indeterminate();

      Uncertain<Comparison_result> cr;

      cr = CGAL::compare( p0.weight(), p1.weight() );
      if ( is_indeterminate(cr) ) { return uncertain_bool; }

      if ( cr == SMALLER ) { return true; }
      if ( cr == LARGER ) { return false; }

      cr = CGAL::compare(p0.x(), p1.x());
      if ( is_indeterminate(cr) ) { return uncertain_bool; }

      if ( cr == SMALLER ) { return true; }
      if ( cr == LARGER ) { return false; }

      cr = CGAL::compare(p0.y(), p1.y());
      if ( is_indeterminate(cr) ) { return uncertain_bool; }
      return cr == SMALLER;
    }

    inline
    int max_radius(const Site_2        &p0, const Site_2 &p1,
            const Site_2 &p2, const Site_2 &p3) const
    {
        int i = 0;
        const Site_2 *p = &p0;

        Uncertain<bool> b;
        b = is_less (*p, p1);
        if ( is_indeterminate(b) ) {
          return -1;
        } else if (b) {
          i = 1; p = &p1;
        }

        b = is_less(*p, p2);
        if ( is_indeterminate(b) ) {
          return -1;
        } else if (b) {
          i = 2; p = &p2;
        }

        b = is_less(*p, p3);
        if ( is_indeterminate(b) ) {
          return -1;
        } else if (b) {
          i = 3;
        }

        return i;
    }

    inline
    Uncertain<Sign>
    predicate(const Site_2 &p1, const Site_2 &p2,
              const Site_2 &p3, const Site_2 &q, bool perturb) const
    {
      RT xq = q.x() - p1.x();
      RT yq = q.y() - p1.y();
      RT wq = q.weight() - p1.weight();
      RT aq = CGAL::square(xq) + CGAL::square(yq) - CGAL::square(wq);

      // q is hiding p1
      Uncertain<Sign> s = CGAL::sign(aq);
      if ( is_indeterminate(s) ) { return s; }
      if (s != POSITIVE){
        // I BELIEVE MENELAOS RETURNS -1 in this case even when degernate
        //if (sign (aq) == ZERO && ! perturb) return ZERO;
        //        return NEGATIVE;
        return POSITIVE;
      }

      RT x2 = p2.x() - p1.x();
      RT y2 = p2.y() - p1.y();
      RT w2 = p2.weight() - p1.weight();
      RT a2 = x2 * x2 + y2 * y2 - w2 * w2;

      CGAL_assertion (a2 > 0);

      RT x3 = p3.x() - p1.x();
      RT y3 = p3.y() - p1.y();
      RT w3 = p3.weight() - p1.weight();
      RT a3 = x3 * x3 + y3 * y3 - w3 * w3;

      CGAL_assertion (a3 > 0);

      RT ax3q = a3 * xq - x3 * aq;
      RT ax2q = a2 * xq - x2 * aq;
      RT ax23 = a2 * x3 - x2 * a3;

      RT ay23 = a2 * y3 - y2 * a3;
      RT ay2q = a2 * yq - y2 * aq;
      RT ay3q = a3 * yq - y3 * aq;

      RT axw23q = ax23 * wq - ax2q * w3 + ax3q * w2;
      RT ayw23q = ay23 * wq - ay2q * w3 + ay3q * w2;

      RT axy23q = y2 * ax3q - y3 * ax2q + yq * ax23;

      // orientation
      Uncertain<Sign> orient = CGAL::sign(axy23q);
      if ( is_indeterminate(orient) ) { return orient; }

      // orientation degenerate
      if (orient == ZERO) {
        Uncertain<Sign> orient1 = CGAL::sign (ax23);
        if ( is_indeterminate(orient1) ) { return orient1; }

        Sign power_test;
        if ( orient1 == ZERO ) {
          Uncertain<Sign> s_ay23 = CGAL::sign(ay23);
          if ( is_indeterminate(s_ay23) ) { return s_ay23; }
          Uncertain<Sign> s_ayw23q = CGAL::sign(ayw23q);
          if ( is_indeterminate(s_ayw23q) ) { return s_ayw23q; }

          power_test = s_ay23 * s_ayw23q;
        } else {
          Uncertain<Sign> s_axw23q = CGAL::sign(axw23q);
          if ( is_indeterminate(s_axw23q) ) { return s_axw23q; }

          power_test = orient1 * s_axw23q;
        }

        if (power_test != ZERO || ! perturb) {
          return -power_test;
        }

        int i = max_radius (p1, p2, p3, q);
        if ( i == -1 ) {
          return Uncertain<Sign>::indeterminate();
        }

        if (i == 3) { return NEGATIVE; }

        Uncertain<Sign> o23, o2q, o3q;

        if (orient1 == ZERO) {
          o23 = CGAL::sign (ay23);
          o2q = CGAL::sign (ay2q);
          o3q = CGAL::sign (ay3q);
        } else {
          o23 = CGAL::sign (ax23);
          o2q = CGAL::sign (ax2q);
          o3q = CGAL::sign (ax3q);
        }
        if ( is_indeterminate(o23) ||
             is_indeterminate(o2q) ||
             is_indeterminate(o3q) ) {
          return Uncertain<Sign>::indeterminate();
        }

        if (o23 != o2q) { return i == 2 ? NEGATIVE : POSITIVE; }

        if (o23 == o3q) { return i == 1 ? NEGATIVE : POSITIVE; }

        return i == 0 ? NEGATIVE : POSITIVE;
      } // if (orient == ZERO )

      // radical side
      RT rs23q = ax23 * axw23q + ay23 * ayw23q;
      Uncertain<Sign> radSide = CGAL::sign (rs23q);
      if ( is_indeterminate(radSide) ) { return radSide; }

      if (radSide == ZERO || radSide != orient) { return orient; }

      // radical intersection
      Uncertain<Sign> radInt =
        CGAL::sign(axw23q * axw23q + ayw23q * ayw23q - axy23q * axy23q);
      if ( is_indeterminate(radInt) ) { return radInt; }

      // radical intersection degenerate
      if (radInt == ZERO) {
        Uncertain<Sign> radSideQ = CGAL::sign(ax23 * axw23q + ay23 * ayw23q);
        if ( is_indeterminate(radSideQ) ) { return radSideQ; }

        CGAL_assertion (radSideQ != ZERO);

        if (! perturb) { return (radSideQ == orient) ? Uncertain<Sign>(ZERO) : orient; }

        int i = max_radius (p1, p2, p3, q);
        if ( i == -1 ) { return Uncertain<Sign>::indeterminate(); }

        if (i == 3) {
          radInt = radSideQ;
        } else if (i == 2) {
          radInt = -CGAL::sign(ax2q * axw23q + ay2q * ayw23q);
          if ( is_indeterminate(radInt) ) { return radInt; }
          if (radInt == ZERO) { return NEGATIVE; }
        } else if (i == 1) {
          radInt = CGAL::sign(ax3q * axw23q + ay3q * ayw23q);
          if ( is_indeterminate(radInt) ) { return radInt; }
          if (radInt == ZERO) { return NEGATIVE; }
        } else {
          CGAL_assertion (i == 0);
          Uncertain<Sign> radSide1 =
            -CGAL::sign(ax2q * axw23q + ay2q * ayw23q);
          if ( is_indeterminate(radSide1) ) { return radSide1; }

          if (radSide1 == ZERO) { return NEGATIVE; }

          Uncertain<Sign> radSide2 = CGAL::sign(ax3q * axw23q + ay3q * ayw23q);
          if ( is_indeterminate(radSide2) ) { return radSide2; }

          if (radSide2 == ZERO) { return NEGATIVE; }

          radInt = -Sign( Sign(radSideQ) + Sign(radSide1) + Sign(radSide2) );
        }
      }

      CGAL_assertion (!perturb || radInt != ZERO);

      if (radInt == NEGATIVE) { return orient; }

      return -radSide;
    }


    inline
    Uncertain<Sign>
    predicate(const Site_2 &p1, const Site_2 &p2,
              const Site_2 &q, bool perturb) const
    {
      // NOTE:***************************************
      // * the perturb boolean variable is not used
      // * for consistancy with Menelaos
      // NOTE:***************************************
      RT x2 = p2.x() - p1.x();
      RT y2 = p2.y() - p1.y();
      RT w2 = p2.weight() - p1.weight();
      RT xq =  q.x() - p1.x();
      RT yq =  q.y() - p1.y();
      RT wq =  q.weight() - p1.weight();

      RT xw2q = x2 * wq - xq * w2;
      RT yw2q = y2 * wq - yq * w2;
      RT xy2q = x2 * yq - xq * y2;

      // orientation
      Uncertain<Sign> orient = CGAL::sign(xy2q);
      if ( is_indeterminate(orient) ) { return orient; }

      // orientation degenerate
      if (orient == ZERO) {
        Uncertain<Sign> o12 = CGAL::sign(x2);
        if ( is_indeterminate(o12) ) { return o12; }

        Uncertain<Sign> o1q, o2q;
        Uncertain<Sign> power_test;

        if (o12 != ZERO) {
          Uncertain<Sign> s_xw2q = CGAL::sign(xw2q);
          if ( is_indeterminate(s_xw2q) ) { return s_xw2q; }
          power_test = o12 * s_xw2q;

          // this results is consistant with Menelaos
          if (power_test != ZERO) { return -power_test; }

            // this result is consistant with the perturb on off idea
            //if (power_test != ZERO || ! perturb) return -power_test;

          o1q = CGAL::sign(xq);
          o2q = CGAL::sign(q.x() - p2.x());
        } else {
          o12 = CGAL::sign(y2);
          if ( is_indeterminate(o12) ) { return o12; }
          Uncertain<Sign> s_yw2q = CGAL::sign(yw2q);
          if ( is_indeterminate(s_yw2q) ) { return s_yw2q; }
          power_test = o12 * s_yw2q;

          // this results is consistant with Menelaos
          if (power_test != ZERO) { return -power_test; }

          // this result is consistant with the perturb on off idea
          //if (power_test != ZERO || ! perturb) return -power_test;

          o1q = CGAL::sign(yq);
          o2q = CGAL::sign(q.y() - p2.y());
        }

        if ( is_indeterminate(o1q) || is_indeterminate(o2q) ) {
          return Uncertain<Sign>::indeterminate();
        }

        if (o1q != o12) { return POSITIVE; }
        if (o2q == o12) { return POSITIVE; }

        return NEGATIVE;
      }

      // radical side
      RT rs12q = x2 * xw2q + y2 * yw2q;
      Uncertain<Sign> radSide = CGAL::sign(rs12q);
      if ( is_indeterminate(radSide) ) { return radSide; }

      if (radSide == ZERO || radSide == orient) {
        return -orient;
      }

      // radical intersection
      Uncertain<Sign> radInt =
        CGAL::sign(CGAL::square(xw2q) + CGAL::square(yw2q)
                   - CGAL::square(xy2q));
      if ( is_indeterminate(radInt) ) { return radInt; }

      // radical intersection degerate
      if (radInt == ZERO) {
        CGAL_assertion (radSide != ZERO);

        // this result is consistant with the perturb on off idea
        //if (! perturb) return (radSide == orient) ? ZERO : orient;

        RT rs2q1 = (p2.x() - q.x()) * xw2q + (p2.y() - q.y()) * yw2q;
        Uncertain<Sign> radSide1 = CGAL::sign(rs2q1);
        if ( is_indeterminate(radSide1) ) { return radSide1; }
        if (radSide1 == ZERO) { return NEGATIVE; }

        RT rsq12 = xq * xw2q + yq * yw2q;
        Uncertain<Sign> radSide2 = CGAL::sign(rsq12);
        if ( is_indeterminate(radSide2) ) { return radSide2; }
        if (radSide2 == ZERO) { return NEGATIVE; }

        return -(radSide1 * radSide2);
      }

      CGAL_assertion (!perturb || radInt != ZERO);

      if (radInt == POSITIVE) { return orient; }
      return radSide;
    }

public:
    inline
    Uncertain<Sign> operator()(const Site_2 &p1, const Site_2 &p2,
                               const Site_2 &p3, const Site_2 &q,
                               bool perturb = true) const
    {
      Uncertain<Sign> newPred = predicate (p1, p2, p3, q, perturb);
      if ( is_indeterminate(newPred) ) { return newPred; }
      CGAL_assertion (! perturb || newPred != ZERO);
      return newPred;
    }

    inline
    Uncertain<Sign> operator()(const Site_2 &p1, const Site_2 &p2,
                               const Site_2 &q, bool perturb = true) const
    {
      Uncertain<Sign> newPred = predicate (p1, p2, q, perturb);
      if ( is_indeterminate(newPred) ) { return newPred; }
      CGAL_assertion (! perturb || newPred != ZERO);
      return newPred;
    }
};


} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_UNCERTAIN_VERTEX_CONFLICT_2_H
