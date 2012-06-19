// Copyright (c) 2006 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H


#include <CGAL/Apollonius_graph_2/basic.h>


namespace CGAL {

namespace ApolloniusGraph_2 {

//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag > 
class Vertex_conflict_new_2 
{
public:
    typedef typename K::Site_2                Site_2;
    typedef typename K::RT                    RT;	
    typedef Sign                              result_type;

private:
    
    inline
    bool is_less (const Site_2 &p0, const Site_2 &p1) const
    {
        if (p0.weight() < p1.weight()) return true;
        if (p0.weight() > p1.weight()) return false;

        if (p0.x() < p1.x()) return true;
        if (p0.x() > p1.x()) return false;

        return p0.y() < p1.y();
    }

    inline
    int max_radius(const Site_2	&p0, const Site_2 &p1,
            const Site_2 &p2, const Site_2 &p3) const
    {
        int i = 0;
        const Site_2 *p = &p0;

        if (is_less (*p, p1)) { i = 1; p = &p1; }
        if (is_less (*p, p2)) { i = 2; p = &p2; }
        if (is_less (*p, p3)) { i = 3; }

        return i;
    }

    inline
    Sign predicate (const Site_2 &p1, const Site_2 &p2,
            const Site_2 &p3, const Site_2 &q, bool perturb) const
    {	
        RT xq = q.x() - p1.x();
	RT yq = q.y() - p1.y();
	RT wq = q.weight() - p1.weight();
        RT aq = CGAL::square(xq) + CGAL::square(yq) - CGAL::square(wq);

        // q is hiding p1
        if (CGAL::sign(aq) != POSITIVE){
            // I BELIEVE MENELAOS RETURNS -1 in this case even when degernate 
            //if (sign (aq) == ZERO && ! perturb) return ZERO;

	  //return NEGATIVE;
	  return POSITIVE;
        }

        RT x2 = p2.x() - p1.x();
	RT y2 = p2.y() - p1.y();
	RT w2 = p2.weight() - p1.weight();
        RT a2 = CGAL::square(x2) + CGAL::square(y2) - CGAL::square(w2);

        CGAL_assertion (a2 > 0);

        RT x3 = p3.x() - p1.x();
	RT y3 = p3.y() - p1.y();
	RT w3 = p3.weight() - p1.weight();
        RT a3 = CGAL::square(x3) + CGAL::square(y3) - CGAL::square(w3);

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
        Sign orient = CGAL::sign(axy23q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign orient1 = CGAL::sign(ax23);

            Sign power_test =
	      ( orient1 == ZERO ?
		(CGAL::sign(ay23) * CGAL::sign(ayw23q)) :
		(orient1 * CGAL::sign(axw23q))
		);

            if (power_test != ZERO || !perturb) {
	      return -power_test;
	    }

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) { return NEGATIVE; }

            Sign o23, o2q, o3q;

            if (orient1 == ZERO) {
                o23 = CGAL::sign(ay23);
                o2q = CGAL::sign(ay2q);
                o3q = CGAL::sign(ay3q);
            } else {
                o23 = CGAL::sign(ax23);
                o2q = CGAL::sign(ax2q);
                o3q = CGAL::sign(ax3q);
            }

            if (o23 != o2q) { return i == 2 ? NEGATIVE : POSITIVE; }

            if (o23 == o3q) { return i == 1 ? NEGATIVE : POSITIVE; }

            return i == 0 ? NEGATIVE : POSITIVE;
        }
 

        // radical side 
        RT rs23q = ax23 * axw23q + ay23 * ayw23q;
        Sign radSide = CGAL::sign(rs23q);

        if (radSide == ZERO || radSide != orient) { return orient; }
       
        // radical intersection
        Sign radInt =
	  CGAL::sign(CGAL::square(axw23q) + CGAL::square(ayw23q)
		     - CGAL::square( axy23q));

        // radical intersection degenerate
        if (radInt == ZERO) {
            Sign radSideQ = CGAL::sign(ax23 * axw23q + ay23 * ayw23q);
            
            CGAL_assertion (radSideQ != ZERO);

            if (!perturb) { return (radSideQ == orient) ? ZERO : orient; }

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) { 
                radInt = radSideQ;
            } else if (i == 2) {
                radInt = -CGAL::sign(ax2q * axw23q + ay2q * ayw23q);
                if (radInt == ZERO) { return NEGATIVE; }
            } else if (i == 1) {
                radInt = CGAL::sign(ax3q * axw23q + ay3q * ayw23q);
                if (radInt == ZERO) { return NEGATIVE; }
            } else {
                CGAL_assertion (i == 0);
                Sign radSide1 = -CGAL::sign(ax2q * axw23q + ay2q * ayw23q);
                if (radSide1 == ZERO) { return NEGATIVE; }

                Sign radSide2 = CGAL::sign(ax3q * axw23q + ay3q * ayw23q);
                if (radSide2 == ZERO) { return NEGATIVE; }

                radInt = Sign (-(radSideQ + radSide1 + radSide2));
            }
        }
        
        CGAL_assertion (!perturb || radInt != ZERO);

        if (radInt == NEGATIVE) { return orient; }
        
        return -radSide;
    }
    

    inline
    Sign predicate(const Site_2 &p1, const Site_2 &p2, 
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
        Sign orient = CGAL::sign(xy2q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign o12 = CGAL::sign(x2);
            Sign o1q, o2q;

            Sign power_test;
            if (o12 != ZERO) {
                power_test = o12 * CGAL::sign(xw2q);
                 
                // this results is consistant with Menelaos
                if (power_test != ZERO) { return -power_test; }

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return -power_test;
                o1q = CGAL::sign(xq);
                o2q = CGAL::sign(q.x() - p2.x());
            } else {
                o12 = CGAL::sign(y2);
                power_test = o12 * CGAL::sign(yw2q);

                // this results is consistant with Menelaos
                if (power_test != ZERO) { return -power_test; }

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return -power_test;
                o1q = CGAL::sign(yq);
                o2q = CGAL::sign(q.y() - p2.y());
            }

            if (o1q != o12) { return POSITIVE; }
            if (o2q == o12) { return POSITIVE; }

            return NEGATIVE;
        }

        // radical side 
        RT rs12q = x2 * xw2q + y2 * yw2q;
        Sign radSide = CGAL::sign(rs12q);

        if (radSide == ZERO || radSide == orient) {
            return -orient;
        }

        // radical intersection
        Sign radInt =
	  CGAL::sign(CGAL::square(xw2q) + CGAL::square(yw2q)
		     - CGAL::square(xy2q));

        // radical intersection degerate
        if (radInt == ZERO) {
            CGAL_assertion (radSide != ZERO);
            
            // this result is consistant with the perturb on off idea
            //if (! perturb) return (radSide == orient) ? ZERO : orient;

            RT rs2q1 = (p2.x() - q.x()) * xw2q + (p2.y() - q.y()) * yw2q;
            Sign radSide1 = CGAL::sign(rs2q1);
            if (radSide1 == ZERO) { return NEGATIVE; }
            
            RT rsq12 = xq * xw2q + yq * yw2q;
            Sign radSide2 = CGAL::sign(rsq12);
            if (radSide2 == ZERO) { return NEGATIVE; }
 
            return -(radSide1 * radSide2);
        }

        CGAL_assertion (!perturb || radInt != ZERO);

        if (radInt == POSITIVE) { return orient; }
        return radSide;
    }

public:
    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &p3, const Site_2 &q,
		    bool perturb = true) const
    {	
        Sign newPred = predicate(p1, p2, p3, q, perturb);	
        CGAL_assertion (!perturb || newPred != ZERO);
        return newPred;
    }

    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &q, bool perturb = true) const
    {
        Sign newPred = predicate(p1, p2, q, perturb);	
        CGAL_assertion (!perturb || newPred != ZERO);
        return newPred;
    }
};

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H
