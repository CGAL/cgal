// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//                 Christophe Delage <Christophe.Delage@sophia.inria.fr>
//                 David Millman <dlm336@cs.nyu.edu>

#ifndef CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H


// FIXME: We include the old traits class file for now to get the functors.
#include <CGAL/Apollonius_graph_traits_2.h>


CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------

template < class K, class Method_tag > 
class Vertex_conflict_new_2 
{
public:
    typedef typename K::Site_2                Site_2;
    typedef typename K::FT                    FT;	
    typedef Sign                              result_type;
    struct Arity {};

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
        FT xq = q.x() - p1.x(), yq = q.y() - p1.y(), wq = q.weight() - p1.weight();
        FT aq = xq * xq + yq * yq - wq * wq;
        
        // q is hiding p1
        if (sign (aq) != POSITIVE){
            // I BELIEVE MENELAOS RETURNS -1 in this case even when degernate 
            //if (sign (aq) == ZERO && ! perturb) return ZERO;
            return NEGATIVE;
        }

        FT x2 = p2.x() - p1.x(), y2 = p2.y() - p1.y(), w2 = p2.weight() - p1.weight();
        FT a2 = x2 * x2 + y2 * y2 - w2 * w2;

        CGAL_assertion (a2 > 0);

        FT x3 = p3.x() - p1.x(), y3 = p3.y() - p1.y(), w3 = p3.weight() - p1.weight();
        FT a3 = x3 * x3 + y3 * y3 - w3 * w3;

        CGAL_assertion (a3 > 0);

        FT ax3q = a3 * xq - x3 * aq; 
        FT ax2q = a2 * xq - x2 * aq;
        FT ax23 = a2 * x3 - x2 * a3;

        FT ay23 = a2 * y3 - y2 * a3;
        FT ay2q = a2 * yq - y2 * aq;
        FT ay3q = a3 * yq - y3 * aq;

        FT axw23q = ax23 * wq - ax2q * w3 + ax3q * w2;
        FT ayw23q = ay23 * wq - ay2q * w3 + ay3q * w2;

        FT axy23q = y2 * ax3q - y3 * ax2q + yq * ax23;

        // orientation
        Sign orient = CGAL::sign(axy23q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign orient1 = sign (ax23);

            Sign power_test = (orient1 == ZERO ?
                    Sign (sign (ay23) * sign (ayw23q)) :
                    Sign (orient1 * sign (axw23q)));

            if (power_test != ZERO || ! perturb) return Sign (- power_test);

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) return NEGATIVE;

            Sign o23, o2q, o3q;

            if (orient1 == ZERO) {
                o23 = sign (ay23);
                o2q = sign (ay2q);
                o3q = sign (ay3q);
            } else {
                o23 = sign (ax23);
                o2q = sign (ax2q);
                o3q = sign (ax3q);
            }

            if (o23 != o2q) return i == 2 ? NEGATIVE : POSITIVE;

            if (o23 == o3q) return i == 1 ? NEGATIVE : POSITIVE;

            return i == 0 ? NEGATIVE : POSITIVE;
        }
 
        // radical side 
        FT rs23q = ax23 * axw23q + ay23 * ayw23q;
        Sign radSide = CGAL::sign (rs23q);

        if (radSide == ZERO || radSide != orient) return orient;
       
        // radical intersection
        Sign radInt = enum_cast<Sign> (CGAL::compare (axw23q * axw23q + ayw23q * ayw23q, axy23q * axy23q));

        // radical intersection degenerate
        if (radInt == ZERO) {
            Sign radSideQ = sign (ax23 * axw23q + ay23 * ayw23q);
            
            CGAL_assertion (radSideQ != ZERO);

            if (! perturb) return (radSideQ == orient) ? ZERO : orient;

            int i = max_radius (p1, p2, p3, q);

            if (i == 3) { 
                radInt = radSideQ;
            } else if (i == 2) {
                radInt = Sign (- sign (ax2q * axw23q + ay2q * ayw23q));
                if (radInt == ZERO) return NEGATIVE;
            } else if (i == 1) {
                radInt = sign (ax3q * axw23q + ay3q * ayw23q);
                if (radInt == ZERO) return NEGATIVE;
            } else {
                CGAL_assertion (i == 0);
                Sign radSide1 = Sign (- sign (ax2q * axw23q + ay2q * ayw23q));      
                if (radSide1 == ZERO) return NEGATIVE;

                Sign radSide2 = sign (ax3q * axw23q + ay3q * ayw23q);	
                if (radSide2 == ZERO) return NEGATIVE;

                radInt = Sign (- (radSideQ + radSide1 + radSide2));
            }
        }
        
        CGAL_assertion (! perturb || radInt != ZERO);

        if (radInt == NEGATIVE) return orient;
        
        return Sign (- radSide);
    }
    

    inline
    Sign predicate(const Site_2 &p1, const Site_2 &p2, 
            const Site_2 &q, bool perturb) const
    {
        // NOTE:***************************************
        // * the perturb boolean variable is not used 
        // * for consistancy with Menelaos
        // NOTE:***************************************
        FT x2 = p2.x() - p1.x(), y2 = p2.y() - p1.y(), w2 = p2.weight() - p1.weight();
        FT xq =  q.x() - p1.x(), yq =  q.y() - p1.y(), wq =  q.weight() - p1.weight();

        FT xw2q = x2 * wq - xq * w2;
        FT yw2q = y2 * wq - yq * w2;
        FT xy2q = x2 * yq - xq * y2;
        
        // orientation
        Sign orient = CGAL::sign(xy2q);

        // orientation degenerate
        if (orient == ZERO) {
            Sign o12 = sign (x2);
            Sign o1q, o2q;

            Sign power_test;
            if (o12 != ZERO) {
                power_test = Sign (o12 * sign (xw2q));
                 
                // this results is consistant with Menelaos
                if (power_test != ZERO) return Sign(- power_test);

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return Sign(- power_test);

                o1q = sign (xq);
                o2q = sign (q.x() - p2.x());
            } else {
                o12 = sign (y2);
                power_test = Sign (o12 * sign (yw2q));

                // this results is consistant with Menelaos
                if (power_test != ZERO) return Sign(- power_test);

                // this result is consistant with the perturb on off idea
                //if (power_test != ZERO || ! perturb) return Sign(- power_test);

                o1q = sign (yq);
                o2q = sign (q.y() - p2.y());
            }

            if (o1q != o12) return POSITIVE;
            if (o2q == o12) return POSITIVE;

            return NEGATIVE;
        }

        // radical side 
        FT rs12q = x2 * xw2q + y2 * yw2q;
        Sign radSide = CGAL::sign (rs12q);

        if (radSide == ZERO || radSide == orient) {
            return Sign(- orient);
        }

        // radical intersection
        Sign radInt = enum_cast<Sign> (CGAL::compare (xw2q * xw2q + yw2q * yw2q, xy2q * xy2q));

        // radical intersection degerate
        if (radInt == ZERO) {
            CGAL_assertion (radSide != ZERO);
            
            // this result is consistant with the perturb on off idea
            //if (! perturb) return (radSide == orient) ? ZERO : orient;

            FT rs2q1 = (p2.x() - q.x()) * xw2q + (p2.y() - q.y()) * yw2q;
            Sign radSide1 = sign (rs2q1);
            if (radSide1 == ZERO) return NEGATIVE;
            
            FT rsq12 = xq * xw2q + yq * yw2q;
            Sign radSide2 = sign (rsq12);
            if (radSide2 == ZERO) return NEGATIVE;
 
            return Sign (- (radSide1 * radSide2));
        }

        CGAL_assertion (! perturb || radInt != ZERO);

        if (radInt == POSITIVE) return orient;
        return radSide;
    }

public:
    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &p3, const Site_2 &q, bool perturb = true) const
    {	
        Sign newPred = predicate (p1, p2, p3, q, perturb);	
        CGAL_assertion (! perturb || newPred != ZERO);
        return newPred;
    }

    inline
    Sign operator()(const Site_2 &p1, const Site_2 &p2,
                    const Site_2 &q, bool perturb = true) const
    {
        Sign newPred = predicate (p1, p2, q, perturb);	
        CGAL_assertion (! perturb || newPred != ZERO);
        return newPred;
    }
};


CGAL_END_NAMESPACE

#endif // CGAL_APOLLONIUS_GRAPH_2_VERTEX_CONFLICT_2_H
