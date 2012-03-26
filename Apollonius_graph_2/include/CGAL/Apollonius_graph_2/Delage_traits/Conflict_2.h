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

#ifndef CGAL_APOLLONIUS_GRAPH_2_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_CONFLICT_2_H 1

#include <CGAL/Apollonius_graph_2/basic.h>
#include <CGAL/Apollonius_graph_2/Predicate_constructions_C2.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

//-----------------------------------------------------------------------
//                        Conflict Base
//-----------------------------------------------------------------------
template < class K, class Method_tag >
class Conflict_2
{
public:
  typedef Inverted_weighted_point_2<K>  Inverted_weighted_point;
  typedef typename K::RT                RT;	
  typedef typename K::Sign              Sign;
  typedef Sign                          result_type;

protected:

    Sign orientation(const Inverted_weighted_point &p1, 
                     const Inverted_weighted_point &p2,
                     const Inverted_weighted_point &p3) const
    {
        return sign_of_determinant( p1.p(), p2.p(), p3.p(),
				       p1.x(), p2.x(), p3.x(),
				       p1.y(), p2.y(), p3.y() );
    }

    Sign radical_intersection(const Inverted_weighted_point &p1,
                              const Inverted_weighted_point &p2,
                              const Inverted_weighted_point &p3, int i) const
    {
        Sign s =  CGAL::sign(
		CGAL::square(determinant<RT>(
                        p1.p(), p1.weight(), p1.y(),
                        p2.p(), p2.weight(), p2.y(),
                        p3.p(), p3.weight(), p3.y()))
                + CGAL::square(determinant<RT>(
                        p1.p(), p1.x(), p1.weight(),
                        p2.p(), p2.x(), p2.weight(),
                        p3.p(), p3.x(), p3.weight()))
                - CGAL::square(determinant<RT>(
                        p1.p(), p1.x(), p1.y(),
                        p2.p(), p2.x(), p2.y(),
                        p3.p(), p3.x(), p3.y())));

        if (s != ZERO || i < 0) { return s; }

        // perturbation
        switch (i) {
            case (1) : return radical_side(p3, p2, p1, 3);
            case (2) : return radical_side(p1, p3, p2, 3);
            case (3) : return radical_side(p1, p2, p3, 3);
            default :
                CGAL_error_msg( "Case should not hapen");
                return ZERO;
        }
    }	

    Sign radical_side(const Inverted_weighted_point &p1,
                      const Inverted_weighted_point &p2,
                      const Inverted_weighted_point &p3, int i) const
    {
        CGAL_assertion(i == -1 || i == 1 || i == 2 || i == 3);

        Sign s =  -CGAL::sign(
                  determinant<RT>(
                      p1.p(), p1.x(), 
                      p2.p(), p2.x()) 
                  * determinant<RT>(
                      p1.p(), p1.weight(), p1.x(),
                      p2.p(), p2.weight(), p2.x(),
                      p3.p(), p3.weight(), p3.x())
                  + determinant<RT>(
                      p1.p(), p1.y(),
                      p2.p(), p2.y())
                  * determinant<RT>(
                      p1.p(), p1.weight(), p1.y(),
                      p2.p(), p2.weight(), p2.y(),
                      p3.p(), p3.weight(), p3.y()));

        if (s != ZERO || i < 1) { return s; }

        // perturbation
        return POSITIVE;
    }

    Sign power_test(const Inverted_weighted_point &p1,
		    const Inverted_weighted_point &p2,
		    const Inverted_weighted_point &p3, int i) const
    {
        CGAL_assertion(i == -1 || i == 1 || i == 2 || i == 3);
        Sign s;

	Sign s_xTest = sign_of_determinant(p2.p(), p2.x(), p1.p(), p1.x());
	if ( s_xTest != ZERO ) {
	  s = s_xTest *
	    sign_of_determinant( p1.p(), p1.x(), p1.weight(),
				    p2.p(), p2.x(), p2.weight(),
				    p3.p(), p3.x(), p3.weight() );
	} else {
            s = sign_determinant2x2( p2.p(), p2.y(), p1.p(), p1.y() ) *
	      sign_of_determinant( p1.p(), p1.y(), p1.weight(),
				      p2.p(), p2.y(), p2.weight(),
				      p3.p(), p3.y(), p3.weight() );
	}
        if (s != ZERO || i < 1) { return s; }

        switch (i) {
	case 1:
	  return ordered_on_line (p1, p2, p3) ? NEGATIVE : POSITIVE;
	case 2:
	  return ordered_on_line (p2, p1, p3) ? NEGATIVE : POSITIVE;
	case 3:
	  return NEGATIVE;
	default:
	  CGAL_error_msg( "this should not happen.");
        }

        // perturbation
        bool ool;
        if (i == 1) {
            ool = ordered_on_line(p2, p1, p3);
        } else if (i == 2) {
            ool = ordered_on_line(p1, p2, p3);
        } else if (i == 3) {
            ool = ordered_on_line(p1, p3, p2); 
        } else 	{		
            CGAL_error_msg( "this does not happen.");
            ool = false;
        }

        if (ool) { return NEGATIVE; }
            
        return POSITIVE;
    }

    Sign ordered_on_line_test(const Inverted_weighted_point &p1,
                              const Inverted_weighted_point &p2) const
    {
      Sign s_det = sign_of_determinant( p1.p(), p1.x(), p2.p(), p2.x() );
      if ( s_det != ZERO ) { return s_det; }
      return sign_of_determinant( p1.p(), p1.y(), p2.p(), p2.y() );
    }

    bool ordered_on_line(const Inverted_weighted_point &p1,
                         const Inverted_weighted_point &p2,
                         const Inverted_weighted_point &p3) const
    {
      if (ordered_on_line_test(p1, p2) == POSITIVE) {
	return ordered_on_line_test(p2, p3) == POSITIVE;
      }
      return ordered_on_line_test(p3, p2) == POSITIVE;
    }
};


} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_CONFLICT_2_H
