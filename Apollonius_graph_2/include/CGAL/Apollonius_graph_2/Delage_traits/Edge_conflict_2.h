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

#ifndef CGAL_APOLLONIUS_GRAPH_2_EDGE_CONFLICT_2_H
#define CGAL_APOLLONIUS_GRAPH_2_EDGE_CONFLICT_2_H 1


#include <CGAL/Apollonius_graph_2/Delage_traits/Conflict_2.h>

namespace CGAL {

namespace ApolloniusGraph_2 {

//-----------------------------------------------------------------------
//                     Edge Conflict Base
//-----------------------------------------------------------------------

template < class K, class Method_tag >
class Edge_conflict_2 : public Conflict_2<K, Method_tag>
{
private:
  typedef Conflict_2<K, Method_tag>       Base;
public:
  typedef typename Base::Inverted_weighted_point   Inverted_weighted_point;
  typedef bool                                     result_type;
  typedef typename Base::Sign                      Sign;

protected:

    bool edge_conflict_test(const Inverted_weighted_point &p2,
                            const Inverted_weighted_point &p3,
                            const Inverted_weighted_point &p4,
                            const Inverted_weighted_point &q,
                            bool b, int /*i23Q*/, int /*i24Q*/) const
    {

        // orientations
        Sign orient23Q = this->orientation(p2, p3, q);
        Sign orient42Q = this->orientation(p4, p2, q);
        Sign orient234 = this->orientation(p2, p3, p4);

        // radical intersections
        Sign radInt23Q = this->radical_intersection(p2, p3, q, -1);
        Sign radInt24Q = this->radical_intersection(p2, p4, q, -1);

        // radical side
        Sign radSid2Q3 = this->radical_side(p2, q, p3, -1);
        Sign radSid2Q4 = this->radical_side(p2, q, p4, -1);

        // order of a line
        bool oolQ24 = this->ordered_on_line(q, p2, p4);
        bool oolQ23 = this->ordered_on_line(q, p2, p3); 

        if ( b )
        {
	  if ( CGAL::sign(q.p()) != POSITIVE ) { return true; }
	  // degenerate case
	  if (orient234 == ZERO && orient23Q == ZERO && orient42Q == ZERO) {
	    return (oolQ23 || oolQ24);
	  } else if (! ((radInt23Q != NEGATIVE && radSid2Q3 == NEGATIVE) && 
			(radInt24Q != NEGATIVE && radSid2Q4 == NEGATIVE))) {
	    // non degenerate case
	    return true;
	  }  else if (orient234 != NEGATIVE) {
	    return orient23Q != POSITIVE && orient42Q != POSITIVE;
	  } else {
	    return orient23Q != POSITIVE || orient42Q != POSITIVE;
	  }
        }
        else
        {
	  CGAL_assertion ( CGAL::sign(q.p()) == POSITIVE );
	  // degenerate case
	  if (orient234 == ZERO && orient23Q == ZERO && orient42Q == ZERO) {
	    return (oolQ23 && oolQ24);
	  } else if (! ((radInt23Q != NEGATIVE && radSid2Q3 == NEGATIVE) && 
                        (radInt24Q != NEGATIVE && radSid2Q4 == NEGATIVE))) {
            // non degenerate case	
	    return false;
	  } else if (orient234 != NEGATIVE) {
	    return orient23Q != POSITIVE || orient42Q != POSITIVE;
	  } else {
	    return orient23Q != POSITIVE && orient42Q != POSITIVE;
	  }
        }
    }
};

} //namespace ApolloniusGraph_2

} //namespace CGAL

#endif // CGAL_APOLLONIUS_GRAPH_2_EDGE_CONFLICT_2_H
