#ifndef CGAL_POLYGON_2_POLYGON_2_INTERSECTION_H
#define CGAL_POLYGON_2_POLYGON_2_INTERSECTION_H

#include <CGAL/basic.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------
template < class Traits , class Container >
bool do_intersect( const Polygon_2<Traits,Container>& P,
		   const Polygon_2<Traits,Container>& Q )
{
    // -------------------------------------------------------------------------
    // This function assumes the polygons are convex.
    CGAL_precondition(P.is_convex() && Q.is_convex());
    // -------------------------------------------------------------------------
    // If P or Q is a point we use the CGAL bounded_side methods
    if (P.size() == 1) 
	return (Q.bounded_side(*P.vertices_begin()) != CGAL::ON_UNBOUNDED_SIDE);
    if (Q.size() == 1) 
	return (P.bounded_side(*Q.vertices_begin()) != CGAL::ON_UNBOUNDED_SIDE);
    // -------------------------------------------------------------------------
    // two booleans to avoid cycling more than twice on a Polygon.
    bool a_has_cycled = false;
    bool b_has_cycled = false;
    // -------------------------------------------------------------------------
    typedef Polygon_2<Traits,Container> Polygon_2;
    typedef typename Polygon_2::Vertex_const_iterator Vertex_const_iterator;
    Vertex_const_iterator a = P.vertices_begin();
    Vertex_const_iterator b = Q.vertices_begin();
    // -------------------------------------------------------------------------
    do {
	// ---------------------------------------------------------------------
	if (a == P.vertices_end()) {
	    a = P.vertices_begin(); a_has_cycled = true; 
	}
	if (b == Q.vertices_end()) {
	    b = Q.vertices_begin(); b_has_cycled = true; 
	}
	// ---------------------------------------------------------------------
	// Computations of key variables. 
	Vertex_const_iterator apred , bpred;
	apred = (a == P.vertices_begin()) ? P.vertices_end() : a; --apred;
	bpred = (b == Q.vertices_begin()) ? Q.vertices_end() : b; --bpred;
	// ---------------------------------------------------------------------
	// If A & B intersect, return true.
	typedef typename Polygon_2::Segment_2 Segment_2;
	/*
	cout << "apred->Ptr() = " << long(apred->Ptr()) << endl;
	cout << "bpred->Ptr() = " << long(bpred->Ptr()) << endl;
	cout << "a->Ptr() = " << long(a->Ptr()) << endl;
	cout << "b->Ptr() = " << long(b->Ptr()) << endl;
	*/
	Segment_2 s1(*apred,*a);
	Segment_2 s2(*bpred,*b);
	if ( do_intersect(s1,s2) ) return true;
	// ---------------------------------------------------------------------
	// Else Advance 
	// ---------------------------------------------------------------------
	Orientation chi2  = orientation( *apred , *a , *b + (*apred - *bpred) );
	Orientation chi1b = orientation( *bpred, *b, *a );
	Orientation chi1a = orientation( *apred, *a, *b );
	// ---------------------------------------------------------------------
	if ( chi2 == COLLINEAR ) {
	    if ( chi1b == RIGHTTURN && chi1a == RIGHTTURN ) return false;
	    if ( chi1b == COLLINEAR && chi1a != RIGHTTURN ) ++a;
	    else ++b;
	}
	else if ( chi2 == LEFTTURN ) {
	    if ( chi1a == LEFTTURN ) ++a; 
	    else                     ++b;
	}
	else {
	    if ( chi1b == LEFTTURN ) ++b;
	    else                     ++a;
	}
	// ---------------------------------------------------------------------
    } while ( (a != P.vertices_end() || !a_has_cycled) &&
	      (b != Q.vertices_end() || !b_has_cycled) );
    // -------------------------------------------------------------------------
    return ((P.bounded_side(*Q.vertices_begin()) != CGAL::ON_UNBOUNDED_SIDE) ||
	    (Q.bounded_side(*P.vertices_begin()) != CGAL::ON_UNBOUNDED_SIDE));
    // -------------------------------------------------------------------------
    //return false;
    // -------------------------------------------------------------------------
}
// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif
