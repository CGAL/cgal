#ifndef VISIBILITY_COMPLEX_FUNCTION_OBJECTS
#define VISIBILITY_COMPLEX_FUNCTION_OBJECTS

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class _Gtr >
struct Less_bitangent {
    typedef typename _Gtr::Bitangent_2 Bitangent_2;
    bool operator()(const Bitangent_2& a, const Bitangent_2& b) const {
	// ---------------------------------------------------------------------
	if (a == b) return false;
	// ---------------------------------------------------------------------
	// If this and b have the same source and target object we sometimes
	// don't need to compute a determinant. This solves some degenerate
	// cases.
	if (a.source_object() == b.source_object() && 
	    a.target_object() == b.target_object()) {
	    if (a.is_right_left() || b.is_left_right()) return true;
	    if (a.is_left_right() || b.is_right_left()) return false;
	    // Case where the left-left and right-right bitangents are
	    // geometrically equal. Arbitrary choice ll < rr.
	    typename _Gtr::Equal_as_segments equal_as_segments;
	    if (equal_as_segments(a,b)) return (a.is_left_left());
	}
	// ---------------------------------------------------------------------
	// Chi2 predicate from geometric traits class
	typename _Gtr::Orientation_object orientation;
	Orientation chi2ab = orientation(a,b);
	// ---------------------------------------------------------------------
	// Non degenerate cases.
	if (chi2ab == RIGHTTURN){ return false;}
	if (chi2ab == LEFTTURN) { return true; }
	// ---------------------------------------------------------------------
	// Treatment of degenerate case when two bitangents are parallel.
	// ---------------------------------------------------------------------
	typename _Gtr::Is_upward_directed is_upward_directed;
	bool upa = is_upward_directed(a); 
	bool upb = is_upward_directed(b);
	// ---------------------------------------------------------------------
	// Treating the case where a and b do not have the same angle.
	// Their angle difference is \pi
	if (!upa && upb)       return false;
	else if (upa && !upb)  return true;
	// Now a and b have the same angle
	// ---------------------------------------------------------------------
	// a and b have distinct types
	if (a.type() != b.type()) {
	    if (a.is_right_left() || b.is_left_right())   return true;
	    if (a.is_left_right() || b.is_right_left())   return false;
	    return (a.is_left_left());
	}
	// ---------------------------------------------------------------------
	// Final case : a and b have the same angle and same type
	typename _Gtr::Orientation_infinite chi3;
	Orientation op = chi3(a,b);
	if (op == RIGHTTURN)     return false;
	else if (op == LEFTTURN) return true;
	else {
	    if (a.target_object() == b.source_object() &&
		a.is_xx_left()    == b.is_left_xx())    return a.is_xx_left();
	    if (a.source_object() == b.target_object() &&
		a.is_left_xx()    == b.is_xx_left())    return a.is_right_xx();
	    typename _Gtr::Compare_extreme_yx compare;
	    if (a.target_object() == b.target_object()) // a or b is not free
		return (compare(a.is_left_xx(),*a.source_object(),
			        b.is_left_xx(),*b.source_object()) == LARGER);
	    if (a.source_object() == b.source_object()) // a or b is not free
		return (compare(a.is_left_xx(),*a.source_object(),
				b.is_left_xx(),*b.source_object()) == SMALLER);
	    Comparison_result cmp;
	    if (upa) cmp = compare(a.is_xx_left(),*a.target_object(),          
				   b.is_xx_left(),*b.target_object());
	    else     cmp = compare(!b.is_left_xx(),*b.source_object(),
				   !a.is_left_xx(),*a.source_object());
	    CGAL_precondition(cmp != EQUAL);
	    return (a.is_xx_right()) ? (cmp == LARGER) : (cmp == SMALLER);
	}
	// ---------------------------------------------------------------------
    }
};

template < class _Gtr >
struct Greater_bitangent {
    typedef typename _Gtr::Bitangent_2 Bitangent_2;
    bool operator()(const Bitangent_2& a, const Bitangent_2& b) {
	if (a == b) return false;
	return Less_bitangent<_Gtr>()(b,a);
    }
};

// -----------------------------------------------------------------------------

template <class _Gtr>
class Less_object {
private:
    // -------------------------------------------------------------------------
    typedef typename _Gtr::Disk             Disk;
    typedef typename _Gtr::Bitangent_2          Bitangent_2;
  typedef Bitangent_2 BT;
    typedef typename _Gtr::Is_upward_directed   Is_upward_directed;
    typedef typename _Gtr::Compare_extreme_yx   Compare;
    typedef typename _Gtr::Orientation_infinite Orientation_infinite;
    // -------------------------------------------------------------------------
public:
    // -------------------------------------------------------------------------
    bool operator() (const Disk& o1, const Disk& o2) const {
	// ---------------------------------------------------------------------
	if (&o1 == &o2) return false;
	BT rl(BT::RL,&o1,&o2);
	return !Is_upward_directed()(rl);
	// ---------------------------------------------------------------------
    }                      
    // -------------------------------------------------------------------------
    bool operator() (const Disk& o, const Bitangent_2& b) const {
	// ---------------------------------------------------------------------
	Is_upward_directed is_upward_directed;
	CGAL_precondition(is_upward_directed(b));
	// ---------------------------------------------------------------------
	// o is the source or target object of b.
	if (b.source_object() == &o) return b.is_left_xx();
	if (b.target_object() == &o) return b.is_xx_left();
	// ---------------------------------------------------------------------
	// Otherwise o and b involve three convex objects.
	Comparison_result cmp = Compare()(false,o,false,b);
	CGAL_precondition(cmp != EQUAL);
	Orientation_infinite chi3;
	if (cmp == SMALLER) return (chi3(b,o) == LEFTTURN);
	else                return (chi3(o,b) == RIGHTTURN);
	// Find somtething better than the above not using chi3...
	// ---------------------------------------------------------------------
    }                      
    // -------------------------------------------------------------------------
    bool operator() (const Bitangent_2& b, const Disk& o) const 
    { return !operator()(o,b); }                      
    // -------------------------------------------------------------------------
    bool operator() (const Bitangent_2& b1, const Bitangent_2& b2) const {
	// ---------------------------------------------------------------------
	if (b1 == b2) return false;
	// ---------------------------------------------------------------------
	Is_upward_directed is_upward_directed;
	CGAL_precondition(is_upward_directed(b1));
	CGAL_precondition(is_upward_directed(b2));
	// ---------------------------------------------------------------------
	// b1 and b2 involve two objects
	if ((b1.source_object() == b2.source_object() &&
	     b1.target_object() == b2.target_object()) ||
	    (b1.source_object() == b2.target_object() &&
	     b1.target_object() == b2.source_object()))
	    return (b1.is_right_right() || b2.is_left_left());
	// ---------------------------------------------------------------------
	// Some cases when b1 and b2 involve three objects
	//if (b1.target_object() == b2.source_object()) return b2.is_left_xx();
	//if (b1.source_object() == b2.target_object()) return b1.is_right_xx();
	Less_bitangent<_Gtr> chi1;
	if (b1.target_object() == b2.target_object()) {
	    if (b1.is_xx_left() == b2.is_xx_left()) return chi1(b1,b2);
	    else                                    return b1.is_xx_right();
	}
	if (b1.source_object() == b2.source_object()) {
	    if (b1.is_left_xx() == b2.is_left_xx()) return chi1(b2,b1);
	    else                                    return b1.is_right_xx();
	}
	// ---------------------------------------------------------------------
	// The remaining cases.
	Comparison_result cmp = Compare()(false,b1,false,b2);
	CGAL_precondition(cmp != EQUAL);
	Orientation_infinite chi3;
	if (cmp == SMALLER) return (chi3(b2,b1) == LEFTTURN);
	else                return (chi3(b1,b2) == RIGHTTURN);
	// Find somtething better than the above not using chi3...
	// ---------------------------------------------------------------------
    }
    // -------------------------------------------------------------------------
};

// -----------------------------------------------------------------------------

template < class _Vc >
struct Less_face_handle {
    template < class Face_handle > 
    bool operator() (const Face_handle& e1, const Face_handle& e2) const {
	if (e1 == e2 || e1->front_view() == 0) return false;
	if (e2->front_view() == 0) return true;

	typedef typename _Vc::Gt _Gtr; 
	typedef typename _Vc::Edge_handle Edge_handle;
	Edge_handle e = e1->front_view();
	Edge_handle f = e2->front_view();
	if (e->object() == 0 && f->object() == 0) 
	    return Less_object<_Gtr>()(*e->sup(),*f->sup());
	else if (e->object() == 0 && f->object() != 0) 
	    return Less_object<_Gtr>()(*e->sup(),*f->object());
	else if (e->object() != 0 && f->object() == 0)
	    return Less_object<_Gtr>()(*e->object(),*f->sup());
	return Less_object<_Gtr>()(*e->object(), *f->object());
    }                      
};

// -----------------------------------------------------------------------------

template < class _Vc >
struct Do_intersect {
    template < class Edge_handle >
    bool operator()(const Edge_handle& e, const Edge_handle& f) const {
	if (e == 0 || f == 0) return false;
	typedef typename _Vc::Gt Gt;
	typename Gt::Do_intersect do_intersect;

	if (e->object() == 0 && f->object() == 0) 
	    return do_intersect(*e->sup(),*f->sup());
	if (e->object() != 0 && f->object() == 0) 
	    return do_intersect(*e->object(),*f->sup());
	if (e->object() == 0 && f->object() != 0) 
	    return do_intersect(*e->sup(),*f->object());
	return do_intersect(*e->object(),*f->object());
    }
};

// -----------------------------------------------------------------------------

template < class _Gtr >
struct Less_edge_handle {
    typedef bool result_type;
    template < class Edge_handle > 
    bool operator () (const Edge_handle& e1, const Edge_handle& e2) const 
    { 
	// ---------------------------------------------------------------------
	// Trivial cases not needing geometry
	if (e1 == e2) {  return false; }
	if (e1->object() != 0 && e1->object() == e2->object()) 
	    return (e1->sign()  && !e2->sign());
	// ---------------------------------------------------------------------
	// Otherwise use predicate from geometric traits class
	// to compare lexicographically the base points
	typename _Gtr::Compare_extreme_yx compare;
	Comparison_result cmp;
	if (e1->object() != 0) {
	    if (e2->object() != 0) cmp = compare(e1->sign(),*e1->object(),
						 e2->sign(),*e2->object());
	    else cmp = compare(e1->sign(),*e1->object(),
			       e2->sign(),*e2->sup());
	}
	else {
	    if (e2->object() != 0) cmp = compare(e1->sign(),*e1->sup(),
						 e2->sign(),*e2->object());
	    else cmp = compare(e1->sign(),*e1->sup(),
			       e2->sign(),*e2->sup());
	}
	// ---------------------------------------------------------------------
	// distinct base points
	if (cmp != EQUAL) return (cmp == SMALLER) ? true : false;
	// ---------------------------------------------------------------------
	// equal base points and e1 or e2 is a constraint edge
	if (e1->object() != 0 && e2->object() == 0) return e1->sign();
	if (e1->object() == 0 && e2->object() != 0) return !e2->sign();
	if (e1->object() == 0 && e2->object() == 0) {
	    bool e1_is_left = (( e1->sign() && e1->sup()->is_left_xx()) ||
			       (!e1->sign() && e1->sup()->is_xx_left()));
	    bool e2_is_left = (( e2->sign() && e2->sup()->is_left_xx()) ||
			       (!e2->sign() && e2->sup()->is_xx_left()));
	    if (e1_is_left  && !e2_is_left) return true;
	    if (!e1_is_left &&  e2_is_left) return false;
	    Less_bitangent<_Gtr> chi2;
	    if (e1_is_left  &&  e2_is_left) {
		if (*e1->sup() == *e2->sup()) 
		    return (long(e1->sup()) < long(e2->sup()));
		return chi2(*e1->sup(),*e2->sup());
	    }
	    if (!e1_is_left && !e2_is_left) {
		if (*e1->sup() == *e2->sup()) 
		    return (long(e2->sup()) < long(e1->sup()));
		return chi2(*e2->sup(),*e1->sup());
	    }
	}
	// ---------------------------------------------------------------------
	// equal base points and e1, e2 are regular edges
	return Less_object<_Gtr>()(*e1->object(),*e2->object());
	// ---------------------------------------------------------------------
    }
};

// -----------------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_FUNCTION_OBJECTS
