#ifndef VISIBILITY_COMPLEX_CIRCLE_TRAITS_H
#define VISIBILITY_COMPLEX_CIRCLE_TRAITS_H

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif

#ifndef CGAL_CIRCLE_BY_RADIUS_2_H
#include <CEP/Visibility_complex/Circle_by_radius_2.h>
#endif

#ifndef CGAL_CONVEX_ARC_2
#include <CEP/Visibility_complex/Arc_2.h>
#endif

#include <CEP/Visibility_complex/Circle_2_Circle_2_intersection.h>
#include <CEP/Visibility_complex/Circle_2_Bitangent_2_intersection.h>
#include <CEP/Visibility_complex/predicates/Visibility_complex_ftC2.h>
#include <CGAL/predicates/kernel_ftC2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class R_ >
struct Visibility_complex_circle_traits
{
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Segment_2         Segment_2;
    typedef Circle_by_radius_2<R>         Disk;
    typedef Arc_2<Disk>        Arc_2;
    typedef Bitangent_2<Disk>         Bitangent_2;
    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    typedef Bitangent_2::Disk_handle Disk_handle;
	    Disk_handle sa(a.source_object()),ta(a.target_object()),
			   sb(b.source_object()),tb(b.target_object());
	    FT ssa = (a.is_left_xx()) ? 1 : -1;
	    FT sta = (a.is_xx_left()) ? 1 : -1;
	    FT ssb = (b.is_left_xx()) ? 1 : -1;
	    FT stb = (b.is_xx_left()) ? 1 : -1;
	    Sign sgn = chi2_testC2(ta->center().x() - sa->center().x(),
				   ta->center().y() - sa->center().y(),
				   sta * ta->radius() - ssa * sa->radius(),
				   tb->center().x() - sb->center().x(),
				   tb->center().y() - sb->center().y(),
				   stb * tb->radius() - ssb * sb->radius());
	    if (sgn == POSITIVE) return LEFTTURN;
	    else if (sgn == NEGATIVE) return RIGHTTURN;
	    return COLLINEAR;
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Bitangent_2& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Bitangent_2& a,
				      bool sb , const Bitangent_2& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Bitangent_2& a,
				      bool sb , const Disk& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Disk& b) const { 
	    FT ar = (sa) ? -a.radius() : a.radius();
	    FT br = (sb) ? -b.radius() : b.radius();
	    return compare_lexicographically_xyC2(a.center().y() + ar,
						  a.center().x(),
						  b.center().y() + br,
						  b.center().x());
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	    Comparison_result comp = 
		compare_lexicographically_xyC2(b.source().y(),b.source().x(),
					       b.target().y(),b.target().x());
	    return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    struct Orientation_infinite {
	// FIXME - not implemented
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const{ return COLLINEAR; }
	// FIXME - not implemented
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const{ return COLLINEAR; } 
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return orientation(a.source(),a.target(),b.target()); } 
    };
    // -------------------------------------------------------------------------
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    if (a == b) return true;
	    if (a.source_object() != b.source_object() ||
		a.target_object() != b.target_object()) return false;
	    if (a.source_object()->radius() == 0 &&
		a.target_object()->radius() == 0) return true;
	    if (a.source_object()->radius() == 0)
		return (a.is_xx_left() == b.is_xx_left());
	    if (a.target_object()->radius() == 0)
		return (a.is_left_xx() == b.is_left_xx());
	    return false;
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const 
	{ return (c.radius() == 0); }
    };
    // -------------------------------------------------------------------------
    // Intersection test. Optional
    typedef Tag_true supports_intersection;
    struct Do_intersect {
	bool operator()(const Disk& o1, const Disk& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& o1, const Disk& o2) {
	    return do_intersect(o2,o1);
	}
	bool operator()(const Disk& o1, const Bitangent_2& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& b1, const Bitangent_2& b2) {
	    // FIXME !!! - implement this
	    return false;
	}
    };
    // -------------------------------------------------------------------------
};


// -----------------------------------------------------------------------------

template < class R_ >
class Visibility_complex_circle_expensive_traits
{
public:
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Segment_2         Segment_2;
    typedef Circle_by_radius_2<R>         Disk;
    typedef Arc_2<Disk>        Arc_2;
    typedef Bitangent_2<Disk>         Bitangent_2;
    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    /*
	    return orientation(a.source() , a.target() ,
			       a.source() + (b.target() - b.source()));
	    */
	    typedef Bitangent_2::Disk_handle Disk_handle;
	    Disk_handle sa(a.source_object()),ta(a.target_object()),
			   sb(b.source_object()),tb(b.target_object());
	    FT ssa = (a.is_left_xx()) ? 1 : -1;
	    FT sta = (a.is_xx_left()) ? 1 : -1;
	    FT ssb = (b.is_left_xx()) ? 1 : -1;
	    FT stb = (b.is_xx_left()) ? 1 : -1;
	    Sign sgn = 
		chi2_test_expensiveC2(ta->center().x() - sa->center().x(),
				      ta->center().y() - sa->center().y(),
				      sta * ta->radius() - ssa * sa->radius(),
				      tb->center().x() - sb->center().x(),
				      tb->center().y() - sb->center().y(),
				      stb * tb->radius() - ssb * sb->radius());
	    if (sgn == POSITIVE) return LEFTTURN;
	    else if (sgn == NEGATIVE) return RIGHTTURN;
	    return COLLINEAR;
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Bitangent_2& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Bitangent_2& a,
				      bool sb , const Bitangent_2& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Bitangent_2& a,
				      bool sb , const Disk& b) const 
	{ return EQUAL; } // FIXME - not implemented
	Comparison_result operator() (bool sa , const Disk& a,
				      bool sb , const Disk& b) const { 
	    FT ar = (sa) ? a.radius() : -a.radius();
	    FT br = (sb) ? b.radius() : -b.radius();
	    return compare_lexicographically_xyC2(a.center().y() + ra,
						  a.center().x(),
						  b.center().y() + rb,
						  b.center().x);
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	    Comparison_result comp = 
		compare_lexicographically_xyC2(b.source().y(),b.source().x(),
					       b.target().y(),b.target().x());
	    return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    struct Orientation_infinite {
	// FIXME - not implemented
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const{ return COLLINEAR; }
	// FIXME - not implemented
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const{ return COLLINEAR; } 
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return orientation(a.source(),a.target(),b.target()); } 
    };
    // -------------------------------------------------------------------------
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    if (a == b) return true;
	    if (a.source_object() != b.source_object() ||
		a.target_object() != b.target_object()) return false;
	    if (a.source_object()->radius() == 0 &&
		a.target_object()->radius() == 0) return true;
	    if (a.source_object()->radius() == 0)
		return (a.is_xx_left() == b.is_xx_left());
	    if (a.target_object()->radius() == 0)
		return (a.is_left_xx() == b.is_left_xx());
	    return false;
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const { return (c.radius() == 0); }
    };
    // -------------------------------------------------------------------------
    // Intersection test. Optional
    typedef Tag_true supports_intersection;
    struct Do_intersect {
	bool operator()(const Disk& o1, const Disk& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& o1, const Disk& o2) {
	    return do_intersect(o2,o1);
	}
	bool operator()(const Disk& o1, const Bitangent_2& o2) {
	    return do_intersect(o1,o2);
	}
	bool operator()(const Bitangent_2& b1, const Bitangent_2& b2) {
	    // FIXME !!! - not implemented
	    return false;
	}
    };
    // -------------------------------------------------------------------------
};

// ----------------------------------------------------------------------------- 

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_CIRCLE_TRAITS_H
