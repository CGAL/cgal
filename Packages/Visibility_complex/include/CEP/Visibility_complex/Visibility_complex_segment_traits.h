#ifndef VISIBILITY_COMPLEX_SEGMENT_TRAITS_H
#define VISIBILITY_COMPLEX_SEGMENT_TRAITS_H


#include <CGAL/basic.h>
#include <CEP/Visibility_complex/Bitangent_2.h>
#include <CEP/Visibility_complex/Arc_2.h>

CGAL_BEGIN_NAMESPACE

// -----------------------------------------------------------------------------

template < class R_ >
class Visibility_complex_segment_traits
{
public:
    // -------------------------------------------------------------------------
    typedef R_                            R;
    typedef typename R::FT                FT;
    typedef typename R::Point_2           Point_2;
    typedef typename R::Segment_2         Segment_2;
    typedef Segment_2                     Disk;
    typedef Arc_2<Disk>               Arc_2;
    typedef Bitangent_2<Disk>         Bitangent_2;
    // -------------------------------------------------------------------------
    // The chi2 predicate
    struct Orientation_object {
	Orientation operator()(const Bitangent_2& a,const Bitangent_2& b) const{ 
	    return orientation(a.source() , a.target() ,
			       a.source() + (b.target() - b.source()));
	}	
    };
    // -------------------------------------------------------------------------
    // The two follwing give the chi2 predicate with a point at infinity
    struct Compare_extreme_yx {
	const Point_2& extreme_point(bool b, const Disk& c) const { 
	    
	    Comparison_result comp = compare_y(c.source(), c.target());
	    comp = (comp != EQUAL) ? comp : compare_x(c.source(), c.target());

	    return ((b == true  && comp == SMALLER) || 
		    (b == false && comp == LARGER)) ?  c.source() : c.target() ;
	}
	const Point_2& extreme_point(bool b, const Bitangent_2& c) const 
	{ return (b) ? c.source() : c.target(); }
	template < class C_ , class D_ >
	Comparison_result operator() (bool sa , const C_& a,
				      bool sb , const D_& b) const { 
	    const Point_2& ap = extreme_point(sa,a);
	    const Point_2& bp = extreme_point(sb,b);

	    Comparison_result cr = compare_y(ap,bp);
	    cr = (cr != EQUAL) ? cr : compare_x(ap,bp);
	    return cr;
	    
	}
    };
    // -------------------------------------------------------------------------
    struct Is_upward_directed {
	bool operator()(const Bitangent_2& b) const {
	    Comparison_result comp = compare_y(b.source(), b.target());
	  comp = (comp != EQUAL) ? comp : compare_x(b.source(), b.target());
	    return (comp != LARGER);
	}
    };
    // -------------------------------------------------------------------------
    // The chi3 predicate
    struct Orientation_infinite {
	const Point_2& top(const Disk& c) const {

	  Comparison_result comp = compare_y(c.source(), c.target());
	  comp = (comp != EQUAL) ? comp : compare_x(c.source(), c.target());

	    return (comp == SMALLER) ? c.target() : c.source();
	}
	const Point_2& bot(const Disk& o) const 
	{ return (top(o) == o.source()) ? o.target() : o.source(); }
	Orientation operator() (const Bitangent_2& a, 
				const Disk& o) const
	{ return orientation(a.source(),a.target(),top(o)); } 
	Orientation operator() (const Disk& o, 
				const Bitangent_2& b) const
	{ return orientation(bot(o),top(o),b.target()); }
	Orientation operator() (const Bitangent_2& a, 
				const Bitangent_2& b) const
	{ return orientation(a.source(),a.target(),b.target()); } 
    };
    // -------------------------------------------------------------------------
    // Detection of degenerate cases
    struct Equal_as_segments {
	bool operator() (const Bitangent_2& a, const Bitangent_2& b) const {
	    return (a.source() == b.source() && a.target() == b.target());
	}
    };
    struct Is_point {
	bool operator() (const Disk& c) const 
	{ return c.source() == c.target(); }
    };
    // -------------------------------------------------------------------------
    struct Do_intersect {
	bool operator()(const Disk& o1, const Disk& o2) {
	    return (o1 == o2);
	}
	bool operator()(const Bitangent_2& o1, const Disk& o2) {
	    return false;		
	}
	bool operator()(const Disk& o1, const Bitangent_2& o2) {
	    return false;
	}
	bool operator()(const Bitangent_2& b1, const Bitangent_2& b2) {
	    return false;
	}
    };
    // -------------------------------------------------------------------------
};

// ----------------------------------------------------------------------------- 

CGAL_END_NAMESPACE

#endif // VISIBILITY_COMPLEX_SEGMENT_TRAITS_H

