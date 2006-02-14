#ifndef CGAL_BBOX_BOX_TRAITS_H
#define CGAL_BBOX_BOX_TRAITS_H

#include <algorithm>
#include <functional>

// BoxAdapter has to provide following static members:
// NumberType get_lo( Box, int dim )
// NumberType get_hi( Box, int dim )
// unsigned int get_num( Box )
// Box may be of type immediate, reference, or pointer

template< class BoxAdapter >
struct Default_Box_Traits : public BoxAdapter {
    typedef typename BoxAdapter::Box         Box;
    typedef typename BoxAdapter::NumberType  NumberType;
    typedef Default_Box_Traits< BoxAdapter > Traits;

    static unsigned int cutoff;

    class Compare : public std::binary_function< Box, Box, bool > {
        unsigned int dim;
    public:
        Compare( unsigned int dim ) : dim( dim ) {}
        bool operator()( const Box& a, const Box& b ) const
        { return Traits::is_lo_less_lo( a, b, dim );  }
    };

    class Lo_Less : public std::unary_function< Box, bool > {
        NumberType value;
        unsigned int dim;
    public:
        Lo_Less( NumberType value, unsigned int dim )
            : value( value ), dim( dim ) {}
        bool operator() ( const Box& box ) const
        { return get_lo( box, dim ) < value; }
    };

    class Hi_Greater : public std::unary_function< Box, bool > {
        NumberType value;
        unsigned int dim;
    public:
        Hi_Greater( NumberType value, unsigned int dim )
            : value( value ), dim( dim ) {}
        bool operator() ( const Box& box ) const
        { return get_hi( box, dim ) > value; }
    };

    // lambda( box ).(get_lo(box,dim) < lo && get_hi(box,dim)  > hi ) )
    class Interval_Spanning_Predicate : public std::unary_function<Box,bool> {
        NumberType lo, hi;
        unsigned int dim;
    public:
        Interval_Spanning_Predicate( NumberType lo, NumberType hi,
                                     unsigned int dim )
            : lo( lo ), hi( hi ), dim( dim ) {}
        // returns true <=> box spans [lo,hi) in dimension dim
        bool operator() ( const Box& box ) const
        { return get_lo( box, dim ) < lo && get_hi( box, dim ) > hi; }
    };

    static bool is_lo_less_lo( const Box& a, const Box& b, unsigned int dim ) {
        return get_lo(a,dim)  < get_lo(b,dim) ||
               get_lo(a,dim) == get_lo(b,dim) && get_num(a) < get_num(b);
    }

    static bool is_lo_less_hi( const Box& a, const Box& b, unsigned int dim )
    { return get_lo(a,dim ) <  get_hi(b,dim); }

    static bool does_intersect ( const Box& a, const Box& b, unsigned int dim )
    { return get_hi(a,dim) > get_lo(b,dim) && get_hi(b,dim) > get_lo(a,dim); }

    static bool contains_lo_point(const Box& a, const Box& b, unsigned int dim)
    { return !is_lo_less_lo( b, a, dim ) && is_lo_less_hi( b, a, dim );  }

    static unsigned int get_cutoff()
    { return cutoff; }
};

template< class BoxAdapter >
unsigned int
Default_Box_Traits<BoxAdapter>::cutoff = 3000;

#endif
