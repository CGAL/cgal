#ifndef CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H

#include <CGAL/basic.h>

#include <algorithm>
#include <functional>
#include <my_limits>

CGAL_BEGIN_NAMESPACE

// trivial default implementations
class UniqueNumbers {
    static unsigned int n;
    unsigned int i;
public:
    UniqueNumbers() : i( n++ ) {}
    unsigned int num() const { return i; }
};
unsigned int UniqueNumbers::n = 0;

template< class T, unsigned int DIM >
struct Default_Bbox_d : public UniqueNumbers {
    typedef T NumberType;
    Default_Bbox_d() { init(); }
    Default_Bbox_d( T l[DIM], T h[dim] ) {
        std::copy( l, l + DIM, _lo );
        std::copy( h, h + DIM, _hi );
    }

    void init ( bool complete = false ) {
        T inf = workaround::numeric_limits< T >::inf();
        T sup = workaround::numeric_limits< T >::sup();
        if( !complete )
                std::swap( inf, sup );
        std::fill( _lo, _lo + DIM, inf );
        std::fill( _hi, _lo + DIM, sup );
    }

    void extend( T p[ DIM ] ) {
        for( unsigned int dim = 0; dim < DIM; ++dim ) {
            _lo [ DIM ] = std::min( _lo[ DIM ], p[ DIM ] );
            _hi [ DIM ] = std::max( _hi[ DIM ], p[ DIM ] );
        }
    }

    T lo( unsigned int dim ) const { return _lo[dim]; }
    T hi( unsigned int dim ) const { return _hi[dim]; }
    static unsigned int get_dim() { return DIM; }
protected:
    T _lo[ DIM ], _hi[ DIM ];
};

template< class _Box >
struct Default_Bbox_d_Adapter {
    typedef _Box Box;
    typedef typename _Box::NumberType NumberType;

    static NumberType get_lo( const Box& b, unsigned int dim )
    { return b.lo( dim ); }

    static NumberType get_hi( const Box& b, unsigned int dim )
    { return b.hi( dim ); }

    static unsigned int get_num( const Box& b )
    { return b.num();     }

    static unsigned int get_dim() { return Box::get_dim(); }
};


// BoxAdapter has to provide following static members:
// NumberType get_lo( Box, int dim )
// NumberType get_hi( Box, int dim )
// unsigned int get_num( Box )
// Box may be of type immediate, reference, or pointer

template< class BoxAdapter, bool closed >
struct Default_Box_Traits : public BoxAdapter {
    typedef typename BoxAdapter::Box         Box;
    typedef typename BoxAdapter::NumberType  NumberType;
    typedef Default_Box_Traits< BoxAdapter, closed > Traits;

    static unsigned int cutoff;

    static bool hi_greater ( NumberType hi, NumberType val )
    { return closed ? hi >= val : hi > val; }

    // cmp dim = \a b -> islolesslo a b dim
    class Compare : public std::binary_function< Box, Box, bool > {
        unsigned int dim;
    public:
        Compare( unsigned int dim ) : dim( dim ) {}
        bool operator()( const Box& a, const Box& b ) const
        { return Traits::is_lo_less_lo( a, b, dim );  }
    };

    // loless val dim = \box -> getlo box dim < val
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
        { return hi_greater( get_hi( box, dim ), value); }
    };

    // p lo hi dim = \box -> getlo box dim < lo && gethi box dim > hi
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
    { return hi_greater( get_hi(b,dim), get_lo(a,dim )); }

    static bool does_intersect ( const Box& a, const Box& b, unsigned int dim )
    { return is_lo_less_hi( a, b, dim ) && is_lo_less_hi( b, a, dim ); }

    static bool contains_lo_point(const Box& a, const Box& b, unsigned int dim)
    { return !is_lo_less_lo( b, a, dim ) && is_lo_less_hi( b, a, dim );  }

    static unsigned int get_cutoff()
    { return cutoff; }
};

template< class BoxAdapter, bool closed >
unsigned int
Default_Box_Traits<BoxAdapter,closed>::cutoff = 3000;

CGAL_END_NAMESPACE

#endif
