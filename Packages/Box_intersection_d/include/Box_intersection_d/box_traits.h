#ifndef CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H

#include <algorithm>
#include <functional>
#include <limits>


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
        NumberType l, h;
        if ( std::numeric_limits< T >::has_infinity ) {
            h = std::numeric_limits< T >::infinity();
            l = -h;
        } else {
            l = std::numeric_limits< T >::min();
            h = std::numeric_limits< T >::max();
        }
        if( !complete )
            std::swap( l, h );
        std::fill( _lo, _lo + DIM, l );
        std::fill( _hi, _lo + DIM, h );
    }

    void extend( T p[ DIM ] ) {
        for( unsigned int dim = 0; dim < DIM; ++dim ) {
            _lo [ DIM ] = std::min( _lo[ DIM ], p[ DIM ] );
            _hi [ DIM ] = std::max( _hi[ DIM ], p[ DIM ] );
        }
    }

    T lo( unsigned int dim ) const { return _lo[dim]; }
    T hi( unsigned int dim ) const { return _hi[dim]; }
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
};


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
        { return get_hi( box, dim ) > value; }
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
