#ifndef CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H

#include <CGAL/basic.h>

#include <algorithm>
#include <functional>
#include <my_limits>

CGAL_BEGIN_NAMESPACE

struct Box_intersection_d {
    enum Setting { COMPLETE, BIPARTITE };
    enum Topology { HALF_OPEN, CLOSED };


    struct Unique_numbers {
        Unique_numbers() : i(n++) {}
        unsigned int get_num() const { return i; }
    private:
        static unsigned int n;
        unsigned int i;
    };
    /*
    template< class Box >
    struct Select_box_traits {
        typedef char True;
        typedef struct { char a[2]; } False;

        template< class T >
        True is_derived_from_base( typename
          T::Unique_default_box_d_type_tag_which_does_not_appear_elsewhere=0);
        template< class T >
        False is_derived_from_base(...);

        template< bool, typename A, typename B >
        struct IF                { typedef A R; };
        template<typename A, typename B>
        struct IF< false, A, B > { typedef B R; };

        typedef typename IF
           < sizeof(is_derived_from_base<Box>()) == sizeof(True),
             Default_box_d_traits< Box >,
             void
           >::R Box_traits;
    };*/

};

unsigned int Box_intersection_d::Unique_numbers::n = 0;

template< class NT, unsigned int N >
struct Default_box_d : public Box_intersection_d::Unique_numbers {
    typedef int Unique_default_box_d_type_tag_which_does_not_appear_elsewhere;
    typedef NT Number_type;
    Default_box_d() { init(); }
    Default_box_d( NT l[N], NT h[N] ) {
        std::copy( l, l + N, lo );
        std::copy( h, h + N, hi );
    }

    void init ( bool complete = false ) {
        NT inf = workaround::numeric_limits< NT >::inf();
        NT sup = workaround::numeric_limits< NT >::sup();
        if( !complete )
            std::swap( inf, sup );
        std::fill( lo, lo + N, inf );
        std::fill( hi, hi + N, sup );
    }

    void extend( NT p[ N ] ) {
        for( unsigned int dim = 0; dim < N; ++dim ) {
            lo[dim] = std::min( lo[dim], p[dim] );
            hi[dim] = std::max( hi[dim], p[dim] );
        }
    }

    NT get_lo( unsigned int dim ) const { return lo[dim]; }
    NT get_hi( unsigned int dim ) const { return hi[dim]; }
    static unsigned int get_dim() { return N; }
protected:
    NT lo[N], hi[N];
};


template< class Box_ >
struct Default_box_d_traits {
    typedef Box_ Box;
    typedef typename Box::Number_type Number_type;

    static Number_type get_lo( const Box& b, unsigned int dim )
    { return b.get_lo( dim ); }

    static Number_type get_hi( const Box& b, unsigned int dim )
    { return b.get_hi( dim ); }

    static unsigned int get_num( const Box& b )
    { return b.get_num();     }

    static unsigned int get_dim() { return Box::get_dim(); }
};

template< class BoxTraits, bool closed >
struct Default_box_predicate_traits : public BoxTraits {
    typedef typename BoxTraits::Box             Box;
    typedef typename BoxTraits::Number_type     Number_type;

    template<bool b> struct type_from_bool {};

    static bool hi_greater ( Number_type hi, Number_type val, type_from_bool<true> b )
    { return hi >= val; }
    static bool hi_greater ( Number_type hi, Number_type val, type_from_bool<false> b )
    { return hi >  val; }
    static bool hi_greater ( Number_type hi, Number_type val )
    { return hi_greater(hi, val, type_from_bool<closed>()); }

    // cmp dim a b = islolesslo a b dim
    class Compare : public std::binary_function< Box, Box, bool > {
        unsigned int dim;
    public:
        Compare( unsigned int dim ) : dim( dim ) {}
        bool operator()( const Box& a, const Box& b ) const
        { return is_lo_less_lo( a, b, dim );  }
    };

    // loless val dim box = getlo box dim < val
    class Lo_Less : public std::unary_function< Box, bool > {
        Number_type value;
        unsigned int dim;
    public:
        Lo_Less( Number_type value, unsigned int dim )
            : value( value ), dim( dim ) {}
        bool operator() ( const Box& box ) const
        { return get_lo( box, dim ) < value; }
    };

    class Hi_Greater : public std::unary_function< Box, bool > {
        Number_type value;
        unsigned int dim;
    public:
        Hi_Greater( Number_type value, unsigned int dim )
            : value( value ), dim( dim ) {}
        bool operator() ( const Box& box ) const
        { return hi_greater( get_hi( box, dim ), value); }
    };

    // p lo hi dim box = getlo box dim < lo && gethi box dim > hi
    class Interval_Spanning_Predicate : public std::unary_function<Box,bool> {
        Number_type lo, hi;
        unsigned int dim;
    public:
        Interval_Spanning_Predicate( Number_type lo, Number_type hi,
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
    { return is_lo_less_lo( a, b, dim ) && is_lo_less_hi( b, a, dim );  }
};

CGAL_END_NAMESPACE

#endif
