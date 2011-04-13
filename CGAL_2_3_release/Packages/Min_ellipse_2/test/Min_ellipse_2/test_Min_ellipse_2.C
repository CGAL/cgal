// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : test/Min_ellipse_2/test_Min_ellipse_2.C
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Min_ellipse_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>, Bernd Gärtner
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: test program for 2D Smallest Enclosing Ellipse
// ============================================================================

#include <CGAL/Cartesian.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Min_ellipse_2.h>
#include <CGAL/Min_ellipse_2_traits_2.h>
#include <CGAL/Min_ellipse_2_adapterC2.h>
#include <CGAL/Min_ellipse_2_adapterH2.h>
#include <CGAL/IO/Verbose_ostream.h>
#include <cassert>
#include <cstring>
#include <fstream>

#ifdef CGAL_USE_LEDA_FOR_OPTIMISATION_TEST
#  include <CGAL/leda_integer.h>
   typedef  leda_integer                       Rt;
   typedef  CGAL::Quotient< leda_integer >     Ft;
#else
#  include <CGAL/Gmpz.h>
   typedef  CGAL::Gmpz                         Rt;
   typedef  CGAL::Quotient< CGAL::Gmpz >       Ft;
#endif

typedef  CGAL::Cartesian< Ft >                 KerC;
typedef  CGAL::Homogeneous< Rt >               KerH;
typedef  CGAL::Min_ellipse_2_traits_2< KerC >  TraitsC;
typedef  CGAL::Min_ellipse_2_traits_2< KerH >  TraitsH;

// code coverage test function
// ---------------------------
template < class Traits, class RT >
void
cover_Min_ellipse_2( bool verbose, const Traits&, const RT&)
{
    CGAL_USING_NAMESPACE_STD

    typedef  CGAL::Min_ellipse_2< Traits >  Min_ellipse;
    typedef  typename Min_ellipse::Point    Point;
    typedef  typename Min_ellipse::Ellipse  Ellipse;

    CGAL::Verbose_ostream verr( verbose);

    // generate `n' points at random
    const int     n = 20;
    CGAL::Random  random_x, random_y;
    Point         random_points[ n];
    int           i;
    verr << n << " random points from [0,128)^2:" << endl;
    for ( i = 0; i < n; ++i) {
        random_points[ i] = Point( RT( random_x( 128)),
                                   RT( random_y( 128)));
        verr << i << ": " << random_points[ i] << endl;
    }

    // cover code
    verr << endl << "default constructor...";
    {
        Min_ellipse  me;
        bool  is_valid = me.is_valid( verbose);
        bool  is_empty = me.is_empty();
        assert( is_valid);
        assert( is_empty);
    }

    verr << endl << "one point constructor...";
    {
        Min_ellipse  me( random_points[ 0]);
        bool  is_valid      = me.is_valid( verbose);
        bool  is_degenerate = me.is_degenerate();
        assert( is_valid);
        assert( is_degenerate);
    }

    verr << endl << "two points constructor...";
    {
        Min_ellipse  me( random_points[ 1],
                         random_points[ 2]);
        bool  is_valid      = me.is_valid( verbose);
        bool  is_degenerate = me.is_degenerate();
        assert( is_valid);
        assert( is_degenerate);
    }

    verr << endl << "three points constructor...";
    {
        Min_ellipse  me( random_points[ 3],
                         random_points[ 4],
                         random_points[ 5]);
        bool  is_valid = me.is_valid( verbose);
        int   num_pts  = me.number_of_points();
        assert( is_valid);
        assert( num_pts == 3);
    }

    verr << endl << "four points constructor...";
    {
        Min_ellipse  me( random_points[ 6],
                         random_points[ 7],
                         random_points[ 8],
                         random_points[ 9]);
        bool  is_valid = me.is_valid( verbose);
        int   num_pts  = me.number_of_points();
        assert( is_valid);
        assert( num_pts == 4);
    }

    verr << endl << "five points constructor...";
    {
        Min_ellipse  me( random_points[ 10],
                         random_points[ 11],
                         random_points[ 12],
                         random_points[ 13],
                         random_points[ 14]);
        bool  is_valid = me.is_valid( verbose);
        int   num_pts  = me.number_of_points();
        assert( is_valid);
        assert( num_pts == 5);
    }

    verr << endl << "Point* constructor...";
    Min_ellipse  me( random_points, random_points+9, false);
    {
        Min_ellipse  me2( random_points, random_points+9, true);
        bool  is_valid  = me .is_valid( verbose);
        bool  is_valid2 = me2.is_valid( verbose);
        assert( is_valid);
        assert( is_valid2);
        assert( me .number_of_points() == 9);
        assert( me2.number_of_points() == 9);
        assert( me.ellipse() == me2.ellipse());
    }

    verr << endl << "list<Point>::const_iterator constructor...";
    {
        Min_ellipse  me1( me.points_begin(), me.points_end(), false);
        Min_ellipse  me2( me.points_begin(), me.points_end(), true);
        bool  is_valid1 = me1.is_valid( verbose);
        bool  is_valid2 = me2.is_valid( verbose);
        assert( is_valid1);
        assert( is_valid2);
        assert( me1.number_of_points() == 9);
        assert( me2.number_of_points() == 9);
        assert( me.ellipse() == me1.ellipse());
        assert( me.ellipse() == me2.ellipse());
    }

    verr << endl << "#points already called above.";

    verr << endl << "points access already called above.";

    verr << endl << "support points access...";
    {
        typedef  typename Min_ellipse::Support_point_iterator
                                          Support_point_iterator;
        Point                   support_point;
        Support_point_iterator  iter( me.support_points_begin());
        for ( i = 0; i < me.number_of_support_points(); ++i, ++iter) {
            support_point = me.support_point( i);
            assert( support_point == *iter); }
        Support_point_iterator  end_iter( me.support_points_end());
        assert( iter == end_iter);
    }

    verr << endl << "ellipse access already called above...";

    verr << endl << "in-ellipse predicates...";
    {
        Point               p;
        CGAL::Bounded_side  bounded_side;
        bool                has_on_bounded_side;
        bool                has_on_boundary;
        bool                has_on_unbounded_side;
        for ( i = 0; i < 9; ++i) {
            p = random_points[ i];
            bounded_side          = me.bounded_side( p);
            has_on_bounded_side   = me.has_on_bounded_side( p);
            has_on_boundary       = me.has_on_boundary( p);
            has_on_unbounded_side = me.has_on_unbounded_side( p);
        assert( bounded_side != CGAL::ON_UNBOUNDED_SIDE);
        assert( has_on_bounded_side || has_on_boundary);
        assert( ! has_on_unbounded_side); }
    }

    verr << endl << "is_... predicates already called above.";

    verr << endl << "single point insert...";
    me.insert( random_points[ 9]);
    {
        bool  is_valid = me.is_valid( verbose);
        assert( is_valid);
        assert( me.number_of_points() == 10);
    }

    verr << endl << "Point* insert...";
    me.insert( random_points+10, random_points+n);
    {
        bool  is_valid = me.is_valid( verbose);
        assert( is_valid);
        assert( me.number_of_points() == n);
    }

    verr << endl << "list<Point>::const_iterator insert...";
    {
        Min_ellipse  me2;
        me2.insert( me.points_begin(), me.points_end());
        bool  is_valid = me2.is_valid( verbose);
        assert( is_valid);
        assert( me2.number_of_points() == n);

        verr << endl << "clear...";
        me2.clear();
              is_valid = me2.is_valid( verbose);
        bool  is_empty = me2.is_empty();
        assert( is_valid);
        assert( is_empty);
    }

    verr << endl << "validity check already called several times.";

    verr << endl << "traits class access...";
    {
        Traits  traits( me.traits());
    }

    verr << endl << "I/O...";
    {
        verr << endl << "  writing `test_Min_ellipse_2.ascii'...";
        ofstream os( "test_Min_ellipse_2.ascii");
        CGAL::set_ascii_mode( os);
        os << me;
    }
    {
        verr << endl << "  writing `test_Min_ellipse_2.pretty'...";
        ofstream os( "test_Min_ellipse_2.pretty");
        CGAL::set_pretty_mode( os);
        os << me;
    }
    {
        verr << endl << "  writing `test_Min_ellipse_2.binary'...";
        ofstream os( "test_Min_ellipse_2.binary");
        CGAL::set_binary_mode( os);
        os << me;
    }
    {
        verr << endl << "  reading `test_Min_ellipse_2.ascii'...";
        Min_ellipse me_in;
        ifstream is( "test_Min_ellipse_2.ascii");
        CGAL::set_ascii_mode( is);
        is >> me_in;
        bool    is_valid = me_in.is_valid( verbose);
        assert( is_valid);
        assert( me_in.number_of_points() == n);
        assert( me_in.ellipse() == me.ellipse());
    }
    verr << endl;
}

// point classes for adapters test
// -------------------------------
// 2D Cartesian point class
class MyPointC2;

std::ostream&  operator << ( std::ostream&, const MyPointC2&);
std::istream&  operator >> ( std::istream&,       MyPointC2&);

class MyPointC2 {
  public:
    typedef  ::Ft  FT;
  private:
    FT x_;
    FT y_;
  public:
    MyPointC2( ) { }
    MyPointC2( const FT& x, const FT& y) : x_( x), y_( y) { }

    const FT&  x( ) const { return( x_); }
    const FT&  y( ) const { return( y_); }

    bool
    operator == ( const MyPointC2& p) const
    {
        return( ( x_ == p.x_) && ( y_ == p.y_));
    }

    bool
    operator != ( const MyPointC2& p) const
    {
        return( ( x_ != p.x_) || ( y_ != p.y_));
    }

    friend
    std::ostream&  operator << ( std::ostream& os, const MyPointC2& p);

    friend
    std::istream&  operator >> ( std::istream& is,       MyPointC2& p);
};

std::ostream&
operator << ( std::ostream& os, const MyPointC2& p)
{
    return( os << p.x_ << ' ' << p.y_);
}

std::istream&
operator >> ( std::istream& is, MyPointC2& p)
{
    return( is >> p.x_ >> p.y_);
}

// 2D Cartesian point class data accessor
class MyPointC2DA {
  public:
    typedef  ::Ft  FT;

    MyPointC2DA( ) { }

    const FT&  get_x( const MyPointC2& p) const { return( p.x()); }
    const FT&  get_y( const MyPointC2& p) const { return( p.y()); }

    void
    get( const MyPointC2& p, FT& x, FT& y) const
    {
        x = get_x( p);
        y = get_y( p);
    }

    void
    set( MyPointC2& p, const FT& x, const FT& y) const
    {
        p = MyPointC2( x, y);
    }
};


// 2D homogeneous point class
class MyPointH2;

std::ostream&  operator << ( std::ostream&, const MyPointH2&);
std::istream&  operator >> ( std::istream&,       MyPointH2&);

class MyPointH2 {
  public:
    typedef  ::Rt  RT;
  private:
    RT hx_;
    RT hy_;
    RT hw_;
  public:
    MyPointH2( ) { }
    MyPointH2( const RT& hx, const RT& hy, const RT& hw = RT( 1))
        : hx_( hx), hy_( hy), hw_( hw) { }

    const RT&  hx( ) const { return( hx_); }
    const RT&  hy( ) const { return( hy_); }
    const RT&  hw( ) const { return( hw_); }

    bool
    operator == ( const MyPointH2& p) const
    {
        return( ( hx_*p.hw_ == p.hx_*hw_) && ( hy_*p.hw_ == p.hy_*hw_));
    }

    bool
    operator != ( const MyPointH2& p) const
    {
        return( ( hx_*p.hw_ != p.hx_*hw_) || ( hy_*p.hw_ != p.hy_*hw_));
    }

    friend
    std::ostream&  operator << ( std::ostream& os, const MyPointH2& p);

    friend
    std::istream&  operator >> ( std::istream& is,       MyPointH2& p);
};

std::ostream&
operator << ( std::ostream& os, const MyPointH2& p)
{
    return( os << p.hx_ << ' ' << p.hy_ << ' ' << p.hw_);
}

std::istream&
operator >> ( std::istream& is, MyPointH2& p)
{
    return( is >> p.hx_ >> p.hy_ >> p.hw_);
}

// 2D homogeneous point class data accessor
class MyPointH2DA {
  public:
    typedef  ::Rt  RT;

    MyPointH2DA( ) { }

    const RT&  get_hx( const MyPointH2& p) const { return( p.hx()); }
    const RT&  get_hy( const MyPointH2& p) const { return( p.hy()); }
    const RT&  get_hw( const MyPointH2& p) const { return( p.hw()); }

    void
    get( const MyPointH2& p, RT& hx, RT& hy, RT& hw) const
    {
        hx = get_hx( p);
        hy = get_hy( p);
        hw = get_hw( p);
    }

    void
    set( MyPointH2& p, const RT& hx, const RT& hy, const RT& hw) const
    {
        p = MyPointH2( hx, hy, hw);
    }
};

CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( MyPointC2)
CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( MyPointH2)

// main
// ----
int
main( int argc, char* argv[])
{
    // command line options
    // --------------------
    // option `-verbose'
    bool  verbose = false;
    if ( ( argc > 1) && ( CGAL_CLIB_STD::strcmp( argv[ 1], "-verbose") == 0)) {
        verbose = true;
        --argc;
        ++argv; }

    // code coverage
    // -------------
    cover_Min_ellipse_2( verbose, TraitsC(), Rt());
    cover_Min_ellipse_2( verbose, TraitsH(), Rt());

    // adapters test
    // -------------
    typedef CGAL::Min_ellipse_2_adapterC2< MyPointC2, MyPointC2DA >  AdapterC2;
    typedef CGAL::Min_ellipse_2_adapterH2< MyPointH2, MyPointH2DA >  AdapterH2;
    cover_Min_ellipse_2( verbose, AdapterC2(), Rt());
    cover_Min_ellipse_2( verbose, AdapterH2(), Rt());

    // external test sets
    // -------------------
    while ( argc > 1) {
    
        typedef  CGAL::Min_ellipse_2< TraitsH >  Min_ellipse;
        typedef  Min_ellipse::Point              Point;
        typedef  Min_ellipse::Ellipse            Ellipse;
    
        CGAL::Verbose_ostream verr( verbose);
    
        // read points from file
        verr << std::endl << "input file: `" << argv[ 1] << "'" << std::flush;
    
        std::list<Point>  points;
        int               n, x, y;
        std::ifstream     in( argv[ 1]);
        in >> n;
        assert( in);
        for ( int i = 0; i < n; ++i) {
            in >> x >> y;
            assert( in);
            points.push_back( Point( x, y)); }
    
        // compute and check min_ellipse
        Min_ellipse  me2( points.begin(), points.end(), false);
        bool  is_valid = me2.is_valid( verbose);
        assert( is_valid);
    
        // next file
        --argc;
        ++argv; }

    return( 0);
}

// ===== EOF ==================================================================
