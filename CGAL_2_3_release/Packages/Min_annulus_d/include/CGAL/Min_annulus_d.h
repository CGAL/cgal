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
// file          : include/CGAL/Min_annulus_d.h
// package       : $CGAL_Package: Min_annulus_d $
// chapter       : Geometric Optimisation
//
// source        : web/Min_annulus_d.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Smallest enclosing annulus in arbitrary dimension
// ============================================================================

#ifndef CGAL_MIN_ANNULUS_D_H
#define CGAL_MIN_ANNULUS_D_H

// includes
// --------
#ifndef CGAL_OPTIMISATION_BASIC_H
#  include <CGAL/Optimisation/basic.h>
#endif

#ifndef CGAL_FUNCTION_OBJECTS_H
#  include <CGAL/function_objects.h>
#endif
#ifndef CGAL_IDENTITY_H
#  include <CGAL/_QP_solver/identity.h>
#endif

#ifndef CGAL_FUNCTION_OBJECTS_ACCESS_BY_INDEX_H
#  include <CGAL/_QP_solver/Access_by_index.h>
#endif
#ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
#  include <CGAL/_QP_solver/Join_random_access_iterator.h>
#endif
#ifndef CGAL_QP_SOLVER_H
#  include <CGAL/_QP_solver/QP_solver.h>
#endif
#ifndef CGAL_JOIN_RANDOM_ACCESS_ITERATOR_H
#  include <CGAL/_QP_solver/Join_random_access_iterator.h>
#endif
#ifndef CGAL_PARTIAL_EXACT_PRICING_H
#  include <CGAL/_QP_solver/Partial_exact_pricing.h>
#endif
#ifndef CGAL_PARTIAL_FILTERED_PRICING_H
#  include <CGAL/_QP_solver/Partial_filtered_pricing.h>
#endif

#ifndef CGAL_PROTECT_VECTOR
#  include <vector>
#  define CGAL_PROTECT_VECTOR
#endif
#ifndef CGAL_PROTECT_IOSTREAM
#  include <iostream>
#  define CGAL_PROTECT_IOSTREAM
#endif
#ifndef CGAL_PROTECT_FUNCTIONAL_H
#  include <functional>
#  define CGAL_PROTECT_FUNCTIONAL_H
#endif


CGAL_BEGIN_NAMESPACE

// Class declarations
// ==================
template < class Traits_ >
class Min_annulus_d;

template < class ET_, class NT_, class Point, class PointIterator,
           class Access_coord, class Access_dim >
struct LP_rep_min_annulus_d;

template < class NT >
struct LP_rep_row_of_a {
    typedef  std::vector<NT>                         argument_type;
    typedef  typename argument_type::const_iterator  result_type;

    result_type
    operator ( ) ( const argument_type& v) const { return v.begin(); }
};


// Class interfaces
// ================
template < class Traits_ >
class Min_annulus_d {
  public:
    // self
    typedef  Traits_                    Traits;
    typedef  Min_annulus_d<Traits>      Self;

    // types from the traits class
    typedef  typename Traits::Point_d   Point;

    typedef  typename Traits::Rep_tag   Rep_tag;

    typedef  typename Traits::RT        RT;
    typedef  typename Traits::FT        FT;

    typedef  typename Traits::Access_dimension_d
                                        Access_dimension_d;
    typedef  typename Traits::Access_coordinates_begin_d
                                        Access_coordinates_begin_d;

    typedef  typename Traits::Construct_point_d
                                        Construct_point_d;

    typedef  typename Traits::ET        ET;
    typedef  typename Traits::NT        NT;

  private:
    // LP solver
    typedef  CGAL::LP_rep_min_annulus_d<
                 ET, NT, Point, typename std::vector<Point>::const_iterator,
                 Access_coordinates_begin_d, Access_dimension_d >
                                        LP_rep;
    typedef  CGAL::QP_solver< LP_rep >  Solver;
    typedef  typename Solver::Pricing_strategy
                                        Pricing_strategy;
    // types from the QP solver
    typedef  typename Solver::Basic_variable_index_iterator
                                        Basic_variable_index_iterator;
    
    // private types
    typedef  std::vector<Point>         Point_vector;
    typedef  std::vector<ET>            ET_vector;
    
    typedef  CGAL::Access_by_index<typename std::vector<Point>::const_iterator>
                                        Point_by_index;
    
    typedef  std::binder2nd< std::divides<int> >
                                        Divide;
    
    typedef  std::vector<int>           Index_vector;
    
    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;
    

  public:
    // public types
    typedef  typename Point_vector::const_iterator
                                        Point_iterator;
    
    typedef  CGAL::Join_random_access_iterator_1<
                 Basic_variable_index_iterator,
                 CGAL::Unary_compose_1<Point_by_index,Divide> >
                                        Support_point_iterator;
    
    typedef  CGAL::Join_random_access_iterator_1<
                 typename Index_vector::const_iterator, Point_by_index >
                                        Inner_support_point_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 typename Index_vector::const_iterator, Point_by_index >
                                        Outer_support_point_iterator;
    
    typedef  typename ET_vector::const_iterator
                                        Coordinate_iterator;
    

    // creation
    Min_annulus_d( const Traits&  traits  = Traits(),
                   int            verbose = 0,
                   std::ostream&  stream  = std::cout)
      : tco( traits), d( -1), solver( verbose, stream)
        {
            set_pricing_strategy( NT());
        }
    
    template < class InputIterator >
    Min_annulus_d( InputIterator  first,
                   InputIterator  last,
                   const Traits&  traits = Traits(),
                   int            verbose = 0,
                   std::ostream&  stream  = std::cout)
      : tco( traits), solver( verbose, stream)
        {
            set_pricing_strategy( NT());
            set( first, last);
        }
    
    // access to point set
    int  ambient_dimension( ) const { return d; }
    
    int  number_of_points( ) const { return points.size(); }
    
    Point_iterator  points_begin( ) const { return points.begin(); }
    Point_iterator  points_end  ( ) const { return points.end  (); }
    
    // access to support points
    int
    number_of_support_points( ) const
        { return number_of_points() < 2 ? number_of_points()
                                        : solver.number_of_basic_variables(); }
    
    Support_point_iterator
    support_points_begin() const
        { return Support_point_iterator(
                     solver.basic_variables_index_begin(),
                     CGAL::compose1_1(
                         Point_by_index( points.begin()),
                         std::bind2nd( std::divides<int>(), 2)));}
    
    Support_point_iterator
    support_points_end() const
        { return Support_point_iterator( number_of_points() < 2
                     ? solver.basic_variables_index_begin()
                     : solver.basic_variables_index_end(),
                     CGAL::compose1_1(
                         Point_by_index( points.begin()),
                         std::bind2nd( std::divides<int>(), 2)));}
    
    int  number_of_inner_support_points() const { return inner_indices.size();}
    int  number_of_outer_support_points() const { return outer_indices.size();}
    
    Inner_support_point_iterator
    inner_support_points_begin() const
        { return Inner_support_point_iterator(
                     inner_indices.begin(),
                     Point_by_index( points.begin())); }
    
    Inner_support_point_iterator
    inner_support_points_end() const
        { return Inner_support_point_iterator(
                     inner_indices.end(),
                     Point_by_index( points.begin())); }
    
    Outer_support_point_iterator
    outer_support_points_begin() const
        { return Outer_support_point_iterator(
                     outer_indices.begin(),
                     Point_by_index( points.begin())); }
    
    Outer_support_point_iterator
    outer_support_points_end() const
        { return Outer_support_point_iterator(
                     outer_indices.end(),
                     Point_by_index( points.begin())); }
    
    // access to center (rational representation)
    Coordinate_iterator
    center_coordinates_begin( ) const { return center_coords.begin(); }
    
    Coordinate_iterator
    center_coordinates_end  ( ) const { return center_coords.end  (); }
    
    // access to squared radii (rational representation)
    ET  squared_inner_radius_numerator( ) const { return sqr_i_rad_numer; }
    ET  squared_outer_radius_numerator( ) const { return sqr_o_rad_numer; }
    ET  squared_radii_denominator     ( ) const { return sqr_rad_denom; }
    
    // access to center and squared radii
    // NOTE: an implicit conversion from ET to RT must be available!
    Point
    center( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return tco.construct_point_d_object()( ambient_dimension(),
                                                 center_coordinates_begin(),
                                                 center_coordinates_end()); }
    
    FT
    squared_inner_radius( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_inner_radius_numerator()) /
                 FT( squared_radii_denominator()); }
    
    FT
    squared_outer_radius( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_outer_radius_numerator()) /
                 FT( squared_radii_denominator()); }
    
    // predicates
    CGAL::Bounded_side
    bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return CGAL::Bounded_side(
                       CGAL::NTS::sign( sqr_d - sqr_i_rad_numer)
                     * CGAL::NTS::sign( sqr_o_rad_numer - sqr_d)); }
    
    bool
    has_on_bounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
            is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return ( ( sqr_i_rad_numer < sqr_d) && ( sqr_d < sqr_o_rad_numer)); }
    
    bool
    has_on_boundary( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return (( sqr_d == sqr_i_rad_numer) || ( sqr_d == sqr_o_rad_numer));}
    
    bool
    has_on_unbounded_side( const Point& p) const
        { CGAL_optimisation_precondition(
              is_empty() || tco.access_dimension_d_object()( p) == d);
          ET sqr_d = sqr_dist( p);
          return ( ( sqr_d < sqr_i_rad_numer) || ( sqr_o_rad_numer < sqr_d)); }
    
    bool  is_empty     ( ) const { return number_of_points() == 0; }
    bool  is_degenerate( ) const
        { return ! CGAL_NTS is_positive( sqr_o_rad_numer); }
    
    // modifiers
    template < class InputIterator >
    void
    set( InputIterator first, InputIterator last)
        { if ( points.size() > 0) points.erase( points.begin(), points.end());
          std::copy( first, last, std::back_inserter( points));
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension(),
              "Not all points have the same dimension.");
          compute_min_annulus(); }
    
    void
    insert( const Point& p)
        { if ( is_empty()) d = tco.access_dimension_d_object()( p);
          CGAL_optimisation_precondition(
              tco.access_dimension_d_object()( p) == d);
          points.push_back( p);
          compute_min_annulus(); }
    
    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
        { CGAL_optimisation_precondition_code( int old_n = points.size());
          points.insert( points.end(), first, last);
          set_dimension();
          CGAL_optimisation_precondition_msg( check_dimension( old_n),
              "Not all points have the same dimension.");
          compute_min_annulus(); }
    
    void
    clear( )
        { points.erase( points.begin(), points.end());
          compute_min_annulus(); }
    
    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
    
    // traits class access
    const Traits&  traits( ) const { return tco; }
    

  private:
    
    Traits                   tco;       // traits class object
    
    Point_vector             points;    // input points
    int                      d;         // dimension of input points
    
    ET_vector                center_coords;     // center of small.encl.annulus
    
    ET                       sqr_i_rad_numer;   // squared inner radius of
    ET                       sqr_o_rad_numer;   // ---"--- outer ----"----
    ET                       sqr_rad_denom;     // smallest enclosing annulus
    
    Solver                   solver;    // linear programming solver
    
    Index_vector             inner_indices;
    Index_vector             outer_indices;
    
    NT_matrix                a_matrix;  // matrix `A' of dual LP
    NT_vector                b_vector;  // vector `b' of dual LP
    NT_vector                c_vector;  // vector `c' of dual LP
    
    typename Solver::Pricing_strategy*  // pricing strategy
                             strategyP; // of the QP solver
    

    
    // squared distance to center
    ET
    sqr_dist( const Point& p) const
        { return std::inner_product(
              center_coords.begin(), center_coords.end()-1,
              tco.access_coordinates_begin_d_object()( p), ET( 0),
              std::plus<ET>(),
              CGAL::compose1_2(
                  CGAL::compose2_1( std::multiplies<ET>(),
                      CGAL::identity<ET>(), CGAL::identity<ET>()),
                  CGAL::compose2_2( std::minus<ET>(),
                      CGAL::identity<ET>(),
                      std::bind2nd( std::multiplies<ET>(),
                                    center_coords.back())))); }
    
    // set dimension of input points
    void
    set_dimension( )
        { d = ( points.size() == 0 ? -1 :
                    tco.access_dimension_d_object()( points[ 0])); }
    
    // check dimension of input points
    bool
    check_dimension( unsigned int  offset = 0)
        { return ( std::find_if( points.begin()+offset, points.end(),
                                 CGAL::compose1_1( std::bind2nd(
                                     std::not_equal_to<int>(), d),
                                     tco.access_dimension_d_object()))
                   == points.end()); }
    
    // compute smallest enclosing annulus
    void
    compute_min_annulus( )
    {
        // clear inner and outer support points
        inner_indices.erase( inner_indices.begin(), inner_indices.end());
        outer_indices.erase( outer_indices.begin(), outer_indices.end());
    
        if ( is_empty()) {
            center_coords.resize( 1);
            sqr_i_rad_numer = -ET( 1);
            sqr_o_rad_numer = -ET( 1);
            return;
        }
    
        if ( number_of_points() == 1) {
            inner_indices.push_back( 0);
            outer_indices.push_back( 0);
            center_coords.resize( d+1);
            std::copy( tco.access_coordinates_begin_d_object()( points[ 0]),
                       tco.access_coordinates_begin_d_object()( points[ 0])+d,
                       center_coords.begin());
            center_coords[ d] = ET( 1);
            sqr_i_rad_numer = ET( 0);
            sqr_o_rad_numer = ET( 0);
            sqr_rad_denom   = ET( 1);
            return;
        }
    
        // set up and solve dual LP
        int i, j;
        NT  nt_0 = 0, nt_1 = 1, nt_2 = 2;
        NT  nt_minus_1 = -nt_1, nt_minus_2 = -nt_2;
        
        // vector b
        b_vector.resize( d+2);
        for ( j = 0; j < d; ++j) b_vector[ j] = nt_0;
        b_vector[ d  ] = nt_1;
        b_vector[ d+1] = nt_minus_1;
        
        // matrix A, vector c
        a_matrix.erase( a_matrix.begin(), a_matrix.end());
        a_matrix.insert( a_matrix.end(), 2*points.size(), NT_vector( d+2));
        c_vector.resize( 2*points.size());
        for ( i = 0; i < number_of_points(); ++i) {
            typename Traits::Access_coordinates_begin_d::Coordinate_iterator
                coord_it = tco.access_coordinates_begin_d_object()( points[i]);
            NT  sum = 0;
            for ( j = 0; j < d; ++j) {
                a_matrix[ 2*i  ][ j] = nt_2*coord_it[ j];
                a_matrix[ 2*i+1][ j] = nt_minus_2*coord_it[ j];
                sum += NT( coord_it[ j])*NT( coord_it[ j]);
            }
            a_matrix[ 2*i  ][ d  ] = nt_1;
            a_matrix[ 2*i+1][ d  ] = nt_0;
            a_matrix[ 2*i  ][ d+1] = nt_0;
            a_matrix[ 2*i+1][ d+1] = nt_minus_1;
            c_vector[ 2*i  ] =  sum;
            c_vector[ 2*i+1] = -sum;
        }
        typedef  typename LP_rep::A_iterator  A_it;
        typedef  typename LP_rep::D_iterator  D_it;
        solver.set( 2*points.size(), d+2, d+2,
                    A_it( a_matrix.begin()), b_vector.begin(),
                    c_vector.begin(), D_it());
        solver.init();
        solver.solve();
    
        // compute center and squared radius
        ET sqr_sum = 0;
        center_coords.resize( ambient_dimension()+1);
        for ( i = 0; i < d; ++i) {
            center_coords[ i] = -solver.dual_variable( i);
            sqr_sum += center_coords[ i] * center_coords[ i];
        }
        center_coords[ d] = solver.variables_common_denominator();
        sqr_i_rad_numer = sqr_sum
                          - solver.dual_variable( d  )*center_coords[ d];
        sqr_o_rad_numer = sqr_sum
                          - solver.dual_variable( d+1)*center_coords[ d];
        sqr_rad_denom   = center_coords[ d] * center_coords[ d];
        
        // split up support points
        for ( i = 0; i < solver.number_of_basic_variables(); ++i) {
            int index = solver.basic_variables_index_begin()[ i];
            if ( index % 2 == 0) {
                inner_indices.push_back( index/2);
            } else {
                outer_indices.push_back( index/2);
            }
        }
    }
    
    template < class NT >
    void  set_pricing_strategy( NT)
        { strategyP = new CGAL::Partial_filtered_pricing<LP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    
    #ifndef _MSC_VER
    void  set_pricing_strategy( ET)
        { strategyP = new CGAL::Partial_exact_pricing<LP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    #endif
    
};

template < class ET_, class NT_, class Point, class PointIterator,
           class Access_coord, class Access_dim >
struct LP_rep_min_annulus_d {
    typedef  ET_                    ET;
    typedef  NT_                    NT;

    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;
    
    typedef  CGAL::Join_random_access_iterator_1<
                 CGAL_TYPENAME_MSVC_NULL NT_matrix::const_iterator,
                 LP_rep_row_of_a<NT> >  A_iterator;
    typedef  typename NT_vector::const_iterator
                                        B_iterator;
    typedef  typename NT_vector::const_iterator
                                        C_iterator;
    
    typedef  A_iterator                 D_iterator;     // dummy
    

    typedef  CGAL::Tag_true         Is_lp;
};

// Function declarations
// =====================
// I/O operators
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os, const Min_annulus_d<Traits_>& min_annulus);

template < class Traits_ >
std::istream&
operator >> ( std::istream& is,       Min_annulus_d<Traits_>& min_annulus);

// ============================================================================

// Class implementation
// ====================

// validity check
template < class Traits_ >
bool
Min_annulus_d<Traits_>::
is_valid( bool verbose, int level) const
{
    CGAL_USING_NAMESPACE_STD

    CGAL::Verbose_ostream verr( verbose);
    verr << "CGAL::Min_annulus_d<Traits>::" << endl;
    verr << "is_valid( true, " << level << "):" << endl;
    verr << "  |P| = " << number_of_points()
         << ", |S| = " << number_of_support_points() << endl;

    // containment check (a)
    // ---------------------
    verr << "  (a) containment check..." << flush;
    
    Point_iterator  point_it = points_begin();
    for ( ; point_it != points_end(); ++point_it) {
        if ( has_on_unbounded_side( *point_it))
            return CGAL::_optimisation_is_valid_fail( verr,
                       "annulus does not contain all points");
    }
    
    verr << "passed." << endl;

    // support set check (b)
    // ---------------------
    verr << "  (b) support set check..." << flush;
    
    // all inner support points on inner boundary?
    Inner_support_point_iterator  i_pt_it = inner_support_points_begin();
    for ( ; i_pt_it != inner_support_points_end(); ++i_pt_it) {
        if ( sqr_dist( *i_pt_it) != sqr_i_rad_numer)
            return CGAL::_optimisation_is_valid_fail( verr,
       "annulus does not have all inner support points on its inner boundary");
    }
    
    // all outer support points on outer boundary?
    Outer_support_point_iterator  o_pt_it = outer_support_points_begin();
    for ( ; o_pt_it != outer_support_points_end(); ++o_pt_it) {
        if ( sqr_dist( *o_pt_it) != sqr_o_rad_numer)
            return CGAL::_optimisation_is_valid_fail( verr,
       "annulus does not have all outer support points on its outer boundary");
    }
    /*
    // center strictly in convex hull of support points?
    typename Solver::Basic_variable_numerator_iterator
        num_it = solver.basic_variables_numerator_begin();
    for ( ; num_it != solver.basic_variables_numerator_end(); ++num_it) {
        if ( ! (    CGAL::NTS::is_positive( *num_it)
                 && *num_it <= solver.variables_common_denominator()))
            return CGAL::_optimisation_is_valid_fail( verr,
              "center does not lie strictly in convex hull of support points");
    }
    */
    
    verr << "passed." << endl;

    verr << "  object is valid!" << endl;
    return( true);
}

// output operator
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Min_annulus_d<Traits_>& min_annulus)
{
    CGAL_USING_NAMESPACE_STD

    typedef  Min_annulus_d<Traits_>::Point  Point;
    typedef  ostream_iterator<Point>       Os_it;
    typedef  typename Traits_::ET          ET;
    typedef  ostream_iterator<ET>          Et_it;

    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << "CGAL::Min_annulus_d( |P| = "
           << min_annulus.number_of_points() << ", |S| = "
           << min_annulus.number_of_inner_support_points() << '+'
           << min_annulus.number_of_outer_support_points() << endl;
        os << "  P = {" << endl;
        os << "    ";
        copy( min_annulus.points_begin(), min_annulus.points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S_i = {" << endl;
        os << "    ";
        copy( min_annulus.inner_support_points_begin(),
              min_annulus.inner_support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S_o = {" << endl;
        os << "    ";
        copy( min_annulus.outer_support_points_begin(),
              min_annulus.outer_support_points_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  center = ( ";
        copy( min_annulus.center_coordinates_begin(),
              min_annulus.center_coordinates_end(),
              Et_it( os, " "));
        os << ")" << endl;
        os << "  squared inner radius = "
           << min_annulus.squared_inner_radius_numerator() << " / "
           << min_annulus.squared_radii_denominator() << endl;
        os << "  squared outer radius = "
           << min_annulus.squared_outer_radius_numerator() << " / "
           << min_annulus.squared_radii_denominator() << endl;
        break;

      case CGAL::IO::ASCII:
        copy( min_annulus.points_begin(), min_annulus.points_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        copy( min_annulus.points_begin(), min_annulus.points_end(),
              Os_it( os));
        break;

      default:
        CGAL_optimisation_assertion_msg( false,
                                         "CGAL::get_mode( os) invalid!");
        break; }

    return( os);
}

// input operator
template < class Traits_ >
std::istream&
operator >> ( std::istream& is, CGAL::Min_annulus_d<Traits_>& min_annulus)
{
    CGAL_USING_NAMESPACE_STD

    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        cerr << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef  CGAL::Min_annulus_d<Traits_>::Point  Point;
        typedef  istream_iterator<Point>             Is_it;
        min_annulus.set( Is_it( is), Is_it());
        break;

      default:
        CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }

    return( is);
}

CGAL_END_NAMESPACE

#endif // CGAL_MIN_ANNULUS_D_H

// ===== EOF ==================================================================
