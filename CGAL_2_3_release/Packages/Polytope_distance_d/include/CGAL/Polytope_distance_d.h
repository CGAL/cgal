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
// file          : include/CGAL/Polytope_distance_d.h
// package       : $CGAL_Package: Polytope_distance_d $
// chapter       : Geometric Optimisation
//
// source        : web/Polytope_distance_d.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Distance of convex polytopes in arbitrary dimension
// ============================================================================

#ifndef CGAL_POLYTOPE_DISTANCE_D_H
#define CGAL_POLYTOPE_DISTANCE_D_H

// includes
// --------
#ifndef CGAL_OPTIMISATION_BASIC_H
#  include <CGAL/Optimisation/basic.h>
#endif

#ifndef CGAL_FUNCTION_OBJECTS_H
#  include <CGAL/function_objects.h>
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
#ifndef CGAL_CONST_VALUE_ITERATOR_H
#  include <CGAL/_QP_solver/Const_value_iterator.h>
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


CGAL_BEGIN_NAMESPACE

// Class declarations
// ==================
template < class Traits_ >
class Polytope_distance_d;

template < class ET_, class NT_, class Point, class Point_iterator,
           class Access_coord, class Access_dim >
struct QP_rep_poly_dist_d;

template < class NT >
struct QP_rep_row_of_a {
    typedef  std::vector<NT>                         argument_type;
    typedef  typename argument_type::const_iterator  result_type;

    result_type
    operator ( ) ( const argument_type& v) const { return v.begin(); }
};

template < class Point, class Point_iterator >
struct QP_rep_signed_point_iterator;

template < class NT, class Point,
           class Access_coord, class Access_dim >
class QP_rep_signed_inner_product;

template < class NT, class Point, class Point_iterator,
           class Access_coord, class Access_dim >
struct QP_rep_row_of_d;


// Class interfaces
// ================
template < class Traits_ >
class Polytope_distance_d {
  public:
    // self
    typedef  Traits_                    Traits;
    typedef  Polytope_distance_d<Traits>
                                        Self;

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
    // QP solver
    typedef  CGAL::QP_rep_poly_dist_d<
                 ET, NT, Point, typename std::vector<Point>::const_iterator,
                 Access_coordinates_begin_d, Access_dimension_d >
                                        QP_rep;
    typedef  CGAL::QP_solver< QP_rep >  Solver;
    typedef  typename Solver::Pricing_strategy
                                        Pricing_strategy;
    // private types
    typedef  std::vector<Point>         Point_vector;
    typedef  std::vector<ET>            ET_vector;
    
    typedef  std::vector<int>           Index_vector;
    
    typedef  CGAL::Access_by_index<typename std::vector<Point>::const_iterator>
                                        Point_by_index;
    
    typedef  std::vector<NT>            NT_vector;
    typedef  std::vector<NT_vector>     NT_matrix;
    

  public:
    // public types
    typedef  typename Point_vector::const_iterator
                                        Point_iterator;
    
    typedef  CGAL::Join_random_access_iterator_1<
                 typename Index_vector::const_iterator, Point_by_index >
                                        Support_point_iterator;
    
    typedef  typename ET_vector::const_iterator
                                        Coordinate_iterator;
    

    // creation
    Polytope_distance_d( const Traits&  traits  = Traits(),
                         int            verbose = 0,
                         std::ostream&  stream  = std::cout)
      : tco( traits), d( -1), solver( verbose, stream)
        {
            set_pricing_strategy( NT());
        }
    
    template < class InputIterator1, class InputIterator2 >
    Polytope_distance_d( InputIterator1 p_first,
                         InputIterator1 p_last,
                         InputIterator2 q_first,
                         InputIterator2 q_last,
                         const Traits&  traits = Traits(),
                         int            verbose = 0,
                         std::ostream&  stream  = std::cout)
      : tco( traits), solver( verbose, stream)
        {
            set_pricing_strategy( NT());
            set( p_first, p_last, q_first, q_last);
        }
    
    // access to point sets
    int  ambient_dimension( ) const { return d; }
    
    int  number_of_points( ) const { return p_points.size()+q_points.size();}
    
    int  number_of_points_p( ) const { return p_points.size(); }
    int  number_of_points_q( ) const { return q_points.size(); }
    
    Point_iterator  points_p_begin( ) const { return p_points.begin(); }
    Point_iterator  points_p_end  ( ) const { return p_points.end  (); }
    
    Point_iterator  points_q_begin( ) const { return q_points.begin(); }
    Point_iterator  points_q_end  ( ) const { return q_points.end  (); }
    
    // access to support points
    int
    number_of_support_points( ) const
        { return is_finite() ? solver.number_of_basic_variables() : 0; }
    
    int  number_of_support_points_p() const { return p_support_indices.size();}
    int  number_of_support_points_q() const { return q_support_indices.size();}
    
    Support_point_iterator
    support_points_p_begin() const
        { return Support_point_iterator(
                     p_support_indices.begin(),
                     Point_by_index( p_points.begin())); }
    
    Support_point_iterator
    support_points_p_end() const
        { return Support_point_iterator(
                     is_finite() ? p_support_indices.end()
                                 : p_support_indices.begin(),
                     Point_by_index( p_points.begin())); }
    
    Support_point_iterator
    support_points_q_begin() const
        { return Support_point_iterator(
                     q_support_indices.begin(),
                     Point_by_index( q_points.begin())); }
    
    Support_point_iterator
    support_points_q_end() const
        { return Support_point_iterator(
                     is_finite() ? q_support_indices.end()
                                 : q_support_indices.begin(),
                     Point_by_index( q_points.begin())); }
    
    // access to realizing points (rational representation)
    Coordinate_iterator
    realizing_point_p_coordinates_begin( ) const { return p_coords.begin(); }
    
    Coordinate_iterator
    realizing_point_p_coordinates_end  ( ) const { return p_coords.end  (); }
    
    Coordinate_iterator
    realizing_point_q_coordinates_begin( ) const { return q_coords.begin(); }
    
    Coordinate_iterator
    realizing_point_q_coordinates_end  ( ) const { return q_coords.end  (); }
    
    // access to squared distance (rational representation)
    ET  squared_distance_numerator  ( ) const
        { return solver.solution_numerator(); }
    
    ET  squared_distance_denominator( ) const
        { return solver.solution_denominator(); }
    
    // access to realizing points and squared distance
    // NOTE: an implicit conversion from ET to RT must be available!
    Point
    realizing_point_p( ) const
        { CGAL_optimisation_precondition( is_finite());
          return tco.construct_point_d_object()( ambient_dimension(),
                     realizing_point_p_coordinates_begin(),
                     realizing_point_p_coordinates_end  ()); }
    
    Point
    realizing_point_q( ) const
        { CGAL_optimisation_precondition( is_finite());
          return tco.construct_point_d_object()( ambient_dimension(),
                     realizing_point_q_coordinates_begin(),
                     realizing_point_q_coordinates_end  ()); }
    
    FT
    squared_distance( ) const
        { CGAL_optimisation_precondition( ! is_empty());
          return FT( squared_distance_numerator  ()) /
                 FT( squared_distance_denominator()); }
    
    bool  is_finite( ) const
        { return ( number_of_points_p() > 0) && ( number_of_points_q() > 0); }
    
    bool  is_zero( ) const
        { return CGAL_NTS is_zero( squared_distance_numerator()); }
    
    bool  is_degenerate( ) const { return ( ! is_finite()); }
    
    // modifiers
    template < class InputIterator1, class InputIterator2 >
    void
    set( InputIterator1 p_first, InputIterator1 p_last,
         InputIterator2 q_first, InputIterator2 q_last)
        { if ( p_points.size() > 0)
              p_points.erase( p_points.begin(), p_points.end());
          if ( q_points.size() > 0)
              q_points.erase( q_points.begin(), q_points.end());
          std::copy( p_first, p_last, std::back_inserter( p_points));
          std::copy( q_first, q_last, std::back_inserter( q_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
                 check_dimension( p_points.begin(), p_points.end())
              && check_dimension( q_points.begin(), q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    template < class InputIterator >
    void
    set_p( InputIterator p_first, InputIterator p_last)
        { if ( p_points.size() > 0)
              p_points.erase( p_points.begin(), p_points.end());
          std::copy( p_first, p_last, std::back_inserter( p_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( p_points.begin(), p_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    template < class InputIterator >
    void
    set_q( InputIterator q_first, InputIterator q_last)
        { if ( q_points.size() > 0)
              q_points.erase( q_points.begin(), q_points.end());
          std::copy( q_first, q_last, std::back_inserter( q_points));
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( q_points.begin(), q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    void
    insert_p( const Point& p)
        { CGAL_optimisation_precondition( ( ! is_finite()) ||
              ( tco.access_dimension_d_object()( p) == d));
          p_points.push_back( p);
          compute_distance(); }
    
    void
    insert_q( const Point& q)
        { CGAL_optimisation_precondition( ( ! is_finite()) ||
              ( tco.access_dimension_d_object()( q) == d));
          q_points.push_back( q);
          compute_distance(); }
    
    template < class InputIterator1, class InputIterator2 >
    void
    insert( InputIterator1 p_first, InputIterator1 p_last,
            InputIterator2 q_first, InputIterator2 q_last)
        { CGAL_optimisation_precondition_code( int old_r = p_points.size());
          CGAL_optimisation_precondition_code( int old_s = q_points.size());
          p_points.insert( p_points.end(), p_first, p_last);
          q_points.insert( q_points.end(), q_first, q_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
                 check_dimension( p_points.begin()+old_r, p_points.end())
              && check_dimension( q_points.begin()+old_s, q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    template < class InputIterator >
    void
    insert_p( InputIterator p_first, InputIterator p_last)
        { CGAL_optimisation_precondition_code( int old_r = p_points.size());
          p_points.insert( p_points.end(), p_first, p_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( p_points.begin()+old_r, p_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    template < class InputIterator >
    void
    insert_q( InputIterator q_first, InputIterator q_last)
        { CGAL_optimisation_precondition_code( int old_s = q_points.size());
          q_points.insert( q_points.end(), q_first, q_last);
          set_dimension();
          CGAL_optimisation_precondition_msg(
              check_dimension( q_points.begin()+old_s, q_points.end()),
              "Not all points have the same dimension.");
          compute_distance(); }
    
    void
    clear( )
        { p_points.erase( p_points.begin(), p_points.end());
          q_points.erase( q_points.begin(), q_points.end());
          compute_distance(); }
    
    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;
    
    // traits class access
    const Traits&  traits( ) const { return tco; }
    

  private:
    
    Traits                   tco;       // traits class object
    
    Point_vector             p_points;  // points of P
    Point_vector             q_points;  // points of Q
    int                      d;         // dimension of input points
    
    ET_vector                p_coords;          // realizing point of P
    ET_vector                q_coords;          // realizing point of Q
    
    Solver                   solver;    // quadratic programming solver
    
    Index_vector             p_support_indices;
    Index_vector             q_support_indices;
    
    NT_matrix                a_matrix;  // matrix `A' of QP
    
    typename Solver::Pricing_strategy*  // pricing strategy
                             strategyP; // of the QP solver
    

    
    // set dimension of input points
    void
    set_dimension( )
        { d = ( p_points.size() > 0 ?
                    tco.access_dimension_d_object()( p_points[ 0]) :
                q_points.size() > 0 ?
                    tco.access_dimension_d_object()( q_points[ 0]) :
                -1); }
    
    // check dimension of input points
    template < class InputIterator >
    bool
    check_dimension( InputIterator first, InputIterator last)
        { return ( std::find_if( first, last,
                                 CGAL::compose1_1( std::bind2nd(
                                     std::not_equal_to<int>(), d),
                                     tco.access_dimension_d_object()))
                   == last); }
    
    // compute (squared) distance
    void
    compute_distance( )
    {
        // clear support points
        p_support_indices.erase( p_support_indices.begin(),
                                 p_support_indices.end());
        q_support_indices.erase( q_support_indices.begin(),
                                 q_support_indices.end());
        if ( ( p_points.size() == 0) || ( q_points.size() == 0)) return;
    
        // set up and solve QP
        int i, j;
        NT  nt_0 = 0, nt_1 = 1;
        
        // matrix A
        a_matrix.erase( a_matrix.begin(), a_matrix.end());
        a_matrix.insert( a_matrix.end(),
                         number_of_points(), NT_vector( 2, nt_0));
        for ( j = 0; j < number_of_points_p(); ++j) a_matrix[ j][ 0] = nt_1;
        for (      ; j < number_of_points  (); ++j) a_matrix[ j][ 1] = nt_1;
        
        // set-up
        typedef  QP_rep_signed_point_iterator< Point, Point_iterator >
                                            Signed_point_iterator;
        Signed_point_iterator
            signed_pts_it(p_points.begin(), p_points.size(), q_points.begin());
        QP_rep_row_of_d< NT, Point, Signed_point_iterator,
                         Access_coordinates_begin_d, Access_dimension_d >
            row_of_d( signed_pts_it,
                      tco.access_coordinates_begin_d_object(),
                      tco.access_dimension_d_object());
        
        typedef  typename QP_rep::A_iterator A_it;
        typedef  typename QP_rep::B_iterator B_it;
        typedef  typename QP_rep::C_iterator C_it;
        typedef  typename QP_rep::D_iterator D_it;
        
        solver.set( number_of_points(), 2, d+2,
                    A_it( a_matrix.begin()), B_it( 1), C_it( 0),
                    D_it( signed_pts_it, row_of_d));
        
        // solve
        solver.init();
        solver.solve();
    
        // compute support and realizing points
        ET  et_0 = 0;
        int r    = number_of_points_p();
        p_coords.resize( ambient_dimension()+1);
        q_coords.resize( ambient_dimension()+1);
        std::fill( p_coords.begin(), p_coords.end(), et_0);
        std::fill( q_coords.begin(), q_coords.end(), et_0);
        for ( i = 0; i < solver.number_of_basic_variables(); ++i) {
            ET  value = solver.basic_variables_numerator_begin()[ i];
            int index = solver.basic_variables_index_begin()[ i];
            if ( index < r) {
                for ( int j = 0; j < d; ++j) {
                    p_coords[ j]
                        += value * tco.access_coordinates_begin_d_object()(
                                                       p_points[ index  ])[ j];
                }
                p_support_indices.push_back( index);
            } else {
                for ( int j = 0; j < d; ++j) {
                    q_coords[ j]
                        += value * tco.access_coordinates_begin_d_object()(
                                                       q_points[ index-r])[ j];
                }
                q_support_indices.push_back( index-r);
            }
        }
        p_coords[ d] = q_coords[ d] = solver.variables_common_denominator();
    }
    
    template < class NT >
    void  set_pricing_strategy( NT)
        { strategyP = new CGAL::Partial_filtered_pricing<QP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    
    #ifndef _MSC_VER
    void  set_pricing_strategy( ET)
        { strategyP = new CGAL::Partial_exact_pricing<QP_rep>;
          solver.set_pricing_strategy( *strategyP); }
    #endif
    
};

template < class Point, class PointIterator >
struct QP_rep_signed_point_iterator {
  public:
    typedef  std::pair<Point,CGAL::Sign>      value_type;
    typedef  ptrdiff_t                        difference_type;
    typedef  value_type*                      pointer;
    typedef  value_type&                      reference;
    typedef  std::random_access_iterator_tag  iterator_category;

    typedef  QP_rep_signed_point_iterator<Point,PointIterator>  Self;
    typedef  value_type                                         Val;
    typedef  difference_type                                    Dist;
    typedef  pointer                                            Ptr;

    // forward operations
    QP_rep_signed_point_iterator(
        const PointIterator& it_p = PointIterator(), Dist n_p = 0,
        const PointIterator& it_q = PointIterator())
        : p_it( it_p), q_it( it_q), n( n_p), curr( 0) { }

    bool   operator == ( const Self& it) const { return (curr == it.curr);}
    bool   operator != ( const Self& it) const { return (curr != it.curr);}

    Val    operator *  ( ) const
        { return ( curr < n) ? std::make_pair( *p_it, CGAL::POSITIVE)
                             : std__make_pair( *q_it, CGAL::NEGATIVE); }

    Self&  operator ++ (    )
               { if ( ++curr <= n) ++p_it; else ++q_it; return *this; }
    Self   operator ++ ( int)
               { Self tmp = *this; operator++(); return tmp; }

    // bidirectional operations
    Self&  operator -- (    )
               { if ( --curr <  n) --p_it; else --q_it; return *this; }
    Self   operator -- ( int)
               { Self tmp = *this; operator--(); return tmp; }

    // random access operations
    Self&  operator += ( Dist i)
               {
                   if ( curr+i <= n) {
                       curr += i;
                       p_it += i;
                   } else {
                       if ( curr < n) p_it += n-curr;
                       curr += i;
                       q_it += curr-n;
                   }
                   return *this;
               }

    Self&  operator -= ( Dist i)
               {
                   if ( curr-i < n) {
                       if ( curr > n) q_it -= curr-n;
                       curr -= i;
                       p_it  -= n-curr;
                   } else {
                       curr -= i;
                       q_it  -= i;
                   }
                   return *this;
               }

    Self   operator +  ( Dist i) const { Self tmp = *this; return tmp+=i; }
    Self   operator -  ( Dist i) const { Self tmp = *this; return tmp-=i; }

    Dist   operator -  ( const Self& it) const { return curr - it.curr; }

    Val    operator [] ( int i) const
        { return ( curr+i < n)
              ? std::make_pair( p_it[ i  ], CGAL::POSITIVE)
              : std::make_pair( q_it[ i-n], CGAL::NEGATIVE); }

    bool   operator <  ( const Self&) const { return ( curr <  it.curr); }
    bool   operator >  ( const Self&) const { return ( curr >  it.curr); }
    bool   operator <= ( const Self&) const { return ( curr <= it.curr); }
    bool   operator >= ( const Self&) const { return ( curr >= it.curr); }

  private:
    PointIterator  p_it;
    PointIterator  q_it;
    Dist           n;
    Dist           curr;
};

template < class NT, class Point,
           class Access_coord, class Access_dim >
class QP_rep_signed_inner_product {
    Point         p_i;
    CGAL::Sign    s_i;
    Access_coord  da_coord;
    Access_dim    da_dim;
  public:
    typedef  std::pair<Point,CGAL::Sign>  argument_type;
    typedef  NT                           result_type;

    QP_rep_signed_inner_product( ) { }
    QP_rep_signed_inner_product( const argument_type& p_signed,
                                 const Access_coord&  ac,
                                 const Access_dim&    ad)
        : p_i( p_signed.first), s_i( p_signed.second),
          da_coord( ac), da_dim( ad) { }

    NT  operator( ) ( const argument_type& p_signed) const
        { NT ip = std::inner_product( da_coord( p_i),
                                      da_coord( p_i)+da_dim( p_i),
                                      da_coord( p_signed.first), NT( 0),
                                      std::plus<NT>(),
                                      std::multiplies<NT>());
          return ( s_i*p_signed.second == CGAL::POSITIVE) ? ip : -ip; }
};

template < class NT, class Point, class Signed_point_iterator,
           class Access_coord, class Access_dim >
class QP_rep_row_of_d {
    Signed_point_iterator  signed_pts_it;
    Access_coord           da_coord;
    Access_dim             da_dim;
public:
    typedef  CGAL::QP_rep_signed_inner_product<
                 NT, Point, Access_coord, Access_dim >
                                    Signed_inner_product;
    typedef  CGAL::Join_random_access_iterator_1<
                 Signed_point_iterator, Signed_inner_product >
                                    Row_of_d;

    typedef  std::pair<Point,CGAL::Sign>
                                    argument_type;
    typedef  Row_of_d               result_type;

    QP_rep_row_of_d( ) { }
    QP_rep_row_of_d( const Signed_point_iterator& it,
                     const Access_coord&          ac,
                     const Access_dim&            ad)
        : signed_pts_it( it), da_coord( ac), da_dim( ad) { }

    Row_of_d  operator( ) ( const argument_type& p_signed) const
    { return Row_of_d( signed_pts_it,
                       Signed_inner_product( p_signed, da_coord, da_dim));}
};

template < class ET_, class NT_, class Point, class Point_iterator,
           class Access_coord, class Access_dim >
struct QP_rep_poly_dist_d {
    typedef  ET_                    ET;
    typedef  NT_                    NT;

    typedef  std::vector< std::vector<NT> >
                                        NT_matrix;
    
    typedef  CGAL::Join_random_access_iterator_1<
                 CGAL_TYPENAME_MSVC_NULL NT_matrix::const_iterator,
                 QP_rep_row_of_a<NT> >  A_iterator;
    typedef  CGAL::Const_value_iterator<NT>
                                        B_iterator;
    typedef  CGAL::Const_value_iterator<NT>
                                        C_iterator;
    
    typedef  CGAL::QP_rep_signed_point_iterator< Point, Point_iterator>
                                        Signed_point_iterator;
    typedef  CGAL::Join_random_access_iterator_1<
                 Signed_point_iterator,
                 QP_rep_row_of_d< NT, Point, Signed_point_iterator,
                                  Access_coord, Access_dim > >
                                        D_iterator;
    

    typedef  CGAL::Tag_false        Is_lp;
};

// Function declarations
// =====================
// I/O operators
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Polytope_distance_d<Traits_>& poly_dist);

template < class Traits_ >
std::istream&
operator >> ( std::istream& is,
              Polytope_distance_d<Traits_>& poly_dist);

// ============================================================================

// Class implementation
// ====================

// validity check
template < class Traits_ >
bool
Polytope_distance_d<Traits_>::
is_valid( bool verbose, int level) const
{
    CGAL_USING_NAMESPACE_STD

    CGAL::Verbose_ostream verr( verbose);
    verr << "CGAL::Polytope_distance_d<Traits>::" << endl;
    verr << "is_valid( true, " << level << "):" << endl;
    verr << "  |P+Q| = " << number_of_points_p()
         <<          '+' << number_of_points_q()
         <<   ", |S| = " << number_of_support_points_p()
         <<          '+' << number_of_support_points_q() << endl;

    if ( is_finite()) {

        // compute normal vector
        ET_vector  normal( d), diff( d);
        ET  et_0 = 0, den = solver.variables_common_denominator();
        int i, j;
        for ( j = 0; j < d; ++j) normal[ j] = p_coords[ j] - q_coords[ j];

        // check P
        // -------
        verr << "  checking P..." << flush;
        
        // check point set
        for ( i = 0; i < number_of_points_p(); ++i) {
            for ( j = 0; j < d; ++j) {
                diff[ j] = p_coords[ j] - den
                  * tco.access_coordinates_begin_d_object()( p_points[ i])[ j];
            }
            if ( std::inner_product( diff.begin(), diff.end(),
                                     normal.begin(), et_0) > et_0)
                return CGAL::_optimisation_is_valid_fail( verr,
                           "polytope P is not separated by its hyperplane");
        }
        
        verr << "passed." << endl;

        // check Q
        // -------
        verr << "  checking Q..." << flush;
        
        // check point set
        for ( i = 0; i < number_of_points_q(); ++i) {
            for ( j = 0; j < d; ++j) {
                diff[ j] = q_coords[ j] - den
                  * tco.access_coordinates_begin_d_object()( q_points[ i])[ j];
            }
            if ( std::inner_product( diff.begin(), diff.end(),
                                     normal.begin(), et_0) < et_0)
                return CGAL::_optimisation_is_valid_fail( verr,
                           "polytope Q is not separated by its hyperplane");
        }
        
        verr << "passed." << endl;
    }

    verr << "  object is valid!" << endl;
    return( true);
}

// output operator
template < class Traits_ >
std::ostream&
operator << ( std::ostream& os,
              const Polytope_distance_d<Traits_>& poly_dist)
{
    CGAL_USING_NAMESPACE_STD

    typedef  Polytope_distance_d<Traits_>::Point  Point;
    typedef  ostream_iterator<Point>       Os_it;
    typedef  typename Traits_::ET          ET;
    typedef  ostream_iterator<ET>          Et_it;

    switch ( CGAL::get_mode( os)) {

      case CGAL::IO::PRETTY:
        os << "CGAL::Polytope_distance_d( |P+Q| = "
           << poly_dist.number_of_points_p() << '+'
           << poly_dist.number_of_points_q() << ", |S| = "
           << poly_dist.number_of_support_points_p() << '+'
           << poly_dist.number_of_support_points_q() << endl;
        os << "  P = {" << endl;
        os << "    ";
        copy( poly_dist.points_p_begin(), poly_dist.points_p_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  Q = {" << endl;
        os << "    ";
        copy( poly_dist.points_q_begin(), poly_dist.points_q_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S_P = {" << endl;
        os << "    ";
        copy( poly_dist.support_points_p_begin(),
              poly_dist.support_points_p_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  S_Q = {" << endl;
        os << "    ";
        copy( poly_dist.support_points_q_begin(),
              poly_dist.support_points_q_end(),
              Os_it( os, ",\n    "));
        os << "}" << endl;
        os << "  p = ( ";
        copy( poly_dist.realizing_point_p_coordinates_begin(),
              poly_dist.realizing_point_p_coordinates_end(),
              Et_it( os, " "));
        os << ")" << endl;
        os << "  q = ( ";
        copy( poly_dist.realizing_point_q_coordinates_begin(),
              poly_dist.realizing_point_q_coordinates_end(),
              Et_it( os, " "));
        os << ")" << endl;
        os << "  squared distance = "
           << poly_dist.squared_distance_numerator() << " / "
           << poly_dist.squared_distance_denominator() << endl;
        break;

      case CGAL::IO::ASCII:
        os << poly_dist.number_of_points_p() << endl;
        copy( poly_dist.points_p_begin(),
              poly_dist.points_p_end(),
              Os_it( os, "\n"));
        os << poly_dist.number_of_points_q() << endl;
        copy( poly_dist.points_q_begin(),
              poly_dist.points_q_end(),
              Os_it( os, "\n"));
        break;

      case CGAL::IO::BINARY:
        os << poly_dist.number_of_points_p() << endl;
        copy( poly_dist.points_p_begin(),
              poly_dist.points_p_end(),
              Os_it( os));
        os << poly_dist.number_of_points_q() << endl;
        copy( poly_dist.points_q_begin(),
              poly_dist.points_q_end(),
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
operator >> ( std::istream& is,
              CGAL::Polytope_distance_d<Traits_>& poly_dist)
{
    CGAL_USING_NAMESPACE_STD
    /*
    switch ( CGAL::get_mode( is)) {

      case CGAL::IO::PRETTY:
        cerr << endl;
        cerr << "Stream must be in ascii or binary mode" << endl;
        break;

      case CGAL::IO::ASCII:
      case CGAL::IO::BINARY:
        typedef  CGAL::Polytope_distance_d<Traits_>::Point  Point;
        typedef  istream_iterator<Point>             Is_it;
        poly_dist.set( Is_it( is), Is_it());
        break;

      default:
        CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
        break; }
    */
    return( is);
}

CGAL_END_NAMESPACE

#endif // CGAL_POLYTOPE_DISTANCE_D_H

// ===== EOF ==================================================================
