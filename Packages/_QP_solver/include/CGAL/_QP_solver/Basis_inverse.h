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
// file          : include/CGAL/_QP_solver/Basis_inverse.h
// package       : $CGAL_Package: _QP_solver $
//
// revision      : 0.2
// revision_date : 2000/08/11
//
// author(s)     : Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: Basis inverse used in the solver for Quadratic Programs
// ============================================================================
                                                                               

#ifndef CGAL_BASIS_INVERSE_H
#define CGAL_BASIS_INVERSE_H

// includes
#include <CGAL/Optimisation/basic.h>
#include <vector>
#include <CGAL/IO/Verbose_ostream.h>


CGAL_BEGIN_NAMESPACE
                    

// Class declaration
// =================
template < class ET_, class Is_lp_ >
class Basis_inverse;
                    

template < class ET, class Is_lp >
class Basis_inverse__entry_iterator;


// Class interface
// ===============
template < class ET_, class Is_lp_ >
class Basis_inverse {
  public:
    // self
    typedef  ET_                        ET;
    typedef  Is_lp_                     Is_lp;
    typedef  Basis_inverse<ET,Is_lp>    Self;

  private:
    
    // private types
    typedef  std::vector<ET>            Row;
    typedef  std::vector<Row>           Matrix;
    
    
    typedef  CGAL::Tag_true             Tag_true;
    typedef  CGAL::Tag_false            Tag_false;
    
    
    // friends
    friend  class Basis_inverse__entry_iterator<ET,Is_lp>;
    
    

  public:
    
    // types
    typedef  Basis_inverse__entry_iterator<ET,Is_lp>
                                        Entry_iterator;
                                                       

    
    // creation
    Basis_inverse( CGAL::Verbose_ostream&  verbose_out)
        : et_0( 0), et_1( 1), vout( verbose_out) { }
    
    
    // initialization
    template < class InputIterator1, class InputIterator2 >
    void
    init( unsigned int qp_m, unsigned int l,
          InputIterator1 u_it, InputIterator2 w_it,
          unsigned int max_size = 0)
        {
            CGAL_optimisation_precondition( qp_m > 0);
            m = qp_m;
            d = et_1;
            init( l, u_it, w_it, max_size, Is_lp());
            
            CGAL_optimisation_debug {
                if ( vout.verbose()) {
                    for ( unsigned int row = 0; row < k; ++row) {
                        std::copy( M[ row].begin(), M[ row].end(),
                                   std::ostream_iterator<ET>( vout.out()," "));
                        vout.out() << std::endl;
                    }
                    vout.out() << "denominator = " << d << std::endl;
                }
            }
             
        }
    
    
    // access
    const ET&  denominator( ) const { return d; }
    
    
    Entry_iterator  column_begin( unsigned int j) const
        { return Entry_iterator( M, 0, j); }
    Entry_iterator  column_end  ( unsigned int j) const
        { return Entry_iterator( M, k, j); }
    
    
    // multiplication functions
    template < class ForwardIterator, class OutputIterator >
    void
    multiply( ForwardIterator z_l, ForwardIterator z_x,
               OutputIterator y_l,  OutputIterator y_x) const
        {
            multiply( z_l, z_x, y_l, y_x, Is_lp());
        }
    
    
    // special multiplication functions for LPs
    template < class ForwardIterator, class OutputIterator >
    void
    multiply_l( ForwardIterator z_x, OutputIterator y_l) const
        {
            multiply_l( z_x, y_l, Is_lp());
        }
    
    template < class ForwardIterator, class OutputIterator >
    void
    multiply_x( ForwardIterator z_l, OutputIterator y_x) const
        {
            multiply_x( z_l, y_x, Is_lp());
        }
    
    
    // update functions
    /*
    template < class ForwardIterator > inline
    void  append( ForwardIterator q_l, ForwardIterator q_x, const ET& nu);
    
    
    inline
    void  remove( unsigned int i);
    
    
    template < class RandomAccessIterator > inline
    void  replace( unsigned int i, RandomAccessIterator q_x);
    */
    
    // swap function
    /*
    inline
    void  swap( unsigned int, unsigned int);
    */
    

  private:
    
    // some constants
    const ET                 et_0, et_1;
                                        

    
    // data members
    unsigned int             m;         // number of constraints
    unsigned int             k;         // size of matrix
    Matrix                   M;         // basis inverse, stored row-wise
    ET                       d;         // denominator
    
    
    CGAL::Verbose_ostream&   vout;      // used for verbose output
    
    

    
    // initialization
    /*
    template < class InputIterator1, class InputIterator2 >
    inline  void  init( unsigned int l,
                        InputIterator1 u_it, InputIterator2 w_it,
                        unsigned int max_size, Tag_false);
    template < class InputIterator1, class InputIterator2 >
    inline  void  init( unsigned int l,
                        InputIterator1 u_it, InputIterator2 w_it,
                        unsigned int max_size, Tag_true );
    */
    
    // multiplication functions
    /*
    template < class ForIt, class OutIt > inline                        // QP
    void  multiply( ForIt z_l, ForIt z_x,
                    OutIt y_l, OutIt y_x, Tag_false) const;
    template < class ForIt, class OutIt > inline                        // LP
    void  multiply( ForIt z_l, ForIt z_x,
                    OutIt y_l, OutIt y_x, Tag_true ) const;
    */
    
    // special multiplication functions for LPs
    /*
    template < class ForIt, class OutIt > inline                        // QP
    void  multiply_l( ForIt z_x, OutIt y_l, Tag_false) const;
    template < class ForIt, class OutIt > inline                        // LP
    void  multiply_l( ForIt z_x, OutIt y_l, Tag_true ) const;
    
    template < class ForIt, class OutIt > inline                        // QP
    void  multiply_x( ForIt z_l, OutIt y_x, Tag_false) const;
    template < class ForIt, class OutIt > inline                        // LP
    void  multiply_x( ForIt z_l, OutIt y_x, Tag_true ) const;
    */
  

  

// ============================================================================
                                                                               

// Class Implementation
// ====================

// initialization
// --------------

// initialization (QP)
template < class InputIterator1, class InputIterator2 > inline
void
init( unsigned int l,
      InputIterator1 u_it, InputIterator2 w_it,
      unsigned int max_size, Tag_false)
{
    k = 2*m;
    M.erase( M.begin(), M.end());
    M.reserve( max_size > m ? m+max_size : k+1);
    unsigned int  i;
    for ( i = 0; i < k; ) M.push_back( Row( ++i, et_0));
    for ( i = 0; i < m; ++i, ++u_it, ++w_it) {
        M[ m+i][ l] = *w_it;
        M[ m+i][ i] = *u_it;
    }
}

// initialization (LP)
template < class InputIterator1, class InputIterator2 > inline
void
init( unsigned int l,
      InputIterator1 u_it, InputIterator2 w_it,
      unsigned int, Tag_true)
{
    k = m;
    M = Matrix( m, Row( m, et_0));
    for ( unsigned int i = 0; i < m; ++i, ++u_it, ++w_it) {
        M[ i][ l] = *w_it;
        M[ i][ i] = *u_it;
    }
}


// multiplication functions
// ------------------------

// multiply (QP)
template < class ForIt, class OutIt > inline
void
multiply( ForIt z_l, ForIt z_x, OutIt y_l, OutIt y_x, Tag_false) const
{
    typename Matrix::const_iterator  matrix_it = M.begin();
    typename Row   ::const_iterator     row_it;     // left  of diagonal
    typename Matrix::const_iterator  column_it;     // right of diagonal
    ForIt                                 z_it;

    unsigned int  row, count;
    ET            sum;

    // compute  P z_l + Q^T z_x
    for ( row = 0; row < m; ++row,                                ++y_l) {
        sum = et_0;

        // P: left of diagonal (inclusive)
        for (   row_it =  matrix_it->begin(),                z_it = z_l;
                row_it != matrix_it->end();
              ++row_it,                                    ++z_it      ) {
            sum += *row_it * *z_it;
        }

        // P: right of diagonal (exclusive)
        for (   count = row+1,   column_it = ++matrix_it;
                count < m;
              ++count,         ++column_it,                ++z_it      ) {
            sum += (*column_it)[ row] * *z_it;
        }

        // Q^T:
        for (                                                z_it = z_x;
                count < k;
              ++count,         ++column_it,                ++z_it      ) {
            sum += (*column_it)[ row] * *z_it;
        }

        // store result
        *y_l = sum;
    }

    // compute  Q z_l + R z_x
    for ( ; row < k; ++row,                                       ++y_x) {
        sum = et_0;

        // Q:
        for (   count = 0,   row_it =  matrix_it->begin(),   z_it = z_l;
                count < m;
              ++count,     ++row_it,                       ++z_it      ) {
            sum += *row_it * *z_it;
        }

        // R: left of diagonal (inclusive)
        for (                                                z_it = z_x;
                             row_it != matrix_it->end();
                           ++row_it,                       ++z_it      ) {
            sum += *row_it * *z_it;
        }

        // R: right of diagonal (exclusive)
        for (   count = row+1,   column_it = ++matrix_it;
                count < k;
              ++count,         ++column_it,                ++z_it      ) {
            sum += (*column_it)[ row] * *z_it;
        }

        // store result
        *y_x = sum;
    }
}

// multiply (LP)
template < class ForIt, class OutIt > inline
void
multiply( ForIt z_l, ForIt z_x, OutIt y_l, OutIt y_x, Tag_true) const
{
    multiply_l( z_x, y_l);
    multiply_x( z_l, y_x);
}

// multiply_l (QP)
template < class ForIt, class OutIt > inline
void
multiply_l( ForIt z_x, OutIt y_l, Tag_false) const
{
    typename Matrix::const_iterator  matrix_it = M.begin()+m;
    typename Matrix::const_iterator  column_it;
    ForIt                                 z_it;

    unsigned int  row, count;
    ET            sum;

    // compute  Q^T z_x
    for ( row = 0; row < m; ++row,                                ++y_l) {
        sum = et_0;

        for (   count = 0,   column_it = matrix_it,   z_it = z_x;
                count < m;
              ++count,     ++column_it,             ++z_it      ) {
            sum += (*column_it)[ row] * *z_it;
        }

        *y_l = sum;
    }
}

// multiply_x (QP)
template < class ForIt, class OutIt > inline
void
multiply_x( ForIt z_l, OutIt y_x, Tag_false) const
{
    typename Matrix::const_iterator  matrix_it = M.begin()+m;
    typename Row   ::const_iterator     row_it;
    ForIt                                 z_it;

    unsigned int  row, count;
    ET            sum;

    // compute  Q z_l
    for ( row = 0; row < m; ++row,  ++matrix_it,                  ++y_x) {
        sum = et_0;

        for (   count = 0,   row_it = matrix_it->begin(),   z_it = z_l;
                count < m;
              ++count,     ++row_it,                      ++z_it      ) {
            sum += *row_it * *z_it;
        }

        *y_x = sum;
    }
}

// multiply_l (LP)
template < class ForIt, class OutIt > inline
void
multiply_l( ForIt z_x, OutIt y_l, Tag_true) const
{
    typename Matrix::const_iterator  matrix_it;
    ForIt                                 z_it;

    unsigned int  row;
    ET            sum;

    // compute  Q^T z_x
    for ( row = 0; row < m; ++row,                   ++y_l) {
        sum = et_0;

        for (   matrix_it =  M.begin(),   z_it = z_x;
                matrix_it != M.end();
              ++matrix_it,              ++z_it      ) {
            sum += (*matrix_it)[ row] * *z_it;
        }

        *y_l = sum;
    }
}

// multiply_x (LP)
template < class ForIt, class OutIt > inline
void
multiply_x( ForIt z_l, OutIt y_x, Tag_true) const
{
    typename Matrix::const_iterator  matrix_it = M.begin();
    typename Row   ::const_iterator     row_it;
    ForIt                                 z_it;

    ET  sum;

    // compute  Q z_l
    for ( ; matrix_it != M.end(); ++matrix_it,              ++y_x) {
        sum = et_0;

        for (   row_it =  matrix_it->begin(),   z_it = z_l;
                row_it != matrix_it->end();
              ++row_it,                       ++z_it      ) {
            sum += *row_it * *z_it;
        }

        *y_x = sum;
    }
}


// update functions
// ----------------

public:

// append
template < class ForIt > inline
void
append( ForIt q_l, ForIt q_x, const ET& nu)
{
    // check for QP
    Assert_compile_time_tag( Tag_false(), Is_lp());

    // handle sign of `nu'
    bool  nu_negative = ( nu < et_0);
    if ( nu_negative) d = -d;       // now `d' is `sgn(nu) * d'

    // update matrix in place
    // ----------------------
    typename Matrix::iterator  matrix_it = M.begin();
    typename Row   ::iterator     row_it;
    ForIt                           q_it1, q_it2;
    unsigned int  row, column;

    // rows: 0..m-1
    for (   row = 0,   q_it1 = q_l;
            row < m;
          ++row,     ++q_it1,      ++matrix_it) {

        // columns: 0..row
        for (                 row_it =  matrix_it->begin(),   q_it2 = q_l;
                              row_it != matrix_it->end();
                            ++row_it,                       ++q_it2      ){

            
            *row_it *= nu;
            *row_it -= *q_it1 * *q_it2;
            *row_it /= d;                       // without remainder!
                                                                     
        }
    }

    // rows: m..k-1
    for (              q_it1 = q_x;
            row < k;
          ++row,     ++q_it1,      ++matrix_it) {

        // columns: 0..m-1
        for (   column = 0,   row_it =  matrix_it->begin(),   q_it2 = q_l;
                column < m;
              ++column,     ++row_it,                       ++q_it2      ){

            
            *row_it *= nu;
            *row_it -= *q_it1 * *q_it2;
            *row_it /= d;                       // without remainder!
                                                                     
        }

        // columns: m..k-1
        for (                                                 q_it2 = q_x;
                              row_it != matrix_it->end();
                            ++row_it,                       ++q_it2      ){

            
            *row_it *= nu;
            *row_it -= *q_it1 * *q_it2;
            *row_it /= d;                       // without remainder!
                                                                     
        }
    }

    // store new row
    // -------------
    // allocate new row, if necessary
    // otherwise `matrix_it' already points to first unused row
    ++k;
    if ( k > M.size()) {
        matrix_it = M.insert( M.end(), Row( k, et_0));
    }

    // store entries in new row
    for (   column = 0,   row_it = matrix_it->begin();
            column < m;
          ++column,     ++row_it                     , ++q_l) {
        *row_it = ( nu_negative ? -( *q_l) : *q_l);
    }
    for ( ;
            column < k-1;
          ++column,     ++row_it                     , ++q_x) {
        *row_it = ( nu_negative ? -( *q_x) : *q_x);
    }
    *row_it = -d;

    // store new denominator
    // ---------------------
    d = ( nu_negative ? -nu : nu);
    CGAL_optimisation_postcondition( d > et_0);

    
    CGAL_optimisation_debug {
        if ( vout.verbose()) {
            for ( unsigned int row = 0; row < k; ++row) {
                std::copy( M[ row].begin(), M[ row].end(),
                           std::ostream_iterator<ET>( vout.out(), " "));
                vout.out() << std::endl;
            }
            vout.out() << "denominator = " << d << std::endl;
        }
    }
     
}

// remove
void
remove( unsigned int i)
{
    // check for QP
    Assert_compile_time_tag( Tag_false(), Is_lp());

    // check for last row/column
    CGAL_optimisation_precondition( m+i == k-1);

    // get `q^T = ( q_l^T | q_x^T)^T' and `d' from last row
    // ----------------------------------------------------
    --k;
    typename Row::const_iterator      q = M[ k].begin();
    ET                            new_d = M[ k][ k];

    // handle sign of `nu'
    if ( new_d > et_0) { d = -d; new_d = -new_d; }

    // update matrix in place
    // ----------------------
    typename Matrix::iterator        matrix_it = M.begin();
    typename Row   ::iterator           row_it;
    typename Row   ::const_iterator       q_it1, q_it2;
    unsigned int  row;

    // rows: 0..k-1
    for (   row = 0,   q_it1 = q;
            row < k;
          ++row,     ++q_it1,      ++matrix_it) {

        // columns: 0..row
        for (   row_it =  matrix_it->begin(),   q_it2 = q;
                row_it != matrix_it->end();
              ++row_it,                       ++q_it2      ) {

            *row_it *= new_d;
            *row_it += *q_it1 * *q_it2;
            *row_it /= d;                       // without remainder!
        }
    }

    // store new denominator
    // ---------------------
    d = -new_d;
    CGAL_optimisation_postcondition( d > et_0);

    
    CGAL_optimisation_debug {
        if ( vout.verbose()) {
            for ( unsigned int row = 0; row < k; ++row) {
                std::copy( M[ row].begin(), M[ row].end(),
                           std::ostream_iterator<ET>( vout.out(), " "));
                vout.out() << std::endl;
            }
            vout.out() << "denominator = " << d << std::endl;
        }
    }
     
}

// replace
template < class RanIt > inline
void
replace( unsigned int i, RanIt q_x)
{
    // check for LP
    Assert_compile_time_tag( Tag_true(), Is_lp());

    // check index
    CGAL_optimisation_precondition( i < m);

    // update matrix in place
    // ----------------------
    typename Matrix::      iterator  matrix_it = M.begin();
    typename Row   ::const_iterator   row_i_it;
    typename Row   ::      iterator     row_it;

    ET            new_d = q_x[ i];
    unsigned int  row;

    // handle sign of `new_d'
    bool  new_d_negative = ( new_d < et_0);
    if ( new_d_negative) { d = -d; new_d = -new_d; }

    // rows: 0..i-1
    for (   row = 0;
            row < i;
          ++row,     ++matrix_it, ++q_x) {

        
        // columns: 0..m-1
        for (   row_it =  matrix_it->begin(),   row_i_it = M[ i].begin();
                row_it != matrix_it->end();
              ++row_it                      , ++row_i_it                ) {
        
            *row_it *= new_d;
            *row_it -= *q_x * *row_i_it;
            *row_it /= d;                       // without remainder!
        }
         
    }

    // rows: i (flip signs, in necessary)
    if ( new_d_negative) {
        for (   row_it =  matrix_it->begin();
                row_it != matrix_it->end();
              ++row_it                      ) {

            *row_it = -( *row_it);
        }
    }

    // rows: i+1..m-1
    for ( ++row,     ++matrix_it, ++q_x;
            row < m;
          ++row,     ++matrix_it, ++q_x) {

        
        // columns: 0..m-1
        for (   row_it =  matrix_it->begin(),   row_i_it = M[ i].begin();
                row_it != matrix_it->end();
              ++row_it                      , ++row_i_it                ) {
        
            *row_it *= new_d;
            *row_it -= *q_x * *row_i_it;
            *row_it /= d;                       // without remainder!
        }
         
    }

    // store new denominator
    // ---------------------
    d = new_d;
    CGAL_optimisation_postcondition( d > et_0);

    
    CGAL_optimisation_debug {
        if ( vout.verbose()) {
            for ( unsigned int row = 0; row < k; ++row) {
                std::copy( M[ row].begin(), M[ row].end(),
                           std::ostream_iterator<ET>( vout.out(), " "));
                vout.out() << std::endl;
            }
            vout.out() << "denominator = " << d << std::endl;
        }
    }
     
}


// swap function
// -------------
void
swap( unsigned int i, unsigned int j)
{
    // check for QP
    Assert_compile_time_tag( Tag_false(), Is_lp());

    // check indices
    CGAL_optimisation_precondition( m+i < k);
    CGAL_optimisation_precondition( m+j < k);
    if ( i == j) return;

    // guarantee `i < j'
    if ( i > j) std::swap( i, j);

    // swap rows and columns
    // ---------------------
    i += m;
    j += m;
    typename    Row::iterator   row_i_it = M[ i].begin();
    typename    Row::iterator   row_j_it = M[ j].begin();
    typename Matrix::iterator  matrix_it = M.begin()+(i+1);
    unsigned int  count;

    // swap entries 0..i-1 (row <-> row)
    for (   count = 0;
            count < i;
          ++count,      ++row_i_it, ++row_j_it) {
        std::iter_swap( row_i_it, row_j_it);
    }

    // swap entries i+1..j-1 (column <-> row)
    for ( ++count,                  ++row_j_it;
            count < j;
          ++count,     ++matrix_it, ++row_j_it) {
        std::swap( ( *matrix_it)[ i], *row_j_it);
    }

    // swap entries j+1..k (column <-> column)
    for ( ++count,     ++matrix_it;
            count < j;
          ++count,     ++matrix_it) {
        std::swap( ( *matrix_it)[ i], ( *matrix_it)[ j]);
    }

    // swap entries i,i with j,j (entry <-> entry)
    std::iter_swap( row_i_it, row_j_it);
}
 
};

template < class ET, class Is_lp >
class Basis_inverse__entry_iterator {
  public:
    typedef  ET                               value_type;
    typedef  ptrdiff_t                        difference_type;
    typedef  value_type*                      pointer;
    typedef  value_type&                      reference;
    typedef  std::random_access_iterator_tag  iterator_category;

    typedef  Basis_inverse__entry_iterator<ET,Is_lp>  Self;
    typedef  value_type                               Val;
    typedef  difference_type                          Dist;

    typedef  typename Basis_inverse<ET,Is_lp>::Matrix          Matrix;

    // creation
    Basis_inverse__entry_iterator( const Matrix& M, Dist i, Dist j)
        : matrix( M), row( i), col( j) { }

    // compare operations
    bool   operator == ( const Self& it) const { return ( row == it.row); }
    bool   operator != ( const Self& it) const { return ( row != it.row); }

    // access
    Val    operator *  ( ) const
        { if ( ( ! CGAL::check_tag( Is_lp())) && ( col > row))
              return matrix[ col][ row];
          return matrix[ row][ col]; }

    // forward operations
    Self&  operator ++ ( )    {                   ++row; return *this; }
    Self   operator ++ ( int) { Self tmp = *this; ++row; return tmp;   }

    // bidirectional operations
    Self&  operator -- (    ) {                   --row; return *this; }
    Self   operator -- ( int) { Self tmp = *this; --row; return tmp;   }

    // random access operations
    Self&  operator += ( Dist n) { row += n; return *this; }
    Self&  operator -= ( Dist n) { row -= n; return *this; }

    Self   operator +  ( Dist n) const { Self tmp=*this; return tmp+=n; }
    Self   operator -  ( Dist n) const { Self tmp=*this; return tmp-=n; }

    Dist   operator - ( const Self& it) const { return row - it.row; }

    Val    operator [] ( Dist i) const
        { if ( ( ! CGAL::check_tag( Is_lp())) && ( col > row+i))
              return matrix[ col][ row+i];
          return matrix[ row+i][ col]; }

    bool   operator <  ( const Self&) const { return ( row <  it.row); }
    bool   operator >  ( const Self&) const { return ( row >  it.row); }
    bool   operator <= ( const Self&) const { return ( row <= it.row); }
    bool   operator >= ( const Self&) const { return ( row >= it.row); }

  private:
    const Matrix&  matrix;
    Dist  row;
    Dist  col;
};

CGAL_END_NAMESPACE
                  

#endif // CGAL_BASIS_INVERSE_H

// ===== EOF ==================================================================
