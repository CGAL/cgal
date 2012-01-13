// Copyright (c) 1997  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Sariel Har-Peled (sariel@math.tau.ac.il)
//                 Eyal Flato (flato@math.tau.ac.il)

#ifndef  CGAL_KDTREE_D_H
#define  CGAL_KDTREE_D_H

#include <CGAL/basic.h>
#include <cstdlib>
#include <cstring>
#include <list>
using std::list; // to avoid compiler crash on MSVC++

namespace CGAL {

/*=======================================================================
 * Kdtree_interface -
 *    This is the default interface of point. It assume that PT (the
 * point type have the following properties:
 *    default constructor
 *    int   dimension()
 *    const   coord_type  & operator[]( int ) const
 *    ...                   operator=( const Pt  & ) - copy operator
\*=======================================================================*/

template <class  PT>
class  Kdtree_interface
{
public:
    typedef  PT    Point;

    static  int   dimension( const PT  & pnt )
    {
      //        return  pnt.dimensions();
      return pnt.dimension();
    }

    static  int    compare( int  d, const PT   & a, const PT    & b )
    {
        if  ( a[ d ] < b[ d ] )
            return  -1;
        if  ( a[ d ] > b[ d ] )
            return  1;

        return  0;
    }
    static  void   copy_coord( int  d, PT   & a, const PT    & b )
    {
        a[ d ] = b[ d ];
    }
};



template <class  PT>
class  Kdtree_interface_2d
{
public:
    typedef  PT    Point;

    static  int   dimension( const PT  & pnt )
    {
      //        return  pnt.dimensions();
      return pnt.dimension();
    }

    static  int    compare( int  d, const PT   & a, const PT    & b )
    {
        if  ( a[ d ] < b[ d ] )
            return  -1;
        if  ( a[ d ] > b[ d ] )
            return  1;

        return  0;
    }
    static  void   copy_coord( int  d, PT   & a, const PT    & b )
    {
      if  ( d == 0 )
        a = PT( b[ 0 ], a[1] );
      else
        if  ( d == 1 )
          a = PT( a[ 0 ], b[1] );
        else {
          CGAL_error();
        }
    }
};


template <class  PT>
class  Kdtree_interface_3d
{
public:
    typedef  PT    Point;

    static  int   dimension( const PT  & pnt )
    {
      //        return  pnt.dimensions();
      return pnt.dimension();
    }

    static  int    compare( int  d, const PT   & a, const PT    & b )
    {
        if  ( a[ d ] < b[ d ] )
            return  -1;
        if  ( a[ d ] > b[ d ] )
            return  1;

        return  0;
    }
    static  void   copy_coord( int  d, PT   & a, const PT    & b )
    {
      if  ( d == 0 )
        a = PT( b[ 0 ], a[1], a[ 2 ] );
      else
        if  ( d == 1 )
          a = PT( a[ 0 ], b[1], a[ 2 ] );
        else
          if  ( d == 2 )
            a = PT( a[ 0 ], a[1], b[ 2 ] );
          else {
            CGAL_error();
          }
    }
};


/*=========================================================================
 * kdtree_d -
 *     A the kdtree class.
 *
 * Remark: The kd-trees allocates all the memory it needs in advance,
 *   This results in a rather efficient memory management, and a fast
 *   iplementation.
\*=========================================================================*/
template <class Traits>
class Kdtree_d
{
  class Plane;

public:
    typedef typename Traits::Point      Point;
    typedef list<Point>   List_points;


    //-------------------------------------------------------------------------
    // Extended_Point_d -
    //    A class for representing an extended d-dimal point: with
    //  ability to support +/- infinity. Templated by the kd-treee interface
    //  type.
    //-------------------------------------------------------------------------
    class  ExtPoint
    {
    private:
        class   coordinate_type
        {
        public:
            const Point     * p_pnt;
            int      type;
            //signed char   type;
        };

        coordinate_type      * p_arr;
        int  dim;
        Point  def_pnt;

        void   init( int  _dim )
        {
            CGAL_precondition( _dim > 0 );
            dim = _dim;

            p_arr = (coordinate_type *)std::malloc( sizeof( coordinate_type )
                                               * dim );
            //printf( "p_arr(new): %p\n", (void *)p_arr );
            CGAL_assertion( p_arr != NULL );

			std::memset( p_arr, 0, sizeof( coordinate_type ) * dim );
        }

    public:
        enum { MINUS_INFINITY = -1, FINITE = 0, PLUS_INFINITY = 1 };

        ExtPoint( int  _type, int  _dim )
        {
            CGAL_precondition( _type == MINUS_INFINITY  ||  _type == PLUS_INFINITY );
            init( _dim );

            for  ( int  ind = 0; ind < dim; ind++ )
                p_arr[ ind ].type = _type;
        }

        ExtPoint()
        {
            p_arr = NULL;
            dim = -1;
        }

        ExtPoint( const ExtPoint  & p )
        {
            init( p.dim );

            def_pnt = p.def_pnt;
            for  ( int  ind = 0; ind < dim; ind++ ) {
                p_arr[ ind ] = p.p_arr[ ind ];
                if  ( p.p_arr[ ind ].p_pnt == &p.def_pnt )
                    p_arr[ ind ].p_pnt = &def_pnt;
            }
        }

        ExtPoint & operator=( const ExtPoint  & p )
        {
            term();

            init( p.dim );

            def_pnt = p.def_pnt;
            for  ( int  ind = 0; ind < dim; ind++ ) {
                p_arr[ ind ] = p.p_arr[ ind ];
                if  ( p.p_arr[ ind ].p_pnt == &p.def_pnt )
                    p_arr[ ind ].p_pnt = &def_pnt;
            }

            return  *this;
        }

        ExtPoint( const Point   & point, int  _dim )
        {
            init( _dim );

            def_pnt = point;
            for  ( int  ind = 0; ind < _dim; ind++ ) {
                p_arr[ ind ].p_pnt = &def_pnt;
                p_arr[ ind ].type = FINITE;
            }
        }

        void  term()
        {
            if  ( p_arr != NULL ) {
                //printf( "a: %p\n", p_arr );
                std::free( p_arr );
                //printf( "_a\n" );
                p_arr = NULL;
            }
            dim = 0;
        }

        ~ExtPoint()
        {
            term();
        }

        void   set_coord( int  k, Point  & point )
        {
            CGAL_precondition( 0 <= k  &&  k < dim );

            p_arr[ k ].type = FINITE;
            //p_arr[ k ].p_pnt = &point;

            Traits::copy_coord( k, def_pnt, point );
            p_arr[ k ].p_pnt =  &def_pnt;
        }

        void   set_coord( int  k, const ExtPoint  & point )
        {
            CGAL_precondition( 0 <= k  &&  k < dim );
            CGAL_precondition( 0 <= k  &&  k < point.dim );

            p_arr[ k ] = point.p_arr[ k ];
            if  ( p_arr[ k ].type == FINITE ) {
                Traits::copy_coord( k, def_pnt, *(p_arr[ k ].p_pnt) );
                p_arr[ k ].p_pnt =  &def_pnt;
            }
        }

        int    compare( int   k, const ExtPoint  & point )  const
        {
            CGAL_precondition( 0 <= k  &&  k < dim );
	    // the following does not compile on msvc++...
	    //            coordinate_type  & a( p_arr[ k ] ),
	    //  & b( point.p_arr[ k ] );
            coordinate_type & a = p_arr[ k ];
            coordinate_type & b = point.p_arr[ k ];

            if  ( a.type != FINITE ) {
                if  ( b.type != FINITE ) {
                    return  a.type - b.type;
                } else {
                    return  a.type;
                }
            } else {
                if  ( b.type != FINITE )
                    return  -b.type;
                else
                    return  Traits::compare( k, *a.p_pnt, *b.p_pnt );
            }
        }


        int    compare_vector( const ExtPoint  & point )  const
        {
            int  ind, res;

            for  ( ind = 0; ind < dim; ind++ ) {
                res = compare( ind, point );
                if  ( res != 0 )
                    return  res;
            }

            return  0;
        }

        int    compare( int   k, const  Point  & point )  const
        {
            CGAL_precondition( 0 <= k  &&  k < dim );

	    //            coordinate_type  & a( p_arr[ k ] );
            coordinate_type  & a = p_arr[ k ];

            if  ( a.type != FINITE )
                return  a.type;
            else
                return  Traits::compare( k, *a.p_pnt, point );
        }

        int  dimension() const
        {
            return  dim;
        }

        int   get_coord_status( int   d ) const
        {
            CGAL_precondition( 0 <= d  &&  d < dim );

            return  p_arr[ d ].type;
        }


        const Point  * get_coord_point( int   d ) const
        {
            CGAL_precondition( 0 <= d  &&  d < dim );

            return  p_arr[ d ].p_pnt;
        }
    };


    // Box - represents an axis parallel box.
    class Box
    {
    public:
        int  dim;
        ExtPoint  left, right;

    private:
        friend  class Plane;

        ExtPoint  & get_vertex( bool  f_left )
        {
            return  f_left? left : right;
        }

    public:
        // constructors

      //CHECK
        Box()
        {
        }

        Box( const Box  & box )
        {
            dim = box.dim;
            left = box.left;
            right = box.right;
        }


        Box( const Point  &l, const Point  &r, int  _dim )
        {
            dim = _dim;
            left = ExtPoint( l, dim );
            right = ExtPoint( r, dim );
        }

        Box( int  _dim ) : left( ExtPoint::MINUS_INFINITY, _dim ),
            right( ExtPoint::PLUS_INFINITY, _dim )
        {
            dim = _dim;
        }


        // data access
        void set_left( Point &l)
        {
            left = ExtPoint( l, dim );
        };

        void set_right( Point &r)
        {
            right = ExtPoint( r, dim );
        }

        const ExtPoint &get_left() const
        {
            return left;
        }

        const ExtPoint &get_right() const
        {
            return right;
        }


        void   set_coord_left( int  k, Point  & p )
        {
            left.set_coord( k, p );
        }


        void   set_coord_right( int  k, Point  & p )
        {
            right.set_coord( k, p );
        }


        // operations
        bool   is_in( const Box &o) const
            // checks if o is completely inside <this>
        {
            int dim = left.dimension();

            for (int i = 0; i < dim; i++)
            {
                if ( (left.compare(i, o.get_left()) > 0 )
                     ||  (right.compare(i, o.get_right()) < 0 ) )
                    return false;
            }

            return true;
        }


        bool   is_in( const Point   & o ) const
            // checks if o is completely inside <this>
        {
            int _dim = left.dimension();
            for (int i = 0; i < _dim; i++)
            {
                if ( (left.compare( i, o ) > 0 ) ||
                     (right.compare( i, o ) <= 0 ) )
                    return false;
            }
            return true;
        }


        bool   is_coord_in_range( int  k, const Point   & o ) const
            // checks if o is completely inside <this>
        {
            return  ( ! ( (left.compare( k, o ) > 0 )
                          ||  (right.compare( k, o ) <= 0 ) ) );
        }


        bool is_intersect( const Box &o) const
            // checks if there is an intersection between o and this
        {
            int dim = left.dimension();
            for (int i = 0; i < dim; i++)
            {
                if ( (left.compare(i, o.get_right()) >= 0) ||
                     (right.compare(i, o.get_left()) <= 0) )
                    return false;
            }
            return true;
        }


        // checks if there is an intersection between o and this box
        // only in a specific coordinate...
        bool is_intersect_in_dim( int d, const Box  & o ) const
        {
            return  (! ( (left.compare( d, o.get_right() ) >= 0 )
                         ||  (right.compare( d, o.get_left() ) <= 0) ));
        }


        // checks if there is an intersection between o and this box
        // only in a specific coordinate...
        bool is_intersect_in_dim_closed( int d, const Box  & o ) const
        {
            return  (! ( (left.compare( d, o.get_right() ) > 0 )
                         ||  (right.compare( d, o.get_left() ) < 0) ));
        }


        bool intersect(Box &o)
            // intersects this with o. the intersection will be in this
            // returns false if intersection is empty
        {
            int dim = left.dimension();
            for (int i = 0; i < dim; i++)
            {
                // left is the maximal of the lefts
                if (left.compare(i, o.get_left()) == -1)
                    left.set_coord(i, o.get_left());

                // right is the minimal of the rights
                if (right.compare(i, o.get_right()) == 1)
                    right.set_coord(i, o.get_right());
            }
            return !(is_empty());
        }

        bool is_empty() const
            // return true if this is not an interval (left[k] > right[k])
        {
            int dim = left.dimension();
            for (int i = 0; i < dim; i++)
            {
                if (left.compare(i, right) == 1)
                    return true;
            }
            return false;
        }

        bool is_empty_open() const
            // return true if this is not an interval (left[k] > right[k])
        {
            int dim = left.dimension();
            for (int i = 0; i < dim; i++)
            {
                if ( left.compare(i, right) >= 0 )
                    return true;
            }
            return false;
        }

        int  comp( const Box     & o ) const
        {
            int  res;

            res = left.compare_vector( o.left );
            if  ( res != 0 )
                return  res;

            return  right.compare_vector( o.right );
        }
	// destructor - needed in ...recursive ... DVP

	~Box()
	  {
	    left.term();
	    right.term();
	    dim = 0;
	  }

    };


private:
    class Plane
    {
    private:
        int  coord;
        Point  * normal;
        bool  f_plus;  // orientation of half space
        //    is (0, 0, ... , +inifinity, 0, ..., 0) inside plane

    public:
        Plane()
        {
            normal = NULL;
            coord = 0;
        }

        Plane( int  k, Point  & p )
        {
            normal = &p;
            coord = k;
        }

        Plane( const Plane  & p )
        {
            coord = p.coord;
            normal = p.normal;
        }

        void  dump( void )
        {
            std::cout << "(" << coord << ": " << *normal << ")";
        }

        bool   is_in( const Point    & p ) const
        {
            int  cmp;

            cmp = Traits::compare( coord, p, *normal );

            if  ( ! f_plus )
                cmp = -cmp;

            return  cmp >= 0;
        }

        void set_plane(int k, Point &p)
        {
            coord = k;
            normal = &p;
            //normal->copy( coord, p );
        }

        void    split( Box   & region, bool  f_neg )
        {
            ExtPoint  * p_p = &(region.get_vertex( ! f_neg ));

            if  ( f_neg ) {
                if  ( p_p->compare( coord, *normal ) > 0 )
                    p_p->set_coord( coord, *normal );
            } else
                if  ( p_p->compare( coord, *normal ) < 0 )
                    p_p->set_coord( coord, *normal );
        }

        void  orient_half_space( bool   f_neg_side )
        {
            f_plus = ! f_neg_side;
        }

        int   get_coord() const
        {
            return  coord;
        }
    };

private:
    class Node
    {
    public:
        Plane  plane;
        Point  * pnt;

        Node  * left, * right;

        enum { LEFT, RIGHT };

        const Plane    & get_hs( int  side ) const
        {

            ((Plane *)&plane)->orient_half_space( side == LEFT );

            return   plane;
        }


        bool  is_points_in_hs( const  Plane    & pl ) const
        {
            if  ( is_point() )
                return  pl.is_in( *pnt );

            if  ( left != NULL  &&  ( ! left->is_points_in_hs( pl ) ) )
                return  false;
            if  ( right != NULL  &&  ( ! right->is_points_in_hs( pl ) ) )
                return  false;

            return  true;
        }


        bool  is_valid()  const
        {
            if  ( is_point() )
                return  true;

            if  ( left != NULL )
                if  ( ! left->is_points_in_hs( get_hs( LEFT ) ) )
                    return  false;
            if  ( right != NULL )
                if  ( ! right->is_points_in_hs( get_hs( RIGHT ) ) )
                    return  false;

            return  true;
        }


        void  dump( int  depth )
        {
            int  ind;

            for  ( ind = 0; ind < depth; ind++ )
                std::cout << " ";

            if  ( is_point() ) {
                std::cout << *pnt << "\n";
                return;
            }

            plane.dump();
            std::cout << "\n";
            left->dump( depth + 1 );
            for  ( ind = 0; ind < depth; ind++ )
                std::cout << " ";

            std::cout << "!!!!!!!!!!!!\n";
            right->dump( depth + 1 );
        }

        bool is_point()  const
        {
            return   ((left == NULL)  &&  (right == NULL));
        }

        typedef std::back_insert_iterator<List_points>  back_iter;

        Node() : plane()
        {
            left = right = NULL;
        }

        void  copy_subtree_points( back_iter  & result,
                                   const Box  & rect )
        {
            if  ( is_point() ) {
                if  ( rect.is_in( *pnt ) )
                    (*result++) = *pnt;
                return;
            }
            if  ( left != NULL )
                left->copy_subtree_points( result, rect );
            if  ( right != NULL )
                right->copy_subtree_points( result, rect );
        }

        static void search_recursive( back_iter  & result,
                                      Node  * node,
                                      const Box  & rect,
                                      Box  & _region,
                                      Plane  & plane,
                                      bool  f_split_plus )
        {
            //printf( "search_recusrive\n" );
            Box  * p_r = new Box( _region );

            //printf( "z" );
            //fflush( stdout );

            plane.split( *p_r, f_split_plus );

            //printf( "c" );
            //fflush( stdout );
            CGAL_precondition( node != NULL );

            //printf( "b" );
            //fflush( stdout );
            if  ( rect.is_in( *p_r ) )
            {
                //printf( "5" );
                //fflush( stdout );
                node->copy_subtree_points( result, rect );
                //printf( "\tsearch_recursive done...\n" );
                //printf( "6" );
                //fflush( stdout );
		delete p_r;
                return;
            }

            //printf( "v" );
            //fflush( stdout );

            if  ( rect.is_intersect_in_dim_closed( plane.get_coord(), *p_r ) )
                node->search( result, rect, *p_r );
            //printf( "x" );
            //fflush( stdout );
	    delete p_r;
        }

        void search( std::back_insert_iterator<List_points>  result,
                     const  Box  &rect, Box &region )
        {
            if (is_point()) {
                if  ( rect.is_in( *pnt ) )
                    (*result++) = *pnt;
                return;
            }

            //this is not a point so it is a hypeplane
            if  ( left != NULL )
                search_recursive( result, left, rect,
                                  region, plane, true );
            if  ( right != NULL )
                search_recursive( result, right, rect,
                                  region, plane, false );
        }
    };

    typedef Point  * Point_ptr;

    int  size;
    Point  * p_arr_pt;
    Node *root;
    int dim;

    Node  * p_node_arr;
    int  node_count;

    Node   * malloc_node( void )
    {
        Node  * p_ret;

        p_ret = &(p_node_arr[ node_count ]);
        node_count++;

        CGAL_assertion( node_count <= ( 2 * size ));

        *p_ret = Node();

        return  p_ret;
    }

    static  int         comp( const Point  & a, const Point  & b, int dim )
    {
        return  Traits::compare( dim, a, b );
    }

    static int  partition( Point_ptr  * arr, int left, int right,
                           Point      * p_pivot, int  dim )
    {
        int  i, j;
        Point_ptr  tmp;

        if  ( left >= right )
            return  left;

        i = left;
        j = right;

        while  ( i < j )  {
            if  ( comp( *(arr[ i ]), *(arr[ j ]), dim ) > 0 ) {
                tmp = arr[ i ];
                arr[ i ] = arr[ j ];
                arr[ j ] = tmp;
            }
            if  ( comp( *(arr[ i ]), *p_pivot, dim ) < 0 ) {
                i++;
            } else
                if  ( comp( *p_pivot, *(arr[ j ]), dim ) <= 0 )
                    j--;
        }

        return  (i > left)? i - 1 : left;
    }


    /* split the array into two sub-arrays, such that all the elements
     * from left to pos_mid are smaller than the elements from pos+1 to
     * right.
     */
    static  void   split_arr( Point_ptr      * arr,
                              int           left,
                              int           right,
                              int           pos_mid,
                              int           dim )
    {
        int  pos;

        if  ( left >= right )
            return;

        pos = partition( arr, left, right, arr[ (left + right ) / 2 ],
                         dim );
        if  ( pos == pos_mid )
            return;

        if  ( pos < pos_mid )
            split_arr( arr, pos+1, right, pos_mid, dim );
        else
            split_arr( arr, left, pos, pos_mid, dim );
    }


    static Point_ptr   get_max_element( Point_ptr   * arr,
                                        int  left, int  right, int  d )
    {
        int  max_pos = left;
        Point  mx = *(arr[ max_pos ]);

        for  ( int  ind = left + 1; ind <= right; ind++ )
            if  ( comp( mx, *(arr[ ind ]), d ) < 0 ) {
                mx = *(arr[ ind ]);
                max_pos = ind;
            }

        return  arr[ max_pos ];
    }

    Node *build_r( Point_ptr      * arr, int  left, int   right,
                   int  d )
    {
        int  num, pos, next_d;
        Node  * n;

        num = right - left + 1;

        if ( num < 1)
            return NULL;

        // if the list contains only one point,
        // construct a leaf for this node
        if  ( num == 1) {
            //n = new node;
            n = malloc_node();
            n->pnt = arr[ left ];

            return  n;
        }

        // else divide space into two regions in
        // dim-dim and cotinue recursively
        pos = (left + right) / 2;
        split_arr( arr, left, right, pos, d );

        Point  * p_median = get_max_element( arr, left, pos, d );

        // create division plane;
        Plane  plane( d, *p_median );
        //n = new node;
        // CGAL_assertion( n != NULL );
        n = malloc_node();

        n->plane = plane;

        next_d = d + 1;
        if  ( next_d >= dim )
            next_d = 0;

        // build left sub-tree
        n->left = build_r( arr, left, pos, next_d );

        // build right sub-tree
        n->right = build_r( arr, pos + 1, right, next_d );

        return n;
    }

public:
    typedef list<Point> list_points;

    Kdtree_d(int k = 2)
    {
        dim = k;
        root = NULL;
        p_arr_pt = NULL;
        p_node_arr = NULL;
    }

    ~Kdtree_d()
    {
        delete_all();
    }


    bool  is_valid( bool  verbose = false, int  level = 0 ) const
    {
        (void)verbose;
        (void)level;

        if  ( root == NULL )
            return  true;

        return  root->is_valid();
    }


    void  dump()
    {
        root->dump( 0);
    }

    void delete_all()
    {
        root = NULL;

        if  ( p_arr_pt != NULL )
            delete[] p_arr_pt;
        p_arr_pt = NULL;
        if  ( p_node_arr != NULL )
            delete[]  p_node_arr;
        p_node_arr = NULL;
    }

    void search( std::back_insert_iterator<list_points>   result,
                 Box  & rect )
    {
        if (root == NULL)
            return; // it is an empty tree - nothing to search in

        Box region = Box( dim );

        root->search( result, rect, region );
    }

    void build(list<Point> &l)
    {
        int  i;
        Point_ptr  * p_arr;

        size = l.size();
        p_arr_pt = new  Point[ size ];
        CGAL_assertion( p_arr_pt != NULL );

        p_arr = new  Point_ptr[ size ];
        CGAL_assertion( p_arr != NULL );

        p_node_arr = new Node[ 2 * size ];
        CGAL_assertion( p_node_arr != NULL );

        node_count = 0;

        /* fill the array */
        i = 0;
        for ( typename std::list<Point>::iterator j = l.begin();
              j != l.end();
              j++ )
        {
            p_arr_pt[ i ] = (*j);
            p_arr[ i ] = p_arr_pt + i;
            i++;
        }

        // recursively build the tree from the sorted list
        // starting to divide it in coordinate 0
        root = build_r( p_arr, 0, size - 1, 0 );

        //printf( "b\n" );
        delete[]  p_arr;
    }
};


} //namespace CGAL

#endif  /* CGAL_KDTREE_D_H */
