// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Michael Seel, Kurt Mehlhorn

#ifndef CGAL_FUNCTION_OBJECTSCD_H
#define CGAL_FUNCTION_OBJECTSCD_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Referenced_argument.h>

#undef CGAL_KD_TRACE
#undef CGAL_KD_TRACEN
#undef CGAL_KD_TRACEV
#define CGAL_KD_TRACE(t)  std::cerr << t
#define CGAL_KD_TRACEN(t) std::cerr << t << std::endl
#define CGAL_KD_TRACEV(t) std::cerr << #t << " = " << (t) << std::endl
 
namespace CGAL {

template <typename K>
class Compute_coordinateCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef FT                         result_type;
  result_type 
    operator()(const Point_d& p, int i) const
  {
    return p.cartesian(i);
  }
};

template <typename K>
class Point_dimensionCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef int                       result_type;
  result_type 
    operator()(const Point_d& p) const
  {
    return p.dimension();
  }
};

template <typename K>
class Less_coordinateCd {
  typedef typename K::FT             FT;
  typedef typename K::Point_d        Point_d;
  public:
  typedef bool                       result_type;
  result_type 
  operator()(const Point_d& p, const Point_d& q, int i) const
  {
    return p.cartesian(i)<q.cartesian(i);
  }
};

template <class R>
class Lift_to_paraboloidCd
{
    typedef typename R::Point_d Point;
    typedef typename R::FT FT;
    typedef typename R::LA LA;
public:
    typedef Point result_type;

    result_type operator()(const Point & p) const
    { 
        int d = p.dimension();
        typename LA::Vector h(d+1);
        FT sum = 0;
        for (int i = 0; i<d; i++) {
            h[i] = p.cartesian(i);
            sum += h[i]*h[i];
        }
        h[d] = sum;
        return Point(d+1,h.begin(),h.end());
    }
};

template <class R>
class Project_along_d_axisCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::FT FT;
public:
    typedef Point_d result_type;

    result_type operator()(const Point_d & p) const
    {
        return Point_d(p.dimension()-1,
               p.cartesian_begin(), p.cartesian_end()-1);
    }
};

template <class R>
class MidpointCd
{
    typedef typename R::Point_d Point_d;
public:
    typedef Point_d result_type;

    result_type operator()(const Point_d & p, const Point_d & q) const
    {
        return Point_d(p + (q-p)/2);
    }
};

template <class R>
class Center_of_sphereCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::FT FT;
    typedef typename R::LA LA;
    typedef typename LA::Vector Vector;
    typedef typename LA::Matrix Matrix;
public:
    typedef Point_d result_type;

    template <class Forward_iterator>
    result_type operator()(Forward_iterator start, Forward_iterator end) const
    {
        CGAL_assertion(start!=end);
        int d = start->dimension();
        Matrix M(d);
        Vector b(d);
        Point_d pd = *start++;
        for (int i = 0; i < d; ++i) { 
            // we set up the equation for p_i
            Point_d pi = *start++;
            b[i] = 0;
            for (int j = 0; j < d; ++j) {
                M(i,j) = FT(2)*(pi.cartesian(j) - pd.cartesian(j));
                b[i] += (pi.cartesian(j) - pd.cartesian(j)) *
                    (pi.cartesian(j) + pd.cartesian(j));
            }
        }
        FT D;
        Vector x;
        LA::linear_solver(M,b,x,D);
        return Point_d(d, x.begin(), x.end());
    }
}; // Center_of_sphereCd

template <class R>
class Squared_distanceCd
{
    typedef typename R::Point_d     Point;
    typedef typename R::Vector_d    Vector;
    typedef typename R::FT          FT;
public:
    typedef FT result_type;

    result_type operator()(const Point & p, const Point & q) const
    {
        Vector v = p - q;
        return v.squared_length();
    }
};

template <class R>
class Position_on_lineCd
{
    typedef typename R::Point_d Point;
    typedef typename R::LA LA;
    typedef typename R::FT FT;
public:
    typedef typename R::Boolean result_type;

    result_type operator()(const Point & p, const Point & s, const Point & t, 
        FT & l) const
    {
        int d = p.dimension(); 
        CGAL_assertion_msg((d==s.dimension())&&(d==t.dimension()&& d>0), 
                "position_along_line: argument dimensions disagree.");
        CGAL_assertion_msg((s!=t), 
                "Position_on_line_d: line defining points are equal.");
        FT lnum = (p.cartesian(0) - s.cartesian(0)); 
        FT lden = (t.cartesian(0) - s.cartesian(0)); 
        FT num(lnum), den(lden), lnum_i, lden_i;
        for (int i = 1; i < d; i++) {  
            lnum_i = (p.cartesian(i) - s.cartesian(i)); 
            lden_i = (t.cartesian(i) - s.cartesian(i)); 
            if (lnum*lden_i != lnum_i*lden)
                return false; 
            if (lden_i != FT(0)) {
                den = lden_i;
                num = lnum_i;
            }
        }
        l = num / den; return true; 
    }
};

template <class R>
class Barycentric_coordinatesCd
{
    typedef typename R::Point_d Point;
    typedef typename R::LA LA;
    typedef typename R::FT FT;
public:

    template <class ForwardIterator, class OutputIterator>
    OutputIterator operator()(ForwardIterator first, ForwardIterator last, 
        const Point & p, OutputIterator result)
    {
        TUPLE_DIM_CHECK(first,last,Barycentric_coordinates_d);
        //int n = std::distance(first,last); //unused variable
        int d = p.dimension();
        typename R::Affine_rank_d affine_rank;
        CGAL_assertion(affine_rank(first,last)==d);
        std::vector< Point > V(first,last);
        typename LA::Matrix M(d+1,V.size());
        typename LA::Vector b(d+1), x;
        int i;
        for (i=0; i<d; ++i) {
            for (int j=0; j<V.size(); ++j) 
                M(i,j)=V[j].cartesian(i);
            b[i] = p.cartesian(i);
        }
        for (int j=0; j<V.size(); ++j) 
            M(d,j) = 1;
        b[d] = 1;
        FT D;
        LA::linear_solver(M,b,x,D);
        for (i=0; i < x.dimension(); ++result, ++i) {
            *result = x[i];
        }
        return result;
    }
};

template <class R>
class OrientationCd
{
    typedef typename R::Point_d     Point;
    typedef typename R::LA          LA;
    typedef typename R::Orientation Orientation;
public:
    typedef Orientation result_type;

    template <class ForwardIterator>
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        TUPLE_DIM_CHECK(first, last, Orientation_d);
        int d = static_cast<int>(std::distance(first,last)) - 1;
        // range contains d+1 points of dimension d
        CGAL_assertion_msg(first->dimension() == d,
                "Orientation_d: needs first->dimension() + 1 many points.");
        typename LA::Matrix M(d);
        ForwardIterator s = first;
        ++s;
        for( int j = 0; j < d; ++s, ++j )
            for( int i = 0; i < d; ++i )
                M(i,j) = s->cartesian(i) - first->cartesian(i);
        return result_type(LA::sign_of_determinant(M));
    }
};

/* This predicates tests the orientation of (k+1) points that span a
 * k-dimensional affine subspace of the ambiant d-dimensional space. We
 * greedily search for an orthogonal projection on a k-dim axis aligned
 * subspace on which the (full k-dim) predicates answers POSITIVE or NEGATIVE.
 * If no such subspace is found, return COPLANAR.
 * IMPORTANT TODO: Current implementation is VERY bad with filters: if one
 * determinant fails in the filtering step, then all the subsequent ones wil be
 * in exact arithmetic :-(
 * TODO: store the axis-aligned subspace that was found in order to avoid
 * re-searching for it for subsequent calls to operator()
 */
template <class R>
class Coaffine_orientationCd
{ 
    typedef typename R::Point_d     Point_d;
    typedef typename R::LA          LA;
    typedef typename R::Orientation Orientation;
public:
    typedef Orientation result_type;

    // typedef internal::stateful_predicate_tag predicate_category;
    typedef std::vector<int>	Axes;
    struct State
    {
        Axes axes_;
        bool axes_found_;
        State(bool b) : axes_(), axes_found_(b) {}
    };
    mutable State state_;

    Coaffine_orientationCd() : state_(false) {}

    State & state()    {        return state_;    }
    const State & state() const   {        return state_;    }

    template < class ForwardIterator >
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        TUPLE_DIM_CHECK(first,last,Coaffine_orientation_d);
        // |k| is the dimension of the affine subspace
        const int k = std::distance(first,last) - 1;
        // |d| is the dimension of the ambiant space
        const int d = first->dimension();
        CGAL_assertion_msg(k <= d, "Coaffine_orientation_d: needs less that (first->dimension() + 1) points.");
        if( false == state_.axes_found_ )
        {
			state_.axes_.resize(d + 1);
			// We start by choosing the first |k| axes to define a plane of projection
            int i = 0;
            for(; i < k;     ++i) state_.axes_[i] = i;
            for(; i < d + 1; ++i) state_.axes_[i] = -1;
        }
        const typename ForwardIterator::value_type & l(*first);
        typename LA::Matrix M(k); // quadratic
        while( true )
        {
            ForwardIterator s = first;
            ++s;
            int j(0);
            while( s != last )
            {
                const typename ForwardIterator::value_type & point(*s);
                for( int i = 0; i < k; ++i )
                    M(i,j) = point.cartesian(state_.axes_[i]) - l.cartesian(state_.axes_[i]);
                ++s;
                ++j;
            }
            Orientation o = Orientation(LA::sign_of_determinant(M));
            if( ( o != COPLANAR ) || state_.axes_found_ )
            {
                state_.axes_found_ = true;
                return o;
            }
            // for generating all possible unordered k-uple in the range
            // [0 .. d-1]... we go to the next unordered k-uple:
            int index = k - 1;
            while( (index >= 0) && (state_.axes_[index] == d - k + index) )
                --index;
            if( index < 0 )
                break;
            ++state_.axes_[index];
            for( int i = 1; i < k - index; ++i )
                state_.axes_[index + i] = state_.axes_[index] + i;
        }
        return COPLANAR;
    }
};

template <class R>
class Side_of_oriented_sphereCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::LA LA;
    typedef typename R::FT FT;
    typedef typename R::Oriented_side Oriented_side;
public:
    typedef Oriented_side result_type;

    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last, 
         const Point_d& x) const
    { 
        TUPLE_DIM_CHECK(first,last,Side_of_oriented_sphere_d);
        int d = static_cast<int>(std::distance(first,last)); // |A| contains |d| points
        CGAL_assertion_msg((d-1 == first->dimension()), 
                "Side_of_oriented_sphere_d: needs first->dimension()+1 many input points.");
        typename LA::Matrix M(d + 1); 
        for (int i = 0; i < d; ++first, ++i) { 
            FT Sum = 0;
            M(i,0) = 1;
            for (int j = 0; j < d-1; j++) { 
                FT cj = first->cartesian(j);
                M(i,j + 1) = cj; Sum += cj*cj; 
            }
            M(i,d) = Sum; 
        }
        FT Sum = 0; 
        M(d,0) = 1; 
        for (int j = 0; j < d-1; j++) { 
            FT hj = x.cartesian(j);
            M(d,j + 1) = hj; Sum += hj*hj; 
        }
        M(d,d) = Sum;
        return result_type( - LA::sign_of_determinant(M));
    }
};


/* This predicates takes k+1 points defining a k-sphere in d-dim space, and a
 * point |x| (assumed to lie in the same affine subspace spanned by the
 * k-sphere). It tests wether the point |x| lies in the positive or negative
 * side of the k-sphere.
 * The parameter |axis| contains the indices of k axis of the canonical base of
 * R^d, on which the affine subspace projects homeomorphically. We can thus
 * "complete" the k+1 points with d-k other points along the "non-used" axes
 * and then call the usual Side_of_oriented_sphereCd predicate.
 */
template < class R >
class Side_of_oriented_subsphereCd
{
	typedef typename R::Point_d			Point;
	typedef typename R::LA				LA;
	typedef typename R::FT				FT;
	typedef typename R::Orientation     Orientation;
	typedef typename R::Oriented_side   Oriented_side;
	typedef typename R::Side_of_oriented_sphere_d	Side_of_oriented_sphere;
	typedef typename R::Coaffine_orientation_d		Coaffine_orientation;
	typedef typename LA::Matrix	        Matrix;	
	typedef typename Coaffine_orientation::Axes		Axes;
	// DATA MEMBERS
	mutable Coaffine_orientation ori_;
	mutable unsigned int adjust_sign_;
    // a square matrix of size (D+1)x(D+1) where D is the ambient dimension 
	mutable typename LA::Matrix M;
public:
	typedef Oriented_side   result_type;
    // typedef internal::stateless_predicate_tag predicate_category;

    // constructor
	Side_of_oriented_subsphereCd(const R & r = R())
	: ori_(r.coaffine_orientation_d_object()), M(), adjust_sign_(0) { }

	template < class ForwardIterator >
	result_type operator()(ForwardIterator first, ForwardIterator last, const Point & q) const
	{
		const int d = first->dimension();
		const int k = std::distance(first, last) - 1; // dimension of affine subspace
		CGAL_assertion_msg( k <= d, "too much points in range.");
		if( k == d )
		{
			Side_of_oriented_sphere sos;
			return sos(first, last, q); // perhaps slap user on the back of the head here?
		}
        if( M.row_dimension() < d+1 )
            M = Matrix(d+1);
		if( ! ori_.state().axes_found_ )
		{
            // the call to ori_(...) will compute a set of axes to complement our base.
			Orientation o = ori_(first, last);
			if( COPLANAR == o )
            {
                std::cerr << "\nAffine base is flat (it should have positive orientation) !!";
                //return ON_ORIENTED_BOUNDARY;
            }
			CGAL_assertion( o == POSITIVE );
			// Now we can setup the fixed part of the matrix:
			int a(0);
			int j(k);
			typename Axes::iterator axis = ori_.state().axes_.begin();
			while( j < d )
			{
				while( a  == *axis )
				{
					++a; ++axis;
				}
	            adjust_sign_ = ( adjust_sign_ + j + a ) % 2;
				int i(0);
				for( ; i < a; ++i )
					M(i, j) = FT(0);
				M(i++, j) = FT(1);  // i.e.: M(a, j) = 1
				for( ; i < d; ++i )
					M(i, j) = FT(0);
				++j;
				++a;
			}
		}
		typename ForwardIterator::value_type p1 = *first;
		FT SumFirst(0); // squared length of first subsphere point, seen as vector.
		for( int i = 0; i < d; ++i )
		{
			FT ci = p1.cartesian(i);
			SumFirst += ci * ci;
		}
		int j(0); // iterates overs columns/subsphere points
		++first;
		while( first != last )
		{
			typename ForwardIterator::value_type v = *first;
			FT Sum = FT(0);
			for( int i = 0; i < d; ++i )
			{
				FT ci = v.cartesian(i);
				M(i, j) = ci - p1.cartesian(i);
				Sum += ci * ci;
			}
			M(d, j) = Sum - SumFirst;
			++first;
			++j;
		}
		int a(0);
		typename Axes::iterator axis = ori_.state().axes_.begin();
		while( j < d )
		{
			while( a  == *axis )
			{
				++a; ++axis;
			}
			M(d, j) = FT(1) + FT(2) * p1.cartesian(a);
			++j;
			++a;
		}
		FT Sum = FT(0);
		for( int i = 0; i < d; ++i )
		{
			FT ci = q.cartesian(i);
			M(i, d) = ci - p1.cartesian(i);
			Sum += ci * ci;
		}
		M(d, d) = Sum - SumFirst;
		if( 0 == ( adjust_sign_ % 2 ) )
			return result_type( - LA::sign_of_determinant( M ) );
		else
			return result_type(   LA::sign_of_determinant( M ) );
	}
};

template <class R>
class Side_of_bounded_sphereCd
{
    typedef typename R::Point_d         Point_d;
    typedef typename R::Orientation_d   Orientation_d;
    typedef typename R::Side_of_oriented_sphere_d Side_of_oriented_sphere_d;
    typedef typename R::Orientation     Orientation;
    typedef typename R::Oriented_side   Oriented_side;
    typedef typename R::Bounded_side    Bounded_side;
public:
    typedef Bounded_side    result_type;

    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last, 
            const Point_d& p) const
    {
        TUPLE_DIM_CHECK(first,last,region_of_sphere);
        Orientation_d _orientation;
        Orientation o = _orientation(first,last);
        CGAL_assertion_msg((o != 0), "Side_of_bounded_sphere_d: \
                A must be full dimensional.");
        Side_of_oriented_sphere_d _side_of_oriented_sphere;
        Oriented_side oside = _side_of_oriented_sphere(first,last,p);
        if (o == POSITIVE) {
            switch (oside) {
                case ON_POSITIVE_SIDE    :   return ON_BOUNDED_SIDE;
                case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
                case ON_NEGATIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
            }       
        } else {
            switch (oside) {
                case ON_POSITIVE_SIDE    :   return ON_UNBOUNDED_SIDE;
                case ON_ORIENTED_BOUNDARY:   return ON_BOUNDARY;
                case ON_NEGATIVE_SIDE    :   return ON_BOUNDED_SIDE;
            }     
        }
        return ON_BOUNDARY; // never reached
    }
};


template <class R>
class Contained_in_simplexCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::FT FT;
    typedef typename R::LA LA;
    typedef typename LA::Vector Vector;
    typedef typename LA::Matrix Matrix;
public:
    typedef typename R::Boolean result_type;

    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last,
            const Point_d& p) const
    {
        TUPLE_DIM_CHECK(first,last,Contained_in_simplex_d);
        int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
        int d = first->dimension(); 
        CGAL_assertion_code( typename R::Affinely_independent_d check_independence; )
        CGAL_assertion_msg(check_independence(first,last),
                "Contained_in_simplex_d: A not affinely independent.");
        CGAL_assertion(d==p.dimension());

        Matrix M(d + 1,k); 
        Vector b(d +1);
        for (int j = 0; j < k; ++first, ++j) {
            for (int i = 0; i < d; ++i) 
                M(i,j) = first->cartesian(i);
            M(d,j) = 1;
        }
        for (int i = 0; i < d; ++i) 
            b[i] = p.cartesian(i);
        b[d] = 1;

        FT D; 
        Vector lambda; 
        if ( LA::linear_solver(M,b,lambda,D) ) {
            for (int j = 0; j < k; j++) { 
                if (lambda[j] < FT(0)) return false;
            }
            return true;
        }
        return false; 
    }
};

template <class R>
class Contained_in_affine_hullCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::LA LA;
    typedef typename R::Boolean Boolean;
public:
    typedef Boolean result_type;

    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last,
                const Point_d& p) const
    {
         TUPLE_DIM_CHECK(first,last,Contained_in_affine_hullCd);
         int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
         int d = first->dimension(); 
         typename LA::Matrix M(d + 1,k); 
         typename LA::Vector b(d + 1); 
         for (int j = 0; j < k; ++first, ++j) {
             for (int i = 0; i < d; ++i) 
                 M(i,j) = first->cartesian(i);
             M(d,j) = 1;
         }
         for (int i = 0; i < d; ++i)
             b[i] = p.cartesian(i);
         b[d] = 1;
         return LA::is_solvable(M,b);
    }
};


template <class R>
class Affine_rankCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::Vector_d Vector_d;
    typedef typename R::LA LA;
public:
    typedef int result_type;
    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        TUPLE_DIM_CHECK(first,last,Affine_rank_d);
        int k = static_cast<int>(std::distance(first,last)); // |A| contains |k| points
        if (k == 0) return -1;
        if (k == 1) return 0; 
        int d = first->dimension();
        typename LA::Matrix M(d,--k);
        Point_d p0 = *first; ++first; // first points to second
        for (int j = 0; j < k; ++first, ++j) {
            Vector_d v = *first - p0;
            for (int i = 0; i < d; i++) 
                M(i,j) = v.cartesian(i); 
        }
        return LA::rank(M);
    }
};

template <class R>
class Affinely_independentCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::LA LA;
    typedef typename R::Affine_rank_d Affine_rank_d;
public:
    typedef typename R::Boolean result_type;

    template <class ForwardIterator> 
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        Affine_rank_d rank;
        int n = static_cast<int>(std::distance(first,last));
        return rank(first,last) == n-1;
    }
};


template <class R>
class Compare_lexicographicallyCd
{
    typedef typename R::Point_d Point_d;
    typedef typename R::Point_d PointD; //MSVC hack
    typedef typename R::Comparison_result Comparison_result;
public:
    typedef Comparison_result result_type;

    result_type operator()(const Point_d & p1, const Point_d & p2) const
    {
        return PointD::cmp(p1,p2);
    }
};

template <class R>
class Contained_in_linear_hullCd
{
    typedef typename R::LA LA;
    typedef typename R::FT FT;
    typedef typename R::Vector_d Vector_d;
public:
    typedef typename R::Boolean result_type;

    template<class ForwardIterator>
    result_type operator()(
        ForwardIterator first, ForwardIterator last, const Vector_d& x) const
    {
        TUPLE_DIM_CHECK(first,last,Contained_in_linear_hull_d);
        int k = static_cast<int>(std::distance(first,last));
        // |A| contains |k| vectors
        int d = first->dimension();
        typename LA::Matrix M(d,k);
        typename LA::Vector b(d); 
        for (int i = 0; i < d; i++) { 
            b[i] = x.cartesian(i); 
            for (int j = 0; j < k; j++) 
                M(i,j) = (first+j)->cartesian(i); 
        }
        return LA::is_solvable(M,b); 
    }
};

template <class R>
class Linear_rankCd
{
    typedef typename R::LA LA;
public:
    typedef int result_type;

    template <class ForwardIterator>
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        TUPLE_DIM_CHECK(first,last,linear_rank);
        int k = static_cast<int>(std::distance(first,last)); // k vectors
        int d = first->dimension(); 
        typename LA::Matrix M(d,k);
        for (int i = 0; i < d  ; i++)
            for (int j = 0; j < k; j++)  
                M(i,j) = (first + j)->cartesian(i);
        return LA::rank(M);
    }
};

template <class R>
class Linearly_independentCd
{
    typedef typename R::LA LA;
public:
    typedef typename R::Boolean result_type;

    template <class ForwardIterator>
    result_type operator()(ForwardIterator first, ForwardIterator last) const
    {
        typename R::Linear_rank_d rank;
        return rank(first,last) == static_cast<int>(std::distance(first,last));
    }
};

template <class R>
class Linear_baseCd
{
    typedef typename R::LA LA;
    typedef typename R::FT FT;
    typedef typename R::Vector_d Vector_d;
public:
    template <class ForwardIterator, class OutputIterator>
    OutputIterator operator()(ForwardIterator first, ForwardIterator last,
        OutputIterator result) const
    {
        TUPLE_DIM_CHECK(first,last,linear_base);
        int k = static_cast<int>(std::distance(first,last)); // k vectors
        int d = first->dimension();
        typename LA::Matrix M(d,k); 
        for (int j = 0; j < k; ++first, ++j)
            for (int i = 0; i < d; i++)
                M(i,j) = first->cartesian(i);

        std::vector<int> indcols;
        int r = LA::independent_columns(M,indcols);

        for (int l=0; l < r; l++) {
            typename LA::Vector v = M.column(indcols[l]);
            *result++ = Vector_d(d,v.begin(),v.end());
        }
        return result;
    }
};

} //namespace CGAL
#endif //CGAL_FUNCTION_OBJECTSCD_H
