// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//                 Sebastian Limbach <slimbach@mpi-inf.mpg.de>
//                 Michael Kerber    <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_D_1_H
#define CGAL_ALGEBRAIC_KERNEL_D_1_H

#include <CGAL/disable_warnings.h>

#ifndef CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#define CGAL_AK_ENABLE_DEPRECATED_INTERFACE 0
#endif

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/Polynomial.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_real_d_1.h>
#include <CGAL/Algebraic_kernel_d/Descartes.h>
#include <CGAL/Algebraic_kernel_d/Real_roots.h>
#include <CGAL/Algebraic_kernel_d/refine_zero_against.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_1.h>
#include <CGAL/Algebraic_kernel_d/bound_between_1.h>
#include <CGAL/ipower.h>

namespace CGAL {

namespace internal {
template< class AlgebraicReal1, class Isolator_ >
class Algebraic_kernel_d_1_base {

public:
  typedef AlgebraicReal1                              Algebraic_real_1;
  typedef Isolator_                                   Isolator;

  typedef typename Algebraic_real_1::Coefficient      Coefficient;
  typedef typename Algebraic_real_1::Bound            Bound;
  typedef typename Algebraic_real_1::Polynomial_1     Polynomial_1;

  // TODO: Other choice?
  typedef int size_type;
  typedef int Multiplicity_type;

private:
  typedef CGAL::Polynomial_traits_d< Polynomial_1 >   PT_1;


protected:

  // Some functors used for STL calls
  template<typename A,typename B>
    struct Pair_first : public CGAL::cpp98::unary_function<std::pair<A,B>,A> {
      A operator() (std::pair<A,B> pair) const { return pair.first; }
  };

  template<typename A,typename B>
    struct Pair_second : public CGAL::cpp98::unary_function<std::pair<A,B>,B> {
      B operator() (std::pair<A,B> pair) const { return pair.second; }
  };

public:
  class Algebraic_real_traits {
  public:
    typedef Algebraic_real_1                      Type;

    struct Bound_between
      : public CGAL::cpp98::binary_function< Type, Type, Bound > {
      Bound operator()( const Type& t1,
          const Type& t2 ) const {
#if CGAL_AK_DONT_USE_SIMPLE_BOUND_BETWEEN
#warning uses deprecated bound_between_1 functor
        return t1.rational_between( t2 );
#else
        return internal::simple_bound_between(t1,t2);
#endif
      }
    };

    struct Lower_bound
      : public CGAL::cpp98::unary_function< Type, Bound > {
      Bound operator()( const Type& t ) const {
        return t.low();
      }
    };

    struct Upper_bound
      : public CGAL::cpp98::unary_function< Type, Bound > {
      Bound operator()( const Type& t ) const {
        return t.high();
      }
    };

    struct Refine
      : public CGAL::cpp98::unary_function< Type, void > {
      void operator()( const Type& t ) const {
        t.refine();
      }

      void operator()( Type& t, int rel_prec ) const {
        // If t is zero, we can refine the interval to
        //  infinite precission
        if( CGAL::is_zero( t ) ) {
          t = Type(0);
        } else {
          // Refine until both boundaries have the same sign
          while( CGAL::sign( t.high() ) !=
              CGAL::sign( t.low() ) )
            t.refine();

          CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
              CGAL::sign( t.low() ) != CGAL::ZERO );

          // Calculate the needed precision
          Bound prec = Bound(1) /
            CGAL::ipower( Bound(2), rel_prec );

          // Refine until precision is reached
          while( CGAL::abs( t.high() - t.low() ) /
                 (CGAL::max)( CGAL::abs( t.high() ),
                              CGAL::abs( t.low() ) ) > prec ) {
            t.refine();

            CGAL_assertion( CGAL::sign( t.high() ) != CGAL::ZERO &&
                CGAL::sign( t.low() ) != CGAL::ZERO );

          }
        }
      }
    };

    struct Approximate_absolute_1:
      public CGAL::cpp98::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
      std::pair<Bound,Bound>
      operator()(const Algebraic_real_1& x, int prec) const {
        Lower_bound lower;
        Upper_bound upper;
        Refine refine;
        Bound l = lower(x);
        Bound u = upper(x);
        Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
        while((prec>0)?((u-l)*error>Bound(1)):((u-l)>error)){
          refine(x);
          u = upper(x);
          l = lower(x);
        }
        return std::make_pair(l,u);
      }
    };

    struct Approximate_relative_1:
      public CGAL::cpp98::binary_function<Algebraic_real_1,int,std::pair<Bound,Bound> >{
      std::pair<Bound,Bound>
      operator()(const Algebraic_real_1& x, int prec) const {

        if(CGAL::is_zero(x)) return std::make_pair(Bound(0),Bound(0));

        Lower_bound lower;
        Upper_bound upper;
        Refine refine;
        Bound l = lower(x);
        Bound u = upper(x);
        Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
        Bound min_b = (CGAL::min)(CGAL::abs(u),CGAL::abs(l));
        while((prec>0)?((u-l)*error>min_b):((u-l)>error*min_b)){
          refine(x);
          u = upper(x);
          l = lower(x);
          min_b = (CGAL::min)(CGAL::abs(u),CGAL::abs(l));
        }
        return std::make_pair(l,u);
      }
    };

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    typedef Lower_bound Lower_boundary;
    typedef Upper_bound Upper_boundary;
    typedef Bound_between Boundary_between;
#endif


  }; // class Algebraic_real_traits

  struct Construct_algebraic_real_1;

  // Functors of Algebraic_kernel_d_1
  struct Solve_1 {
  public:
    template <class OutputIterator>
    OutputIterator
    operator()(const Polynomial_1& p, OutputIterator oi) const {
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#else
      CGAL_precondition(!CGAL::is_zero(p));
#endif
      internal::Real_roots< Algebraic_real_1, Isolator > real_roots;
      std::list< int > mults;
      std::list< Algebraic_real_1 > roots;
      real_roots( p, std::back_inserter(roots), std::back_inserter( mults ) );
      CGAL_assertion(roots.size()==mults.size());
      std::list<int>::iterator mit =mults.begin();
      typename std::list< Algebraic_real_1 >::iterator rit = roots.begin();
      while(rit != roots.end()) {
        //*oi++ = std::make_pair(*rit, (unsigned int)(*mit));
        *oi++ = std::make_pair(*rit, *mit);
        rit++;
        mit++;
      }
      return oi;
    }

#if 1 || CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    template< class OutputIterator >
    OutputIterator operator()(
        const Polynomial_1& p,
        OutputIterator oi ,
        bool known_to_be_square_free) const {
      return this->operator()(p,known_to_be_square_free,oi);
    }
#endif

    template< class OutputIterator >
    OutputIterator operator()(
        const Polynomial_1& p,
        bool known_to_be_square_free,
        OutputIterator oi) const {

      internal::Real_roots< Algebraic_real_1, Isolator > real_roots;
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
#else
      CGAL_precondition(!CGAL::is_zero(p));
#endif
      std::list<Algebraic_real_1> roots;
      if( known_to_be_square_free ){
        real_roots(p,std::back_inserter(roots));
      }else{
        std::list<int> dummy;
        real_roots(p,std::back_inserter(roots),std::back_inserter(dummy));
      }
      return std::copy(roots.begin(),roots.end(),oi);
    }

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    template< class OutputIteratorRoots , class OutputIteratorMults >
    std::pair<OutputIteratorRoots,OutputIteratorMults>
    operator()(
        const Polynomial_1& p,
        OutputIteratorRoots roi,
        OutputIteratorMults moi) const {

      internal::Real_roots< Algebraic_real_1, Isolator > real_roots;
      real_roots(p,roi,moi);
      return std::make_pair(roi,moi);
    }
#endif

    protected:

    /*
    // TODO: Can we avoid to use this?
    struct Greater_compare :
      public CGAL::cpp98::binary_function<Algebraic_real_1,Algebraic_real_1,bool> {

      bool operator() (const Algebraic_real_1& a, const Algebraic_real_1& b)
        const {
        return a>b;
      }

    };
    */

    public:

    template< class OutputIterator >
    OutputIterator operator()(const Polynomial_1& p, Bound l, Bound u,
                              OutputIterator res) const {

      std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > roots;
      this->operator() (p,std::back_inserter(roots));
      Algebraic_real_1 alg_l=Construct_algebraic_real_1()(l);
      Algebraic_real_1 alg_u=Construct_algebraic_real_1()(u);
      typedef typename
        std::vector<std::pair<Algebraic_real_1,Multiplicity_type> >::iterator
        Iterator;
      Pair_first<Algebraic_real_1,Multiplicity_type> pair_first;
      Iterator it_start=std::lower_bound
        (::boost::make_transform_iterator(roots.begin(),pair_first),
         ::boost::make_transform_iterator(roots.end(),pair_first),
         alg_l).base();
      Iterator it_end=std::upper_bound
        (::boost::make_transform_iterator(it_start,pair_first),
         ::boost::make_transform_iterator(roots.end(),pair_first),
         alg_u).base();
      std::copy(it_start,it_end,res);
      return res;
    }

    template< class OutputIterator >
    OutputIterator operator()(const Polynomial_1& p,
                              bool known_to_be_square_free,
                              Bound l, Bound u,
                              OutputIterator res) const {

      std::vector<Algebraic_real_1 > roots;
      this->operator() (p,known_to_be_square_free,std::back_inserter(roots));
      Algebraic_real_1 alg_l=Construct_algebraic_real_1()(l);
      Algebraic_real_1 alg_u=Construct_algebraic_real_1()(u);
      typedef typename
        std::vector<Algebraic_real_1>::iterator
        Iterator;
      Iterator it_start=std::lower_bound(roots.begin(),roots.end(),alg_l);
      Iterator it_end=std::upper_bound(it_start,roots.end(),alg_u);
      std::copy(it_start,it_end,res);
      return res;
    }

  };

  class Number_of_solutions_1
      : public CGAL::cpp98::unary_function<Polynomial_1,size_type> {

    public:

      size_type operator()
        (const Polynomial_1& p) const {

        std::vector<std::pair<Algebraic_real_1,Multiplicity_type> > roots;
        Solve_1()(p,std::back_inserter(roots));
        return static_cast<size_type>(roots.size());
      }

  };


  struct Sign_at_1
    : public CGAL::cpp98::binary_function< Polynomial_1, Algebraic_real_1, CGAL::Sign > {
    CGAL::Sign operator()( const Polynomial_1& p, const Algebraic_real_1& ar ) const {
      if(CGAL::is_zero(p)) return ZERO;
      if(CGAL::degree(p)==0) return p.sign_at(0);
      if( ar.low() == ar.high() ) return p.sign_at( ar.low() );

      if (p == ar.polynomial()) {
        return ZERO;
      }

      Polynomial_1 g = gcd_utcf(p,ar.polynomial());
      if (g.sign_at(ar.low()) != g.sign_at(ar.high())) return ZERO;

      while(internal::descartes(p,ar.low(),ar.high()) > 0) ar.refine();
      while( p.sign_at(ar.low())  == ZERO )  ar.refine();
      while( p.sign_at(ar.high()) == ZERO )  ar.refine();

      CGAL::Sign result = p.sign_at(ar.low());
      CGAL_assertion(result == p.sign_at(ar.high()));
      return result;
    }
  };
  struct Is_zero_at_1
    : public CGAL::cpp98::binary_function< Polynomial_1, Algebraic_real_1, bool > {
    bool operator()( const Polynomial_1& p, const Algebraic_real_1& ar ) const {
      if(CGAL::is_zero(p)) return true;
      if( ar.low() == ar.high() ) return p.sign_at( ar.low() ) == ZERO;
      Polynomial_1 g = gcd_utcf(p,ar.polynomial());
      return g.sign_at(ar.low()) != g.sign_at(ar.high());
    }
  };

  struct Is_square_free_1
    : public CGAL::cpp98::unary_function< Polynomial_1, bool > {
    bool operator()( const Polynomial_1& p ) const {
      typename CGAL::Polynomial_traits_d< Polynomial_1 >::Is_square_free isf;
      return isf(p);
    }
  };

  struct Is_coprime_1
    : public CGAL::cpp98::binary_function< Polynomial_1, Polynomial_1, bool > {
    bool operator()( const Polynomial_1& p1, const Polynomial_1& p2 ) const {
      typename CGAL::Polynomial_traits_d< Polynomial_1 >::Total_degree total_degree;

      // TODO: Is GCD already filtered?
      return( total_degree( gcd_utcf( p1, p2 ) ) == 0 );
    }
  };

  struct Make_square_free_1
    : public CGAL::cpp98::unary_function< Polynomial_1, Polynomial_1 > {
    Polynomial_1 operator()( const Polynomial_1& p ) const {
      return typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()( p );
    }
  };

  struct Make_coprime_1 {
    typedef bool         result_type;
    typedef Polynomial_1 first_argument_type;
    typedef Polynomial_1 second_argument_type;
    typedef Polynomial_1 third_argument_type;
    typedef Polynomial_1 fourth_argument_type;
    typedef Polynomial_1 fifth_argument_type;

    bool operator()( const Polynomial_1& p1,
        const Polynomial_1& p2,
        Polynomial_1& g, // ggT utcf
        Polynomial_1& q1, // Rest utcf
        Polynomial_1& q2 ) const {
      g = typename CGAL::Polynomial_traits_d< Polynomial_1 >::Gcd_up_to_constant_factor()( p1, p2 );
      q1 = p1 / g;
      q2 = p2 / g;
      return CGAL::is_one(g);
    }
  };


  struct Square_free_factorize_1 {
    template< class OutputIterator>
    OutputIterator operator()( const Polynomial_1& p, OutputIterator it) const {
      typename PT_1::Square_free_factorize_up_to_constant_factor sqff;
      return sqff(p,it);
    }
  };

  struct Compute_polynomial_1 : public CGAL::cpp98::unary_function<Algebraic_real_1,
                                                           Polynomial_1> {
    Polynomial_1 operator()(const Algebraic_real_1& x) const {
      return x.polynomial();
    }
  };

  struct Construct_algebraic_real_1 {

    public:

      typedef Algebraic_real_1 result_type;

      result_type operator() (int a) const {
        return Algebraic_real_1(a);
      }

      result_type operator() (Bound a) const {
        return Algebraic_real_1(a);
      }

      result_type operator()
      (typename CGAL::First_if_different<Coefficient,Bound>::Type a) const {
        Coefficient coeffs[2] = {a,Coefficient(-1)};
        Polynomial_1 p = typename PT_1::Construct_polynomial()
          (coeffs,coeffs+2);
        std::vector<Algebraic_real_1 > roots;
        Solve_1()(p,true,std::back_inserter(roots));
        CGAL_assertion(roots.size() == size_type(1));
        return roots[0];
      }


      result_type operator() (Polynomial_1 p,size_type i)
        const {
        std::vector<Algebraic_real_1 > roots;
        Solve_1()(p,true,std::back_inserter(roots));
        CGAL_assertion( size_type(roots.size()) > i);
        return roots[i];
      }

      result_type operator() (Polynomial_1 p,
                              Bound l, Bound u) const {
        CGAL_precondition(l<u);
        return Algebraic_real_1(p,l,u);
      }

  };

  struct Compare_1
    : public CGAL::cpp98::binary_function<Algebraic_real_1,
                                  Algebraic_real_1,
                                  CGAL::Comparison_result>{

    typedef CGAL::Comparison_result result_type;

    result_type operator() (Algebraic_real_1 a,Algebraic_real_1 b) const {
      return typename Real_embeddable_traits<Algebraic_real_1>
                        ::Compare() (a,b);
    }

    result_type operator() (Algebraic_real_1 a,int b) const {
      return this->operator()(a,Construct_algebraic_real_1()(b));
    }

    result_type operator() (Algebraic_real_1 a,Bound b) const {
      return this->operator()(a,Construct_algebraic_real_1()(b));
    }


    result_type operator()
      (Algebraic_real_1 a,
       typename CGAL::First_if_different<Coefficient,Bound>::Type b) const {
      return this->operator()(a,Construct_algebraic_real_1()(b));
    }

    result_type operator() (int a, Algebraic_real_1 b) const {
      return this->operator()(Construct_algebraic_real_1()(a),b);
    }

    result_type operator() (Bound a,Algebraic_real_1 b) const {
      return this->operator()(Construct_algebraic_real_1()(a),b);
    }


    result_type operator()
      (typename CGAL::First_if_different<Coefficient,Bound>::Type a,
       Algebraic_real_1 b) const {
      return this->operator()(Construct_algebraic_real_1()(a),b);
    }

  };

 public:

  struct Isolate_1 : public CGAL::cpp98::binary_function
    < Algebraic_real_1,Polynomial_1,std::pair<Bound,Bound> > {

    public:

    std::pair<Bound,Bound> operator() (const Algebraic_real_1 a,
                                       const Polynomial_1 p) const {

      if(p == a.polynomial()) return std::make_pair(a.low(),a.high());

      std::vector<Algebraic_real_1> roots;
      // First isolate p...
      Solve_1()(p,false,std::back_inserter(roots));
      typedef typename std::vector<Algebraic_real_1>::iterator Iterator;
      // Binary search on the root to find a place where a could be inserted
      std::pair<Iterator,Iterator> it_pair
        = equal_range(roots.begin(),roots.end(),a);
      CGAL_assertion(std::distance(it_pair.first,it_pair.second)==0 ||
                     std::distance(it_pair.first,it_pair.second)==1);
      // If we can insert a in two places, it must have been in roots already
      bool a_in_roots = std::distance(it_pair.first,it_pair.second)==1;
      if(a_in_roots) {
        // TODO: can we rely on the property that the isolating intervals
        // of the roots in p are isolating from each other. What
        // if p was factorized during isolation? Is that still
        // guaranteed? To be sure, we do it this way:
        if(it_pair.first!=roots.begin()) {
          it_pair.first->strong_refine(*(it_pair.first-1));
        }
        if(it_pair.second!=roots.end()) {
          it_pair.first->strong_refine(*(it_pair.second));
        }
        return std::make_pair(it_pair.first->low(),it_pair.first->high());
      } else {
        // Refine a until disjoint from neighbors
        // This is probably not even necessary since the isolating
        // interval of a isolates against all roots of p thanks to the
        // comparisons. But to be sure...
        if(it_pair.first!=roots.begin()) {
          a.strong_refine(*(it_pair.first-1));
        }
        if(it_pair.first!=roots.end()) {
          a.strong_refine(*(it_pair.first));
        }
        return std::make_pair(a.low(),a.high());
      }
    }

  };

  typedef typename Algebraic_real_traits::Bound_between Bound_between_1;
  typedef typename Algebraic_real_traits::Approximate_absolute_1 Approximate_absolute_1;
  typedef typename Algebraic_real_traits::Approximate_relative_1 Approximate_relative_1;




#define CGAL_ALGEBRAIC_KERNEL_1_PRED(Y,Z) Y Z() const { return Y(); }
#define CGAL_ALGEBRAIC_KERNEL_1_PRED_WITH_KERNEL \
            Y Z() const { return Y((const Algebraic_kernel_d_1*)this); }


  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_square_free_1,
      is_square_free_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_square_free_1,
      make_square_free_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Square_free_factorize_1,
      square_free_factorize_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_coprime_1,
      is_coprime_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Make_coprime_1,
      make_coprime_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Solve_1,
      solve_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Number_of_solutions_1,
      number_of_solutions_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Construct_algebraic_real_1,
      construct_algebraic_real_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Sign_at_1,
      sign_at_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Is_zero_at_1,
      is_zero_at_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Compare_1,compare_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Bound_between_1,
      bound_between_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_absolute_1,
      approximate_absolute_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Approximate_relative_1,
      approximate_relative_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Compute_polynomial_1,
      compute_polynomial_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Isolate_1,
      isolate_1_object);

  // Deprecated
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
  typedef Bound Boundary;
  typedef typename Algebraic_real_traits::Refine Refine_1;
  typedef typename Algebraic_real_traits::Lower_bound Lower_bound_1;
  typedef typename Algebraic_real_traits::Upper_bound Upper_bound_1;
  typedef typename Algebraic_real_traits::Lower_bound Lower_boundary_1;
  typedef typename Algebraic_real_traits::Upper_bound Upper_boundary_1;
  typedef Bound_between_1 Boundary_between_1;

  CGAL_ALGEBRAIC_KERNEL_1_PRED(Refine_1,      refine_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Lower_bound_1, lower_bound_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Upper_bound_1, upper_bound_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Lower_boundary_1, lower_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Upper_boundary_1, upper_boundary_1_object);
  CGAL_ALGEBRAIC_KERNEL_1_PRED(Boundary_between_1, boundary_between_1_object);
#endif

#undef CGAL_ALGEBRAIC_KERNEL_1_PRED

};
} // namespace internal


template< class Coefficient,
          class Bound = typename CGAL::Get_arithmetic_kernel< Coefficient >::Arithmetic_kernel::Rational,
          class RepClass = internal::Algebraic_real_rep< Coefficient, Bound >,
          class Isolator = internal::Descartes< typename CGAL::Polynomial_type_generator<Coefficient,1>::Type, Bound > >
class Algebraic_kernel_d_1
  : public internal::Algebraic_kernel_d_1_base<

    // Template argument #1 (AlgebraicReal1)
        internal::Algebraic_real_d_1<
            Coefficient,
            Bound,
            ::CGAL::Handle_policy_no_union,
            RepClass >,

    // Template argument #2 (Isolator_)
        Isolator >

{};


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ALGEBRAIC_KERNEL_D_1_H
