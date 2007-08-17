// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_2_H

#include <CGAL/basic.h>
#include <CGAL/Algebraic_kernel_1.h>

#include <CGAL/Algebraic_kernel_d/Real_embeddable_extension.h>

#include <CGAL/Algebraic_curve_kernel_2/Xy_coordinate_2.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_vertical_line_1.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_analysis_2.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_pair_vertical_line_1.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_pair_analysis_2.h>

#include <CGAL/Algebraic_curve_kernel_2/LRU_hashed_map.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template < class AlgebraicCurvePair_2, class AlgebraicKernel_1 >
class Algebraic_curve_kernel_2 {

// for each functor defines a member function returning an instance of this
// functor
#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y(); }

private:
    //! \name wrapping types
    //!@{

    //! type of an internal curve pair
    typedef AlgebraicCurvePair_2 Internal_curve_pair_2;

    //! type of an internal curve
    typedef typename AlgebraicCurvePair_2::Algebraic_curve_2 Internal_curve_2;

    //! type of internal x_coordinate
    typedef typename Internal_curve_2::X_coordinate Internal_x_coordinate;
    
    //! type of internal coefficient
    typedef typename Internal_curve_2::Coefficient Internal_coefficient;

    //! type of internal polynomial
    typedef typename Internal_curve_2::Poly_d Internal_polynomial_2;
    
    typedef typename NiX::Polynomial_traits<Internal_polynomial_2>::
        Innermost_coefficient Innermost_coefficient;
   
    //!@}
public:
    //! \name types and functors for \c Curved_kernel_2< >
    //!@{
    
    //! type of 1D algebraic kernel 
    typedef AlgebraicKernel_1 Algebraic_kernel_1;

    //! myself
    typedef Algebraic_curve_kernel_2<AlgebraicCurvePair_2, AlgebraicKernel_1>
        Self;
    
    //! type of coefficient
    typedef Internal_coefficient Coefficient;

    //! type of curve pair
    typedef Internal_curve_pair_2 Curve_pair_2;

    //! type of single curve
    typedef Internal_curve_2 Curve_2;

    //! type of x-coordinate
    typedef Internal_x_coordinate X_coordinate_1;
        
    //! new CGAL univariate polynomial type
    typedef ::CGAL::Polynomial<Innermost_coefficient> Polynomial_1;
    //! new CGAL bivariate polynomial type
    typedef ::CGAL::Polynomial<Polynomial_1> Polynomial_2;
    //! bivariate polynomial traits
    typedef ::CGAL::Polynomial_traits_d< Polynomial_2 > Polynomial_traits_2;
    
    //!@}
private:
    //! \name private functors
    //!@{
    
    //! temporary functor providing conversion from \c Poly_in type to
    //! \c Poly_out type, required for NumeriX <-> CGAL polynomial type
    //! conversion
    template <class Poly_2_from, class Poly_2_to>
    struct Polynomial_converter 
    {
        typedef typename Poly_2_from::NT Poly_1_from;
        typedef typename Poly_2_to::NT Poly_1_to;
        // needed for make_transform_iterator
        typedef Poly_1_to result_type;
        
        Poly_2_to operator()(const Poly_2_from& p) const
        {
            return Poly_2_to(
                ::boost::make_transform_iterator(p.begin(), *this),
                ::boost::make_transform_iterator(p.end(), *this));
        }
        Poly_1_to operator()(const Poly_1_from& p) const
        {
            return Poly_1_to(p.begin(), p.end());
        }
    };
    
    //! polynomial canonicalizer: temporarily we use NiX functors since
    //! \c Poly is NiX-type polynomial
    template <class Poly> 
    struct Poly_canonicalizer : public Unary_function< Poly, Poly >
    {
        Poly operator()(Poly p) 
        {
            typedef NiX::Scalar_factor_traits<Poly> Sf_traits;
            typedef typename Sf_traits::Scalar Scalar;
            typename Sf_traits::Scalar_factor scalar_factor;
            typename Sf_traits::Scalar_div scalar_div;
            Scalar g = scalar_factor(p);
            CGAL_assertion(g != Scalar(0));
            if(g != Scalar(1)) 
                scalar_div(p,g);
            if(p.lcoeff().lcoeff() < 0) 
                scalar_div(p,Scalar(-1));
            return p;        
        }
           
    };
    
    //! \brief a simple bivariate polynomial hasher
    //!
    //! picks up at most 12 non-zero low-degree coefficients, sums or
    //! multiplies them and compute floor_log2_abs out of the result
    template <class Poly_2> 
    struct Poly_hasher_2 : public Unary_function< Poly_2, size_t >
    {
        typedef typename Poly_2::NT Poly_1;
        typedef typename Poly_1::NT NT;
        
        size_t operator()(const Poly_2& p) const
        {
            int n_count = 12;  
            NT res(1), sum;
            typename Poly_2::const_iterator it_2 = p.begin();
            while(it_2 != p.end() && n_count > 0) {
                typename Poly_1::const_iterator it_1 = it_2->begin();
                sum = NT(0);
                while(it_1 != it_2->end() && n_count > 0) {
                    NT tmp = *it_1++;
                    if(tmp == NT(0))
                        continue;
                    sum += tmp;
                    n_count--;
                }
                if(sum != 0)
                    res *= sum;
                it_2++;   
            }
            return static_cast<size_t>(CGALi::floor_log2_abs<NT>(res));
        }
    };
    
    //! \brief a simple poly-pair hasher
    template <class Poly_2> 
    struct Poly_pair_hasher_2 
    {
        typedef std::pair<Poly_2, Poly_2> Poly_pair;
        typedef Poly_pair argument_type;
        typedef size_t result_type;
        
        size_t operator()(const Poly_pair& p) const
        {
            Poly_hasher_2<Poly_2> hasher;
            return hasher(p.first);
        }  
    };
    
    //! polynomial pair canonicalizer
    template <class Poly> 
    struct Poly_pair_canonicalizer  
    {
        typedef std::pair<Poly, Poly> Poly_pair;
        typedef Poly_pair argument_type;
        typedef Poly_pair result_type;
        
        Poly_pair operator()(const Poly_pair& pair) const
        {
            if(pair.first > pair.second) 
                return std::make_pair(pair.second,pair.first);
            return pair;
        }
    };
    
    //! polynomial pair gcd creator
    template <class Poly> 
    struct Poly_pair_gcd_creator
    {
        typedef std::pair<Poly, Poly> Poly_pair;
        typedef Poly_pair argument_type;
        typedef Poly result_type;
            
        Poly operator()(const Poly_pair& p) const
        {
            return NiX::gcd(p.first,p.second);
        }
    };     
    
    template <class Poly> 
    struct Poly_pair_creator 
    {
        typedef std::pair<Poly, Poly> Poly_pair;
        typedef Poly_pair argument_type;
        typedef Curve_pair_2 result_type;
        
        Curve_pair_2 operator()(const Poly_pair& p) const 
        {
            return Curve_pair_2(p.first, p.second);
        }
    };

    typedef CGAL::Pair_lexicographical_less_than<Internal_polynomial_2,
        Internal_polynomial_2> Poly_pair_compare;
    
    //! type of curve cache
    typedef CGALi::LRU_hashed_map<Internal_polynomial_2, Curve_2,
        Poly_canonicalizer<Internal_polynomial_2>,
        Poly_hasher_2<Internal_polynomial_2> > Curve_cache;
        
    typedef std::pair<Internal_polynomial_2, Internal_polynomial_2>
        Poly_pair_2;
        
    //! type of curve pair cache 
    typedef CGALi::LRU_hashed_map<Poly_pair_2, Internal_polynomial_2,
        Poly_pair_canonicalizer<Internal_polynomial_2>, 
        Poly_pair_hasher_2<Internal_polynomial_2>, 
        Poly_pair_creator<Internal_polynomial_2>,
        Poly_pair_compare> Curve_pair_cache;
    
    //! an instance of curve cache
    mutable Curve_cache _m_curve_cache;
    
    //! an instance of curve pair cache
    mutable Curve_pair_cache _m_curve_pair_cache;
    //!@}
public:
    //! \name public functors and predicates
    //!@{
       
    //! NumeriX to CGAL polynomial type conversion
    typedef Polynomial_converter<Internal_polynomial_2, Polynomial_2>
                NiX2CGAL_converter;
    //! CGAL to NumeriX polynomial type conversion
    typedef Polynomial_converter<Polynomial_2, Internal_polynomial_2>
                CGAL2NiX_converter;
                
    //! \brief default constructor
    Algebraic_curve_kernel_2() //: _m_curve_cache() 
    {  }
    
    //! \brief constructs \c Curve_2 object, uses caching if appropriate
    struct Construct_curve_2 :
            public Unary_function< Internal_polynomial_2, Curve_2 >
    {
        //! \brief constructs an object from \c Algebraic_curve_kernel_2 type
        ////! no default constructor provided
        Construct_curve_2(Self *pkernel_2) :
             _m_pkernel_2(pkernel_2)
         {  }
     
        Curve_2 operator()(const Internal_polynomial_2& f) const
        {
            Curve_2 cc = _m_pkernel_2->get_curve_cache()(f);
              //Self::get_curve_cache()(f);
            /*std::cout << ";hasher: " << 
                Poly_hasher_2<Internal_polynomial_2>()(f) << "\n";  */
             return cc;
        }
        Curve_2 operator()(const Polynomial_2& f) const
        {
            CGAL2NiX_converter cvt;
            return _m_pkernel_2->get_curve_cache()(cvt(f));
            //Self::get_curve_cache()(cvt(f));
        }
        
    private:
        //! \c Algebraic_curve_kernel_2 pointer 
        Self *_m_pkernel_2; 
    };
        
    Construct_curve_2 construct_curve_2_object() 
    {
        return Construct_curve_2(this);
    }
    
    //! access to curve cache
    Curve_cache& get_curve_cache() const
    {
        //static Curve_cache _m_curve_cache;
        return _m_curve_cache;
    }
    
    //! access to curve pair cache
    Curve_cache& get_curve_pair_cache() const
    {
        //static Curve_cache _m_curve_cache;
        return _m_curve_pair_cache;
    }
    
    //! type of a curve point 
    typedef CGALi::Xy_coordinate_2<Self> Xy_coordinate_2;
    
    //! comparison of x-coordinates 
    struct Compare_x_2 :
         public Binary_function<X_coordinate_1, X_coordinate_1, 
                Comparison_result >
    {
        Comparison_result operator()(const X_coordinate_1& x1, 
                                         const X_coordinate_1& x2) const {
        // not yet implemented in Algebraic_kernel_1
//             Algebraic_kernel_1::C ak;
//             return (ak.compare_x_2_object()(x1, x2));
            return x1.compare(x2);
        }
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                         const Xy_coordinate_2& xy2) const {
            return ((*this)(xy1.x(), xy2.x()));
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_x_2, compare_x_2_object);

    //! comparison of y-coordinates
    struct Compare_y_2 {
        // later!
    };
    CGAL_Algebraic_Kernel_pred(Compare_y_2, compare_y_2_object);
    
    //! lexicographical comparison of two objects of type \c Xy_coordinate_2
    struct Compare_xy_2 :
          public Binary_function<Xy_coordinate_2, Xy_coordinate_2, 
                Comparison_result > 
    {
        Comparison_result operator()(const Xy_coordinate_2& xy1, 
                                          const Xy_coordinate_2& xy2) const {
             return xy1.compare_xy(xy2);
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);
    
    //! \brief checks whether curve has only finitely many self-intersection 
    //! points, i.e., it has no self-overlapped continuous parts
    //!
    //! for algerbaic curves this means that supporting polynomial is 
    //! square-free
    struct Has_finite_number_of_self_intersections_2 :
            public Unary_function< Curve_2, bool > {
        
        bool operator()(const Curve_2& c) const {
            typename Polynomial_traits_2::Is_square_free is_square_free;
            NiX2CGAL_converter cvt;
            Polynomial_2 res = cvt(c.f());
            return is_square_free(res);
        }
    };
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_self_intersections_2, 
            has_finite_number_of_self_intersections_2_object);
            
    //! \brief checks whether a curve pair has finitely many intersections,
    //! in other words, whether two curves have no continuous common part
    //!
    //! in case of algerbaic curves: checks whether supporting polynomials are
    //! coprime
    struct Has_finite_number_of_intersections_2 :
            public Binary_function< Curve_2, Curve_2, bool > { 
               
        bool operator()(const Curve_2& c1, const Curve_2& c2) {
            // if curve ids are the same - non-decomposable
            if(c1.id() == c2.id()) 
                return true;
            typename Polynomial_traits_2::Gcd_up_to_constant_factor gcd_utcf;
            typename Polynomial_traits_2::Total_degree total_degree;
            NiX2CGAL_converter cvt;
            Polynomial_2 p1 = cvt(c1.f()), p2 = cvt(c2.f());
            return (total_degree(gcd_utcf(p1, p2)) == 0);  
        }
    };
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_intersections_2, 
            has_finite_number_of_intersections_2_object);
    
    //! a set of verious curve and curve pair decomposition functions
    struct Decompose_2 {
        
        //! \brief returns a curve without self-overlapping parts 
        //!
        //! in case of algebraic curves computes square-free part of supporting
        //! polynomial
        Curve_2 operator()(const Curve_2& c) {
            typename Polynomial_traits_2::Make_square_free make_square_free;
            NiX2CGAL_converter cvt;
            CGAL2NiX_converter cvt_back;
            // Construct_curve_2_object must be used !!
            return Curve_2(cvt_back(make_square_free(cvt(c.f()))));
        }
        
        //! \brief computes a square-free factorization of a curve \c c, 
        //! retunrs the number of pairwise coprime square-free factors
        //! 
        //! returns square-free pairwise coprime factors in \c fit and
        //! multiplicities in \c mit. Template argument type of \c fit is
        //! \c Curve_2, and \c mit is \c int
        template< class OutputIterator1, class OutputIterator2 >
        int operator()( const Curve_2& c, OutputIterator1 fit, 
                OutputIterator2 mit ) const {
            typename Polynomial_traits_2::
                Square_free_factorization_up_to_constant_factor factorize;
            NiX2CGAL_converter cvt;
            CGAL2NiX_converter cvt_back;
            std::vector<Polynomial_2> factors;
            int n_factors = factorize(cvt(c.f()), std::back_inserter(factors),
                    mit); 
            // Construct_curve_2_object must be used !!
            for(int i = 0; i < (int)factors.size(); i++) {
                *fit++ = Curve_2(cvt_back(factors[i]));
            }
            return n_factors;
        }
        
        //! \brief computes for a given pair of curves \c c1 and \c c2 their 
        //! common  part \c oib and coprime parts \c oi1 and \c oi2 
        //! respectively; returns \c true if the curves were decomposed
        //!
        //! returns true if \c c1 and \c c2 are coprime. Template argument
        //! type of \c oi* is \c Curve_2
        template < class OutputIterator > 
        bool operator()(const Curve_2& c1, const Curve_2& c2,
            OutputIterator oi1, OutputIterator oi2, OutputIterator oib) {
            
            typedef std::vector<Curve_2> Curves;
            Curves parts_f, parts_g;
            if(Curve_2::decompose(c1, c2, 
                std::back_inserter(parts_f), std::back_inserter(parts_g))) {
                // copy the common part returned through both iterators
                // we need to move the common part from oi1/oi2 to oib
                *oib++ = parts_f[0];
                if(parts_f.size() > 1)
                    std::copy(parts_f.begin() + 1, parts_f.end(), oi1);
                if(parts_g.size() > 1)
                    std::copy(parts_g.begin() + 1, parts_g.end(), oi2);
                return true;
            }
            return false;
        }
    };
    CGAL_Algebraic_Kernel_pred(Decompose_2, decompose_2_object);
    
    //!@}
public:
    //! \name types and functors for \c Curved_kernel_2<Algebraic_kernel_2>
    //!@{
    
    typedef Curve_2 Polynomial_2;
    
    typedef Construct_curve_2 Construct_polynomial_2;

    typedef X_coordinate_1 Algebraic_real_1;
    typedef Xy_coordinate_2 Algebraic_real_2;
    
    typedef Has_finite_number_of_self_intersections_2 Is_square_free_2;
    typedef Has_finite_number_of_intersections_2 Is_coprime_2;

    typedef Decompose_2 Make_square_free_2;
    typedef Decompose_2 Square_free_factorization;
    typedef Decompose_2 Make_coprime_2;
    
    struct Derivative_x_2 {
        // later
    };

    struct Derivative_y_2 {
        // later
    };

    struct X_critical_points_2 {
        // later
    };

    struct Y_critical_points_2 {
        // later
     };

    struct Sign_at_2 {
        // can be implemented using Curve_pair_analysis -> later
    };

    struct Solve_2 {
        // can be implemented using Curve_pair_analysis -> later
    };

    // TYPES AND FUNCTORS for Curved_kernel_2< both >
    
    //! type of 1-curve analysis
    typedef CGALi::Curve_analysis_2<Self> Curve_analysis_2; 

    //! type of a curve pair analysis
    typedef CGALi::Curve_pair_analysis_2<Self> Curve_pair_analysis_2; 
    
    //!@}
      
}; // class Algebraic_curve_kernel_2

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_2_H
