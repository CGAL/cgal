// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

/*! \file NiX/Bitstream_descartes.h
  \brief Defines class NiX::Bitstream_descartes. 
  
  Isolate real roots of polynomials.

  This file provides a class to isolate real roots of polynomials,
  using the algorithm based on the method of Descartes.

  The polynomial has to be a univariat polynomial over any number
  type which is contained in the real numbers.
*/

#ifndef CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_H
#define CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_H

#include <CGAL/basic.h>

#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes_rndl_tree_traits.h>

// TODO: Support this...
/*#ifdef NiX_ISOLATION_TIMER
#include <LiS/Real_timer.h>
#endif*/

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class Polynomial_ , 
           class Rational_, 
           class Tree_ = Bitstream_descartes_rndl_tree
               < Bitstream_descartes_rndl_tree_traits
                   < typename Polynomial_::NT > 
               >
         > 
class Bitstream_descartes{

/*#ifdef NiX_ISOLATION_TIMER
public:
    static LiS::Real_timer timer_1;
    static LiS::Real_timer timer_2;
    static LiS::Real_timer timer_3;
    static LiS::Real_timer timer_4;
#endif*/

public:
    typedef Polynomial_ Polynomial;
    typedef Rational_   Boundary;
    typedef Tree_       Tree;
private:
    typedef typename Polynomial::NT Coefficient; 
 
    typedef typename Get_arithmetic_kernel<Coefficient>::Arithmetic_kernel AT;
    
    typedef typename AT::Integer Integer;
    typedef typename AT::Rational Rational;
    typedef typename AT::Field_with_sqrt FWS;

    //BOOST_STATIC_ASSERT((::boost::is_same<Rational_,Rational>::value));
    
private: // member
    std::vector<Rational> lower_bounds;
    std::vector<Rational> upper_bounds;
    int number_of_real_roots_; 
    Polynomial polynomial_; 

public:
    Bitstream_descartes(): number_of_real_roots_(-1),polynomial_(0){};
    
    Bitstream_descartes(const Polynomial& poly): number_of_real_roots_(-1), polynomial_(poly) {
        if(polynomial_ == Polynomial(0)) return; 
        

        typedef typename Tree::TRAITS Traits;  
        Traits traits(poly);
        typename Traits::Lower_bound_log2_abs lbd = traits.lower_bound_log2_abs_object();
        typename Traits::Upper_bound_log2_abs ubd = traits.upper_bound_log2_abs_object();
    
        Tree tree(
                Fujiwara_root_bound_log(traits.begin(),traits.end(),lbd,ubd)+1,
                traits.begin(), traits.end(),
                typename Tree::Monomial_basis_tag(),
                traits);
          
        typename Tree::Node_iterator it = tree.begin();
        typename Tree::Node_iterator chld_first, chld_beyond; 
        
        // std::cout << " find roots " << std::endl; 
        while (it != tree.end()) {
            if (tree.max_var(it) == 1) {
                ++it;
            } else {
                tree.subdivide(it, chld_first, chld_beyond);
                it = chld_first;
            }
        }
        
        it = tree.begin();
        while (it != tree.end()) {
            lower_bounds.push_back(tree.lower(it));
            upper_bounds.push_back(tree.upper(it));
            // CGAL_postcondition(CGALi::descartes(polynomial_,tree.lower(it),tree.upper(it)) == 1);
            it++;
            }
        number_of_real_roots_ = lower_bounds.size();
    };

public: // functions
    
    /*! \brief returns the defining polynomial*/ 
    Polynomial polynomial() const { return polynomial_; }
    
    //! returns the number of real roots
    int number_of_real_roots() const { return number_of_real_roots_; }

    bool is_exact_root(int i) const { return false; } 
  
    void left_boundary(int i, Integer& numerator_, Integer& denominator_) const {
        typedef CGAL::Fraction_traits<Rational> Fraction_traits; 
        typename Fraction_traits::Decompose decompose;
        decompose(lower_bounds[i],numerator_,denominator_);
    }
    
    void right_boundary(int i,Integer& numerator_, Integer& denominator_) const {
        typedef CGAL::Fraction_traits<Rational> Fraction_traits; 
        typename Fraction_traits::Decompose decompose;
        decompose(upper_bounds[i],numerator_,denominator_);
    }
    
    Rational left_boundary(int i)  const { return lower_bounds[i]; }
    Rational right_boundary(int i) const { return upper_bounds[i]; }    
};

} // namespace CGALi 

CGAL_END_NAMESPACE

/*#ifdef NiX_ISOLATION_TIMER
template <class T, class S> LiS::Real_timer NiX::Bitstream_descartes<T,S>::timer_1;
template <class T, class S> LiS::Real_timer NiX::Bitstream_descartes<T,S>::timer_2;
template <class T, class S> LiS::Real_timer NiX::Bitstream_descartes<T,S>::timer_3;
template <class T, class S> LiS::Real_timer NiX::Bitstream_descartes<T,S>::timer_4;
#endif*/

#endif //  CGAL_ALGEBRAIC_KERNEL_D_BITSTREAM_DESCARTES_H

