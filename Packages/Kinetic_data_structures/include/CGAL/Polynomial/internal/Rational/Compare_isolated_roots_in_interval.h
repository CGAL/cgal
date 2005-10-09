#ifndef CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H
#define CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H

#include <CGAL/Polynomial/basic.h>


CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template<class Kernel, class M = CGAL::Ring_tag>
class Compare_isolated_roots_in_interval
{
protected:
  typedef M                        Method_tag;
  typedef typename Kernel::NT  NT;

  typedef typename 
  Kernel::Sign_Sturm_sequence Sign_Sturm_sequence;

public:
  Compare_isolated_roots_in_interval(const typename Kernel::Function& p,
				     const typename Kernel::Function& q,
				     Kernel k): kernel_(k),q(q), seq(kernel_.sign_Sturm_sequence_object(p, q)) {}

  Compare_isolated_roots_in_interval(){}


  typedef POLYNOMIAL_NS::Comparison_result  result_type;
  typedef NT                       argument_type;
  typedef NT                       first_argument_type;
  typedef NT                       second_argument_type;
  
  result_type operator()(const NT& a, const NT& b) const
  {
    int sgn = seq.sum_of_signs(a, b);

    if ( sgn == POLYNOMIAL_NS::ZERO ) { return POLYNOMIAL_NS::EQUAL; }

    int s_a = static_cast<int>(  kernel_.sign_at_object(q)(a) );

    if ( sgn == s_a ) { return POLYNOMIAL_NS::SMALLER; }
    return POLYNOMIAL_NS::LARGER;
  }


protected:
  Kernel kernel_;
  typename Kernel::Function          q;
  Sign_Sturm_sequence seq;
};






CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE


#endif // CGAL_COMPARE_ISOLATED_ROOTS_IN_INTERVAL_H
