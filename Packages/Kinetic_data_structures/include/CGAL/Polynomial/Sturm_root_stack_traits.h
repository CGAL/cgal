#ifndef CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H
#define CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/internal/Rational/Euclidean_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Monic_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Primitive_part_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Reduced_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Subresultant_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Sturm_isolating_interval.h>

//#include <CGAL/Polynomial/Tools/Hybrid_isolating_interval.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template<class Polynomial, class M = Field_tag>
class Sturm_root_stack_traits
  : public Root_stack_default_traits<Polynomial>
{
private:
  typedef Root_stack_default_traits<Polynomial>  Base;
  typedef Sturm_root_stack_traits<Polynomial>    Self;

public:
  typedef internal::Sturm_root_counter<Self> Root_count;
  typedef typename Base::Function            Function;

  //! The the sturm sequence
  //  typedef internal::Sturm_sequence<Self> Sturm_sequence;
  //  typedef internal::Euclidean_Sturm_sequence<Self> Sturm_sequence;
  typedef internal::Monic_Sturm_sequence<Self> Sturm_sequence;
  //  typedef internal::Primitive_part_Sturm_sequence<Self> Sturm_sequence;
  //  typedef internal::Reduced_Sturm_sequence<Self> Sturm_sequence;
  //  typedef internal::Subresultant_Sturm_sequence<Self> Sturm_sequence;
  Sturm_sequence Sturm_sequence_object(const Function &f,
				       const Function &g) const {
    return Sturm_sequence(f, g, *this);
  }


  //! The the standard sequence
  typedef internal::Standard_sequence<Sturm_sequence> Standard_sequence;
  Standard_sequence standard_sequence_object(const Function &f) const {
    return Standard_sequence(f, *this);
  }

  //! The the sign sturm sequence
  typedef internal::Sign_Sturm_sequence<Sturm_sequence> Sign_Sturm_sequence;
  Sign_Sturm_sequence sign_Sturm_sequence_object(const Function &f, const
						 Function &g) const {
    return Sign_Sturm_sequence(f, g, *this);
  }


  //! A bound on the size of roots
  typedef internal::Root_bound_evaluator<Function,M> Root_bound;
  Root_bound root_bound_object(bool b = true, const M& m = M()) const {
    return Root_bound(b, m);
  }

  //  typedef
  //  POLYNOMIAL_NS::Hybrid_isolating_interval<typename Polynomial::NT>
  //  Isolating_interval;
  typedef CGAL_POLYNOMIAL_NS::internal::Isolating_interval<typename Polynomial::NT>
  Isolating_interval;
};


CGAL_POLYNOMIAL_END_NAMESPACE

#endif // CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H
