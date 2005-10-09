#ifndef CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_KERNEL_H
#define CGAL_POLYNOMIAL_FILTERED_POLYNOMIAL_KERNEL_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Filtered_rational/Filtered_rational_traits.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_sign_at.h>
#include <CGAL/Polynomial/Filtered_kernel/Filtered_root_multiplicity.h>
#include <CGAL/Polynomial/Kernel/Root_container.h>
#include <CGAL/Polynomial/Kernel/Is_even_multiplicity.h>
#include <CGAL/Polynomial/Kernel/Is_rational.h>
#include <CGAL/Polynomial/Kernel/To_rational.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE


//! A filtered polynomial kernel.
/*!  
*/
template < class Traits_t, class RE>
class Filtered_kernel: public internal::Filtered_rational_traits<Traits_t> {
  typedef Filtered_kernel< Traits_t, RE > This;
  typedef internal::Filtered_rational_traits<Traits_t> P;
public:
  typedef RE Root_stack;
  typedef typename RE::Root Root;
  typedef typename P::Function Function;
  typedef typename Function::NT NT;
  typedef typename Root_stack::Traits Root_stack_traits;

  Filtered_kernel(const Root_stack_traits &tr = Root_stack_traits()): ret_(tr){}

  typedef typename internal::Filtered_sign_at<This> Sign_at;
  Sign_at sign_at_object(const Function &p) const {
    return Sign_at(p, *this);
  }
  
  typedef internal::Filtered_root_multiplicity<This> Multiplicity;
  Multiplicity multiplicity_object(const Function &p0) const {
    return Multiplicity(p0, *this);
  }

  typedef internal::Is_even_multiplicity<This> Is_even_multiplicity;
  Is_even_multiplicity is_even_multiplicity_object(const Function &) const {
    return Is_even_multiplicity();
  }

  typedef internal::Is_rational<This> Is_rational;
  Is_rational is_rational_object() const {
    return Is_rational();
  }

  typedef internal::To_rational<This> To_rational;
  To_rational to_rational_object() const {
    return To_rational();
  }

  //! Return a container for roots in an interval
  /*!
    \todo make sure that the iterator has all the right types.
  */
  typedef internal::Root_container<This> Root_container;
  friend class internal::Root_container<This>;
  Root_container root_container_object(const typename P::Function &f, 
				       const Root &lb=-Root::infinity(),
				       const Root &ub= Root::infinity()) const {
    return Root_container(f, lb, ub, *this);
  }

  Root_stack root_stack_object(const typename P::Function &f, 
					const Root &lb=-Root::infinity(),
					const Root &ub= Root::infinity()) const {
    return Root_stack(f, lb, ub, root_stack_traits_object());
  }

  

  Root_stack_traits root_stack_traits_object() const {
    return ret_;
  }
protected:
  Root_stack_traits ret_;
 };

CGAL_POLYNOMIAL_END_NAMESPACE

#endif
