#ifndef CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H
#define CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Root_stack_traits_base.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Polynomial>
class Root_stack_default_traits: public internal::Root_stack_traits_base<Polynomial>{
private:
  typedef internal::Root_stack_traits_base<Polynomial> Base;

public:
};

CGAL_POLYNOMIAL_END_NAMESPACE

#endif
