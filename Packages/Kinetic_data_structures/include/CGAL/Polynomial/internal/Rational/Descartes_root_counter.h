#ifndef CGAL_POLYNOMIAL_DESCARTES_ROOT_COUNTER_H
#define CGAL_POLYNOMIAL_DESCARTES_ROOT_COUNTER_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Alternation_counter.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE;

template <class Kernel>
class Descartes_root_counter {
public:
  Descartes_root_counter(){}
  Descartes_root_counter(const typename Kernel::Function &f,
			 const Kernel &k): map_(k.map_rational_interval_to_positive_object(f)),
					   kernel_(k){
    
  }

  typedef unsigned int result_type;
  typedef typename Kernel::NT first_argument_type;
  typedef typename Kernel::NT second_argument_type;

  //! Note, the result is an upper bound
  template <class NTT>
  result_type operator()(const NTT &lb, const NTT &ub,
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO,
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO) const {
    typename Kernel::Function mf= map_(lb, ub);
    
    typename POLYNOMIAL_NS::Alternation_counter<first_argument_type> ac;
    for (int i=0; i<= mf.degree(); ++i){
      ac.push_back(mf[i]);
    }
    //std::cout << "Num alternations is " << ac.number_of_alternations() << std::endl;
    return ac.number_of_alternations();
  }
protected:
  typename Kernel::Map_rational_interval_to_positive map_;
  //! What are these?
  Kernel kernel_;
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE;

#endif
