#ifndef CGAL_POLYNOMIAL_FILTERED_DESCARTES_HAS_ROOT_H
#define CGAL_POLYNOMIAL_FILTERED_DESCARTES_HAS_ROOT_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Filtered_rational/Filtered_Descartes_root_counter.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE;





template <class Kernel>
class Filtered_Descartes_has_root {
public:
  Filtered_Descartes_has_root(){}

  Filtered_Descartes_has_root(const typename Kernel::Function &fh, Kernel k= Kernel()): h_(fh), kernel_(k) {
  }

  ~Filtered_Descartes_has_root(){
  }

  typedef bool result_type;

  template <class NTT>
  result_type operator()(const NTT &begin, const NTT &end,
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO, 
			 POLYNOMIAL_NS::Sign=POLYNOMIAL_NS::ZERO) const {
    return filtered_Descartes_root_counter(h_, begin, end, false, kernel_)!= 0;
  }

protected:
  typename Kernel::Function h_;
  Kernel kernel_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE;

#endif
