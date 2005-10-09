#ifndef CGAL_KDS_NUMERIC_SOLVER_H
#define CGAL_KDS_NUMERIC_SOLVER_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Numeric_root_stack_core.h>


CGAL_KDS_BEGIN_NAMESPACE



template <class Traits>
struct Derivitive_filter_function_kernel: public Traits {
  
  class Root_stack: public POLYNOMIAL_NS::internal::Numeric_root_stack_core<Traits, true> {
    typedef POLYNOMIAL_NS::internal::Numeric_root_stack_core<Traits, true> Parent;
  public:
    typedef typename Parent::Root Root;
    typedef typename Traits::Function Function;
    Root_stack(const typename Traits::Function &f, 
	       Root lb, Root ub, 
	       const Traits&k): Parent(f, lb, ub, k){
      CGAL_KDS_LOG(LOG_LOTS, "Solved " << f << " from " << lb << " to " << ub << " to get ");
      for (unsigned int i=0; i< Parent::roots_.size(); ++i){
	CGAL_KDS_LOG(LOG_LOTS, Parent::roots_[i] << " ");
      }
      CGAL_KDS_LOG(LOG_LOTS, std::endl);
    }

    Root_stack(){};
  };

  typedef typename Root_stack::Root Root;

  Derivitive_filter_function_kernel(Traits tr): Traits(tr){}
  Derivitive_filter_function_kernel(){}

  Root_stack root_stack_object(const typename Traits::Function &f,
			       const Root &lb,
			       const Root &ub) {
    return Root_stack(f, lb, ub, *this);
  }

};

CGAL_KDS_END_NAMESPACE




#endif // inclusion guard
