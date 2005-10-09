#ifndef CGAL_INSTANTANEOUS_ADAPTOR_H
#define CGAL_INSTANTANEOUS_ADAPTOR_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE;

#define compute_static(k) ik_.static_object(k)

//! An object to help convert between moving objects and their static representations to wrap a predicate.
/*!  A ref counted pointer is stored to the part of the
  Instantaneous_kernel which actually performs the conversions.

  Look at the source to figure out how it works. 
*/
template <class Predicate, class IK_handle, class Object>
class Instantaneous_adaptor{
public:
  Instantaneous_adaptor(const IK_handle &ik,
			Predicate pred=Predicate()): ik_(ik), pred_(pred){
  }


  typedef typename Predicate::result_type result_type;
  typedef Object argument_type;
  typedef argument_type first_argument_type;
  typedef argument_type second_argument_type;
  typedef argument_type third_argument_type;
  typedef argument_type fourth_argument_type;
  typedef argument_type fifth_argument_type;



  result_type operator()(const first_argument_type &arg0) const {
    return pred_(compute_static(arg0));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1) const {
    /*std::cout << "Args " << compute_static(arg0) <<", " 
	      << compute_static(arg1) << " result " << pred_(compute_static(arg0), compute_static(arg1)) 
	      << " antiresult " << pred_(compute_static(arg1), compute_static(arg0)) << std::endl;*/

    return pred_(compute_static(arg0), compute_static(arg1));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2) const {
    return pred_(compute_static(arg0), compute_static(arg1), compute_static(arg2));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3) const {
    return pred_(compute_static(arg0), compute_static(arg1), compute_static(arg2), compute_static(arg3));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3,
			 const fifth_argument_type &arg4) const {
    return pred_(compute_static(arg0), compute_static(arg1), compute_static(arg2),
		 compute_static(arg3), compute_static(arg4));
  }

protected:
  
  const IK_handle ik_;
  Predicate pred_;
};

#undef compute_static

CGAL_KDS_END_NAMESPACE;
#endif

