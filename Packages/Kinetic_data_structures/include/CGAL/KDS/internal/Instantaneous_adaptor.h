#ifndef CGAL_INSTANTANEOUS_ADAPTOR_H
#define CGAL_INSTANTANEOUS_ADAPTOR_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE;


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
    return pred_(ik_.static_object(arg0));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1) const {
    /*std::cout << "Args " << ik_.static_object(arg0) <<", " 
	      << ik_.static_object(arg1) << " result " << pred_(ik_.static_object(arg0), ik_.static_object(arg1)) 
	      << " antiresult " << pred_(ik_.static_object(arg1), ik_.static_object(arg0)) << std::endl;*/

    return pred_(ik_.static_object(arg0), ik_.static_object(arg1));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2) const {
    return pred_(ik_.static_object(arg0), ik_.static_object(arg1), 
		 ik_.static_object(arg2));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3) const {
    return pred_(ik_.static_object(arg0), ik_.static_object(arg1), 
		 ik_.static_object(arg2), ik_.static_object(arg3));
  }

  result_type operator()(const first_argument_type &arg0, 
			 const second_argument_type &arg1,
			 const third_argument_type &arg2,
			 const fourth_argument_type &arg3,
			 const fifth_argument_type &arg4) const {
    return pred_(ik_.static_object(arg0), ik_.static_object(arg1),
		 ik_.static_object(arg2), ik_.static_object(arg3),
		 ik_.static_object(arg4));
  }

protected:
  
  const IK_handle ik_;
  Predicate pred_;
};


CGAL_KDS_END_NAMESPACE;
#endif

