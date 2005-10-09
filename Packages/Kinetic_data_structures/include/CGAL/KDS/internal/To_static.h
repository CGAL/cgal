#ifndef CGAL_KDS_TO_STATIC_H
#define CGAL_KDS_TO_STATIC_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Simple_cartesian.h>

CGAL_KDS_BEGIN_NAMESPACE

//! a functor that turns a moving object into a static object. It needs to pick the static type.
/*!  This is the default implementation which looks for a class named
  Static_traits in the moving object type (the first template
  argument) and uses the Static_type field of that traits when
  instantiated with the static kernel to figure out the return
  type. It then calls the to_static member method of the class.
*/
template <class Arg,
	  class SK= CGAL::Simple_cartesian<typename Arg::Coefficient::NT> >
class To_static{
  typedef typename Arg::template Static_traits<SK> Traits;
public:
  //! The way time is represented.
  typedef typename Arg::Coordinate::NT Time;

  To_static(){}

  //! Construct it with a static kernel
  To_static(const SK &sk): sk_(sk){}

  typedef typename Traits::Static_type result_type;
  typedef Arg argument_type;
  //! Convert an appropriate moving object to a static object
  result_type operator()(const argument_type &arg) const {
    return Traits::to_static(arg, time(), sk_);
  }
  //! What this believes the time to be.
  const Time& time() const {
    return t_;
  }
  //! Set the time.
  void set_time(const Time &t) {
    t_=t;
  }
protected:
  SK sk_;
  Time t_;
};
CGAL_KDS_END_NAMESPACE
#endif
