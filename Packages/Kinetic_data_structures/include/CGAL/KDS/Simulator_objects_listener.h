#ifndef CGAL_KDS_NOTIFICATION_HELPERS_H
#define CGAL_KDS_NOTIFICATION_HELPERS_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE



//! A helper function to handle the simulator reversing time.
/*!  This helper is only useful if you are implementing a variant on a
  moving object table. It makes sure that the direct_of_time fields of
  the MOT and the Simulator agree.

  See CGAL::Listener for a description of what the Simulator_listener
  template paramenter should provide.
*/
template <class Simulator_listener, class MOT>
class Simulator_objects_listener: public Simulator_listener {
  typedef Simulator_listener P;
public:
  //! THe only constructor
  Simulator_objects_listener(typename Simulator_listener::Notifier_pointer sim, 
			      MOT *kds): Simulator_listener(sim), t_(kds){
    CGAL_precondition(kds != NULL);
    if (P::notifier()->direction_of_time() != t_->direction_of_time()){
      t_->reverse_time();
    }
  }
  //! Pass DIRECTION_OF_TIME notifications via the set_direction_of_time method
  void new_notification(typename Simulator_listener::Notification_type t){
    if (t== Simulator_listener::DIRECTION_OF_TIME){
      if (P::notifier()->direction_of_time() != t_->direction_of_time()){
	t_->reverse_time();
	CGAL_postcondition(P::notifier()->direction_of_time() ==  t_->direction_of_time());
      } else {
	std::cerr << "ndir= "<< P::notifier()->direction_of_time() << " dir = " << t_->direction_of_time() << std::endl;
      }
    }
  }
protected:
  MOT *t_;
};


CGAL_KDS_END_NAMESPACE
#endif
