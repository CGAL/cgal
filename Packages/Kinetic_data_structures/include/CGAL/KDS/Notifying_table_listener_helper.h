#ifndef CGAL_KDS_NOTIFICATION_LH2_HELPERS_H
#define CGAL_KDS_NOTIFICATION_LH2_HELPERS_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE


//! A helper for KDSs interested in events from MovingObjectTables.
/*!
  The lister helper listens for EDITING notifications and calls the functions
  -new_object
  -change_object
  -delete_object
  with the appropriate object keys from the moving object table. 

  This helper should be stored as a member variable in the KDS since
  its pointer to the KDS is not protected.

  See CGAL::Listener for a description of what the
  Moving_object_table_listener template paramenter should provide.
*/
template <class Moving_object_table_listener, class KDS>
class Notifying_table_listener_helper: public Moving_object_table_listener {
  typedef typename Moving_object_table_listener::Notifier_pointer::element_type MOT;
  typedef Moving_object_table_listener P;
public:

  //! The constructor
  Notifying_table_listener_helper(typename Moving_object_table_listener::Notifier_pointer h,
				  KDS *kds):
    Moving_object_table_listener(h), t_(kds){
    for (typename Moving_object_table_listener::Notifier::Keys_iterator it= P::notifier()->keys_begin(); 
	 it != P::notifier()->keys_end(); ++it){
      t_->insert(*it);
    }
  }
  //! Pass EDITING notifications
  /*!  When editing changes to false, call new_object, changed_object,
    deleted_object for each new, changed or deleted object in the
    MovingObjectTable.
  */
  virtual void new_notification(typename Moving_object_table_listener::Notification_type et){
    if (et== P::IS_EDITING){
      if (P::notifier()->is_editing()==false){
	//! Note, this order is important
	for (typename MOT::Inserted_iterator it= P::notifier()->inserted_begin(); 
	     it != P::notifier()->inserted_end(); ++it){
	  t_->insert(*it);
	}
	for (typename MOT::Changed_iterator it= P::notifier()->changed_begin();
	     it != P::notifier()->changed_end(); ++it){
	  t_->set(*it);
	}
	for (typename MOT::Erased_iterator it= P::notifier()->erased_begin();
	     it != P::notifier()->erased_end(); ++it){
	  t_->erase(*it);
	}
      }
    }
  }
    
  

protected:
  KDS *t_;
};





CGAL_KDS_END_NAMESPACE
#endif
