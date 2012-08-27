// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_NOTIFICATION_SK_HELPERS_H
#define CGAL_KINETIC_NOTIFICATION_SK_HELPERS_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Multi_listener.h>

namespace CGAL { namespace Kinetic {

//! This is a helper for KDSs which want to audit themselves when appropriate..
/*!  This listener object listens for
  Simulator_listener::HAS_EVENT_FREE_TIME notificatations and calls
  the audit method of the KDS. This object should be stored as a
  member variable in the KDS, since the pointer it contains to the KDS
  is not protected.

  See CGAL::Listener for a description of what the Simulator_listener
  template paramenter should provide.
*/
template <class SL, class KDS>
class Simulator_listener: public SL
{
  typedef SL P;
public:
  CGAL_KINETIC_LISTENER_BASICS(Simulator_listener, KDS)

  //! Pass HAS_AUDIT_TIME notifications via a call to the audit() function
  void new_notification(typename P::Notification_type t) {
    if (t== P::HAS_AUDIT_TIME) {
      if (P::notifier()->has_audit_time()) {
	recipient()->audit();
      }
    }
  }
};


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
class Active_objects_listener: public Moving_object_table_listener
{
  typedef typename Moving_object_table_listener::Notifier_handle::element_type MOT;
  typedef Moving_object_table_listener P;
public:
  CGAL_KINETIC_LISTENER_BASICS(Active_objects_listener, KDS)
public:
  //! Pass EDITING notifications
  /*!  When editing changes to false, call new_object, changed_object,
    deleted_object for each new, changed or deleted object in the
    MovingObjectTable.
  */
  virtual void new_notification(typename Moving_object_table_listener::Notification_type et) {
    if (et== P::IS_EDITING) {
      if (P::notifier()->is_editing()==false) {
	//! Note, this order is important
	for (typename MOT::Inserted_iterator it= P::notifier()->inserted_begin();
	     it != P::notifier()->inserted_end(); ++it) {
	  recipient()->insert(*it);
	}
	for (typename MOT::Changed_iterator it= P::notifier()->changed_begin();
	     it != P::notifier()->changed_end(); ++it) {
	  recipient()->set(*it);
	}
	for (typename MOT::Erased_iterator it= P::notifier()->erased_begin();
	     it != P::notifier()->erased_end(); ++it) {
	  recipient()->erase(*it);
	}
      }
    }
  }

  void catch_up() {
    for (typename Moving_object_table_listener::Notifier::Key_iterator it= P::notifier()->keys_begin();
	 it != P::notifier()->keys_end(); ++it) {
      recipient()->insert(*it);
    }
  }

};


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
class Active_objects_batch_listener: public Moving_object_table_listener
{
  typedef typename Moving_object_table_listener::Notifier_handle::element_type MOT;
  typedef Moving_object_table_listener P;
public:
  CGAL_KINETIC_LISTENER_BASICS(Active_objects_batch_listener, KDS)
  //! Pass EDITING notifications
  /*!  When editing changes to false, call new_object, changed_object,
    deleted_object for each new, changed or deleted object in the
    MovingObjectTable.
  */
public:
  virtual void new_notification(typename Moving_object_table_listener::Notification_type et) {
    if (et== P::IS_EDITING) {
      if (P::notifier()->is_editing()==false) {
	CGAL_precondition(!recipient()->is_batch_editing());
	recipient()->set_is_batch_editing(true);
	//! Note, this order is important
	for (typename MOT::Inserted_iterator it= P::notifier()->inserted_begin();
	     it != P::notifier()->inserted_end(); ++it) {
	  recipient()->insert(*it);
	}
	for (typename MOT::Changed_iterator it= P::notifier()->changed_begin();
	     it != P::notifier()->changed_end(); ++it) {
	  recipient()->set(*it);
	}
	for (typename MOT::Erased_iterator it= P::notifier()->erased_begin();
	     it != P::notifier()->erased_end(); ++it) {
	  recipient()->erase(*it);
	}
	recipient()->set_is_batch_editing(false);
      }
    }
  }

  void catch_up() {
    for (typename Moving_object_table_listener::Notifier::Key_iterator it= P::notifier()->keys_begin();
	 it != P::notifier()->keys_end(); ++it) {
      recipient()->insert(*it);
    }
  }
};

#define CGAL_KINETIC_DECLARE_LISTENERS(S, A) private:			\
   typedef typename CGAL::Kinetic::Simulator_listener<S::Listener, This> Simulator_listener;\
  friend  class CGAL::Kinetic::Simulator_listener<S::Listener, This>;	\
  typedef typename CGAL::Kinetic::Active_objects_listener<A::Listener, This> Moving_point_table_listener; \
  friend class CGAL::Kinetic::Active_objects_listener<A::Listener, This>; \
  Simulator_listener siml_;\
  Moving_point_table_listener motl_;


#define CGAL_KINETIC_DECLARE_AOT_LISTENER(A) private:			\
  typedef typename CGAL::Kinetic::Active_objects_listener<A::Listener, This> Moving_point_table_listener; \
  friend class CGAL::Kinetic::Active_objects_listener<A::Listener, This>; \
  Moving_point_table_listener motl_;


#define CGAL_KINETIC_INITIALIZE_LISTENERS(sh, ph)			\
  siml_.set_recipient(this); siml_.set_notifier(sh);			\
  motl_.set_recipient(this); motl_.set_notifier(ph);			\
  motl_.catch_up()


#define CGAL_KINETIC_INITIALIZE_AOT_LISTENER(ph)			\
  motl_.set_recipient(this); motl_.set_notifier(ph);			\
  motl_.catch_up()


#define CGAL_KINETIC_DECLARE_BATCH_LISTENERS(S, A) private:		\
  typedef typename CGAL::Kinetic::Simulator_listener<S::Listener, This> Simulator_listener; \
  friend  class CGAL::Kinetic::Simulator_listener<S::Listener, This>;\
  typedef typename CGAL::Kinetic::Active_objects_batch_listener<A::Listener, This> Moving_point_table_listener; \
  friend class CGAL::Kinetic::Active_objects_batch_listener<A::Listener, This>;\
  Simulator_listener siml_;						\
  Moving_point_table_listener motl_;

#define CGAL_KINETIC_INITIALIZE_BATCH_LISTENERS(sh, ph, insert) \
  siml_.set_recipient(this); siml_.set_notifier(sh);		\
  motl_.set_recipient(this); motl_.set_notifier(ph);		\
  if (insert) motl_.catch_up()

} } //namespace CGAL::Kinetic
#endif
