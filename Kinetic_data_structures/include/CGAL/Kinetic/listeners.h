// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kinetic_data_structures/include/CGAL/Kinetic/Simulator_kds_listener.h $
// $Id: Simulator_kds_listener.h 29334 2006-03-10 00:00:09Z drussel $
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_NOTIFICATION_SK_HELPERS_H
#define CGAL_KINETIC_NOTIFICATION_SK_HELPERS_H
#include <CGAL/Kinetic/basic.h>

CGAL_KINETIC_BEGIN_NAMESPACE

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
  //! The only constructor
  Simulator_listener(typename P::Notifier_handle sim,
			 KDS *kds): P(sim), t_(kds) {
    CGAL_precondition(kds != NULL);
  }
  Simulator_listener(){}
  //! Pass HAS_AUDIT_TIME notifications via a call to the audit() function
  void new_notification(typename P::Notification_type t) {
    if (t== P::HAS_AUDIT_TIME) {
      if (P::notifier()->has_audit_time()) {
	t_->audit();
      }
    }
  }
protected:
  KDS *t_;
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

  //! The constructor
  Active_objects_listener(typename Moving_object_table_listener::Notifier_handle h,
				 KDS *kds):
    Moving_object_table_listener(h), t_(kds) {
    for (typename Moving_object_table_listener::Notifier::Key_iterator it= P::notifier()->keys_begin();
	 it != P::notifier()->keys_end(); ++it) {
      t_->insert(*it);
    }
  }
  Active_objects_listener(){}
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
	  t_->insert(*it);
	}
	for (typename MOT::Changed_iterator it= P::notifier()->changed_begin();
	     it != P::notifier()->changed_end(); ++it) {
	  t_->set(*it);
	}
	for (typename MOT::Erased_iterator it= P::notifier()->erased_begin();
	     it != P::notifier()->erased_end(); ++it) {
	  t_->erase(*it);
	}
      }
    }
  }

protected:
  KDS *t_;
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

  //! The constructor
  Active_objects_batch_listener(typename Moving_object_table_listener::Notifier_handle h,
				       KDS *kds,
				       bool insert):
    Moving_object_table_listener(h), t_(kds) {
    if (insert) {
      for (typename Moving_object_table_listener::Notifier::Key_iterator it= P::notifier()->keys_begin();
	   it != P::notifier()->keys_end(); ++it) {
	t_->insert(*it);
      }
    }
  }
  Active_objects_batch_listener(){}
  //! Pass EDITING notifications
  /*!  When editing changes to false, call new_object, changed_object,
    deleted_object for each new, changed or deleted object in the
    MovingObjectTable.
  */
  virtual void new_notification(typename Moving_object_table_listener::Notification_type et) {
    if (et== P::IS_EDITING) {
      if (P::notifier()->is_editing()==false) {
	CGAL_precondition(!t_->is_batch_editing());
	t_->set_is_batch_editing(true);
	//! Note, this order is important
	for (typename MOT::Inserted_iterator it= P::notifier()->inserted_begin();
	     it != P::notifier()->inserted_end(); ++it) {
	  t_->insert(*it);
	}
	for (typename MOT::Changed_iterator it= P::notifier()->changed_begin();
	     it != P::notifier()->changed_end(); ++it) {
	  t_->set(*it);
	}
	for (typename MOT::Erased_iterator it= P::notifier()->erased_begin();
	     it != P::notifier()->erased_end(); ++it) {
	  t_->erase(*it);
	}
	t_->set_is_batch_editing(false);
      }
    }
  }

protected:
  KDS *t_;
};

#define CGAL_KINETIC_DECLARE_LISTENERS(S, A) private:			\
   typedef typename CGAL::Kinetic::Simulator_listener<S::Listener, This> Simulator_listener;\
  friend  class CGAL::Kinetic::Simulator_listener<S::Listener, This>;	\
  typedef typename CGAL::Kinetic::Active_objects_listener<A::Listener, This> Moving_point_table_listener; \
  friend class CGAL::Kinetic::Active_objects_listener<A::Listener, This>; \
  Simulator_listener siml_;\
  Moving_point_table_listener motl_; 


#define CGAL_KINETIC_INITIALIZE_LISTENERS(sh, ph)			\
  siml_= Simulator_listener(sh, this);\
  motl_= Moving_point_table_listener(ph, this); \



#define CGAL_KINETIC_DECLARE_BATCH_LISTENERS(S, A) private:		\
  typedef typename CGAL::Kinetic::Simulator_listener<S::Listener, This> Simulator_listener; \
  friend  class CGAL::Kinetic::Simulator_listener<S::Listener, This>;\
  typedef typename CGAL::Kinetic::Active_objects_batch_listener<A::Listener, This> Moving_point_table_listener; \
  friend class CGAL::Kinetic::Active_objects_batch_listener<A::Listener, This>;\
  Simulator_listener siml_;						\
  Moving_point_table_listener motl_; 

#define CGAL_KINETIC_INITIALIZE_BATCH_LISTENERS(sh, ph, insert) \
  siml_= Simulator_listener(sh, this);				\
  motl_= Moving_point_table_listener(ph, this, insert);		\

CGAL_KINETIC_END_NAMESPACE
#endif
