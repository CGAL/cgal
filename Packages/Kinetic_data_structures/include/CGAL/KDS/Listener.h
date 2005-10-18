// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_TOOLS_LISTENER_BASE_H
#define CGAL_TOOLS_LISTENER_BASE_H
#include <CGAL/KDS/basic.h>
CGAL_KDS_BEGIN_NAMESPACE


//! This is the base class for all listener objects.
/*!  An object which wishes to provide runtime notifications to other
  objects should provide a contained class which uses this base. To
  use this base class, implement a class, here called Interface, which
  defines a type Interface::Notification_type and a type
  Interface::Notifier_pointer. 

  The Notification_type is generally an enum with one value for each
  type of notification which can be used.

  The Notifier_pointer is the type of a (ref counted) pointer to the
  object providing the notifications. The ref counter pointer must
  provide a nested type Pointer which is the type of a raw pointer.

  This object maintains a ref counted pointer to the object performing
  notifications. It is registered for notifications on construction
  and unregistered on destruction using the function set_listener on
  the object providing the notifications. The use of ref counted
  pointers means that as long as the notification object exists, the
  object providing the notifications must exist, ensuring that the
  object providing the notifications is not prematurely destroyed.

  These objects cannot be copied since the notifier only support one
  listener.

  Boost provides a similar functionality in the Boost.Signal
  package. However, it is quite a bit more complex (and
  flexible). This complexity add significantly to compile time and
  (although I did not test this directly), I suspect it is much slower
  at runtime due to the overhead of worrying about signal orders and
  not supporting single signals. In addition, it does not get on well
  with Qt due to collisions with the Qt moc keywords.

  There is also the TinyTL library which implements signals. As of
  writing it did not have any easy support for making sure all
  pointers are valid, so it did not seem to offer significant code
  saving over writing my own.
*/
template <class Interface>
class Listener: public Interface {
  typedef Listener<Interface> LB_this;
public:
  typedef typename Interface::Notifier_pointer::element_type Notifier;

  //typedef typename Notifier::Handle Notifier_handle;
  Listener(typename Interface::Notifier_pointer &nh): h_(nh){
    CGAL_precondition(h_->listener()==NULL);
    h_->set_listener(this);
  }

  Listener(Notifier* nh): h_(nh){
    CGAL_precondition(h_->listener()==NULL);
    h_->set_listener(this);
  }

  virtual ~Listener(){
    CGAL_precondition(h_->listener()==this);
    h_->set_listener(NULL);
  }
  //! Access the object providing notifications.
  /*!  This pointer is not intended to be stored, so a bare pointer is
    passed. Wrap it with a Notifier_pointer if you for some reason
    which to store it.
  */
  typename Interface::Notifier_pointer::element_type* notifier(){
    return h_.get();
  }
  //! Constant version. 
  /*!
    See Listener::notifier()
  */
  const typename Interface::Notifier_pointer::element_type* notifier() const {
    return h_.get();
  }

  //! The method called when there is a runtime notification to be made
  /*!  The Notification_type is the type of notification to be made,
    generally an enum with the name of a field of the object providing
    notifications.
  */
  virtual void new_notification(typename Interface::Notification_type nt)=0;

private:
  Listener(const LB_this &o){
    h_->this_is_not_a_function();
    CGAL_assertion(0);
  }
  const LB_this& operator=(const LB_this &o){
    CGAL_assertion(0);
    h_->this_is_not_a_function();
    return *this;
  }
protected:
  typename Interface::Notifier_pointer h_;
};



CGAL_KDS_END_NAMESPACE
#endif
