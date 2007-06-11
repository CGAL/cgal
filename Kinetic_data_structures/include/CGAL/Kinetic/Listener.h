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
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_TOOLS_LISTENER_BASE_H
#define CGAL_TOOLS_LISTENER_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <boost/utility.hpp>
CGAL_KINETIC_BEGIN_NAMESPACE

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
class Listener_base: public Interface
{
  typedef Listener_base<Interface> LB_this;
public:
  typedef typename Interface::Notifier_handle::element_type Notifier;

  //typedef typename Notifier::Handle Notifier_handle;
  Listener_base(typename Interface::Notifier_handle &nh): h_(nh) {
    CGAL_precondition(h_->listener()==NULL);
    h_->set_listener(this);
  }

  Listener_base(Notifier* nh): h_(nh) {
    CGAL_precondition(h_->listener()==NULL);
    h_->set_listener(this);
  }

  Listener_base(){}

  virtual ~Listener_base() {
    CGAL_precondition(h_->listener()==this);
    h_->set_listener(NULL);
  }
  //! Access the object providing notifications.
  /*!  This pointer is not intended to be stored, so a bare pointer is
    passed. Wrap it with a Notifier_pointer if you for some reason
    which to store it.
  */
  typename Interface::Notifier_handle::element_type* notifier() {
    return h_.get();
  }
  //! Constant version.
  /*!
    See Listener::notifier()
  */
  const typename Interface::Notifier_handle::element_type* notifier() const
  {
    return h_.get();
  }

  //! The method called when there is a runtime notification to be made
  /*!  The Notification_type is the type of notification to be made,
    generally an enum with the name of a field of the object providing
    notifications.
  */
  virtual void new_notification(typename Interface::Notification_type nt)=0;

  struct Pointer:public boost::noncopyable {
    Pointer(): p_(NULL){}
    Pointer(LB_this *p):p_(p){}
    LB_this &operator*(){return *p_;}
    const LB_this &operator*() const {return *p_;}
    LB_this *operator->(){return p_;}
    const LB_this *operator->() const {return p_;}
    bool operator==(void* o) const {
      return p_== o;
    }
    bool operator!=(void* o) const {
      return p_!= o;
    }
    void operator=(LB_this *o) {
      p_=o;
    }
    /*CGAL_COPY_CONSTRUCTOR(Pointer);
    void copy_from(const LB_this &o) {
      
    }*/
    LB_this* get() {return p_;}
  private:
    
    LB_this *p_;
  };
  
private:
  Listener_base(const LB_this &o) {
    h_->this_is_not_a_function();
    CGAL_assertion(0);
  }
  const LB_this& operator=(const LB_this &o) {
    CGAL_assertion(0);
    h_->this_is_not_a_function();
    return *this;
  }
protected:
  typename Interface::Notifier_handle h_;
};


#define CGAL_KINETIC_LISTENER(...) private:		  \
  struct Listener_core{\
      typedef typename This::Handle Notifier_handle;\
      typedef enum {__VA_ARGS__} Notification_type;\
  };						   \
public:\
  typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
  friend class Listener_base<Listener_core>;				\
private:							\
  void set_listener(Listener *sk) {				\
    listener_=sk;						\
  }								\
  Listener* listener() const {return listener_.get();}		\
  typename Listener::Pointer listener_;

#define CGAL_KINETIC_SIGNAL(field) if (listener_!= NULL) listener_->new_notification(Listener::field)

#define CGAL_KINETIC_LISTENER_DESTRUCTOR CGAL_assertion(listener_==NULL);

#define CGAL_KINETIC_LISTEN1(Notifier, Signal, func) private:\
  class Notifier##_listener: public Notifier::Listener {\
    This *t_;						\
  public:					\
    Notifier##_listener(Notifier* tm, This &t): Notifier::Listener(tm), t_(t){}	\
    void new_notification(typename Notifier::Listener::Notification_type t) {\
      if (t== Notifier::Listener::Signal) t_->func();			\
    }\
  };\
  friend class Notifier##_listener;\
  Notifier##_listener listener_##Notifier##_;

#define CGAL_KINETIC_INIT_LISTEN(Notifier, ptr) listener_##Notifier##_=Notifier##_listener(ptr, const_cast<This*>(this))
 
CGAL_KINETIC_END_NAMESPACE
#endif
