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

#ifndef CGAL_TOOLS_LISTENER_BASE_H
#define CGAL_TOOLS_LISTENER_BASE_H
#include <CGAL/Kinetic/basic.h>
#include <boost/utility.hpp>
namespace CGAL { namespace Kinetic {

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
class Listener_base: public Interface, public boost::noncopyable
{
  typedef Listener_base<Interface> LB_this;
public:
  typedef typename Interface::Notifier Notifier;  

  Listener_base(){}

  virtual ~Listener_base() {
    if (h_ != NULL) {
      CGAL_precondition(h_->listener()==this);
      h_->set_listener(NULL);
    }
  }
  //! Access the object providing notifications.
  /*!  This pointer is not intended to be stored, so a bare pointer is
    passed. Wrap it with a Notifier_pointer if you for some reason
    which to store it.
  */
  Notifier* notifier() {
    return h_.get();
  }
  //! Constant version.
  /*!
    See Listener::notifier()
  */
  const Notifier* notifier() const
  {
    return h_.get();
  }

  template <class Ptr>
  void set_notifier(Ptr t) {
    if (h_!= NULL) h_->set_listener(NULL);
    h_=t;
    CGAL_precondition(h_->listener() == NULL);
    h_->set_listener(this);
  }

 
  //! The method called when there is a runtime notification to be made
  /*!  The Notification_type is the type of notification to be made,
    generally an enum with the name of a field of the object providing
    notifications.
  */
  virtual void new_notification(typename Interface::Notification_type nt)=0;

  class Handle: public boost::noncopyable {
  public:
    Handle(): p_(NULL){}
    Handle(LB_this *p):p_(p){}
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
    LB_this* get() {return p_;}
  private:
    
    LB_this *p_;
  };
  
protected:
  typename Notifier::Handle h_;
};

/*
#define CGAL_KINETIC_LISTENER(...) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {__VA_ARGS__} Notification_type;		\
  };								\
  public:							\
  typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
  friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
  void set_listener(Listener *sk) {				\
    listener_=sk;						\
  }								\
  Listener* listener() {return listener_.get();}		\
  typename Listener::Handle listener_;


#define CGAL_KINETIC_LISTENERNT(...) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {__VA_ARGS__} Notification_type;		\
  };								\
public:								\
 typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
 friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
 void set_listener(Listener *sk) {				\
   listener_=sk;						\
 }								\
 Listener* listener() {return listener_.get();}			\
 Listener::Handle listener_;
*/
#define CGAL_KINETIC_LISTENER1(A) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {A} Notification_type;				\
  };								\
  public:							\
  typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
  friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
  void set_listener(Listener *sk) {				\
    listener_=sk;						\
  }								\
  Listener* listener() {return listener_.get();}		\
  typename Listener::Handle listener_;


#define CGAL_KINETIC_LISTENERNT1(A) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {A} Notification_type;				\
  };								\
public:								\
 typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
 friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
 void set_listener(Listener *sk) {				\
   listener_=sk;						\
 }								\
 Listener* listener() {return listener_.get();}			\
 Listener::Handle listener_;

#define CGAL_KINETIC_LISTENER2(A,B) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {A,B} Notification_type;			\
  };								\
  public:							\
  typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
  friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
  void set_listener(Listener *sk) {				\
    listener_=sk;						\
  }								\
  Listener* listener() {return listener_.get();}		\
  typename Listener::Handle listener_;


#define CGAL_KINETIC_LISTENERNT2(A,B) private:			\
  struct Listener_core{						\
    typedef This Notifier;					\
    typedef enum {A,B} Notification_type;			\
  };								\
public:								\
 typedef CGAL::Kinetic::Listener_base<Listener_core> Listener;	\
 friend class CGAL::Kinetic::Listener_base<Listener_core>;	\
private:							\
 void set_listener(Listener *sk) {				\
   listener_=sk;						\
 }								\
 Listener* listener() {return listener_.get();}			\
 Listener::Handle listener_;

#define CGAL_KINETIC_NOTIFY(field) if (listener_!= NULL) listener_->new_notification(Listener::field)

#define CGAL_KINETIC_LISTENER_DESTRUCTOR CGAL_assertion(listener_==NULL)


#define CGAL_KINETIC_LISTENER_BASICS(Name, KDS)		\
  public:						\
  Name(): recipient_(NULL){}				\
  typedef KDS Recipient;				\
  Recipient* recipient() const {return recipient_;}	\
  void set_recipient(Recipient *r){recipient_=r;}	\
private:						\
 Recipient* recipient_;




#define CGAL_KINETIC_LISTEN1(Notifier, NOTIF, function)\
  private:								\
  class Notifier##_listener: public Notifier::Listener {		\
    CGAL_KINETIC_LISTENER_BASICS(Notifier##_listener, This)		\
  public:								\
    virtual void new_notification(typename Notifier::Listener::Notification_type t) { \
      if (recipient() != NULL && t== Notifier::Listener::NOTIF) recipient()->function; \
      else {								\
      }									\
    }									\
  };								\
  friend class Notifier##_listener;					\
  Notifier##_listener listener_##Notifier##_;

#define CGAL_KINETIC_LISTEN2(Notifier, NOTIF, function, NOTIF2, function2) \
  private:								\
  class Notifier##_listener: public Notifier::Listener {                \
    CGAL_KINETIC_LISTENER_BASICS(Notifier##_listener, This)		\
  public:								\
    virtual void new_notification(typename Notifier::Listener::Notification_type t) { \
      if (recipient()== NULL) return;					\
      if (t== Notifier::Listener::NOTIF) 	recipient()->function();	\
      else if (t== Notifier::Listener::NOTIF2) 	recipient()->function2; \
      else {								\
      }									\
    }									\
  };                                                                    \
  friend class Notifier##_listener;					\
  Notifier##_listener listener_##Notifier##_;

#define CGAL_KINETIC_NOTIFIER(Notifier) listener_##Notifier##_.notifier()

#define CGAL_KINETIC_INIT_LISTEN(Notifier, ptr)	\
  listener_##Notifier##_.set_recipient(this);	\
  listener_##Notifier##_.set_notifier(ptr);
  

  
} } //namespace CGAL::Kinetic
#endif
