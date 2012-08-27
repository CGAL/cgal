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

#ifndef CGAL_TOOLS_MULTI_LISTENER_BASE_H
#define CGAL_TOOLS_MULTI_LISTENER_BASE_H
#include <CGAL/basic.h>
namespace CGAL { namespace Kinetic {

//! This is a variant of Listener which supports more than one object receiving notifications
/*!
  See Listener for full documentation.

  In contrast to listener, this can be copied.
*/
template <class Interface>
class Multi_listener_base: public Interface
{
  typedef Multi_listener_base<Interface> This;
public:
  typedef typename Interface::Notifier_handle::element_type Notifier;

  Multi_listener_base(){}

  virtual ~Multi_listener_base() {
    if (h_ != NULL) h_->delete_listener(this);
  }

  void set_notifier(Notifier * n) {
    if (h_!= NULL) h_->delete_listener(this);
    h_=n;
    h_->new_listener(this);
  }
  void set_notifier(typename Interface::Notifier_handle t) {
    set_notifier(t.get());
  }


  typename Interface::Notifier_handle::element_type* notifier() {
    return h_.get();
  }

  const typename Interface::Notifier_handle::element_type* notifier() const
  {

    return h_.get();
  }

  virtual void new_notification(typename Interface::Notification_type nt)=0;


  Multi_listener_base(const This &o) {
    h_= o.h_;
    h_->new_listener(this);
  }

  const This& operator=(const This &o) {
    h_= o.h_;
    h_->new_listener(this);
    return *this;
  }
protected:
  typename Interface::Notifier_handle h_;
};

/*
#define CGAL_KINETIC_MULTILISTENER(...) private:		\
  struct Listener_core{						\
    typedef typename This::Handle Notifier_handle;		\
      typedef enum {__VA_ARGS__} Notification_type;		\
  };								\
public:								\
 typedef Multi_listener_base<Listener_core> Listener;		\
 friend class Multi_listener_base<Listener_core>;		\
private:							\
 void new_listener(Listener *sk) {				\
   listeners_.push_back(sk);					\
 }								\
 void delete_listener(Listener *kds) {				\
   for (unsigned int i=0; i< listeners_.size(); ++i){		\
     if (listeners_[i] == kds) {				\
       std::swap(listeners_[i], listeners_.back());		\
       listeners_.pop_back();					\
       return;							\
     }								\
   }								\
 }								\
 std::vector<Listener*> listeners_;


#define CGAL_KINETIC_MULTILISTENERNT(...) private:		\
  struct Listener_core{						\
    typedef This::Handle Notifier_handle;			\
    typedef enum {__VA_ARGS__} Notification_type;		\
  };								\
public:								\
 typedef Multi_listener_base<Listener_core> Listener;		\
 friend class Multi_listener_base<Listener_core>;		\
private:							\
 void new_listener(Listener *sk) {				\
   listeners_.push_back(sk);					\
 }								\
 void delete_listener(Listener *kds) {				\
   for (unsigned int i=0; i< listeners_.size(); ++i){		\
     if (listeners_[i] == kds) {				\
       std::swap(listeners_[i], listeners_.back());		\
       listeners_.pop_back();					\
       return;							\
     }								\
   }								\
 }								\
 std::vector<Listener*> listeners_;
*/



#define CGAL_KINETIC_MULTILISTENER1(A) private:                         \
  struct Listener_core{                                                 \
    typedef typename This::Handle Notifier_handle;                      \
    typedef enum {A} Notification_type;                                 \
  };                                                                    \
public:                                                                 \
 typedef CGAL::Kinetic::Multi_listener_base<Listener_core> Listener;    \
 friend class CGAL::Kinetic::Multi_listener_base<Listener_core>;        \
private:                                                                \
 void new_listener(Listener *sk) {                                      \
   listeners_.push_back(sk);                                            \
 }                                                                      \
 void delete_listener(Listener *kds) {                                  \
   for (unsigned int i=0; i< listeners_.size(); ++i){                   \
     if (listeners_[i] == kds) {                                        \
       std::swap(listeners_[i], listeners_.back());                     \
       listeners_.pop_back();                                           \
       return;                                                          \
     }                                                                  \
   }                                                                    \
 }                                                                      \
 std::vector<Listener*> listeners_;


#define CGAL_KINETIC_MULTILISTENERNT1(A) private:                       \
  struct Listener_core{                                                 \
    typedef This::Handle Notifier_handle;                               \
    typedef enum {A} Notification_type;                                 \
  };                                                                    \
public:                                                                 \
 typedef CGAL::Kinetic::Multi_listener_base<Listener_core> Listener;    \
 friend class CGAL::Kinetic::Multi_listener_base<Listener_core>;        \
private:                                                                \
 void new_listener(Listener *sk) {                                      \
   listeners_.push_back(sk);                                            \
 }                                                                      \
 void delete_listener(Listener *kds) {                                  \
   for (unsigned int i=0; i< listeners_.size(); ++i){                   \
     if (listeners_[i] == kds) {                                        \
       std::swap(listeners_[i], listeners_.back());                     \
       listeners_.pop_back();                                           \
       return;                                                          \
     }                                                                  \
   }                                                                    \
 }                                                                      \
 std::vector<Listener*> listeners_;



#define CGAL_KINETIC_MULTILISTENER2(A,B) private:			\
  struct Listener_core{                                                 \
    typedef typename This::Handle Notifier_handle;                      \
    typedef enum {A,B} Notification_type;                               \
  };                                                                    \
public:                                                                 \
 typedef CGAL::Kinetic::Multi_listener_base<Listener_core> Listener;    \
 friend class CGAL::Kinetic::Multi_listener_base<Listener_core>;        \
private:                                                                \
 void new_listener(Listener *sk) {                                      \
   listeners_.push_back(sk);                                            \
 }                                                                      \
 void delete_listener(Listener *kds) {                                  \
   for (unsigned int i=0; i< listeners_.size(); ++i){                   \
     if (listeners_[i] == kds) {                                        \
       std::swap(listeners_[i], listeners_.back());                     \
       listeners_.pop_back();                                           \
       return;                                                          \
     }                                                                  \
   }                                                                    \
 }                                                                      \
 std::vector<Listener*> listeners_;


#define CGAL_KINETIC_MULTILISTENERNT2(A,B) private:                     \
  struct Listener_core{                                                 \
    typedef This::Handle Notifier_handle;                               \
    typedef enum {A,B} Notification_type;                               \
  };                                                                    \
public:                                                                 \
 typedef CGAL::Kinetic::Multi_listener_base<Listener_core> Listener;    \
 friend class CGAL::Kinetic::Multi_listener_base<Listener_core>;        \
private:                                                                \
 void new_listener(Listener *sk) {                                      \
   listeners_.push_back(sk);                                            \
 }                                                                      \
 void delete_listener(Listener *kds) {                                  \
   for (unsigned int i=0; i< listeners_.size(); ++i){                   \
     if (listeners_[i] == kds) {                                        \
       std::swap(listeners_[i], listeners_.back());                     \
       listeners_.pop_back();                                           \
       return;                                                          \
     }                                                                  \
   }                                                                    \
 }                                                                      \
 std::vector<Listener*> listeners_;

#define CGAL_KINETIC_MULTINOTIFY(field) do {                            \
  for(typename std::vector<Listener*>::iterator it= listeners_.begin(); it != listeners_.end(); ++it){ \
    (*it)->new_notification(Listener::field);				\
  }} while (false)

#define CGAL_KINETIC_MULTILISTENER_DESTRUCTOR CGAL_assertion(listeners_.empty())

} } //namespace CGAL::Kinetic
#endif
