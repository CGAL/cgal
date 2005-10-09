#ifndef CGAL_TOOLS_MULTI_LISTENER_BASE_H
#define CGAL_TOOLS_MULTI_LISTENER_BASE_H
#include <CGAL/basic.h>
CGAL_KDS_BEGIN_NAMESPACE

//! This is a variant of Listener which supports more than one object receiving notifications
/*!
  See Listener for full documentation. 

  In contrast to listener, this can be copied. 
*/
template <class Interface>
class Multi_listener: public Interface {
  typedef Multi_listener<Interface> This;
public:
  typedef typename Interface::Notifier_pointer::element_type Notifier;

  //typedef typename Notifier::Handle Notifier_handle;
  Multi_listener(typename Interface::Notifier_pointer &nh): h_(nh){
    h_->new_listener(this);
  }

  Multi_listener(Notifier* nh): h_(nh){
    h_->new_listener(this);
  }

  virtual ~Multi_listener(){
    h_->delete_listener(this);
  }
  //! The object doing the notifying
  typename Interface::Notifier_pointer::element_type* notifier(){
    return h_.get();
  }
  //! Constant access to the object doing the notifying
  const typename Interface::Notifier_pointer::element_type* notifier() const {
    //const typename Interface::Notifier_pointer p= h_.pointer();
    return h_.get();
  }
  //! This is called when there is a notification.
  virtual void new_notification(typename Interface::Notification_type nt)=0;

  //! Copy and subscribe the new object.
  Multi_listener(const This &o){
    h_= o.h_;
    h_->new_listener(this);
  }
  //! Copy and subscribe the new object
  This operator=(const This &o){
    h_= o.h_;
    h_->new_listener(this);
  }
protected:
  typename Interface::Notifier_pointer h_;
};

CGAL_KDS_END_NAMESPACE
#endif
