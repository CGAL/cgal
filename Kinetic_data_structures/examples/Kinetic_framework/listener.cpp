#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

template <typename H>
struct Listener_interface_impl
  {
  public:
    typedef enum Notification_type {DATA_CHANGED}
      Notification_type;
   typedef H Notifier_handle;
  };


struct Notifier: public CGAL::Kinetic::Ref_counted<Notifier>
{
public:
  Notifier(): data_(0), listener_(NULL){}

  typedef CGAL::Kinetic::Ref_counted<Notifier> Base;
  typedef Listener_interface_impl<Base::Handle> Listener_interface;
  typedef CGAL::Kinetic::Listener<Listener_interface> Listener;


  int data() const {return data_;}
  void set_data(int d) {
    if (d != data_) {
      data_=d;
      if (listener_ != NULL) listener_->new_notification(Listener_interface::DATA_CHANGED);
    }
  }

  //protected:
  void set_listener(Listener *l) {
    listener_= l;
  }
  Listener* listener() const {return listener_;}

  int data_;
  Listener *listener_;
};

struct My_listener: public Notifier::Listener, public CGAL::Kinetic::Ref_counted<My_listener>
{
  typedef Notifier::Listener P;
  typedef P::Notifier_handle PP;
  My_listener(PP p): P(p){}

  void new_notification(P::Notification_type nt) {
    if (nt == P::DATA_CHANGED) {
      std::cout << "Data is now " << P::notifier()->data()<< std::endl;
    }
    else {
      std::cerr << "Unknown notification type: " << nt << std::endl;
    }
  }
};

int main(int, char *[])
{
  {
    Notifier::Handle n= new Notifier;
    My_listener::Handle ml= new My_listener(n);

    n->set_data(4);
    n->set_data(4);
    n->set_data(5);
    ml=NULL;
    n->set_data(6);
  }
  {
    Notifier::Handle n= new Notifier;
    My_listener::Handle ml= new My_listener(n);

    n->set_data(4);
    n->set_data(4);
    n->set_data(5);
    n=NULL;
  }

  return EXIT_SUCCESS;
}
