#include <CGAL/KDS/Listener.h>
#include <CGAL/KDS/Ref_counted.h>

struct Notifier: public CGAL::KDS::Ref_counted<Notifier> {
  Notifier(): data_(0), listener_(NULL){}

private:
  struct Listener_interface {
    typedef enum Notification_type {DATA_CHANGED} Notification_type;
    typedef Notifier::Pointer Notifier_pointer;
  };


public:

  typedef CGAL::KDS::Listener<Listener_interface> Listener;
  friend class CGAL::KDS::Listener<Listener_interface>;

  int data() const {return data_;}
  void set_data(int d) { 
    if (d != data_) {
      data_=d;
      if (listener_ != NULL) listener_->new_notification(Listener_interface::DATA_CHANGED);
    }
  }

protected:
  void set_listener(Listener *l) {
    listener_= l;
  }
  Listener* listener() const {return listener_;}

  int data_;
  Listener *listener_;
};


struct My_listener: public Notifier::Listener, public CGAL::KDS::Ref_counted<My_listener> {
  typedef Notifier::Listener P;
  typedef P::Notifier_pointer PP;
  My_listener(PP p): P(p){}
  
  
  void new_notification(P::Notification_type nt) {
    if (nt == P::DATA_CHANGED) {
      std::cout << "Data is now " << P::notifier()->data()<< std::endl;
    } else {
      std::cerr << "Unknown notification type: " << nt << std::endl;
    }
  }
};

int main(int, char *[]){
  {
    Notifier::Pointer n= new Notifier;
    My_listener::Pointer ml= new My_listener(n);
    
    n->set_data(4);
    n->set_data(4);
    n->set_data(5);
    ml=NULL;
    n->set_data(6);
  }
  {
    Notifier::Pointer n= new Notifier;
    My_listener::Pointer ml= new My_listener(n);
    
    n->set_data(4);
    n->set_data(4);
    n->set_data(5);
    n=NULL;
  }
  return EXIT_SUCCESS;
}
