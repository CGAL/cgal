#include <CGAL/Kinetic/Listener.h>
#include <CGAL/Kinetic/Ref_counted.h>

template <class Data>
struct Notifier: public CGAL::Kinetic::Ref_counted<Notifier<Data> >
{
  typedef Notifier<Data> This;
public:
  Notifier(): data_(0){}

  typedef CGAL::Kinetic::Ref_counted<Notifier> Base;
  CGAL_KINETIC_LISTENER1(DATA_CHANGED)
public:

  int data() const {return data_;}
  void set_data(int d) {
    if (d != data_) {
      data_=d;
      CGAL_KINETIC_NOTIFY(DATA_CHANGED);
    }
  }


  Data data_;
};

template <class Data>
class Receiver: public CGAL::Kinetic::Ref_counted<Receiver<Data> > {
  typedef Receiver<Data> This;
  typedef ::Notifier<Data> Notifier;
  CGAL_KINETIC_LISTEN1(Notifier, DATA_CHANGED, ping())
public:
  Receiver( Notifier* p){
    CGAL_KINETIC_INIT_LISTEN(Notifier, p);
  }
  void ping() const {
    std::cout << "Data changed " << std::endl;
  }
};


int main(int, char *[])
{
  {
    Notifier<int>::Handle n= new Notifier<int>();
    Receiver<int>::Handle r= new Receiver<int>(n.get());

    n->set_data(4);
    n->set_data(4);
    n->set_data(5);
    r=NULL;
    n->set_data(6);
  }

  return EXIT_SUCCESS;
}
