#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Event_base.h>
#include <iostream>
#include <set>
#include <vector>

#ifdef _MSC_VER
#pragma warning(disable:4355) // complaint about using 'this' to
                                // initialize a member
#endif

// This must be external since operator<< has to be defined
template <class Object, class Time, class KDS>
struct Trivial_event: public CGAL::Kinetic::Event_base<KDS*>
{
  typedef CGAL::Kinetic::Event_base<KDS*> P;
  Trivial_event(){}
  template <class It, class Map>
  Trivial_event(It beg, It end, const Map &m, KDS *kds): P(kds) {
    for (; beg != end; ++beg) {
      objects_.push_back(m[*beg]);
    }
  }
  void write(std::ostream &out) const {
    out << "An event";
  }
  void process() const
  {
    std::cout << "The following objects are in the table: ";
    for (typename std::vector<Object>::const_iterator
	   cit= objects_.begin();
	 cit != objects_.end(); ++cit) {
      std::cout << *cit << " ";
    }
    std::cout << std::endl;
    P::kds()->set_processed(true);
  }

  std::vector<Object> objects_;
};



/*!  This is a trivial kinetic data structure which doesn't actually
  do anything, other than print what it should be doing and maintain 1
  certificate in the event queue which has a list of all the moving
  objects in the moving object table.
*/
template <class Traits>
struct Trivial_kds: CGAL::Kinetic::Ref_counted<Trivial_kds<Traits> >
{
  typedef Trivial_kds<Traits> This;
  typedef typename Traits::Active_points_1_table::Data Point;
  typedef typename Traits::Simulator::Time Time;
  typedef typename Traits::Active_points_1_table::Key Point_key;
  typedef typename Traits::Simulator::Event_key Event_key;
  CGAL_KINETIC_DECLARE_LISTENERS(typename Traits::Simulator, typename Traits::Active_points_1_table)
public:
  typedef Trivial_event<Point, Time, This> Event;

  Trivial_kds(Traits tr): has_certificates_(true),
			  tr_(tr){
    CGAL_KINETIC_INITIALIZE_LISTENERS(tr_.simulator_handle(),
				      tr_.active_points_1_table_handle());
  }

  // this method is called with the value true when the event is processed
  void set_processed(bool tf) {
    if (tf== true) {
      event_= Event_key();
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void audit() const
  {
    /* In a real KDS you could use tr_.instantaneous_kernel() to build
       a static version of the kinetic data structure and check it for
       equality with the kinetic version.
    */
    CGAL_assertion(event_.is_valid());
    std::cout << "The structure is trivially correct at time "
	      << tr_.simulator_handle()->audit_time() << std::endl;
  }

  void set_has_certificates(bool tf) {
    typename Traits::Simulator::Handle sp= tr_.simulator_handle();
    if (has_certificates_ != tf) {
      has_certificates_=tf;
      if (has_certificates_) {
	bool ev= event_.is_valid();
	CGAL_assertion(!ev);
	Time t= CGAL::to_interval(sp->current_time()).second+1;
	event_= sp->new_event(t, Event(objects_.begin(),
				       objects_.end(),
				       *tr_.active_points_1_table_handle(),
				       this));
	std::cout << "Created event (" << event_ << ") at time " << t << std::endl;
      } else {
	if (event_.is_valid()) {
	  std::cout << "Deleting event " << event_ << std::endl;
	  sp->delete_event(event_);
	  event_=Event_key();
	}
      }
    }
  }

  bool has_certificates() const
  {
    // you can't use !event_ because the simulation
    // might end before the time you tried to schedule the event.
    return has_certificates_;
  }

  void insert(Point_key k) {
    std::cout << "Updating structure to include new object "
	      << k << "." << std::endl;
    objects_.insert(k);
    if (has_certificates_) {
      std::cout << "Updating all certificates which depend on "
                << k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void set(Point_key k) {
    if (has_certificates_) {
      std::cout << "Updating all certificates which depend on "
                << k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    } else {
      std::cout << "An object changed, but there was no certificate."<< std::endl;
    }
  }

  void erase(Point_key k) {
    std::cout << "An object " << k << " was removed."<< std::endl;
    objects_.erase(k);
    if (has_certificates_) {
      std::cout << "Updating all certificates which depend on "
                << k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

protected:
  bool has_certificates_;
  std::set<Point_key > objects_;
  Event_key event_;
  Traits tr_;
};

int main(int, char *[])
{
  typedef CGAL::Kinetic::Exact_simulation_traits Traits;
  typedef Trivial_kds<Traits> TKDS;

  Traits tr(1,100);
  TKDS::Handle tk= new TKDS(tr);

  Traits::Simulator::Handle sp=tr.simulator_handle();

  Traits::Simulator::Function_kernel::Construct_function cf
    = tr.kinetic_kernel_object().function_kernel_object().construct_function_object();

  tk->set_has_certificates(true);

  typedef Traits::Kinetic_kernel::Point_2 Point;
  Traits::Active_points_1_table::Key k
    = tr.active_points_1_table_handle()->insert(Point(cf(1), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.active_points_1_table_handle()
    ->set(k, Traits::Kinetic_kernel::Point_2(cf(2), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.active_points_1_table_handle()->erase(k);
  sp->set_current_event_number(sp->current_event_number()+10);
  return EXIT_SUCCESS;
}
