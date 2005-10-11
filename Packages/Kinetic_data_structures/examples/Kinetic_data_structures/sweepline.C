#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Cartesian_instantaneous_kernel.h>
#include <algorithm>
#include <map>
#include <list>
#include <iterator>
#include <iostream>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Insert_event.h>
#include <CGAL/KDS/Erase_event.h>
#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <CGAL/KDS/Simulator_kds_listener.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>


template <class Sort, class Id, class Solver> 
class Swap_event;

//! A simple KDS which maintains objects sorted by their x coordinate
/*!  This does not use Simple_kds_base for now irrelevant
  reasons. Ditto for the lack of protection of any of the fields. The
  code is designed to be read, so read it if you want to figure out
  what is going on.
*/
template <class TraitsT> class Planar_arrangement:
// for ref counted pointers
  public CGAL::KDS::Ref_counted<Planar_arrangement<TraitsT> >  {
  typedef TraitsT Traits;
  typedef Planar_arrangement<TraitsT> This;

  typedef std::pair<typename Traits::Moving_point_table::Key, int> Cut_pair;

  // this is used to identify pairs of objects in the list
  typedef typename std::list<Cut_pair>::iterator iterator;

  typedef Swap_event<This,iterator,
		     typename Traits::Simulator::Root_stack> Event;

  
  
  typedef std::pair<int,int> Edge;

  struct Compare_first {
    Compare_first(typename Traits::Instantaneous_kernel::Less_x_1 l): less_(l){}
    typename Traits::Instantaneous_kernel::Less_x_1 less_;
    bool operator()(const Cut_pair &a, const Cut_pair &b) const {
      return less_(a.first, b.first);
    }
  };

  typedef typename CGAL::KDS::
  Simulator_kds_listener<typename Traits::Simulator::Listener, 
			 This> Sim_listener;
  // Redirects the MovingObjectTable notifications to function calls
  typedef typename CGAL::KDS::
  Notifying_table_listener_helper<typename Traits::Moving_point_table::Listener,
				  This> MOT_listener;

public:
  typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 Approximate_point;

  // Register this KDS with the MovingObjectTable and the Simulator
  Planar_arrangement(Traits tr):tr_(tr),
				less_(tr.instantaneous_kernel_object().less_x_1_object()),
				siml_(tr.simulator_pointer(), this),
				motl_(tr.moving_point_table_pointer(), this){}


  Approximate_point vertex(int i) const {
    return approx_coords_[i];
  }

  size_t vertices_size() const {
    return approx_coords_.size();
  }

  typedef std::vector<Edge >::const_iterator Edges_iterator;
  Edges_iterator edges_begin() const {
    return edges_.begin();
  }
  Edges_iterator edges_end() const {
    return edges_.end();
  }

  void write(std::ostream &out) const {
    for (typename std::list<Cut_pair>::const_iterator it
	   = sorted_.begin(); it != sorted_.end(); ++it){
      out << it->first << " "; 
    }
    out << std::endl;
  }


  /* Insert k and update the affected certificates. std::upper_bound
     returns the first place where an item can be inserted in a sorted
     list. Called by the MOT_listener. This adds an vertex for the
     start of the curve.*/
  void insert(typename Traits::Moving_point_table::Key k) {
    //std::cout << "Inserting " << k <<std::endl;
    tr_.instantaneous_kernel_object().set_time(tr_.simulator_pointer()->rational_current_time());
    Cut_pair cp(k, new_point(k));
    iterator it = std::upper_bound(sorted_.begin(), sorted_.end(),
				   cp,less_);

    std::cout << "Curve " << k << " starts at " 
	      << tr_.simulator_pointer()->rational_current_time();
    if (it != sorted_.end()) {
      std::cout << " below curve " << next(it)->first;
    }
    std::cout << std::endl;

    sorted_.insert(it, cp);
    rebuild_certificate(--it); rebuild_certificate(--it);
    
  }



  /* Rebuild the certificate for the pair of points *it and *(++it).
     If there is a previous certificate there, deschedule it.*/
  void rebuild_certificate(const iterator it) {
    typename Traits::Simulator::Pointer sp= tr_.simulator_pointer();
    typename Traits::Moving_point_table::Pointer mp= tr_.moving_point_table_pointer();
    
    if (it == sorted_.end()) return;
    if (events_.find(it->first) != events_.end()) {
      tr_.simulator_pointer()->delete_event(events_[it->first]); events_.erase(it->first);
    }
    assert(find(it->first) == it);
    iterator n= next(it);
    if (n== sorted_.end()) return;     
    typename Traits::Kinetic_kernel::Less_x_1 less= tr_.kinetic_kernel_object().less_x_1_object();
    typename Traits::Simulator::Root_stack s 
      = sp->root_stack_object(less(mp->at(it->first), 
					mp->at(n->first)));
    // the Simulator will detect if the failure time is at infinity
    if (!s.empty()) {
      typename Traits::Simulator::Time t= s.top(); 
      s.pop();
      Event e(it, this,s);
      events_[it->first]= sp->new_event(t, e);
    } else events_[it->first]= sp->null_event();
  }

  /* Swap the pair of objects with *it as the first element.  The old
     solver is used to compute the next root between the two points
     being swapped. This method is called by an Event object. This
     addes a new vertex and two new edges.*/
  void swap(const typename Traits::Simulator::Time &t, iterator it,
	    typename Traits::Simulator::Root_stack &s) {
    std::cout << "Curves " << it->first << " and " << next(it)->first 
	      << " intersect at " << t << std::endl;
    typename Traits::Simulator::Pointer sp= tr_.simulator_pointer();
    assert(find(it->first) == it);
    events_.erase(it->first);
    iterator n= next(it);
    if (*n!= sorted_.back()) sp->delete_event(events_[n->first]);
    events_.erase(next(it)->first);

    int swap_point= new_point(it->first);
    edges_.push_back(Edge(swap_point, it->second));
    edges_.push_back(Edge(swap_point, next(it)->second));
    it->second= swap_point;
    next(it)->second=swap_point;
    
    

    std::swap(*it, *next(it));
    rebuild_certificate(next(it));
    if (!s.empty()){
      typename Traits::Simulator::Time t= s.top(); s.pop(); 
      events_[it->first]= sp->new_event(t, Event(it, this,s));
    } else events_[it->first]= sp->null_event();
    if (it != sorted_.begin()) rebuild_certificate(--it);
  }

  /* Verify the structure by checking that the current coordinates are
     properly sorted for time t. This function is called by the Sim_listener.*/
  void audit() const {
    if (sorted_.size() <2) return;
    typename Traits::Simulator::Const_pointer sp= tr_.simulator_pointer();
    tr_.instantaneous_kernel_object().set_time(sp->rational_current_time());
    //typename Instantaneous_kernel::Less_x_2 less= kernel_i_.less_x_2_object();
    for (typename std::list<Cut_pair>::const_iterator it
	   = sorted_.begin(); it->first != sorted_.back().first; ++it){
      if (!less_(*next(it), *it)){
	write(std::cerr);
      }
    }
  }

  /* Update the certificates adjacent to object k. This method is called by
     the MOT_listener. std::equal_range finds all items equal 
     to a key in a sorted list (there can only be one).*/
  void set(typename Traits::Moving_point_table::Key k) {
    iterator it =  find(k);
    rebuild_certificate(it); rebuild_certificate(--it);
  }

  /* Remove object k and destroy 2 certificates and create one new one.
     This function is called by the MOT_listener. One new vertex and one new edge are added.*/
  void erase(typename Traits::Moving_point_table::Key k) {
    std::cout << "Curve " << k << " ends at " 
	      << tr_.simulator_pointer()->rational_current_time() << std::endl;
    iterator it =  find(k);
    int lastp= it->second;
    iterator p= it; --p;
    tr_.simulator_pointer()->delete_event(events_[k]);
    sorted_.erase(it);
    rebuild_certificate(p);
    events_.erase(k);
       
    edges_.push_back(Edge(lastp, new_point(k)));
  }

  template <class It> static It next(It it){ return ++it;}

 
  iterator find(typename Traits::Moving_point_table::Key k) {
    iterator it = sorted_.begin();
    while (it != sorted_.end()){
      if (it->first ==k) return it;
      ++it;
    }
    assert(0);
    return it;
  }

  int new_point(typename Traits::Moving_point_table::Key k) {
    double tv= CGAL::to_double(tr_.simulator_pointer()->current_time());
    double dv= CGAL::to_double(tr_.moving_point_table_pointer()->at(k).x()(tv));
    approx_coords_.push_back(Approximate_point(tv, dv));
    return approx_coords_.size()-1;
  }

  typedef typename std::list<typename Traits::Moving_point_table::Key>::const_iterator Key_iterator;
  Key_iterator begin() const {
    return sorted_.begin();
  }
  Key_iterator end() const {
    return sorted_.end();
  }
  ~Planar_arrangement() {
    for (typename std::map<typename Traits::Moving_point_table::Key,
	   typename Traits::Simulator::Event_key >::iterator it= events_.begin(); 
	 it != events_.end(); ++it){
      tr_.simulator_pointer()->delete_event(it->second);
    }
  }

  Traits tr_;
  // The points in sorted order
  std::list<Cut_pair > sorted_;
  // events_[k] is the certificates between k and the object after it
  std::map<typename Traits::Moving_point_table::Key,
	   typename Traits::Simulator::Event_key > events_;
  std::vector<Edge > edges_;
  std::vector<Approximate_point > approx_coords_;
  Compare_first less_;
  Sim_listener siml_; 
  MOT_listener motl_; 
};

/* It needs to implement the time() and process() functions and 
   operator<< */
template <class Planar_arrangement, class Id, class Solver> 
class Swap_event {
public:
  Swap_event(Id o, Planar_arrangement* sorter, 
	     const Solver &s): left_object_(o), sorter_(sorter), s_(s){}
  void process(const typename Solver::Root &t){
    sorter_->swap(t, left_object_, s_);
  }
  Id left_object_; Planar_arrangement* sorter_; Solver s_;
};
template <class S, class I, class SS>
std::ostream &operator<<(std::ostream &out,
			 const Swap_event<S,I,SS> &ev){
  return out << "swap " << ev.left_object_->first ;
}

double snap(double v){
  double ival= static_cast<int>((v*1000))/1000.0;
  if (ival > 10) return 20;
  else if (ival <-10) return 0;
  else return ival+10;
}


int main(int, char *[]){
  typedef CGAL::KDS::Exact_simulation_traits_1 Traits;
  typedef Traits::Kinetic_kernel::Point_1 Moving_point;
  typedef Traits::Simulator::Time Time;
  typedef CGAL::KDS::Insert_event<Traits::Moving_point_table> Insert_event;
  typedef CGAL::KDS::Erase_event<Traits::Moving_point_table> Erase_event;
  typedef Planar_arrangement<Traits> Sort;

  Traits tr;
  Sort sort(tr);
  
  typedef Traits::NT NT;
  typedef std::pair<NT, NT> Ival;
  typedef std::pair<Ival,  Traits::Function_kernel::Function> Cv;
  typedef Traits::Function_kernel::Construct_function CF;
  std::vector<Cv> curves;
  CF fc= tr.function_kernel_object().construct_function_object();

  curves.push_back(Cv(Ival(NT(0), NT(9.9)), fc(NT(0), NT(.9), NT(-.05))));
  curves.push_back(Cv(Ival(NT(.1), NT(8)), fc(NT(3),NT(0), NT(.05))));
  curves.push_back(Cv(Ival(NT(.2), NT(7)), fc(NT(2),NT(.8),NT(-.06))));
  curves.push_back(Cv(Ival(NT(.5), NT(8)), fc(NT(6),NT(.2),NT(-.03))));
  curves.push_back(Cv(Ival(NT(1), NT(9.9)), fc(NT(2.3),NT(-.03))));
  curves.push_back(Cv(Ival(NT(3), NT(5)), fc(NT(1),NT(.7),NT( -.04))));
  curves.push_back(Cv(Ival(NT(6), NT(9)), fc(NT(-9),NT(.5),NT(-.013), NT(.0012), NT(-.00040))));

  Traits::Simulator::Pointer sp= tr.simulator_pointer();

  for (unsigned int i=0; i< curves.size(); ++i){
    tr.simulator_pointer()->new_event(Time(curves[i].first.first), 
				      Insert_event(Moving_point(curves[i].second),
						   tr.moving_point_table_pointer()));
    tr.simulator_pointer()->new_event(Time(curves[i].first.second), 
				      Erase_event(Traits::Moving_point_table::Key(i), 
						  tr.moving_point_table_pointer()));
  }

  while (sp->next_event_time() < sp->end_time()) {
    sp->set_current_event_number(sp->current_event_number()+1);
  }

  // output to metapost
  std::ofstream mpostfile("arrangement.mp");

  mpostfile << "beginfig(0)\n";
  mpostfile << "u:=20.0pt;\n";
  mpostfile << "pickup pencircle scaled 0.05u;\n";
  mpostfile << "%" << std::ios::fixed<< std::endl;

  for (Sort::Edges_iterator it = sort.edges_begin(); it != sort.edges_end(); ++it){
    mpostfile << "draw(" << snap(sort.vertex(it->first).x()) << "u, " 
	      << snap(sort.vertex(it->first).y()) << "u)";
    mpostfile << "--(" << snap(sort.vertex(it->second).x()) << "u, " 
	      << snap(sort.vertex(it->second).y()) << "u) withcolor black;\n";
  }

  for (unsigned int i=0; i< curves.size(); ++i){
    mpostfile << "draw(";
    for (double t= CGAL::to_double(curves[i].first.first); t < CGAL::to_double(curves[i].first.second); t += .01){
      double val= snap(CGAL::to_double(curves[i].second(t)));

      mpostfile << snap(t) << "u, " << val << "u)\n--(";
      
    }
    mpostfile << snap(CGAL::to_double(curves[i].first.second)) << "u, "
	      << snap(CGAL::to_double(curves[i].second(curves[i].first.second))) 
	      << "u) withcolor red;\n";
  }
  
  //for (double t =0; t< 
  
  mpostfile << "endfig\n";
  return EXIT_SUCCESS;
}
