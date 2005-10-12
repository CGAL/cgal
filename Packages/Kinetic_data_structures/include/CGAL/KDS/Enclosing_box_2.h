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

#ifndef CGAL_KDS_ENCLOSING_BOX_2_H
#define CGAL_KDS_ENCLOSING_BOX_2_H
#include <CGAL/basic.h>
#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>

CGAL_KDS_BEGIN_NAMESPACE

template <class EB2>
class Enclosing_box_bounce_event_2 {
public:
  Enclosing_box_bounce_event_2(){}
  Enclosing_box_bounce_event_2(EB2* eb,
			       typename EB2::Side s,
			       typename EB2::Point_key k,
			       typename EB2::NT t): eb_(eb),
						      k_(k),
						      t_(t),
						      s_(s){

  }
  void process(const typename EB2::Simulator::Time &){
    eb_->bounce(k_, t_, s_);
  }
  void write(std::ostream &out) const {
    out << "[Bounce " << k_ << " off " << s_ << " at " << t_ <<"]";
  }
  EB2* eb_;
  typename EB2::Point_key k_;
  typename EB2::NT t_;
  typename EB2::Side s_;
};

template <class EB2>
std::ostream &operator<<(std::ostream &out, const Enclosing_box_bounce_event_2<EB2> &eb) {
  eb.write(out);
  return out;
}


template <class Traits>
class  Enclosing_box_2: public Ref_counted<Enclosing_box_2<Traits> > {

  typedef Enclosing_box_2<Traits> This;
  typedef typename Traits::Simulator Simulator;
  typedef typename Traits::Moving_point_table Moving_point_table;

  typedef typename CGAL::KDS::Simulator_kds_listener<typename Simulator::Listener, This> Simulator_listener;
  friend  class CGAL::KDS::Simulator_kds_listener<typename Simulator::Listener, This>;
  typedef typename CGAL::KDS::Notifying_table_listener_helper<typename Moving_point_table::Listener, This> Moving_point_table_listener;
  friend class CGAL::KDS::Notifying_table_listener_helper<typename Moving_point_table::Listener, This>;

  typedef typename Simulator::Event_key Event_key;
  typedef typename Simulator::Time Time;

  typedef Enclosing_box_bounce_event_2<This> Event;
  friend class Enclosing_box_bounce_event_2<This>;
  typedef typename Simulator::Function_kernel::Function Function;
public:
  enum Side {TOP=0, BOTTOM=1, LEFT=2, RIGHT=3};

  typedef typename Moving_point_table::Data Point;
  typedef typename Moving_point_table::Key Point_key;

  typedef typename Traits::Simulator::NT NT;
  //typedef double NT;
  Enclosing_box_2( Traits tr, NT xmin=-10, NT xmax=10, NT ymin=-10, NT ymax=10):traits_(tr), 
								 motl_(tr.moving_point_table_pointer(), this) {
    CGAL_assertion(xmin<xmax);
    CGAL_assertion(ymin<ymax);
    bounds_[LEFT]=xmin;
    bounds_[RIGHT]=xmax;
    bounds_[TOP]=ymin;
    bounds_[BOTTOM]=ymax;
  };
  
  ~Enclosing_box_2() {
    for (typename std::map<Point_key, Event_key>::iterator it= certs_.begin(); it!= certs_.end(); ++it){
      traits_.simulator_pointer()->delete_event(it->second);
    }
  }

  void set(Point_key k){
    erase(k);
    insert(k);
  }

  void insert(Point_key k) {
    double tb=std::numeric_limits<double>::infinity();
    Side bs= try_bound(LEFT, k, bs, tb);
    bs= try_bound(RIGHT, k, bs, tb);
    bs= try_bound(TOP, k, bs, tb);
    bs= try_bound(BOTTOM, k, bs, tb);
    if (tb != std::numeric_limits<double>::infinity()){
      certs_[k]= traits_.simulator_pointer()->new_event(tb, Event(this,bs,k,tb));
      //std::cout << certs_[k] << std::endl;
    }
    /*std::cout << "Scheduled event for point " << k << " with motion " << traits_.moving_point_table_pointer()->at(k) 
	      << " for time " << tb << " on wall " << bs << std::endl;*/
  }

  void erase(Point_key k) {
    if (certs_.find(k) != certs_.end()) {
      traits_.simulator_pointer()->delete_event(certs_[k]);
      certs_.erase(k);
    }
  }

  void bounce(Point_key k, NT time, Side s) {
    certs_.erase(k);
    std::vector<NT> coefs[2];
    if (s==TOP || s== BOTTOM){
      coefs[0].insert(coefs[0].end(), 
		      traits_.moving_point_table_pointer()->at(k).x().begin(),
		      traits_.moving_point_table_pointer()->at(k).x().end());
      compute_bounce(traits_.moving_point_table_pointer()->at(k).y(),time, coefs[1]);
    } else {
      compute_bounce(traits_.moving_point_table_pointer()->at(k).x(),time, coefs[0]);
      coefs[1].insert(coefs[1].end(), 
		      traits_.moving_point_table_pointer()->at(k).y().begin(),
		      traits_.moving_point_table_pointer()->at(k).y().end());
    }
    
    /*typename Traits::Simulator::Function_kernel::Create_function cf
      = traits_.simulator_pointer()->function_kernel().create_function_object();*/

    Point pt(Function(coefs[0].begin(), coefs[0].end()),
	     Function(coefs[1].begin(), coefs[1].end()));
    //std::cout << "Changing motion from " << traits_.moving_point_table_pointer()->at(k) << " to " << pt << std::endl;
    traits_.moving_point_table_pointer()->set(k, pt);
    CGAL_assertion(traits_.moving_point_table_pointer()->at(k) == pt);
  }
  
protected:

  typename Simulator::Function_kernel function_kernel_object() const {
    return traits_.simulator_pointer()->function_kernel_object();
  }

  Side try_bound(Side try_side, Point_key k,Side old_side,  double& old_time) const {
    Function nf;
    NT bound=bounds_[try_side];
    if (try_side== TOP || try_side==BOTTOM){
      nf=traits_.moving_point_table_pointer()->at(k).y()-bound;
    } else {
      nf=traits_.moving_point_table_pointer()->at(k).x()-bound;
    }
    if (try_side == BOTTOM || try_side == RIGHT) {
      nf=-nf;
    }

    typename Simulator::Root_stack re
      = traits_.simulator_pointer()->root_stack_object(nf);

    double dv = std::numeric_limits<double>::infinity();
    if (!re.empty()) {
      typename Simulator::Time rec=re.top();
      dv= CGAL::to_interval(rec).first;
    }
    
    /*while (!re.finished()) {
      CGAL_assertion(!function_kernel_object().is_even_multiplicity_object(nf)(re.current()));
      dv= CGAL::to_interval(re.current()).second;
      if (!re.finished()) {
	re.advance();
	if (!re.finished()){
	  CGAL_assertion(!function_kernel_object().is_even_multiplicity_object(nf)(re.current()));
	  if (dv < CGAL::to_interval(re.current()).first) {
	    break;
	  } else {
	    re.advance();
	  }
	}
      }
      }*/
    if (dv < old_time) {
      old_time=dv;
      return try_side;
    } else {
      return old_side;
    }
  }

  void compute_bounce(const Function& f, NT time, std::vector<NT> &out) {
    // out(time)=f(time)
    // out'(time)= -f'(time)
    typename Simulator::Function_kernel::Differentiate cd
      = function_kernel_object().differentiate_object();
    NT v= -cd(f)(time);
    NT x= f(time);
    if (f.degree() >=2) {
      Function fh= Function(f.begin()+2, f.end());
      NT fhv= fh(time);
      out.push_back(x-v*time-fhv);
      out.push_back(v-fhv);
      out.insert(out.end(), f.begin()+2, f.end());
    } else {
      out.push_back(x-v*time);
      out.push_back(v);
      //out.push_back(x);
    }
    CGAL_exactness_assertion_code(Function ft(out.begin(), out.end()););
    CGAL_exactness_assertion(ft(time) == f(time));
  }
  

  NT bounds_[4];
  Traits traits_;
  std::map<Point_key, Event_key> certs_;
  Moving_point_table_listener motl_;

};

CGAL_KDS_END_NAMESPACE

#endif
