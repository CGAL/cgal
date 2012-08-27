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

#ifndef CGAL_KINETIC_ENCLOSING_BOX_3_H
#define CGAL_KINETIC_ENCLOSING_BOX_3_H
#include <CGAL/basic.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/listeners.h>
#include <CGAL/Kinetic/Enclosing_box_3.h>
#include <CGAL/Kinetic/Event_base.h>

namespace CGAL { namespace Kinetic {

template <class EB3>
class Enclosing_box_bounce_event_3: public Event_base<EB3*>
{
public:
  Enclosing_box_bounce_event_3(){}
  Enclosing_box_bounce_event_3(EB3* eb,
			       typename EB3::Side s,
			       typename EB3::Point_key k,
			       typename EB3::NT t): eb_(eb),
						    k_(k),
						    t_(t),
						    s_(s) {

  }
  void process() {
    eb_->bounce(k_, t_, s_);
  }
  std::ostream& write(std::ostream &out) const
  {
    out << "Bounce " << k_ << " off " << s_ << " at " << t_;
    return out;
  }
  EB3* eb_;
  typename EB3::Point_key k_;
  typename EB3::NT t_;
  typename EB3::Side s_;
};

template <class EB3>
std::ostream &operator<<(std::ostream &out, const Enclosing_box_bounce_event_3<EB3> &eb)
{
  eb.write(out);
  return out;
}


template <class Traits>
class  Enclosing_box_3: public Ref_counted<Enclosing_box_3<Traits> >
{

  typedef Enclosing_box_3<Traits> This;
  typedef typename Traits::Simulator Simulator;
  typedef typename Traits::Kinetic_kernel Kinetic_kernel;
  typedef typename Traits::Active_points_3_table Active_points_3_table;

  CGAL_KINETIC_DECLARE_AOT_LISTENER(typename Active_points_3_table);

  typedef typename Simulator::Event_key Event_key;
  typedef typename Simulator::Time Time;

  typedef Enclosing_box_bounce_event_3<This> Event;
  friend class Enclosing_box_bounce_event_3<This>;
  typedef typename Simulator::Function_kernel::Function Fn;
  typedef Fn Coordinate;

  template <class It>
  Coordinate make_coordinate(It b, It e) const {
    return Coordinate(b,e);
  }
public:
  enum Side {TOP=0, BOTTOM=1, LEFT=2, RIGHT=3, FRONT=4, BACK=5};

  typedef typename Active_points_3_table::Data Point;
  typedef typename Active_points_3_table::Key Point_key;

  typedef typename Fn::NT NT;
  //typedef double NT;
  Enclosing_box_3(Traits tr, NT xmin=-10, NT xmax=10, NT ymin=-10, NT ymax=10, NT zmin=-10, NT zmax=10):traits_(tr) {
    CGAL_assertion(xmin<xmax);
    CGAL_assertion(ymin<ymax);
    CGAL_assertion(zmin<zmax);
    bounds_[LEFT]=xmin;
    bounds_[RIGHT]=xmax;
    bounds_[TOP]=ymax;
    bounds_[BOTTOM]=ymin;
    bounds_[FRONT]=zmax;
    bounds_[BACK]=zmin;
    CGAL_KINETIC_INITIALIZE_AOT_LISTENER(tr.active_points_3_table_handle());
  };

  ~Enclosing_box_3() {
    for (typename std::map<Point_key, Event_key>::iterator it= certs_.begin(); it!= certs_.end(); ++it) {
      traits_.simulator_handle()->delete_event(it->second);
    }
  }

  void set(Point_key k) {
    erase(k);
    insert(k);
  }

  void insert(Point_key k) {
    double tb=std::numeric_limits<double>::infinity();
    Side bs= try_bound(LEFT, k, bs, tb);
    bs= try_bound(RIGHT, k, bs, tb);
    bs= try_bound(TOP, k, bs, tb);
    bs= try_bound(BOTTOM, k, bs, tb);
    bs= try_bound(FRONT, k, bs, tb);
    bs= try_bound(BACK, k, bs, tb);
    if (tb != std::numeric_limits<double>::infinity()) {
      certs_[k]= traits_.simulator_handle()->new_event(tb, Event(this,bs,k,tb));
      //std::cout << certs_[k] << std::endl;
    }
    /*std::cout << "Scheduled event for point " << k << " with motion " << traits_.active_points_3_table_handle()->at(k)
      << " for time " << tb << " on wall " << bs << std::endl;*/
  }

  void erase(Point_key k) {
    if (certs_.find(k) != certs_.end()) {
      traits_.simulator_handle()->delete_event(certs_[k]);
      certs_.erase(k);
    }
  }

  void bounce(Point_key k, NT time, Side s) {
    CGAL_LOG(Log::LOTS, "Bouncing " << k << " off side " << s << std::endl);
    certs_.erase(k);
    std::vector<NT> coefs[3];
    if (s==TOP || s== BOTTOM) {
      coefs[0].insert(coefs[0].end(),
		      traits_.active_points_3_table_handle()->at(k).x().begin(),
		      traits_.active_points_3_table_handle()->at(k).x().end());
      compute_bounce(traits_.active_points_3_table_handle()->at(k).y(),time, coefs[1]);
      coefs[2].insert(coefs[2].end(),
		      traits_.active_points_3_table_handle()->at(k).z().begin(),
		      traits_.active_points_3_table_handle()->at(k).z().end());
    }
    else if (s==LEFT || s == RIGHT) {
      compute_bounce(traits_.active_points_3_table_handle()->at(k).x(),time, coefs[0]);
      coefs[1].insert(coefs[1].end(),
		      traits_.active_points_3_table_handle()->at(k).y().begin(),
		      traits_.active_points_3_table_handle()->at(k).y().end());
      coefs[2].insert(coefs[2].end(),
		      traits_.active_points_3_table_handle()->at(k).z().begin(),
		      traits_.active_points_3_table_handle()->at(k).z().end());
    }
    else {
      coefs[0].insert(coefs[0].end(),
		      traits_.active_points_3_table_handle()->at(k).x().begin(),
		      traits_.active_points_3_table_handle()->at(k).x().end());
      coefs[1].insert(coefs[1].end(),
		      traits_.active_points_3_table_handle()->at(k).y().begin(),
		      traits_.active_points_3_table_handle()->at(k).y().end());
      compute_bounce(traits_.active_points_3_table_handle()->at(k).z(),time, coefs[2]);
    }

    /*typename Traits::Simulator::Function_kernel::Create_function cf
      = traits_.simulator_handle()->function_kernel().create_function_object();*/

    Point pt(make_coordinate(coefs[0].begin(), coefs[0].end()),
	     make_coordinate(coefs[1].begin(), coefs[1].end()),
	     make_coordinate(coefs[2].begin(), coefs[2].end()));
    //std::cout << "Changing motion from " << traits_.active_points_3_table_handle()->at(k) << " to " << pt << std::endl;
    traits_.active_points_3_table_handle()->set(k, pt);
    // CGAL_assertion(traits_.active_points_3_table_handle()->at(k) == pt);
  }

protected:


  typename Simulator::Function_kernel function_kernel_object() const
  {
    return traits_.kinetic_kernel_object().function_kernel_object();
  }

  Side try_bound(Side try_side, Point_key k,Side old_side,  double& old_time) const
  {
    Coordinate nf;
    NT bound=bounds_[try_side];
   typename Kinetic_kernel::Certificate re;
    if (try_side== TOP || try_side==BOTTOM) {
      typename Kinetic_kernel::Compare_y_3 ily = traits_.kinetic_kernel_object().compare_y_3_object();     
      if (try_side== BOTTOM) {
	re= ily(traits_.active_points_3_table_handle()->at(k), bound,
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      } else {
	re= ily(bound, traits_.active_points_3_table_handle()->at(k),
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      }

    } else if (try_side == LEFT || try_side == RIGHT) {
      typename Kinetic_kernel::Compare_x_3 ily = traits_.kinetic_kernel_object().compare_x_3_object();
      if (try_side== LEFT) {
	re= ily(traits_.active_points_3_table_handle()->at(k), bound,
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      } else {
	re= ily(bound, traits_.active_points_3_table_handle()->at(k),
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      }
    } else {
      typename Kinetic_kernel::Compare_z_3 ily = traits_.kinetic_kernel_object().compare_z_3_object();
      if (try_side== BACK) {
	re= ily(traits_.active_points_3_table_handle()->at(k), bound,
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      } else {
	re= ily(bound, traits_.active_points_3_table_handle()->at(k),
		traits_.simulator_handle()->current_time(), traits_.simulator_handle()->end_time());
      }
    }
    if (re.will_fail()) {
      CGAL_LOG(Log::LOTS, "Side fails at " << re.failure_time() << std::endl);
      double dv= CGAL::to_interval(re.failure_time()).first;
      if (dv < old_time) {
	old_time=dv;
	return try_side;
      } else {
	return old_side;
      }
    } else {
      return old_side;
    }
  }

  void compute_bounce(const Coordinate& f, NT t, std::vector<NT> &out) {
    // x is contant
    // v is negative v
    // higher order coefs on constant
    // out(time)=f(time)
    // out'(time)= -f'(time)
    typename Simulator::Function_kernel::Differentiate cd
      = function_kernel_object().differentiate_object();

    if (f.degree() >=2) {
      std::vector<NT> hcoefs(f.begin(), f.end());
      hcoefs[0]=0;
      hcoefs[1]=0;
      Coordinate fh(hcoefs.begin(), hcoefs.end());
      Coordinate dfh= cd(fh);
      out.push_back(f[0]+2*f[1]*t+2*t*dfh(t));
      out.push_back(-f[1]-2*dfh(t));
      out.insert(out.end(), f.begin()+2, f.end());
    }
    else {
      NT v= -cd(f)(t);
      NT x= f(t);
      out.push_back(x-v*t);
      out.push_back(v);
      //out.push_back(x);
    }
    /*{
      Function ft(out.begin(), out.end());
      NT nt= ft(t);
      NT ot= f(t);
      NT nd= cd(ft)(t);
      NT od= cd(f)(t);
      }*/

    CGAL_exactness_assertion_code(Coordinate ft(out.begin(), out.end()););
    CGAL_exactness_assertion(ft(t) == f(t));
    CGAL_exactness_assertion(cd(ft)(t) == -cd(f)(t));
  }

  NT bounds_[6];
  Traits traits_;
  std::map<Point_key, Event_key> certs_;
};

} } //namespace CGAL::Kinetic
#endif
