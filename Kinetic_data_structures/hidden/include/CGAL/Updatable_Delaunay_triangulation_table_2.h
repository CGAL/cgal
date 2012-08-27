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

#ifndef CGAL_KINETIC_UPDATE_DEL_NOTIFYING_TABLE_BASE_3_H
#define CGAL_KINETIC_UPDATE_DEL_NOTIFYING_TABLE_BASE_3_H
#include <CGAL/basic.h>
#include <CGAL/Tools/Label.h>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Tools/Counter.h>
#include <CGAL/Kinetic/Multi_listener.h>
#include <boost/dynamic_bitset.hpp>
#include <iostream>
#include <vector>
#include <sstream>
#include <map>
#include <CGAL/protected_arrays.h>

#ifndef NDEBUG
#include <CGAL/Polynomial/Interval_polynomial.h>
#endif

namespace CGAL {;


bool disable_filter_0_=false;
bool disable_filter_1_=false;
bool disable_filter_2_=false;
bool disable_filter_3_=false;


unsigned int stat_interval_certificate_functions_;
unsigned int stat_interval_certificate_function_evaluations_;
unsigned int stat_interval_certificate_function_interval_evaluations_;
unsigned int stat_deriv_filter_calls_;
unsigned int stat_newton_isolate_calls_;
unsigned int stat_newton_refine_calls_;
unsigned int stat_exact_points_;
unsigned int stat_exact_certificate_functions_;
unsigned int stat_exact_certificate_functions_from_compare_;
unsigned int stat_exact_certificate_functions_from_compare_curt_;
unsigned int stat_exact_certificate_functions_from_advance_;
unsigned int stat_unfailing_exact_certificate_functions_;
unsigned int stat_exact_certificate_function_advances_;
unsigned int stat_interval_predicate_evaluations_;
unsigned int stat_point_predicate_evaluations_;
unsigned int stat_certificate_computations_;
unsigned int stat_certificate_advances_;
unsigned int stat_number_of_interpolations_;
unsigned int stat_number_of_vertices_;
unsigned int stat_number_of_activated_vertices_;
unsigned int stat_number_of_edges_=0;
unsigned int stat_number_of_bad_edges_;
double stat_total_interval_width_;
unsigned int stat_number_of_intervals_;

template <class UI>
struct URefiner {
  typedef typename UI::Cert_tuple Cert_tuple;
  typedef typename UI::Exact_certificate Exact_certificate;
  typedef typename UI::Exact_time Exact_root;
  typedef CGAL::Interval_nt_advanced INT;
  struct Cs: public CGAL::Kinetic::Ref_counted<Cs> {
    Cs(const Exact_certificate &c): cert_(c){}
    Exact_certificate cert_;
  };
  
  // enum State {INTERVAL, EXACT, INVALID};
  
  URefiner(Cert_tuple t,
	   typename UI::Certificate_function_const_pointer certf,
	   typename UI::Certificate_derivitive_const_pointer certpf,
	   typename UI::Handle tbl): tuple_(t),
				     ui_(tbl),
				     certf_(certf)
				   , certpf_(certpf)
  {
  }
  URefiner() {
  }
  bool is_invalid() const {
    return tuple_.is_invalid();
  }
  
  /*This operator-() const {
    return *this;
    }*/
    
  bool refine(std::pair<double,double> &iv) const {
    CGAL_precondition(ui_ != typename UI::Handle());
    if (iv.first == iv.second) return false;
    if (has_exact_root()) return false;

    double dd= iv.second-iv.first;
    if (cert_ == typename Cs::Handle()) {
      CGAL::Protect_FPU_rounding<true> prot;

      INT niv=ui_->Newton_refine(certf_, certpf_, INT(iv));
      bool ret= (niv.inf() != iv.first || niv.sup() != iv.second);
      iv= std::pair<double,double>(niv.inf(), niv.sup());
      return ret;
    } else if (dd != 0) {
      // compute exact;
      CGAL_assertion(!has_exact_root());
      //++comparison_certificates_;
	
      compute_exact_root(iv);
      iv= CGAL::to_interval(exact_root());
      return true;
    } else {
      compute_exact_root(iv);
      return false;
    }
  }

  const Exact_root& 
  exact_root() const {
    return cert_->cert_.failure_time();
  }

  /*const Exact_root& 
    exact_root(std::pair<double,double> iv) const {
    compute_exact_root(iv);
    return exact_root();
    }*/

  void write(std::ostream &out) const {
    if (tuple_ == Cert_tuple()) {
    } else{
      out << tuple_;
      if (has_exact_root()) {
	out << " " << exact_root();
      }
    }
  }

  const Cert_tuple tuple() const {
    return tuple_;
  }

  bool equal_description(const URefiner &o) const {
    if (tuple_ == Cert_tuple()) return false;
    return tuple_== o.tuple_;
  }

  bool has_exact_root() const {
    return cert_!= typename Cs::Handle() ;
  }
    
  CGAL::Sign sign_at(INT t) const {
    return UI::sign_at(certf_, t);
  }

  bool ensure_exact_root(std::pair<double,double> iv) const {
    if (!has_exact_root()) {
      compute_exact_root(iv);
      return true;
    } else return false;
  }

  void compute_exact_root(std::pair<double,double> iv) const {
    CGAL_precondition(!has_exact_root());
    if (ui_ != typename UI::Handle()) {
      ++stat_exact_certificate_functions_;
      ++stat_exact_certificate_functions_from_compare_;
      CGAL_UD_DEBUG("Generating exact with interval " << iv.first << " " << iv.second << std::endl);
      cert_= new Cs(ui_->compute_exact_failure_time(tuple_, iv.first));
      CGAL_UD_DEBUG("Got " << cert_->cert_.failure_time() << std::endl);
      CGAL_assertion(check_.failure_time() == cert_->cert_.failure_time());
      //iv= CGAL::to_interval(exact_root());
    } else {
      CGAL_UD_DEBUG("Faking exact" << std::endl);
      double cs[2];
      cs[0]=iv.first;
      cs[1]=-1.0;
      typename Exact_certificate::Function f(cs, cs+2);
      // hack
      cert_= new Cs(Exact_certificate(f, typename UI::Kinetic_kernel::Function_kernel(), -1, 2));
      CGAL_UD_DEBUG(f << ": " << exact_root()  << std::endl); 
      CGAL_assertion(exact_root() == Exact_root(iv.first));
#ifndef NDEBUG
      check_= cert_->cert_;
#endif
    }
  }

  void advance_exact_root(const Exact_root &et) const {
    CGAL_precondition(has_exact_root());
    ++stat_exact_certificate_functions_;
    while (CGAL::compare(exact_root(), et) == CGAL::SMALLER) {
      CGAL_UD_DEBUG("Advancing exact root " << cert_->cert_ << std::endl);
      cert_->cert_.pop_failure_time();
    }
  }

  const Exact_certificate& exact_certificate() const {
    return cert_->cert_;
  }

  std::pair<double,double> interval_from_exact() const {
    CGAL_precondition(has_exact_root());
    if (exact_root() >= Exact_root(1)) {
      return std::make_pair(1.0, -1.0);
    } else {
      std::pair<double,double> ip= CGAL::to_interval(exact_root());
      return ip;
    }
  }

  void set_exact_certificate(const Exact_certificate& ec) {
    CGAL_precondition(!has_exact_root());
    cert_= new Cs(ec);
    CGAL_postcondition(exact_root() == check_.failure_time());
  }

  typename UI::Certificate_function_const_pointer certificate_function() const {
    return certf_;
  }

  typename UI::Certificate_derivitive_const_pointer certificate_derivitive() const {
    return certpf_;
  }

  Cert_tuple tuple_;
  typename UI::Handle ui_;
  mutable typename Cs::Handle cert_;
  typename UI::Certificate_function certf_;
  typename UI::Certificate_derivitive certpf_;
#ifndef NDEBUG
  mutable Exact_certificate check_;
#endif
};


template <class Point_key>
struct Update_cert_tuple {
  Update_cert_tuple(){}
  Update_cert_tuple(Point_key k[4]) {
    std::copy(k, k+4, k_);
    int swaps=0;
    for (unsigned int i=0; i< 3; ++i){
      for (unsigned int j=i+1; j< 4; ++j) {
	if (k_[i] > k_[j]) {
	  std::swap(k_[i], k_[j]);
	  ++swaps;
	}
      }
    }
    
    if (swaps%2 ==1) {
      std::swap(k_[2], k_[3]);
    }
  }

  /*static std::pair<Cert_tuple, bool> make(Point_key k[4]) {
    Cert_tuple r(k);
    int swaps=0;
    for (unsigned int i=0; i< 3; ++i){
    for (unsigned int j=i+1; j< 4; ++j) {
    if (r.k_[i] > r.k_[j]) {
    std::swap(r.k_[i], r.k_[j]);
    ++swaps;
    }
    }
    }
    return std::make_pair(r, swaps%2==0);
    }*/

  Update_cert_tuple canonicalize() const {
    Update_cert_tuple ret= *this;
    if (ret.k_[2] > ret.k_[3]) {
      std::swap(ret.k_[2], ret.k_[3]);
    }
    return ret;
  }
  Update_cert_tuple opposite() const {
    Update_cert_tuple ret= *this;
    std::swap(ret.k_[2], ret.k_[3]);
    return ret;
  }

  bool operator<(const Update_cert_tuple o) const {
    for (unsigned int i=0; i< 4; ++i){
      if (k_[i] < o.k_[i]) return true;
      else if (k_[i] > o.k_[i]) return false;
    }
    return false;
  }
  bool operator==(const Update_cert_tuple o) const {
    for (unsigned int i=0; i< 4; ++i){
      if (k_[i] != o.k_[i]) return false;
    }
    return true;
  }
  bool is_invalid() const {
    return k_[0] == Point_key();
  }
  Point_key operator[](int i) const {
    return k_[i];
  }
  void write(std::ostream&out) const {
    out << k_[0] << " " << k_[1] << " " << k_[2] << " " << k_[3];
  }
  Point_key k_[4];
};

template <class K>
std::ostream &operator<<(std::ostream &out, 
			 const Update_cert_tuple<K> &ct) {
  out << ct[0] << " " << ct[1] << " " << ct[2] << " " << ct[3];
  return out;
};



template <class UI>
struct Interpolate_event: public CGAL::Kinetic::Free_event_base {
  typedef typename UI::Key Point_key;
  template<class SH>
  Interpolate_event(typename UI::Handle ui, SH sh): ui_(ui) {
    motions_.reserve(20);
    std::pair<double, double> dp=CGAL::to_interval(sh->current_time());
    if (dp.second != 0) {
      time_= nextafter(dp.second,
		       std::numeric_limits<double>::max());
    } else {
      time_= 0.0;
    }
    //}
  }
  Interpolate_event(typename UI::Handle ui, double t): ui_(ui), time_(t) {
  }
  
  void insert(Point_key k) {
    motions_.push_back(k);
    /*if (ic(*c) != fc(*c)) {
      motions_.push_back(*c);
      } else {
      CGAL_UD_DEBUG("Point " << *c << " doesn't move." << std::endl);
      }*/
  }
  
  double time() const {
    return time_;
  }
  
  std::ostream & write(std::ostream&out) const {
    out << "Update ";
    for (unsigned int i=0; i< motions_.size(); ++i){
      out << motions_[i] << " ";
    }
    return out;
  }
  
  void process() {
    std::sort(motions_.begin(), motions_.end());
    motions_.erase(std::unique(motions_.begin(), motions_.end()), motions_.end());
    
    ui_->set_is_activating(true);
    //INT it= CGAL::to_interval(time_);
    ui_->set_is_editing(UI::LOGGED);
    CGAL::Protect_FPU_rounding<true> prot;
    for (unsigned int i=0; i< motions_.size(); ++i) {
      ui_->preactivate(time_, motions_[i]);
      ui_->activate(time_, motions_[i]);
    }
    ui_->set_is_editing(UI::NOT);
    ui_->set_is_activating(false);
  }
  typename UI::Handle ui_;
  typename std::vector<Point_key> motions_;
  double time_;
};



/*!
 */
template <class Indirect_kernel_t, class Kinetic_kernel_t>
class Updatable_delaunay_triangulation_table_2:
  public Kinetic::Ref_counted<Updatable_delaunay_triangulation_table_2<Indirect_kernel_t, Kinetic_kernel_t> >
{
public:
  typedef Kinetic_kernel_t Kinetic_kernel;
  typedef Indirect_kernel_t Indirect_kernel;
  typedef Updatable_delaunay_triangulation_table_2<Indirect_kernel, Kinetic_kernel> This;
  typedef typename Kinetic_kernel::Point_2 Exact_point_2;
  typedef typename Indirect_kernel::Point_2 Key;
  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;
  typedef CGAL::Interval_nt_advanced INT;
  typedef typename Kinetic_kernel::Certificate Exact_certificate;
  typedef typename Exact_certificate::Time Exact_time;
  typedef Update_cert_tuple<Key> Cert_tuple;
  typedef Interpolate_event<This> IE;
 
  
  template <int D>
  struct FN: Protected_array<INT, D+1> {
    typedef Protected_array<INT, D+1> P;
    FN(){}
    FN(Protected_array_const_pointer<INT, D+1> cp): P(cp){}
    FN(Protected_array_const_pointer<INT, D+1> cp, bool){
      for (unsigned int i=0; i<= D; ++i) {
	P::operator[](i)= -cp[i];
      }
    }
    
    INT operator()(INT t) const {
      return evaluate_ipoly(Protected_array_const_pointer<INT, D+1>(*this), t);
    }
    double inf(INT t) const {
      return evaluate_ipoly_inf(Protected_array_const_pointer<INT, D+1>(*this), t);
    }
    FN<D-1> prime() const {
      FN<D-1> ret;
      for (unsigned int i=1; i<=D; ++i) {
	double di= i;
	ret[i-1]= INT(-(di* (-P::operator[](i).inf())), di *P::operator[](i).inf()) ;
      }
      return ret;
    }
    /*Function<D> flip() const {
      Function<D> ret;
      for (unsigned int i=0; i<=D; ++i) {
      ret[i]= -P::operator[](i);
      }
      return ret;
      }*/
  };

  typedef FN<4> Certificate_function;
  typedef FN<3> Certificate_derivitive;
  typedef FN<2> Certificate_acceleration;
  typedef Protected_array_pointer<INT, 5> Certificate_function_pointer;
  typedef Protected_array_pointer<INT, 4> Certificate_derivitive_pointer;
  typedef Protected_array_pointer<INT, 3> Certificate_acceleration_pointer;
  typedef Protected_array_const_pointer<INT, 5> Certificate_function_const_pointer;
  typedef Protected_array_const_pointer<INT, 4> Certificate_derivitive_const_pointer;
  typedef Protected_array_const_pointer<INT, 3> Certificate_acceleration_const_pointer;
 

  typedef URefiner<This> Refiner;
  typedef typename Kinetic::Interval_simulator_traits<Refiner> Simulator_traits;
  typedef typename Kinetic::Two_list_pointer_event_queue<Simulator_traits, false, 2> Queue;
  typedef typename Kinetic::Default_simulator<Simulator_traits, Queue > Simulator;
  typedef typename Simulator::Event_key Event_key;
protected:
  //! Convenience
  struct Listener_core
  {
    typedef enum {IS_EDITING}
      Notification_type;
    typedef typename This::Handle Notifier_handle;
  };
public:
  //! The interface to implement if you want to receive notifications.
  typedef Kinetic::Multi_listener<Listener_core> Listener;

protected:
  // This is evil here.
  typedef std::set<Listener*> Subscribers;

#ifndef MOVE_ALL
  struct Coef_data {
    Coef_data(double st, INT x_0, INT x_1, INT y_0, INT y_1){
      c_[0][0]=x_0;
      c_[0][1]=x_1;
      c_[1][0]=y_0;
      c_[1][1]=y_1;
      start_time_=st;
    }
    Coef_data(double x, double y){
      c_[0][0]=x;
      c_[0][1]=0;
      c_[1][0]=y;
      c_[1][1]=0;
      start_time_=-2;
    }
    Coef_data(){}
    Protected_array_const_pointer<INT,2> operator[](unsigned int i) const {
      CGAL_precondition(i <2);
      return Protected_array_const_pointer<INT,2>(c_[i]);
    }
    double start_time() const {
      return start_time_;
    }
    double start_time_;
    INT c_[2][2];
  };
#endif

  typedef boost::dynamic_bitset<> Active;
public:



  enum Editing_state {LOGGED, UNLOGGED, NOT};

  //! default constructor
  Updatable_delaunay_triangulation_table_2(Indirect_kernel ik,
                                           Kinetic_kernel kk):
    editing_(NOT), 
    notifying_(false),
    soc_(kk.positive_side_of_oriented_circle_2_object()),
    ik_(ik), 
#ifndef MOVE_ALL
    active_(ik.number_of_point_2s(), false),
    coef_cache_(ik.number_of_point_2s()),
    is_activating_(false)
#else
    fk_(ik.clone())
#endif
  {
    reset(false);
    clear_stats();

    max_=0;
    typename Indirect_kernel::Current_coordinates co= ik_.current_coordinates_object();
    for (typename Indirect_kernel::Key_iterator kit= ik_.point_2s_begin(); kit != ik_.point_2s_end(); ++kit) {
      max_= std::max(max_, std::abs(co(*kit)[0]));
      max_= std::max(max_, std::abs(co(*kit)[1]));
    }
  }

  ~Updatable_delaunay_triangulation_table_2(){CGAL_precondition(subscribers_.empty());}


  // backwards compat
  typedef Exact_point_2 Data;
  const Exact_point_2 &at(Key key) const {
    return exact_point(key);
  }

 

  //! access a point
  const Exact_point_2 &exact_point(Key key) const
  {
    CGAL_precondition(key.is_valid());
    if (exact_points_.find(key) == exact_points_.end()) {
      ++stat_exact_points_;
      Static_point_2 ip= initial(key);
      typedef typename Exact_point_2::Coordinate MF;
      typedef typename MF::NT NT;
      MF mf[2];
      CGAL_UD_DEBUG("Computing exact motions for " << key << " from " << ip );
     
#ifdef MOVE_ALL
      Static_point_2 fp= final(key);
      
      for (unsigned int i=0; i< 2; ++i){
        NT c[2];
        c[1]=-(NT(ip[i])-NT(fp[i]));
        c[0]=NT(fp[i])-c[1];
        mf[i]=MF(c, c+2);
      }
#else
      if (is_active(key)) {
        Static_point_2 fp= final(key);
        double start_time= coef_cache_[key.index()].start_time();
        
        for (unsigned int i=0; i< 2; ++i){
          NT c[2];
          c[1]=(NT(ip[i])-NT(fp[i]))/(NT(start_time)-NT(1));
          c[0]=NT(fp[i])-c[1];
          mf[i]=MF(c, c+2);
        }
	CGAL_postcondition(mf[0](NT(start_time)) == NT(ip[0]));
	CGAL_postcondition(mf[1](NT(start_time)) == NT(ip[1]));
      } else {
        mf[0] = MF(NT(ip.x()));
        mf[1] = MF(NT(ip.y()));
      }
#endif

      /*CGAL_UD_DEBUG(" starting at " << start_time 
        << " to " << fp);*/

      Exact_point_2 ep(mf[0], mf[1]);
     
      CGAL_UD_DEBUG(" got " << ep << std::endl);
      exact_points_[key]= ep;
    }

    return exact_points_.find(key)->second;
  }

  
  

  void set_is_editing(Editing_state is_e) {
    if (is_e== editing_) return;
    if (editing_==LOGGED) {
      finish_editing();
    }
    editing_=is_e;
  }

  void set_is_editing(bool tf) {
    if (tf) set_is_editing(LOGGED);
    else set_is_editing(NOT);
  }

  //! return the editing state
  bool is_editing() const
  {
    return editing_;
  }

  Key ith_key(unsigned int i) const {
    //CGAL_precondition(i < storage_.size());
    return Key(i);
  }

 
  /*void set(Key key, const Data &new_value) {
    CGAL_precondition(key.is_valid());
    CGAL_assertion(storage_.size() > static_cast<unsigned int>(key.index()));
    storage_[key.index()]=new_value;
    if (editing_== LOGGED) changed_objects_.push_back(key);
    }*/


  /*Key insert(const Data &ob) {
  //CGAL_precondition(editing_);
  storage_.push_back(ob);
  return Key(storage_.size()-1);
  }*/

  typedef CGAL::Counter<int, Key> Key_iterator;
  Key_iterator keys_begin() const {
    return Key_iterator(0);
  }
  Key_iterator keys_end() const {
    return Key_iterator(size());
  }

#ifdef MOVE_ALL

  typedef Key_iterator Changed_iterator;
  Changed_iterator changed_begin() const
  {
    return keys_begin();
  }
  Changed_iterator changed_end() const
  {
    return keys_end();
  }


  typedef Key_iterator Inserted_iterator;
  Changed_iterator inserted_begin() const
  {
    return keys_end();
  }
  Changed_iterator inserted_end() const
  {
    return keys_end();
  }
  typedef Key_iterator Erased_iterator;
  Erased_iterator erased_begin() const
  {
    return keys_end();
  }
  Erased_iterator erased_end() const
  {
    return keys_end();
  }

#else
  
  typedef typename std::vector<Key>::const_iterator Changed_iterator;
  Changed_iterator changed_begin() const
  {
    return changed_objects_.begin();
  }
  Changed_iterator changed_end() const
  {
    return changed_objects_.end();
  }


  typedef typename std::vector<Key>::const_iterator Inserted_iterator;
  Changed_iterator inserted_begin() const
  {
    return changed_objects_.end();
  }
  Changed_iterator inserted_end() const
  {
    return changed_objects_.end();
  }
  typedef typename std::vector<Key>::const_iterator Erased_iterator;
  Changed_iterator erased_begin() const
  {
    return changed_objects_.end();
  }
  Changed_iterator erased_end() const
  {
    return changed_objects_.end();
  }
#endif



  unsigned int size() const
  {
    return ik_.number_of_point_2s();
  }

  std::ostream &write(std::ostream &out) const {
    /*for (unsigned int i=0; i< storage_.size(); ++i){
      out << storage_[i].second << std::endl;
      }*/
    return out;
  }
  
  std::istream &read(std::istream &in) {
    /*if (!storage_.empty()) {
      set_is_editing(true);
      for (Key_iterator kit= keys_begin(); kit != keys_end(); ++kit){
      erase(*kit);
      }
      set_is_editing(false);
      storage_.clear();
      }
      set_is_editing(true);
      do {
      char buf[10000];
      in.getline(buf, 10000);
      if (!in || buf[0]=='\0' || buf[0]=='#') break;
      std::istringstream iss(buf);
      Data d; 
      iss >> d;
      if (!iss) {
      std::cerr << "ERROR reading object from line " << buf << std::endl;
      internal::fail__=true;
      } else {
      insert(d);
      }
      } while (true);
      set_is_editing(false);*/
    return in;
  }

  bool is_set(Key k) const {
    CGAL_precondition(notifying_);
#ifdef MOVE_ALL
    return true;
#else
    return std::binary_search(k, changed_objects_.begin(), changed_objects_.end());
#endif
  }



  bool is_initializing() const {
    return is_init_;
  }

  void set_is_initializing(bool ft) const {
    is_init_=ft;
  }

  /*std::vector<Key>& activate_list() {
    return activate_list_;
    }*/

  void clear_stats() {
    stat_interval_certificate_functions_=0;
    stat_interval_certificate_function_evaluations_=0;
    stat_interval_certificate_function_interval_evaluations_=0;
    stat_deriv_filter_calls_=0;
    stat_newton_isolate_calls_=0;
    stat_newton_refine_calls_=0;
    stat_exact_points_=0;
    stat_exact_certificate_functions_=0;
    stat_exact_certificate_functions_from_compare_=0;
    stat_exact_certificate_functions_from_advance_=0;
    stat_unfailing_exact_certificate_functions_=0;
    stat_exact_certificate_function_advances_=0;
    stat_exact_certificate_functions_from_compare_curt_ =0;
    stat_interval_predicate_evaluations_=0;
    stat_point_predicate_evaluations_=0;
    stat_certificate_computations_=0;
    stat_certificate_advances_=0;
    stat_number_of_interpolations_=0;
    stat_number_of_edges_=0;
    stat_number_of_activated_vertices_=0;
    stat_number_of_bad_edges_=0;
    stat_total_interval_width_=0;
    stat_number_of_intervals_=0;
  }


  INT interval(Key a, int i) const {
#ifndef MOVE_ALL
    if (!is_active(a)) return CGAL::to_interval(initial(a)[i]);
#endif
    INT ii= CGAL::to_interval(initial(a)[i]);
    INT fi= CGAL::to_interval(final(a)[i]);
    return INT(std::min(fi.inf(), ii.inf()),
	       std::max(fi.sup(), ii.sup()));
  }


  // need start time for each point
  INT cur_interval(Key a, INT ct, int i) const {
    CGAL_precondition(ct.inf()>=0 && ct.sup() <=1);
#ifdef MOVE_ALL
    INT c[2];
    coef(a,i,c);
    return c[0] + ct*c[1];
#else

    if (is_active(a)) {
      return coef_cache_[a.index()][i][0] + ct*coef_cache_[a.index()][i][1];
    } else {
      return CGAL::to_interval(initial(a)[i]);
    }

#endif
  }

  /*		    
    INT cur_interval(Key a, Key base, const std::set<Key> &pks, int i) const {
    double bi= initial_coordinates_object()(base)[i];
    double bf;
    if (pks.find(base) == pks.end()) {
    bf=bi;
    } else {
    bf= final_coordinates_object()(base)[i];
    }
    double ai= initial_coordinates_object()(a)[i];
    double af;
    if (pks.find(a) == pks.end()) {
    af=ai;
    } else {
    af= final_coordinates_object()(a)[i];
    }
    INT ii= INT(ai)- INT(bi);
    INT fi= INT(af)- INT(bf);
    return INT(std::min(ii.inf(), fi.inf()), std::max(ii.sup(), fi.sup()));
    }
  */

  INT cur_interval(Key a, Key base, INT ct, int i) const {
    CGAL_precondition(ct.inf()>=0 && ct.sup() <=1);
#ifdef MOVE_ALL
    Protected_array<INT,2> c, bc;
    coef(a,i,c);
    coef(base, i, bc);
    return c[0]- bc[0] + ct*(c[1]-bc[1]);
#else
    if (is_active(a) && is_active(base)) {
      return coef_cache_[a.index()][i][0]
	- coef_cache_[base.index()][i][0] 
	+ ct*(coef_cache_[a.index()][i][1]- coef_cache_[base.index()][i][1]);
    } else if (is_active(a)) {
      return coef_cache_[a.index()][i][0] - INT(initial(base)[i]) + ct*coef_cache_[a.index()][i][1];
    } else if (is_active(base)) {
      return INT(initial(a)[i])- coef_cache_[base.index()][i][0] - ct*coef_cache_[base.index()][i][1];
    } else {
      return INT(initial(a)[i]) - INT(initial(base)[i]);
    }

#endif
  }


  void write_stats(std::ostream &out) {
    out << "Triangulation had " << stat_number_of_edges_ << " edges and " << stat_number_of_vertices_ << " vertices"
	      << " and " << stat_number_of_bad_edges_ << " bad edges." << std::endl;
    out << "Points were interpolated " << stat_number_of_interpolations_ << " times resulting in " 
	      << stat_number_of_activated_vertices_ <<" activate" << std::endl;
    out << "Certificates were requested " << stat_certificate_computations_ << " times "
	      << " and advanced " << stat_certificate_advances_ << " times" << std::endl;
    out << "Point predicate evals " << stat_point_predicate_evaluations_ << std::endl;
    out << "Interval predicate evals " << stat_interval_predicate_evaluations_ << std::endl;
    out << "Interval certificate functions " << stat_interval_certificate_functions_ << std::endl;
    out << "\t deriv filtered " << stat_deriv_filter_calls_ << std::endl;
    out << "\t newton isolated " << stat_newton_isolate_calls_ << std::endl;
    out << "\t newton refined " << stat_newton_refine_calls_ << std::endl;
    out << "\t evaluations " << stat_interval_certificate_function_evaluations_ << std::endl;
    out << "\t interval evaluations " << stat_interval_certificate_function_interval_evaluations_ << std::endl;
    out << "\t interval average final width " << stat_total_interval_width_/double(stat_number_of_intervals_) << std::endl;
    out << "Exact points " << stat_exact_points_ << std::endl;
    out << "Exact certificate functions " << stat_exact_certificate_functions_ << std::endl;
    out << "\t from compare " << stat_exact_certificate_functions_from_compare_ << std::endl;
    out << "\t from compare curt " << stat_exact_certificate_functions_from_compare_curt_ << std::endl;
    out << "\t from advance " << stat_exact_certificate_functions_from_advance_ << std::endl;
    out << "\t unfailing " << stat_unfailing_exact_certificate_functions_ << std::endl;
    out << "\t advanced " << stat_exact_certificate_function_advances_ << std::endl;
  }

#ifdef MOVE_ALL
  void coef(Key a, int i, Protected_array_pointer<INT,2> r) const {
    r[0]= CGAL::to_interval(initial(a)[i]);
    r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];
  }
#endif

#ifndef MOVE_ALL



  bool is_active(Key k) const {
    return active_[k.index()];
  }

  void preactivate(double t, Key k) {
    //INT t= next_activation_;
    Protected_array<INT, 2> x, y;

    interpolate(INT(CGAL::to_interval(initial(k)[0])),
		INT(CGAL::to_interval(final(k)[0])),
		INT(t), Protected_array_pointer<INT, 2>(x));
    interpolate(INT(CGAL::to_interval(initial(k)[1])),
		INT(CGAL::to_interval(final(k)[1])),
		INT(t), Protected_array_pointer<INT, 2>(y));
  
    coef_cache_[k.index()]= Coef_data(t, x[0], x[1], y[0], y[1]);
  }

  /*bool is_activating(Key k) const {
    return coef_cache_[k.index()].start_time() == next_activation_;
    }*/

  void activate(double, Key k) {

    ++stat_number_of_activated_vertices_;
    //CGAL_precondition(t.inf() == t.sup());
    //CGAL_precondition(active_[k.index()]==false);
    if (!active_[k.index()]) {
      active_[k.index()]=true;
      changed_objects_.push_back(k);
      exact_points_.erase(k);
      //aot_->set(k, kp);
    }
  }

 
#endif

  void set_final_kernel(Indirect_kernel &fk){
    fk_=fk;
#ifndef NDEBUG
    vmax_=0;
    CGAL::Protect_FPU_rounding<true> prot;
    typename Indirect_kernel::Current_coordinates ico= ik_.current_coordinates_object();
    typename Indirect_kernel::Current_coordinates fco= fk_.current_coordinates_object();
    for (typename Indirect_kernel::Key_iterator kit= fk_.point_2s_begin(); kit != fk_.point_2s_end(); ++kit) {
      max_= std::max(max_, std::abs(fco(*kit)[0]));
      max_= std::max(max_, std::abs(fco(*kit)[1]));
      INT vx= INT(ico(*kit)[0]) - INT(fco(*kit)[0]);
      INT vy= INT(ico(*kit)[1]) - INT(fco(*kit)[1]);
      vmax_= std::max(vmax_, std::max(std::abs(vx.inf()), std::abs(vx.sup())));
      vmax_= std::max(vmax_, std::max(std::abs(vy.inf()), std::abs(vy.sup())));
    }
#endif
  }

  void reset(bool clear=true) {
    //activate_time_=-1;
#ifndef MOVE_ALL
    //next_activation_=1.0;
#endif
    start_time_=-1;
    is_init_=false;
    exact_points_.clear();
    typename Indirect_kernel::Current_coordinates cc= ik_.current_coordinates_object();
#ifndef MOVE_ALL
    for (unsigned int i=0; i< coef_cache_.size(); ++i){
      coef_cache_[i]= Coef_data(cc(Key(i))[0],cc(Key(i))[1]) ;
    }
    if (clear) active_.reset();
#endif
  }
    
  const Static_point_2& initial(Key pk) const {
    return ik_.current_coordinates_object()(pk);
  }
  const Static_point_2& final(Key pk) const {
    return fk_.current_coordinates_object()(pk);
  }
    
  /* 96 ops for poly
     25 for regular
  */

  // assumes a is 0
  // 16 ops
  template <class CNT>
  CNT  eval_incircle(CNT bx, CNT by, 
		     CNT cx, CNT cy, CNT dx, CNT dy) const {
    
    CNT qpx = bx;
    CNT qpy = by;
    CNT rpx = cx;
    CNT rpy = cy;
    CNT tpx = dx;
    CNT tpy = dy;
    CNT det=CGAL::determinant(qpx*tpy - qpy*tpx, 
				    tpx*(dx - bx) + tpy*(dy - by),
				    qpx*rpy - qpy*rpx, 
				    rpx*(cx - bx) + rpy*(cy - by));
    return det;
  }

  template <class CNT>
  CNT  eval_incircle(CNT ax, CNT ay, CNT bx, CNT by, 
		     CNT cx, CNT cy, CNT dx, CNT dy) const {
    
    CNT qpx = bx - ax;
    CNT qpy = by - ay;
    CNT rpx = cx - ax;
    CNT rpy = cy - ay;
    CNT tpx = dx - ax;
    CNT tpy = dy - ay;
    CNT det=CGAL::determinant(qpx*tpy - qpy*tpx,
				    tpx*(dx - bx) + tpy*(dy - by),
				    qpx*rpy - qpy*rpx,
				    rpx*(cx - bx) + rpy*(cy - by));
    return det;
  }

  template <class NT>
  void  incircle_p(Protected_array_const_pointer<NT,2> ax, Protected_array_const_pointer<NT,2> ay,
		   Protected_array_const_pointer<NT,2> bx, Protected_array_const_pointer<NT,2> by, 
		   Protected_array_const_pointer<NT,2> cx, Protected_array_const_pointer<NT,2> cy,
		   Protected_array_const_pointer<NT,2> dx, Protected_array_const_pointer<NT,2> dy,
		   Protected_array_pointer<NT,5> ret) const {
    Protected_array<NT,2> qx, qy, rx, ry, tx, ty;
    qx[0]= bx[0]-ax[0];
    qy[0]= by[0]-ay[0];
    rx[0]= cx[0]-ax[0];
    ry[0]= cy[0]-ay[0];
    tx[0]= dx[0]-ax[0];
    ty[0]= dy[0]-ay[0];

    qx[1]= bx[1]-ax[1];
    qy[1]= by[1]-ay[1];
    rx[1]= cx[1]-ax[1];
    ry[1]= cy[1]-ay[1];
    tx[1]= dx[1]-ax[1];
    ty[1]= dy[1]-ay[1];
    
    NT t0 = (-bx[0]+cx[0]);
    NT t1 = (-by[0]+cy[0]);
    NT t2 = (-bx[1]+cx[1]);
    NT t3 = (-by[1]+cy[1]);
    NT t4 = (-bx[0]+dx[0]);
    NT t5 = (-bx[1]+dx[1]);
    NT t6 = (-by[1]+dy[1]);
    NT t7 = (-by[0]+dy[0]);
    NT ca0 = qx[0]*ty[0]-qy[0]*tx[0];
    NT ca1 = qx[0]*ty[1]+qx[1]*ty[0]-qy[0]*tx[1]-qy[1]*tx[0];
    NT ca2 = qx[1]*ty[1]-qy[1]*tx[1]; //33
    
    NT cb0 = rx[0]*t0 + ry[0]*t1;
    NT cb1 = rx[0]*t2 + rx[1]*t0 + ry[0]*t3 + ry[1]*t1;
    NT cb2 = rx[1]*t2 + ry[1]*t3;
    
    NT cc0 = tx[0]*t4+ty[0]*t7;
    NT cc1 = tx[0]*t5+tx[1]*t4+ty[0]*t6+ty[1]*t7;
    NT cc2 = tx[1]*t5+ty[1]*t6;
    
    NT cd0 = qx[0]*ry[0]-qy[0]*rx[0];
    NT cd1 = qx[0]*ry[1]+qx[1]*ry[0]-qy[0]*rx[1]-qy[1]*rx[0];
    NT cd2 = qx[1]*ry[1]-qy[1]*rx[1];
    
    NT oa0 = ca0*cb0;
    NT oa1 = ca1*cb0 + ca0*cb1;
    NT oa2 = ca2*cb0 + ca1*cb1 + ca0*cb2;
    NT oa3 = ca2*cb1 + ca1*cb2;
    NT oa4 = ca2*cb2;
    NT ob0 = cc0*cd0;
    NT ob1 = cc1*cd0 + cc0*cd1;
    NT ob2 = cc2*cd0 + cc1*cd1 + cc0*cd2;
    NT ob3 = cc2*cd1 + cc1*cd2;
    NT ob4 = cc2*cd2;
    ret[0] = oa0 - ob0;
    ret[1] = oa1 - ob1;
    ret[2] = oa2 - ob2;
    ret[3] = oa3 - ob3;
    ret[4] = oa4 - ob4; // 51 +/- 50 *
  }

  template <class NT> 
  void interpolate(NT x0, NT x1, NT t, Protected_array_pointer<NT, 2> r) {
    NT stm1=(t-1);

    r[0]= (t*x1 - x0)/stm1;
    r[1]= (x0-x1)/stm1;
  }

#if 0
  // these are before and after rather than coefficients. 
  template <class Pt, class NT>
  void  static_incircle_p(const Pt &a0, const Pt &a1, double ta,
			  const Pt &b0, const Pt &b1, double tb,
			  const Pt &c0, const Pt &c1, double tc, 
			  const Pt &d0, const Pt &d1, double td,
			  Protected_array_pointer<NT,5> rets) const {
    
    
    Protected_array<NT,2> ax; Protected_array<NT,2> ay;
    Protected_array<NT,2> bx; Protected_array<NT,2> by;
    Protected_array<NT,2> cx; Protected_array<NT,2> cy;
    Protected_array<NT,2> dx; Protected_array<NT,2> dy;
    NT t1a, t1b, t1c, t1d;
    t1a=1-ta;
    t1b=1-tb;
    t1c=1-tc;
    t1d=1-td;
    ax[0]= (ta*a1[0] - a0[0])*t1b*t1c*t1d;
    ax[1]= (a0[0]-a1[0])*t1b*t1c*t1d;
    ay[0]= (ta*a1[1] - a0[1])*t1b*t1c*t1d;
    ay[1]= (a0[1]-a1[1])*t1b*t1c*t1d;

    bx[0]= (tb*b1[0] - b0[0])*t1a*t1c*t1d;
    bx[1]= (b0[0]-b1[0])*t1a*t1c*t1d;
    by[0]= (tb*b1[1] - b0[1])*t1a*t1c*t1d;
    by[1]= (b0[1]-b1[1])*t1a*t1c*t1d;

    cx[0]= (tc*c1[0] - c0[0])*t1b*t1a*t1d;
    cx[1]= (c0[0]-c1[0])*t1b*t1a*t1d;
    cy[0]= (tc*c1[1] - c0[1])*t1b*t1a*t1d;
    cy[1]= (c0[1]-c1[1])*t1b*t1a*t1d;

    dx[0]= (td*d1[0] - d0[0])*t1b*t1c*t1a;
    dx[1]= (d0[0]-d1[0])*t1b*t1c*t1a;
    dy[0]= (td*d1[1] - d0[1])*t1b*t1c*t1a;
    dy[1]= (d0[1]-d1[1])*t1b*t1c*t1a;

    Protected_array<NT,2> qx, qy, rx, ry, tx, ty;
    qx[0]= bx[0]-ax[0];
    qy[0]= by[0]-ay[0];
    rx[0]= cx[0]-ax[0];
    ry[0]= cy[0]-ay[0];
    tx[0]= dx[0]-ax[0];
    ty[0]= dy[0]-ay[0];

    qx[1]= bx[1]-ax[1];
    qy[1]= by[1]-ay[1];
    rx[1]= cx[1]-ax[1];
    ry[1]= cy[1]-ay[1];
    tx[1]= dx[1]-ax[1];
    ty[1]= dy[1]-ay[1];
    
    NT t0 = (-bx[0]+cx[0]);
    NT t1 = (-by[0]+cy[0]);
    NT t2 = (-bx[1]+cx[1]);
    NT t3 = (-by[1]+cy[1]);
    NT t4 = (-bx[0]+dx[0]);
    NT t5 = (-bx[1]+dx[1]);
    NT t6 = (-by[1]+dy[1]);
    NT t7 = (-by[0]+dy[0]);
    NT ca0 = qx[0]*ty[0]-qy[0]*tx[0];
    NT ca1 = qx[0]*ty[1]+qx[1]*ty[0]-qy[0]*tx[1]-qy[1]*tx[0];
    NT ca2 = qx[1]*ty[1]-qy[1]*tx[1]; //33
    
    NT cb0 = rx[0]*t0 +  ry[0]*t1;
    NT cb1 = rx[0]*t2 + rx[1]*t0 + ry[0]*t3 + ry[1]*t1;
    NT cb2 = rx[1]*t2 + ry[1]*t3;
    
    NT cc0 = tx[0]*t4+ty[0]*t7;
    NT cc1 = tx[0]*t5+tx[1]*t4+ty[0]*t6+ty[1]*t7;
    NT cc2 = tx[1]*t5+ty[1]*t6;
    
    NT cd0 = qx[0]*ry[0]-qy[0]*rx[0];
    NT cd1 = qx[0]*ry[1]+qx[1]*ry[0]-qy[0]*rx[1]-qy[1]*rx[0];
    NT cd2 = qx[1]*ry[1]-qy[1]*rx[1];
    
    NT oa0 = ca0*cb0;
    NT oa1 = ca1*cb0 + ca0*cb1;
    NT oa2 = ca2*cb0 + ca1*cb1 + ca0*cb2;
    NT oa3 = ca2*cb1 + ca1*cb2;
    NT oa4 = ca2*cb2;
    NT ob0 = cc0*cd0;
    NT ob1 = cc1*cd0 + cc0*cd1;
    NT ob2 = cc2*cd0 + cc1*cd1 + cc0*cd2;
    NT ob3 = cc2*cd1 + cc1*cd2;
    NT ob4 = cc2*cd2;
    ret[0] = oa0 - ob0;
    ret[1] = oa1 - ob1;
    ret[2] = oa2 - ob2;
    ret[3] = oa3 - ob3;
    ret[4] = oa4 - ob4; // 51 +/- 50 *
  }
#endif
  // t is 1-t upperbound
  template <class NT> 
  void incircle_equiv_p(NT a0, NT a1, NT t, Protected_array_pointer<NT, 5> ret) const {
    //NT tm1= 1-t;
    NT tm13= t*t*t;
    NT q0= /*tm13*/(a0+a0);
    NT q1= /*tm13*/(a1+a1);
    
    NT q02= q0*q0; // 6
    NT q12= q1*q1; // 6
    NT sq02= q02+q02; // 6,1
    NT sq12= q12+q12; // 6,1
    NT q0q1= q0*q1; // 6
    NT sq0q1= q0q1+q0q1; // 6,1
    NT ssq0q1= sq0q1+sq0q1; // 6,2
    
    NT ssq0q1sq12= ssq0q1*sq12;
    NT ssq0q1sq02= ssq0q1*sq02;
    NT sq02sq12= sq12*sq02;
    
    NT oa0 = sq02*sq02;
    NT oa1 = ssq0q1sq02 + ssq0q1sq02;
    NT oa2 = sq02sq12 + ssq0q1*ssq0q1 + sq02sq12;
    NT oa3 = ssq0q1sq12 + ssq0q1sq12;
    NT oa4 = sq12*sq12;
  
    ret[0] = oa0 + oa0;
    ret[1] = oa1 + oa1;
    ret[2] = oa2 + oa2;
    ret[3] = oa3 + oa3;
    ret[4] = oa4 + oa4; // 51 +/- 50 *
  }

  template <class NT>
  void static_function_bounds(Protected_array_const_pointer<NT, 5> ret) const {
    NT t(1);
    NT vb= ret[0] + ret[1]*t+ ret[2]*t*t + ret[3]*t*t*t+ ret[4]*t*t*t*t;
    NT db= ret[1]+ NT(2.0)*ret[2]*t + NT(3.0)*ret[3]*t*t+ NT(4.0)*ret[4]*t*t*t;
    std::cout << "Function error is " << vb.error() << " and deriv error is " << db.error() << std::endl;
  }

  void incircle_static(Key a, Key b, Key c, Key d) const {
    typedef CGAL::Static_filter_error NT;
    typedef Protected_array_pointer<NT, 5> RP;
    typedef Protected_array_const_pointer<NT,2> P;
    typedef Protected_array_pointer<double, 5> DRP;
    typedef Protected_array_const_pointer<double,2> DP;

    Protected_array<double, 2> a0; a0[0]=initial(a)[0]; a0[1]=final(a)[0]-initial(a)[0];
    Protected_array<double, 2> a1; a1[0]=initial(a)[1]; a1[1]=final(a)[1]-initial(a)[1];
    Protected_array<double, 2> b0; b0[0]=initial(b)[0]; b0[1]=final(b)[0]-initial(b)[0];
    Protected_array<double, 2> b1; b1[0]=initial(b)[1]; b1[1]=final(b)[1]-initial(b)[1];
    Protected_array<double, 2> c0; c0[0]=initial(c)[0]; c0[1]=final(c)[0]-initial(c)[0];
    Protected_array<double, 2> c1; c1[0]=initial(c)[1]; c1[1]=final(c)[1]-initial(c)[1];
    Protected_array<double, 2> d0; d0[0]=initial(d)[0]; d0[1]=final(d)[0]-initial(d)[0];
    Protected_array<double, 2> d1; d1[0]=initial(d)[1]; d1[1]=final(d)[1]-initial(d)[1];    
    /*r[0]= CGAL::to_interval(initial(a)[i]);
      r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];*/
    Protected_array<double, 5> dret;
    incircle_p(DP(a0), DP(a1), DP(b0), DP(b1), DP(c0), DP(c1), DP(d0), DP(d1),
	       DRP(dret));
    std::cout << "Approx is " << dret[0]
	      << " + " << dret[1] << "*t + " 
	      << dret[2] << "*t^2 + " 
	      << dret[3] << "*t^3 + " 
	      << dret[4] << "*t^4" << std::endl; 
    {
      
      Protected_array<NT, 2> a0; a0[0]=max_; a0[1]=vmax_;
      Protected_array<NT, 2> a1; a1[0]=max_; a1[1]=vmax_;
      Protected_array<NT, 2> b0; b0[0]=max_; b0[1]=vmax_;
      Protected_array<NT, 2> b1; b1[0]=max_; b1[1]=vmax_;
      Protected_array<NT, 2> c0; c0[0]=max_; c0[1]=vmax_;
      Protected_array<NT, 2> c1; c1[0]=max_; c1[1]=vmax_;
      Protected_array<NT, 2> d0; d0[0]=max_; d0[1]=vmax_;
      Protected_array<NT, 2> d1; d1[0]=max_; d1[1]=vmax_;
      /*r[0]= CGAL::to_interval(initial(a)[i]);
	r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];*/
      Protected_array<NT, 5> ret;

   
      incircle_p(P(a0), P(a1), P(b0), P(b1), P(c0), P(c1), P(d0), P(d1),
		 RP(ret));
      Protected_array<NT, 5> eret;
      incircle_equiv_p(a0[0], a0[1], NT(.5), RP(eret));
      static_function_bounds(Protected_array_const_pointer<NT,5>(eret));
      std::cout << "Static bounds are " << ret[0].bound() << "+/- " << ret[0].error() 
		<< ", " << ret[1].bound() << "+/- " << ret[1].error() 
		<< ", " << ret[2].bound() << "+/- " << ret[2].error() 
		<< ", " << ret[3].bound() << "+/- " << ret[3].error() 
		<< ", " << ret[4].bound() << "+/- " << ret[4].error() << std::endl;

      std::cout << "Other errors are "<< eret[0].error() 
		<< ", " << eret[1].error() 
		<< ", " << eret[2].error() 
		<< ", " << eret[3].error() 
		<< ", " << eret[4].error() << std::endl;
      std::cout << "Max is " << max_ << " and vmax is " << vmax_ << std::endl;
      std::cout << "Static functions lb " << dret[0]- ret[0].error() 
		<< " + " << dret[1] - ret[1].error() << "*t + " 
		<< dret[2] - ret[2].error() << "*t^2 + " 
		<< dret[3] - ret[3].error() << "*t^3 + " 
		<< dret[4] - ret[4].error() << "*t^4" << std::endl; 
      std::cout << "Static functions ub " << dret[0]+ ret[0].error() 
		<< " + " << dret[1] + ret[1].error() << "*t + " 
		<< dret[2] + ret[2].error() << "*t^2 + " 
		<< dret[3] + ret[3].error() << "*t^3 + " 
		<< dret[4] + ret[4].error() << "*t^4" << std::endl;
    }
    if (0) {
      Protected_array<NT, 2> a0; a0[0]=initial(a)[0]; a0[1]=NT(final(a)[0])-NT(initial(a)[0]);
      Protected_array<NT, 2> a1; a1[0]=initial(a)[1]; a1[1]=NT(final(a)[1])-NT(initial(a)[1]);
      Protected_array<NT, 2> b0; b0[0]=initial(b)[0]; b0[1]=NT(final(b)[0])-NT(initial(b)[0]);
      Protected_array<NT, 2> b1; b1[0]=initial(b)[1]; b1[1]=NT(final(b)[1])-NT(initial(b)[1]);
      Protected_array<NT, 2> c0; c0[0]=initial(c)[0]; c0[1]=NT(final(c)[0])-NT(initial(c)[0]);
      Protected_array<NT, 2> c1; c1[0]=initial(c)[1]; c1[1]=NT(final(c)[1])-NT(initial(c)[1]);
      Protected_array<NT, 2> d0; d0[0]=initial(d)[0]; d0[1]=NT(final(d)[0])-NT(initial(d)[0]);
      Protected_array<NT, 2> d1; d1[0]=initial(d)[1]; d1[1]=NT(final(d)[1])-NT(initial(d)[1]);
    
      /*r[0]= CGAL::to_interval(initial(a)[i]);
	r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];*/
      Protected_array<NT, 5> ret;
      incircle_p(P(a0), P(a1), P(b0), P(b1), P(c0), P(c1), P(d0), P(d1),
		 RP(ret));
      std::cout << "Pseudo static bounds are " << ret[0].bound() << "+/- " << ret[0].error() 
		<< ", " << ret[1].bound() << "+/- " << ret[1].error() 
		<< ", " << ret[2].bound() << "+/- " << ret[2].error() 
		<< ", " << ret[3].bound() << "+/- " << ret[3].error() 
		<< ", " << ret[4].bound() << "+/- " << ret[4].error() << std::endl;
      std::cout << "Pseudo static functions lb " << dret[0]- ret[0].error() 
		<< " + " << dret[1] - ret[1].error() << "*t + " 
		<< dret[2] - ret[2].error() << "*t^2 + " 
		<< dret[3] - ret[3].error() << "*t^3 + " 
		<< dret[4] - ret[4].error() << "*t^4" << std::endl; 
      std::cout << "Pseudo static functions ub " << dret[0]+ ret[0].error() 
		<< " + " << dret[1] + ret[1].error() << "*t + " 
		<< dret[2] + ret[2].error() << "*t^2 + " 
		<< dret[3] + ret[3].error() << "*t^3 + " 
		<< dret[4] + ret[4].error() << "*t^4" << std::endl; 
		
    }
   
  }
		       

  void certificate_function(Key a, Key b, Key c, Key d, Protected_array_pointer<INT, 5> ret) const {
    ++stat_interval_certificate_functions_;
#ifndef NDEBUG
    incircle_static(a,b,c,d);
#endif

#ifdef MOVE_ALL
    Protected_array<INT,2> a0; coef(a, 0, a0);
    Protected_array<INT,2>  a1; coef(a, 1, a1);
    Protected_array<INT,2>  b0; coef(b, 0, b0);
    Protected_array<INT,2>  b1; coef(b, 1, b1);
    Protected_array<INT,2>  c0; coef(c, 0, c0);
    Protected_array<INT,2>  c1; coef(c, 1, c1);
    Protected_array<INT,2>  d0; coef(d, 0, d0);
    Protected_array<INT,2>  d1; coef(d, 1, d1);
    typedef Protected_array_const_pointer<INT,2> P;
    incircle_p(P(a0), P(a1), P(b0), P(b1), P(c0), P(c1), P(d0), P(d1),
	       ret);
    
#else
    incircle_p(coef_cache_[a.index()][0],
	       coef_cache_[a.index()][1],
	       coef_cache_[b.index()][0],
	       coef_cache_[b.index()][1],
	       coef_cache_[c.index()][0],
	       coef_cache_[c.index()][1],
	       coef_cache_[d.index()][0],
	       coef_cache_[d.index()][1],
	       ret);
#endif

#ifndef NDEBUG
    Exact_certificate ec= soc_(exact_point(a),
			       exact_point(b),
			       exact_point(c), 
			       exact_point(d),
			       0, 1);
#endif
    CGAL_UD_DEBUG("Certificate func is " << ret[0] << " + "
		  << ret[1] << "*t + " << ret[2] << "*t^2 + "
		  << ret[3] << "*t^3 + " << ret[4] << "*t^4" 
		  << std::endl );
  }

   
  /*void point_changed(Key ){
    }*/

  Exact_certificate compute_exact_failure_time(Cert_tuple ct, 
					       Exact_time et) const {
    CGAL_UD_DEBUG("Computing exact time for " << ct 
		  << " from " << et  << std::endl);
    Exact_certificate ec= soc_(exact_point(ct[0]),
			       exact_point(ct[1]),
			       exact_point(ct[2]), 
			       exact_point(ct[3]),
			       et, 1);
    CGAL_UD_DEBUG("Got " << ec.failure_time() << std::endl);
    return ec;
  }

  void update_exact_failure_time(Cert_tuple ct, Exact_time et,
				 Cert_tuple ot, 
				 Exact_certificate& oc) const {
    //CGAL_precondition(check_.find(ct) != check_.end());
    CGAL_UD_DEBUG("Updating exact time for  " << ct 
		  << " starting at " << et << std::endl);

    if (ot != ct) {
      CGAL_precondition(ot.opposite() == ct);
      CGAL_UD_DEBUG("Flipping exact time for " << ct 
		    << " from " << oc.failure_time()
		    << std::endl);
      oc.pop_failure_time();
    }
      
    while (oc.failure_time() < et) {
      oc.pop_failure_time();
      oc.pop_failure_time();
      CGAL_UD_DEBUG("Advancing exact time for " << ct 
		    << " from " << oc.failure_time() 
		    << std::endl);
    }
    CGAL_UD_DEBUG("Got " << oc.failure_time() << std::endl);
      
  }
  
   

  std::pair<double,double> join(std::pair<double,double> a, INT b) const {
    return std::make_pair(std::min(a.first, b.inf()),
			  std::max(a.second, b.sup()));
  }


  CGAL::Sign sign_at(Key a, Key b, 
		     Key c, Key d,
		     INT ct) const {
    INT det= eval_incircle(//cur_interval(a,ct,0),
			   //cur_interval(a,ct,1),
			   cur_interval(b,a,ct,0),
			   cur_interval(b,a,ct,1),
			   cur_interval(c,a,ct,0),
			   cur_interval(c,a,ct,1),
			   cur_interval(d,a,ct,0),
			   cur_interval(d,a,ct,1));
    
    if (det.sup() < 0) return CGAL::NEGATIVE;
    else if (det.inf() > 0) return CGAL::POSITIVE;
    else return CGAL::ZERO;
  }

  CGAL::Sign sign_at(Key a, Key b, 
		     Key c, Key d,
		     const std::set<Key> &pks) const {
    INT det= eval_incircle(//cur_interval(a,ct,0),
			   //cur_interval(a,ct,1),
			   cur_interval(b,a,pks,0),
			   cur_interval(b,a,pks,1),
			   cur_interval(c,a,pks,0),
			   cur_interval(c,a,pks,1),
			   cur_interval(d,a,pks,0),
			   cur_interval(d,a,pks,1));
    
    if (det.sup() < 0) return CGAL::NEGATIVE;
    else if (det.inf() > 0) return CGAL::POSITIVE;
    else return CGAL::ZERO;
  }
#ifndef MOVE_ALL
  bool ok_at_end(Key a, Key b, Key c, Key d) const {
    typedef typename Indirect_kernel::Side_of_oriented_circle_2 OC;
    OC oc= ik_.side_of_oriented_circle_2_object();
    return oc(a,b,c,d) == CGAL::ON_POSITIVE_SIDE;
  }

  bool mixed_ok_at_end(Key a, Key b, Key c, Key d) const {
    typedef typename Indirect_kernel::Direct_kernel::Side_of_oriented_circle_2 OC;
    OC oc= ik_.direct_kernel_object().side_of_oriented_circle_2_object();
    typename Indirect_kernel::Direct_kernel::Point_2 pa, pb, pc, pd;
    typename Indirect_kernel::Current_coordinates fko= final_coordinates_object(), iko= initial_coordinates_object();
    if (is_active(a)) {
      pa= fko(a);
    } else {
      pa = iko(a);
    }
    if (is_active(b)) {
      pb= fko(b);
    } else {
      pb = iko(b);
    }
    if (is_active(c)) {
      pc= fko(c);
    } else {
      pc = iko(c);
    }
    if (is_active(d)) {
      pd= fko(d);
    } else {
      pd = iko(d);
    }
    CGAL::Oriented_side os = oc(pa, pb, pc, pd); 
    return (os == CGAL::ON_POSITIVE_SIDE);
  }
#endif

  void set(Key k, Exact_point_2 ep) {
    exact_points_[k]=ep;
#ifndef MOVE_ALL
    changed_objects_.push_back(k);
#endif
  }

  

  enum Isolate_result {NO_FAILURE=-1, POSSIBLE_FAILURE=0, CERTAIN_FAILURE=1};

  /*struct Static_evaluator {
    Static_evaluator(Cert_tuple ct, const This *t): tuple_(ct), ui_(t) {}
    
    CGAL::Sign operator()(INT t) const {
    return ui_->sign_at(tuple_[0], tuple_[1], tuple_[2], tuple_[3], t);
    }
    Cert_tuple tuple_;
    const This *ui_;
    };*/

  template <unsigned int D>
  static INT evaluate_ipoly(Protected_array_const_pointer<INT, D> coefs, const INT t) {
    CGAL_precondition(t.inf() >=0);
    CGAL_precondition(t.sup() <=1);
    CGAL_precondition(t.inf() <=1);
    {
      INT tt(t.inf());
      CGAL_assertion(tt.inf()==tt.sup());
    }
    ++stat_interval_certificate_function_evaluations_;
    if (t.inf() != t.sup()) {
      ++stat_interval_certificate_function_interval_evaluations_;
    }
    if (D==-1) return INT(0);
    INT cum=coefs[0];
    CGAL_assertion_code(INT debug_cum= coefs[0]);
    if (D==0) {
      CGAL_postcondition(cum.inf() == debug_cum.inf() && cum.sup() == debug_cum.sup());
      return cum;
    }
    if (coefs[1].inf() > 0) {
      cum += INT(-(coefs[1].inf()*(-t.inf())), coefs[1].sup()*t.sup());
    } else if (coefs[1].sup() < 0) {
      cum += INT(-(coefs[1].inf()*(-t.sup())), coefs[1].sup()*t.inf());
    } else {
      cum += INT(-(coefs[1].inf()*(-t.sup())), coefs[1].sup()*t.sup());
    }
    CGAL_assertion_code(debug_cum += coefs[1]*t);
    CGAL_postcondition(cum.inf() == debug_cum.inf() && cum.sup() == debug_cum.sup());
    if (D==1) {
      CGAL_postcondition(cum.inf() == debug_cum.inf() && cum.sup() == debug_cum.sup());
      return cum;
    }
    INT tm=t;
    for (unsigned int i=2; i<=D-1; ++i) {
      tm = INT(-(tm.inf()*(-t.inf())), tm.sup()*t.sup());
      if (coefs[i].inf() > 0) {
	cum += INT(-(coefs[i].inf()*(-tm.inf())), coefs[i].sup()*tm.sup());
      } else if (coefs[i].sup() < 0) {
	cum += INT(-(coefs[i].inf()*(-tm.sup())), coefs[i].sup()*tm.inf());
      } else {
	cum += INT(-(coefs[i].inf()*(-tm.sup())), coefs[i].sup()*tm.sup());
      }
      CGAL_assertion_code(debug_cum += coefs[i]*tm);
      CGAL_postcondition(cum.inf() == debug_cum.inf() && cum.sup() == debug_cum.sup());
 
      //cum += INT(-(coefs[i].inf()*(-tm.inf())), coefs[i].sup()*tm.sup());
    }
    return cum;
  }


 template <unsigned int D>
  static double evaluate_ipoly_inf(Protected_array_const_pointer<INT, D> coefs, const INT t) {
    CGAL_precondition(t.inf() >=0);
    CGAL_precondition(t.sup() <=1);
    CGAL_precondition(t.inf() <=1);
    ++stat_interval_certificate_function_evaluations_;
    if (t.inf() != t.sup()) {
      ++stat_interval_certificate_function_interval_evaluations_;
    }
    if (D==-1) return 0.0;
    double cum=-coefs[0].inf();
    if (D==0) {
      return -cum;
    }
    if (coefs[1].inf() > 0) {
      cum += coefs[1].inf()*(-t.inf());
    } else {
      cum += coefs[1].inf()*(-t.sup());
    } 
    if (D==1) {
      return -cum;
    }
    INT tm=t;
    for (unsigned int i=2; i<=D-1; ++i) {
      tm = INT(-(tm.inf()*(-t.inf())), tm.sup()*t.sup());
      if (coefs[i].inf() > 0) {
	cum += coefs[i].inf()*(-tm.inf());
      } else {
	cum += coefs[i].inf()*(-tm.sup());
      } 
    }
    return -cum;
  }

  bool can_fail(Cert_tuple ct, double cur_time
#ifndef MOVE_ALL
		, double end_time
#endif
		) const {
    CGAL::Protect_FPU_rounding<true> prot;
#ifdef MOVE_ALL
    INT rct(cur_time, 1);
#else
    INT rct(cur_time, end_time);
#endif
    bool ret= (sign_at(ct[0],ct[1],ct[2],ct[3], rct) != CGAL::POSITIVE);
    ++stat_interval_predicate_evaluations_;
   
    return ret;
  }

  template <int N>
  static CGAL::Sign sign_at(Protected_array_const_pointer<INT, N> coefs, const INT t) {
    INT cum= evaluate_ipoly(coefs, t);
    if (cum.inf() > 0) return CGAL::POSITIVE;
    else if (cum.sup() < 0) return CGAL::NEGATIVE;
    else return CGAL::ZERO;
  }

  template <int N>
  struct Certificate_sign_at {
    Certificate_sign_at(Protected_array_const_pointer<INT,N> coefs): coefs_(coefs) {
    }
    
    CGAL::Sign operator()(INT t) const {
      return sign_at(coefs_, t);
    }
    Protected_array_const_pointer<INT, N> coefs_;
  };

 

  /*
    I think I can only handle the case where I am positive at the beginning anyway. So use this assumption. 

    I can prove if a root is not there by having the derivitive not be able to hit the interval.

    I can prove that a root is there by having two evaluations of opposite signs.

    - One solution: isolate all roots of the function and the derivitive to some scale. 


    Newton_refine takes an interval, a poly and a deriv and returns a smaller interval (the interval must contain a root)
    Even this is a bit hard. 
    Evaluate sign in middle.
    if it is negative then take left half.
    If is is positive, and sign of deriv on interval is negative, take right (or if I can't hit left).
    What if it is positive, and deriv has positive and negative parts then call Newton_search on left and use its return.

    Newton_search takes an interval which is the same at both end points and f and fp and fpp and sign at ends.
    - it tries to determine if there is a root in the interval. 
    - To do this, search for minimum/max points in the interval-eval deriv at endpoints, call Newton_search if both positive/negative Newton_refine otherwise.
    - If there is an isolated min, split on it. 
    

    Newton_search takes the function and two derivitives.
    the goal is to find the closest positive to the left and negative to the right values.
    eval function at left of interval (possiibly update positive) and right (and update negative)
    find leftmost root of derivitive using newton (when do I stop?--ideally when the signs of f on both ends are the same)
    evaluate there and see if I have a positive and a negative, is so done, otherwise repeat
  
  */

  template <class F, class FP>
  double derivitive_filter(F f, FP fp, double lb, double ub) {
    ++stat_deriv_filter_calls_;
    double clb=lb;
    double dvp= fp.inf(INT(lb,ub));
#ifndef NDEBUG
    /*{
      INT ovp = fp(INT(lb, ub));
      CGAL_UD_DEBUG("New lb for deriv is " << dvp  << " old lb is " << ovp << std::endl);
      }*/
#endif
    double step=0;
    int nsteps=0;
    // step is negative to get rounding right
    if (dvp < 0) {
      INT vp = INT(-dvp);
      do {
	double v= f.inf(clb);
	if (v <= 0) {
	  CGAL_UD_DEBUG("Deriv return " << clb << " from " << lb << " with " << nsteps << std::endl);
	  return clb;
	}
	step = (-v)/vp.sup();
	clb -= step;
	++nsteps;
      } while (clb < ub && step < .01);
      
    } else {
      clb=ub;
    }
    CGAL_UD_DEBUG("Deriv return " << clb << " from " << lb << " with " << nsteps << std::endl);
    return clb;
  }



  template <class F, class FP, bool OPT>
  bool internal_Newton_refine(F f, FP fp, double max_size, std::vector<INT> &stack) const {
    CGAL_precondition(!stack.empty());
    CGAL_UD_DEBUG("Newton refine " <<  stack.back() << std::endl);
    // can skip second if mp is negative
    bool front_is_sure= (stack.size()==1);
    do {
      CGAL_precondition(!stack.empty());
      if (stack.back().sup() - stack.back().inf() < max_size) {
	CGAL_UD_DEBUG("Returning small interval " << stack.back() << std::endl);
	return front_is_sure;
      }
      INT c= stack.back();
      stack.pop_back();
      CGAL_UD_DEBUG("Cur is " << c << std::endl);
      INT m(.5*(c.inf()+ c.sup()));
      CGAL_UD_DEBUG("m is " << m << std::endl);
      CGAL_assertion(m.inf() >= c.inf() && m.sup() <= c.sup());
      INT fm= f(m);
      CGAL_UD_DEBUG("fm is " << fm << std::endl);
      {
	if (CGAL::sign(fm.inf()) != CGAL::sign(fm.sup())) {
	  // put in random fallback
	  m= (.25*c.inf()+ .75*c.sup());
	  CGAL_UD_DEBUG("m is " << m << std::endl);
	  fm= f(m);
	  CGAL_UD_DEBUG("fm is " << fm << std::endl);
	  if (CGAL::sign(fm.inf()) != CGAL::sign(fm.sup())) {
	    // put in random fallback
	    m= (.75*c.inf()+ .25*c.sup());
	    CGAL_UD_DEBUG("m is " << m << std::endl);
	    fm= f(m);
	    CGAL_UD_DEBUG("fm is " << fm << std::endl);
	    if (CGAL::sign(fm.inf()) != CGAL::sign(fm.sup())) {
	      stack.push_back(c); 
	      CGAL_UD_DEBUG("Returning giving up " << stack.back() << std::endl);
	      return front_is_sure;
	    }
	  }
	}

      }

      CGAL_assertion(sign(fm.inf()) == sign(fm.sup()));
      INT fpc= fp(c);
      CGAL_UD_DEBUG("fpc is " << fpc << std::endl);
    
      if (fpc.inf() > 0 || fpc.sup() < 0){
	if (OPT && fm.sup() < 0 && fpc.inf() > 0) front_is_sure = true;
	INT ni= m - fm/fpc;
	double mn= std::max(c.inf(), ni.inf());
	double mx= std::min(c.sup(), ni.sup());
	CGAL_UD_DEBUG("ni is " << ni << std::endl);
	if (mn <= mx) stack.push_back(INT(mn, mx));
      } else {
	INT fpcp= INT(0, fpc.sup());
	INT fpcn= INT(fpc.inf(), -0.0);
	if (fm.inf() > 0) {
	  if (OPT) front_is_sure=false;
	  // could skip the second one
	  //process_newton(fm, fpcp, m, c, stack);
	  {
	    INT ni= m - fm/INT(fpc.inf());
	    CGAL_UD_DEBUG("ni is " << ni << std::endl);
	    if (ni.inf() <= c.sup()) {
	      stack.push_back(INT(ni.inf(), c.sup()));
	    }
	  }
	  {
	    INT ni= m - fm/INT(fpc.sup());
	    CGAL_UD_DEBUG("ni is " << ni << std::endl);
	    if (ni.sup() >= c.inf()) {
	      stack.push_back(INT(c.inf(), ni.sup()));
	    }
	  }
	} else {
	  if (OPT) front_is_sure=true;
	  if (!OPT) {
	    INT ni= m - fm/INT(fpc.sup());
	    CGAL_UD_DEBUG("ni is " << ni << std::endl);
	    if (ni.inf() <= c.sup()) {
	      stack.push_back(INT(ni.inf(), c.sup()));
	    }
	  }
	  {
	    INT ni= m - fm/INT(fpc.inf());
	    CGAL_UD_DEBUG("ni is " << ni << std::endl);
	    if (ni.sup() >= c.inf()) {
	      stack.push_back(INT(c.inf(), ni.sup()));
	    }
	  }
	}
      }
    } while (!stack.empty());
    CGAL_assertion(!OPT);
    return front_is_sure;
  }

  
  template <class F, class FP> 
  INT Newton_refine(F f, FP fp, INT ii) const {
    ++stat_newton_refine_calls_;
    std::vector<INT> stack;
    stack.push_back(ii);
    do {
      CGAL_assertion(!stack.empty());
      bool fis= internal_Newton_refine<F, FP, true>(f, fp, .0625*(ii.sup() - ii.inf()), stack);
      if (!fis && stack.size() >1) {
	// we don't know if the front interval actually contains a root. 
	//CGAL_assertion(stack.size() >1);
	for (unsigned int i=2; i<= stack.size(); ++i) {
	  INT mp= .5*(stack[stack.size()-i+1].inf()+ INT(stack[stack.size()-i].sup()));
	  // check that it is in interval
	  if (f(mp).sup() < 0) {
	    CGAL_UD_DEBUG( "Found split at " <<  mp << std::endl);
	    return INT(stack.back().inf(), stack[stack.size()-i].sup());
	  }
	}
	CGAL_UD_DEBUG( "Failed to find split for " << ii << " returning "
		       << INT(stack.back().inf(), stack.front().sup()) << std::endl);
	return INT(stack.back().inf(), stack.front().sup());
      } else {
	CGAL_precondition(!stack.empty());
	CGAL_UD_DEBUG( "Returning " << stack.back() << std::endl);
	return stack.back();
      }
    } while (true);
  }

  /* void process_newton(INT fm, INT fmp, INT m, INT c, std::vector<INT> &stack) const {
     INT ni= m - fm/fmp;
     double mn= std::max(c.inf(), ni.inf());
     double mx= std::min(c.sup(), ni.sup());
     if (mn <= mx) stack.push_back(INT(mn, mx));
     }*/



  template <class F, class FP, class FPP> 
  Isolate_result Newton_isolate(F f, FP fp, FPP fpp, bool lb_is_positive,
				double lb, double ub, INT &ret) const {
    ++stat_newton_isolate_calls_;
    CGAL_UD_DEBUG("Newton isolate " <<  lb <<  " ... " << ub << std::endl);
    if (lb==ub) {
      CGAL_UD_DEBUG("Empty interval " << std::endl);
      return NO_FAILURE;
    }
    if (!lb_is_positive) {
      INT flb = f(lb);
      if (flb.inf() <= 0) {
	CGAL_UD_DEBUG( "Returning possible failure at start " << std::endl);
	return POSSIBLE_FAILURE;
      } else {
	lb= nextafter(lb, std::numeric_limits<double>::max());
      }
    }
    bool possible=false;
    std::vector<INT> stack;
    stack.push_back(INT(lb,ub));
    do {
      double wid=.0625;
      do {
	CGAL_precondition(!stack.empty());
	internal_Newton_refine<FP,FPP,false>(fp, fpp, wid*(ub-lb), stack);
	if (stack.empty()) break;
	CGAL_UD_DEBUG("Trying " << stack.back() << std::endl);
	INT fe= f(stack.back());
	CGAL_UD_DEBUG("Got fe= " << fe << std::endl);
	if (fe.sup() < 0) {
	  ret= INT(lb, stack.back().sup());
	  return CERTAIN_FAILURE;
	} else if (fe.inf() > 0) {
	  if (!possible) {
	    lb= stack.back().sup();
	  }
	  stack.pop_back();
	  break;
	} else if (wid < .0000001) {
	  CGAL_UD_DEBUG("Failed to resolve sign on interval " << stack.back() << std::endl);
	  possible=true;
	  stack.pop_back();
	  break;
	} else {
	  wid=.25*wid;
	}
      } while (true);
    } while (!stack.empty());

    CGAL_UD_DEBUG("Trying ub " << ub << std::endl);
    INT fub= f(ub);
    CGAL_UD_DEBUG( "Got fub= " << fub << std::endl);
    if (fub.inf() < 0) {
      ret= INT(lb,ub);
      return CERTAIN_FAILURE;
    } else if (fub.sup() > 0 && !possible) {
      return NO_FAILURE;
    } else {
      CGAL_UD_DEBUG( "Returning possible failure at end " << std::endl);
      return POSSIBLE_FAILURE;
    }
  }
  
  double start_time() const {
    return start_time_;
  }
  void set_start_time(double f) {
    start_time_=f;
  }

  const Indirect_kernel& final_kernel_object() const {
    return fk_;
  }

  const Indirect_kernel& initial_kernel_object() const {
    return ik_;
  }

  typename Indirect_kernel::Current_coordinates initial_coordinates_object() const {
    return ik_.current_coordinates_object();
  }
  typename Indirect_kernel::Current_coordinates final_coordinates_object() const {
    return fk_.current_coordinates_object();
  }

  const typename Kinetic_kernel::Positive_side_of_oriented_circle_2& in_circle_object() const {
    return soc_;
  }

  //template <class SH>
#ifndef MOVE_ALL
  void schedule_activation(Key k, typename Simulator::Handle sh) {
    if (activation_event_== Event_key()) {
      if (is_activating()) {
	CGAL_precondition(CGAL::to_interval(sh->current_time()).first == (CGAL::to_interval(sh->current_time()).second));
	IE ev(this, CGAL::to_interval(sh->current_time()).first);
	activation_event_= sh->new_event(ev.time(), ev);
      } else {
	IE ev(this, sh);
	activation_event_= sh->new_event(ev.time(), ev);
      }
      ++stat_number_of_interpolations_;
    }
    sh->template event<IE>(activation_event_).insert(k);
  }

  //template <class SH>
  double activation_time(typename Simulator::Handle sh) const {
    if (activation_event_ != Event_key()) {
      return sh->template event<IE>(activation_event_).time();
    } else {
      return 1;
    }
  }
  bool is_activating() const {
    return is_activating_;
  }
  void set_is_activating(bool tf) {
    if (tf) {
      activation_event_= Event_key();
    }
    is_activating_= tf;
  }
#endif

private:
  friend class Kinetic::Multi_listener<Listener_core>;
  //! listen for changes
  /*!
    This method alerts the subscribe to all exising objects.
  */
  void new_listener(Listener *sub) const
  {
    subscribers_.insert(sub);
  }

  //! end listening for changes
  void delete_listener(Listener *sub) const
  {
    subscribers_.erase(sub);
  }

  void finish_editing() {
#ifdef MOVE_ALL
    notifying_=true;
    for (typename Subscribers::iterator it= subscribers_.begin(); it != subscribers_.end(); ++it) {
      (*it)->new_notification(Listener::IS_EDITING);
    }
    notifying_=false;
#else
    if (!changed_objects_.empty()) {
      notifying_=true;
      std::sort(changed_objects_.begin(), changed_objects_.end());
      for (typename Subscribers::iterator it= subscribers_.begin();
	   it != subscribers_.end(); ++it) {
	(*it)->new_notification(Listener::IS_EDITING);
      }
      changed_objects_.clear();
      notifying_=false;
    }
#endif
  }




protected:
  mutable std::map<Key, Exact_point_2> exact_points_;
  //std::vector<Coef_data > coef_cache_;
  bool is_init_;
  //std::vector<Key> activate_list_;
  double start_time_;
  
  mutable Subscribers subscribers_;
  Editing_state editing_;
  bool notifying_;


  typename Kinetic_kernel::Positive_side_of_oriented_circle_2 soc_;
  Indirect_kernel ik_, fk_;

  double max_, vmax_;

#ifndef MOVE_ALL
  Active active_;
  std::vector<Coef_data> coef_cache_;
  std::vector<Key> changed_objects_;
  Event_key activation_event_;
  bool is_activating_;
#endif


  

};

/*template <class V, class K >
  inline std::ostream &operator<<(std::ostream &out, const Active_objects_update_vector<V, K> &v) {
  return v.write(out);
  }


  template <class V, class K >
  inline std::istream &operator>>(std::istream &in, Active_objects_update_vector<V, K> &v) {
  return v.read(in);
  }*/
//typename Updatable_Delaunay_triangulation_table_2<IK, KK>::Protected_array_pointer<T,n>
  
} //namespace CGAL;
#endif
