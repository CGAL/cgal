// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kinetic_data_structures/include/CGAL/Kinetic/Active_objects_vector.h $
// $Id: Active_objects_vector.h 31436 2006-06-05 12:44:12Z drussel $
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

CGAL_BEGIN_NAMESPACE;






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
unsigned int stat_number_of_edges_;
unsigned int stat_number_of_bad_edges_;
double stat_total_interval_width_;
unsigned int stat_number_of_intervals_;


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
}

/*!
 */
template <class Indirect_kernel, class Kinetic_kernel >
class Updatable_delaunay_triangulation_table_2:
  public Kinetic::Ref_counted<Updatable_delaunay_triangulation_table_2<Indirect_kernel, Kinetic_kernel> >
{
public:
  typedef typename Kinetic_kernel::Point_2 Exact_point_2;
  typedef typename Indirect_kernel::Point_2 Key;
  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;
  typedef CGAL::Interval_nt_advanced INT;
  typedef typename Kinetic_kernel::Certificate Exact_certificate;
  typedef typename Exact_certificate::Time Exact_time;
  typedef Update_cert_tuple<Key> Cert_tuple;






 
protected:
  typedef Updatable_delaunay_triangulation_table_2<Indirect_kernel, Kinetic_kernel> This;
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
    coef_cache_(ik.number_of_point_2s())
#else
    fk_(ik.clone())
#endif
 {
    reset(false);
    clear_stats();
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
        double start_time= coef_cache_[key.to_index()].start_time();
        
        for (unsigned int i=0; i< 2; ++i){
          NT c[2];
          c[1]=(NT(ip[i])-NT(fp[i]))/(NT(start_time)-NT(1));
          c[0]=NT(fp[i])-c[1];
          mf[i]=MF(c, c+2);
        }
	CGAL_postcondition(mf[0](start_time) == NT(ip[0]));
	CGAL_postcondition(mf[1](start_time) == NT(ip[1]));
      } else {
        mf[0] = MF(ip.x());
        mf[1] = MF(ip.y());
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
    CGAL_assertion(storage_.size() > static_cast<unsigned int>(key.to_index()));
    storage_[key.to_index()]=new_value;
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
      return coef_cache_[a.to_index()][i][0] + ct*coef_cache_[a.to_index()][i][1];
    } else {
      return CGAL::to_interval(initial(a)[i]);
    }

#endif
  }

  INT cur_interval(Key a, Key base, INT ct, int i) const {
    CGAL_precondition(ct.inf()>=0 && ct.sup() <=1);
#ifdef MOVE_ALL
    Protected_array<INT,2> c, bc;
    coef(a,i,c);
    coef(base, i, bc);
    return c[0]- bc[0] + ct*(c[1]-bc[1]);
#else
    if (is_active(a) && is_active(base)) {
      return coef_cache_[a.to_index()][i][0]
	- coef_cache_[base.to_index()][i][0] 
	+ ct*(coef_cache_[a.to_index()][i][1]- coef_cache_[base.to_index()][i][1]);
    } else if (is_active(a)) {
      return coef_cache_[a.to_index()][i][0] - INT(initial(base)[i]) + ct*coef_cache_[a.to_index()][i][1];
    } else if (is_active(base)) {
      return INT(initial(a)[i])- coef_cache_[base.to_index()][i][0] - ct*coef_cache_[base.to_index()][i][1];
    } else {
      return INT(initial(a)[i]) - INT(initial(base)[i]);
    }

#endif
  }


  void write_stats(std::ostream &out) {
    std::cout << "Triangulation had " << stat_number_of_edges_ << " edges "
	      << " and " << stat_number_of_bad_edges_ << " bad edges." << std::endl;
    std::cout << "Points were interpolated " << stat_number_of_interpolations_ << " times" << std::endl;
    std::cout << "Certificates were requested " << stat_certificate_computations_ << " times "
	      << " and advanced " << stat_certificate_advances_ << " times" << std::endl;
    std::cout << "Point predicate evals " << stat_point_predicate_evaluations_ << std::endl;
    std::cout << "Interval predicate evals " << stat_interval_predicate_evaluations_ << std::endl;
    std::cout << "Interval certificate functions " << stat_interval_certificate_functions_ << std::endl;
    std::cout << "\t deriv filtered " << stat_deriv_filter_calls_ << std::endl;
    std::cout << "\t newton isolated " << stat_newton_isolate_calls_ << std::endl;
    std::cout << "\t newton refined " << stat_newton_refine_calls_ << std::endl;
    std::cout << "\t evaluations " << stat_interval_certificate_function_evaluations_ << std::endl;
    std::cout << "\t interval evaluations " << stat_interval_certificate_function_interval_evaluations_ << std::endl;
    std::cout << "\t interval average final width " << stat_total_interval_width_/double(stat_number_of_intervals_) << std::endl;
    std::cout << "Exact points " << stat_exact_points_ << std::endl;
    std::cout << "Exact certificate functions " << stat_exact_certificate_functions_ << std::endl;
    std::cout << "\t from compare " << stat_exact_certificate_functions_from_compare_ << std::endl;
    std::cout << "\t from compare curt " << stat_exact_certificate_functions_from_compare_curt_ << std::endl;
    std::cout << "\t from advance " << stat_exact_certificate_functions_from_advance_ << std::endl;
    std::cout << "\t unfailing " << stat_unfailing_exact_certificate_functions_ << std::endl;
    std::cout << "\t advanced " << stat_exact_certificate_function_advances_ << std::endl;
  }

#ifdef MOVE_ALL
  void coef(Key a, int i, Protected_array_pointer<INT,2> r) const {
    r[0]= CGAL::to_interval(initial(a)[i]);
    r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];
  }
#endif

#ifndef MOVE_ALL
  bool is_active(Key k) const {
    return active_[k.to_index()];
  }

  void preactivate( Key k) {
    INT t= next_activation_;
    INT stm1=(t-1);
    
    INT xii= CGAL::to_interval(initial(k)[0]);
    INT xfi= CGAL::to_interval(final(k)[0]);
    INT xc0= (t*xfi - xii)/stm1;
    INT xc1=(xii-xfi)/stm1;
    
    INT yii= CGAL::to_interval(initial(k)[1]);
    INT yfi= CGAL::to_interval(final(k)[1]);
    INT yc0= (t*yfi - yii)/stm1;
    INT yc1=(yii-yfi)/stm1;
    coef_cache_[k.to_index()]= Coef_data(next_activation_, xc0, xc1, yc0, yc1);
  }

  bool is_activating(Key k) const {
    return coef_cache_[k.to_index()].start_time() == next_activation_;
  }

  void activate(double t, Key k) {
    //CGAL_precondition(t.inf() == t.sup());
    //CGAL_precondition(active_[k.to_index()]==false);
    if (!active_[k.to_index()]) {
      active_[k.to_index()]=true;
      changed_objects_.push_back(k);
      exact_points_.erase(k);
      //aot_->set(k, kp);
    }
  }
#endif

  void set_final_kernel(Indirect_kernel &fk){
    fk_=fk;
  }

  void reset(bool clear=true) {
    //activate_time_=-1;
#ifndef MOVE_ALL
    next_activation_=1.0;
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
    CNT det=CGAL::det2x2_by_formula(qpx*tpy - qpy*tpx, 
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
    CNT det=CGAL::det2x2_by_formula(qpx*tpy - qpy*tpx,
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
#if 1
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
#else
    NT xb0=bx[0]-ax[0], xb1= bx[1]-ax[1];
    NT yb0=by[0]-ay[0], yb1= by[1]-ay[1];
    NT xc0=cx[0]-ax[0], xc1=cx[1]-ax[1];
    NT yc0=cy[0]-ay[0], yc1=cy[1]-ay[1];
    NT xd0=dx[0]-ax[0], xd1=dx[1]-ax[1];
    NT yd0=dy[0]-ay[0], yd1=dy[1]-ay[1];
    NT t37 = 2*yc0;		
    NT t38 = xb0*t37;			
    NT t51 = 2*xb1;		
    NT t50 = 2*yb1;		
    NT t44 = -CGAL::square(yd0)-CGAL::square(xd0);
    NT t49 = -2*xc0;
    NT t48 = xc0*xb0;
    NT t47 = yb0*yc0;
    NT t46 = xb0*yc1;
    NT t45 = -CGAL::square(xc0)-CGAL::square(yc0);
    NT t43 = CGAL::square( xc1)+CGAL::square(yc1);
    NT t42 = CGAL::square(xd1)+CGAL::square(yd1);
    NT t41 = -CGAL::square(xb0)-CGAL::square(yb0);
    NT t40 = -CGAL::square(yb1)-CGAL::square(xb1);
    NT t39 = yb0*t49;
    NT t36 = 2*yc1;
    NT t35 = xb1*t36;
    NT t34 = t42*yb1;
    NT t33 = t42*xb0;
    NT t32 = t42*yb0;
    NT t30 = (-yd1+yc1)*t37;
    NT t29 = (-yc0+yb0)*t36;
    NT t28 = 2*(-yb1+yd1)*yb0;
    NT t27 = (xc0-xb0)*t51;
    NT t26 = -t40*yc1-t43*yb1;
    NT t25 = -t43*yb0-t40*yc0;
    NT t24 = -t41*yc0+t45*yb0;
    NT t23 = xc1*yb0+yb1*xc0-yc0*xb1-t46;
    ret[0] = (t44*yc0-t45*yd0)*xb0+(t41*yd0-t44*yb0)*xc0+t24*xd0;
    ret[1] = (-t45*xb0+t41*xc0)*yd1+(-t41*yc1+(-xd1+xb1)*t38+t45*yb1+((yb1-yc1)*t37+2*(xd1-xc1)*xc0)*yb0)*xd0+t24*xd1+(-t45*xb1+xb0*t30+t41*xc1+(2*(-xb1+xc1)*xb0+t28)*xc0)*yd0-t44*t23;
    ret[2] = -yc0*t33+(xb1*t38+xc1*t39+(-2*t47-t41)*yc1+(2*t47+t45)*yb1)*xd1+xc0*t32+(t46*t51+(xc1*t49+t29)*yb1+2*t23*xd1+t25)*xd0+(yb1*t39+yc1*t38+(-2*t48-t45)*xb1+(2*t48+t41)*xc1)*yd1+(xb1*t30+(-2*yd1*yc1+t43)*xb0+(yd1*t50+t40)*xc0+(t27+t28)*xc1)*yd0-t44*(xc1*yb1-xb1*yc1);
    ret[3] =
      -yc1*t33+xc0*t34+(-t42*yc0+t43*yd0)*xb1+(t40*yd0+t32)*xc1+t26*xd0+(t43*xb0+(-yd0+yc0)*t35+t40*xc0+(t27+(yd0-yb0)*t50)*xc1)*yd1+((-xd0+xb0)*t35+(t29+2*(-xc0+xd0)*xc1)*yb1+t25)*xd1;
    ret[4] = (-t42*yc1+t43*yd1)*xb1+(t40*yd1+t34)*xc1+t26*xd1;
    //CGAL_UD_DEBUG("New function is " << ret[0] << " + " << ret[1] << "*t+ " << ret[2] << "*t^2+ " << ret[3] << "*t^3+ " << ret[4] << "*t^4" << std::endl);
    //CGAL_UD_DEBUG("Old function is " << oret[0] << " + " << oret[1] << "*t+ " << oret[2] << "*t^2+ " << oret[3] << "*t^3+ " << oret[4] << "*t^4" << std::endl);
#endif
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
      Protected_array<NT, 2> a0; a0[0]=1; a0[1]=1;
      Protected_array<NT, 2> a1; a1[0]=1; a1[1]=1;
      Protected_array<NT, 2> b0; b0[0]=1; b0[1]=1;
      Protected_array<NT, 2> b1; b1[0]=1; b1[1]=1;
      Protected_array<NT, 2> c0; c0[0]=1; c0[1]=1;
      Protected_array<NT, 2> c1; c1[0]=1; c1[1]=1;
      Protected_array<NT, 2> d0; d0[0]=1; d0[1]=1;
      Protected_array<NT, 2> d1; d1[0]=1; d1[1]=1;
      /*r[0]= CGAL::to_interval(initial(a)[i]);
	r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];*/
      Protected_array<NT, 5> ret;
   
      incircle_p(P(a0), P(a1), P(b0), P(b1), P(c0), P(c1), P(d0), P(d1),
		 RP(ret));
      std::cout << "Static bounds are " << ret[0].bound() << "+/- " << ret[0].error() 
		<< ", " << ret[1].bound() << "+/- " << ret[1].error() 
		<< ", " << ret[2].bound() << "+/- " << ret[2].error() 
		<< ", " << ret[3].bound() << "+/- " << ret[3].error() 
		<< ", " << ret[4].bound() << "+/- " << ret[4].error() << std::endl;
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
    {
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
#ifdef MOVE_ALL
#ifndef NDEBUG
    incircle_static(a,b,c,d);
#endif
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
    incircle_p(coef_cache_[a.to_index()][0],
               coef_cache_[a.to_index()][1],
               coef_cache_[b.to_index()][0],
	       coef_cache_[b.to_index()][1],
	       coef_cache_[c.to_index()][0],
	       coef_cache_[c.to_index()][1],
               coef_cache_[d.to_index()][0],
               coef_cache_[d.to_index()][1],
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
    if (D==0) return cum;
    cum += coefs[1]*t;
    if (D==1) return cum;
    INT tm=t;
    for (unsigned int i=2; i<=D-1; ++i) {
      tm *= t;
      cum += coefs[i]*tm;
    }
    return cum;
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
    FN<D-1> prime() const {
      FN<D-1> ret;
      for (unsigned int i=1; i<=D; ++i) {
	ret[i-1]= INT(static_cast<double>(i))*P::operator[](i);
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
    INT vp= fp(INT(lb,ub));
    INT step=0;
    int nsteps=0;
    if (vp.inf() < 0) {
      vp = INT(-vp.inf());
      do {
	INT v= f(clb);
	if (v.inf() <= 0) {
	  CGAL_UD_DEBUG("Deriv return " << clb << " from " << lb << " with " << nsteps << std::endl);
	  return clb;
	}
	step = v/vp;
	clb += step.inf();
      } while (clb < ub && step.inf() > .01);
      
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
  Isolate_result Newton_isolate(F f, FP fp, FPP fpp, double lb, double ub, INT &ret) const {
    ++stat_newton_isolate_calls_;
    CGAL_UD_DEBUG("Newton isolate " <<  lb <<  " ... " << ub << std::endl);
    INT flb = f(lb);
    if (flb.inf() <= 0) {
      CGAL_UD_DEBUG( "Returning possible failure at start " << std::endl);
      return POSSIBLE_FAILURE;
    } else if (lb==ub) {
      CGAL_UD_DEBUG("Empty interval " << std::endl);
      return NO_FAILURE;
    } else {
      lb= nextafter(lb, std::numeric_limits<double>::max());
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


#ifndef MOVE_ALL  
  double next_activation() const {
    return next_activation_;
  }
  void set_next_activation(double f) {
    next_activation_=f;
  }
#endif
  
  double start_time() const {
    return start_time_;
  }
  void set_start_time(double f) {
    start_time_=f;
  }

  const Indirect_kernel& final_kernel_object() const {
    return fk_;
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

#ifndef MOVE_ALL
  Active active_;
  std::vector<Coef_data> coef_cache_;
  double next_activation_;
  std::vector<Key> changed_objects_;
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
  
CGAL_END_NAMESPACE;
#endif
