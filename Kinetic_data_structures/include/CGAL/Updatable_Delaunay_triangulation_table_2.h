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

#ifndef NDEBUG
#include <CGAL/Polynomial/Interval_polynomial.h>
#endif

CGAL_BEGIN_NAMESPACE;


unsigned int kinetic_certificates_;
unsigned int unfailing_kinetic_certificates_;
unsigned int certificate_advances_;
unsigned int static_certificates_;
//unsigned int num_events_=0;
unsigned int stat_bad_edges_;
unsigned int stat_num_edges_;
unsigned int num_interpolations_;
unsigned int constant_filtered_;
unsigned int interval_filtered_;
unsigned int interval2_filtered_;
unsigned int interval3_filtered_;
unsigned int uncertain_exact_computations_;
unsigned int exact_current_time_certificates_;
unsigned int comparison_certificates_;
unsigned int ipolynomials_evaluated_;
unsigned int ipolynomials_interval_evaluated_;
unsigned int icertificates_generated_;
unsigned int deriv_filtered_;

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

#ifndef ALL_MOVE
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
      start_time_=-1;
    }
    Coef_data(){}
    const INT* operator[](int i) const {
      return c_[i];
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
#ifndef ALL_MOVE
    active_(ik.number_of_point_2s(), false),
    coef_cache_(ik.number_of_point_2s()),
#endif
    editing_(NOT), 
    notifying_(false),
    soc_(kk.positive_side_of_oriented_circle_2_object()),
    ik_(ik)
#ifdef ALL_MOVE
    , fk_(ik.clone())
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
      Static_point_2 ip= initial(key);
      typedef typename Exact_point_2::Coordinate MF;
      typedef typename MF::NT NT;
      MF mf[2];
      CGAL_UD_DEBUG("Computing exact motions for " << key << " from " << ip );
     
#ifdef ALL_MOVE
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
          c[1]=(NT(ip[i])-NT(fp[i]))/NT(start_time-1);
          c[0]=NT(fp[i])-c[1];
          mf[i]=MF(c, c+2);
        }
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
    static_certificates_=0;
    //unsigned int num_events_=0;
    stat_bad_edges_=0;
    num_interpolations_=0;
    constant_filtered_=0;
    interval_filtered_=0;
    interval2_filtered_=0;
    interval3_filtered_=0;
    kinetic_certificates_=0;
    unfailing_kinetic_certificates_=0;
    certificate_advances_=0;
    uncertain_exact_computations_=0;
    exact_current_time_certificates_=0;
    comparison_certificates_=0;
    icertificates_generated_=0;
    ipolynomials_evaluated_=0;
    ipolynomials_interval_evaluated_=0;
    deriv_filtered_=0;
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


  void write_stats(std::ostream &out) {
    out << "The triangululation had " << stat_bad_edges_ << " bad edges" << std::endl;
    out << "Points were interpolated on " << num_interpolations_ << " occasions" << std::endl;
    out << "There were " << static_certificates_ 
	<< " static computations for " << stat_num_edges_
	<< " edges" << std::endl;
    out << "Interval polynomials were generated " << icertificates_generated_ 
	<< " and evaluated " << ipolynomials_evaluated_ << " times " 
	<< ipolynomials_interval_evaluated_ << " were not on point intervals" << std::endl;
    out << "Whole interval filtering removed " << interval2_filtered_
	<< " derivitive filtering remove " << deriv_filtered_ 
	<< " isolation remove " << interval3_filtered_ << " certificates "
	<< std::endl;
    out << "There remained " << kinetic_certificates_ << " certificates that were "
	<< "computed which were advanced " << certificate_advances_ 
	<< " times and " << unfailing_kinetic_certificates_ 
	<< " of which could have been filtered" << std::endl;
    out << "Of these " << uncertain_exact_computations_ << " were caused "
	<< " by inability to determine if the root was there and " 
	<< exact_current_time_certificates_ 
	<< " were caused to make the current time exact " 
	<< " and " << comparison_certificates_ 
        << " were caused in order to compare." << std::endl;
  }

#ifdef MOVE_ALL
  void coef(Key a, int i, INT r[2]) const {
    r[0]= CGAL::to_interval(initial(a)[i]);
    r[1]= INT(CGAL::to_interval(final(a)[i]))-r[0];
  }
#endif

#ifndef MOVE_ALL
  bool is_active(Key k) const {
    return active_[k.to_index()];
  }


  void activate(double t, Key k) {
    //CGAL_precondition(t.inf() == t.sup());
    //CGAL_precondition(active_[k.to_index()]==false);
    if (!active_[k.to_index()]) {
      active_[k.to_index()]=true;
      //start_time_[k.to_index()]= t.inf();
      INT stm1=(t-1);

      INT xii= CGAL::to_interval(initial(k)[0]);
      INT xfi= CGAL::to_interval(final(k)[0]);
      INT xc0= (t*xfi - xii)/stm1;
      INT xc1=(xii-xfi)/stm1;

      INT yii= CGAL::to_interval(initial(k)[1]);
      INT yfi= CGAL::to_interval(final(k)[1]);
      INT yc0= (t*yfi - yii)/stm1;
      INT yc1=(yii-yfi)/stm1;
      coef_cache_[k.to_index()]= Coef_data(t, xc0, xc1, yc0, yc1);
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
  template <class CNT>
  CNT  incircle(CNT ax, CNT ay, CNT bx, CNT by, 
		CNT cx, CNT cy, CNT dx, CNT dy,
		bool stat=false) const {
    if (stat) ++static_certificates_;
    CNT qpx = bx - ax;
    CNT qpy = by - ay;
    CNT rpx = cx - ax;
    CNT rpy = cy - ay;
    CNT tpx = dx - ax;
    CNT tpy = dy - ay;
    CNT det=CGAL::det2x2_by_formula(qpx*tpy - qpy*tpx, tpx*(dx - bx) 
				    + tpy*(dy - by),
				    qpx*rpy - qpy*rpx, rpx*(cx - bx) 
				    + rpy*(cy - by));
    return det;
  }

  void  incircle_p(const INT ax[], const INT ay[], const INT bx[], const INT by[], 
		   const INT cx[], const INT cy[], const INT dx[], const INT dy[],
		   INT *ret) const {
    ++icertificates_generated_;
    INT xb0=bx[0]-ax[0], xb1= bx[1]-ax[1];
    INT yb0=by[0]-ay[0], yb1= by[1]-ay[1];
    INT xc0=cx[0]-ax[0], xc1=cx[1]-ax[1];
    INT yc0=cy[0]-ay[0], yc1=cy[1]-ay[1];
    INT xd0=dx[0]-ax[0], xd1=dx[1]-ax[1];
    INT yd0=dy[0]-ay[0], yd1=dy[1]-ay[1];
    INT t37 = 2*yc0;		
    INT t38 = xb0*t37;			
    INT t51 = 2*xb1;		
    INT t50 = 2*yb1;		
    INT t44 = -CGAL::square(yd0)-CGAL::square(xd0);
    INT t49 = -2*xc0;
    INT t48 = xc0*xb0;
    INT t47 = yb0*yc0;
    INT t46 = xb0*yc1;
    INT t45 = -CGAL::square(xc0)-CGAL::square(yc0);
    INT t43 =CGAL::square( xc1)+CGAL::square(yc1);
    INT t42 = CGAL::square(xd1)+CGAL::square(yd1);
    INT t41 = -CGAL::square(xb0)-CGAL::square(yb0);
    INT t40 = -CGAL::square(yb1)-CGAL::square(xb1);
    INT t39 = yb0*t49;
    INT t36 = 2*yc1;
    INT t35 = xb1*t36;
    INT t34 = t42*yb1;
    INT t33 = t42*xb0;
    INT t32 = t42*yb0;
    INT t30 = (-yd1+yc1)*t37;
    INT t29 = (-yc0+yb0)*t36;
    INT t28 = 2*(-yb1+yd1)*yb0;
    INT t27 = (xc0-xb0)*t51;
    INT t26 = -t40*yc1-t43*yb1;
    INT t25 = -t43*yb0-t40*yc0;
    INT t24 = -t41*yc0+t45*yb0;
    INT t23 = xc1*yb0+yb1*xc0-yc0*xb1-t46;
    ret[0] = (t44*yc0-t45*yd0)*xb0+(t41*yd0-t44*yb0)*xc0+t24*xd0;
    ret[1] = (-t45*xb0+t41*xc0)*yd1+(-t41*yc1+(-xd1+xb1)*t38+t45*yb1+((yb1-yc1)*t37+2*(xd1-xc1)*xc0)*yb0)*xd0+t24*xd1+(-t45*xb1+xb0*t30+t41*xc1+(2*(-xb1+xc1)*xb0+t28)*xc0)*yd0-t44*t23;
    ret[2] = -yc0*t33+(xb1*t38+xc1*t39+(-2*t47-t41)*yc1+(2*t47+t45)*yb1)*xd1+xc0*t32+(t46*t51+(xc1*t49+t29)*yb1+2*t23*xd1+t25)*xd0+(yb1*t39+yc1*t38+(-2*t48-t45)*xb1+(2*t48+t41)*xc1)*yd1+(xb1*t30+(-2*yd1*yc1+t43)*xb0+(yd1*t50+t40)*xc0+(t27+t28)*xc1)*yd0-t44*(xc1*yb1-xb1*yc1);
    ret[3] =
      -yc1*t33+xc0*t34+(-t42*yc0+t43*yd0)*xb1+(t40*yd0+t32)*xc1+t26*xd0+(t43*xb0+(-yd0+yc0)*t35+t40*xc0+(t27+(yd0-yb0)*t50)*xc1)*yd1+((-xd0+xb0)*t35+(t29+2*(-xc0+xd0)*xc1)*yb1+t25)*xd1;
    ret[4] = (-t42*yc1+t43*yd1)*xb1+(t40*yd1+t34)*xc1+t26*xd1;
  }

  void certificate_function(Key a, Key b, Key c, Key d, INT *ret) const {
#ifdef MOVE_ALL
    INT a0[2]; coef(a, 0, a0);
    INT a1[2]; coef(a, 0, a1);
    INT b0[2]; coef(b, 0, b0);
    INT b1[2]; coef(b, 0, b1);
    INT c0[2]; coef(c, 0, c0);
    INT c1[2]; coef(c, 0, c1);
    INT d0[2]; coef(d, 0, d0);
    INT d1[2]; coef(d, 0, d1);
   
    incircle_p(a0, a1, b0, b1, c0, c1, d0, d1,
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
    /*
      typedef POLYNOMIAL::Interval_polynomial IP;
      IP ax(coef_cache_[a.to_index()][0], coef_cache_[a.to_index()][0]+2);
      IP ay(coef_cache_[a.to_index()][1], coef_cache_[a.to_index()][1]+2);
      IP bx(coef_cache_[b.to_index()][0], coef_cache_[b.to_index()][0]+2);
      IP by(coef_cache_[b.to_index()][1], coef_cache_[b.to_index()][1]+2);
      IP cx(coef_cache_[c.to_index()][0], coef_cache_[c.to_index()][0]+2);
      IP cy(coef_cache_[c.to_index()][1], coef_cache_[c.to_index()][1]+2);
      IP dx(coef_cache_[d.to_index()][0], coef_cache_[d.to_index()][0]+2);
      IP dy(coef_cache_[d.to_index()][1], coef_cache_[d.to_index()][1]+2);
      IP ic= incircle(ax, ay, bx, by, cx, cy, dx, dy, false);
      CGAL_UD_DEBUG(ic << std::endl);*/
    
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
    CGAL_UD_DEBUG("Exact func is " << ec
		  << std::endl );
    /*ret[0]=ic[0];
      ret[1]=ic[1];
      ret[2]=ic[2];
      ret[3]=ic[3];*/
  }

   
  void point_changed(Key ){
  }

  Exact_certificate compute_exact_failure_time(Cert_tuple ct, 
					       Exact_time et) const {
    CGAL_UD_DEBUG("Computing exact time for " << ct 
		  << " from " << et  << std::endl);
    ++kinetic_certificates_;
    Exact_certificate ec= soc_(exact_point(ct[0]),
			       exact_point(ct[1]),
			       exact_point(ct[2]), 
			       exact_point(ct[3]),
			       et, 1);
    if (ec.failure_time() > 1) {
      ++unfailing_kinetic_certificates_;
    }
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
      ++certificate_advances_;
      oc.pop_failure_time();
    }
      
    while (oc.failure_time() < et) {
      oc.pop_failure_time();
      oc.pop_failure_time();
      CGAL_UD_DEBUG("Advancing exact time for " << ct 
		    << " from " << oc.failure_time() 
		    << std::endl);
      ++certificate_advances_;
      ++certificate_advances_;
      if (oc.failure_time() > 1) {
	++unfailing_kinetic_certificates_;
      }
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
    INT det= incircle(cur_interval(a,ct,0),
		      cur_interval(a,ct,1),
		      cur_interval(b,ct,0),
		      cur_interval(b,ct,1),
		      cur_interval(c,ct,0),
		      cur_interval(c,ct,1),
		      cur_interval(d,ct,0),
		      cur_interval(d,ct,1),
		      true);
      
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

  static INT evaluate_ipoly(const INT *coefs, const INT t) {
    CGAL_precondition(t.inf() >=0);
    CGAL_precondition(t.sup() <=1);
    CGAL_precondition(t.inf() <=1);
    {
      INT tt(t.inf());
      CGAL_assertion(tt.inf()==tt.sup());
    }
    ++ipolynomials_evaluated_;
    if (t.inf() != t.sup()) {
      ++ipolynomials_interval_evaluated_;
    }
    INT cum=coefs[0];
    cum += coefs[1]*t;
    INT tm=CGAL::square(t);
    cum += coefs[2]*tm;
    tm *= t;
    cum += coefs[3]*tm;
    tm *= t;
    cum += coefs[4]*tm;
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
    
    if (!ret) {
      ++interval2_filtered_;
    } 
    return ret;
  }


  static CGAL::Sign sign_at(const INT *coefs, const INT t) {
    INT cum= evaluate_ipoly(coefs, t);
    if (cum.inf() > 0) return CGAL::POSITIVE;
    else if (cum.sup() < 0) return CGAL::NEGATIVE;
    else return CGAL::ZERO;
  }

  struct Certificate_evaluator {
    Certificate_evaluator(const INT* coefs): coefs_(coefs) {
    }
    
    CGAL::Sign operator()(INT t) const {
      return sign_at(coefs_, t);
    }
    const INT *coefs_;
  };

  template <class F>
  Isolate_result isolate_failure(F f, double lb, double ub, 
				 bool starts_positive, int rem_depth,
				 int NS,
				 std::pair<double,double> &ret) const {
    //const int NS=10;
    double ld= lb;
    //const double growth=1.5;
    double step= (ub-lb)/(NS-1); //std::pow(growth,NS-1);
    //if (step <0) std::cerr << "Step is " << step << std::endl;
    std::pair<double,double> failures(std::numeric_limits<double>::infinity(),
				      -std::numeric_limits<double>::infinity());
    bool has_zero=false;
    //bool has_negative=false;
    bool has_positive = starts_positive;
    Isolate_result certain=POSSIBLE_FAILURE;
      
    /*if (!has_positive) {
      CGAL::Sign sn= sign_at(ct[0],ct[1],ct[2],ct[3], INT(lb, lb));
      if (sn== CGAL::POSITIVE) has_positive=true;
      }*/


    CGAL_UD_DEBUG( "Interval is " << lb << " to " << ub 
		   << " with initial step " << step << "(" << rem_depth << ")" 
		   << std::endl);

    for (int i=0; i< NS; ++i) {
      double nld= std::min(ld+step, ub);
      if (nld == ld) break;

      //step *= growth;

      CGAL_assertion(i != NS-1 || nld == ub);

      INT ci(ld, nld);
      ld= nld;
      CGAL::Sign csn= f(ci);
      CGAL_UD_DEBUG("Sign on " << ci << " is " << csn << "(" << i << ")" 
		    << std::endl);
      if (csn == CGAL::NEGATIVE) {
	//CGAL_assertion(has_zero); 
	// could start negative in the (large) initial interval
	//has_negative=true;
	if (has_positive) {
	  certain=CERTAIN_FAILURE;
	  break;
	}
      } else if (csn == CGAL::ZERO) {
	has_zero=true;
	//if (has_zero){
	failures = join(failures, ci);
	INT ci(nld, nld);
	CGAL::Sign sn= f(ci);
	CGAL_UD_DEBUG("Sign at " << ci << " is " << sn << "(" << i << ")" 
		      << std::endl);
	if (sn == CGAL::NEGATIVE) {
	  CGAL_assertion(csn != CGAL::POSITIVE);
	  if (has_positive) {
	    certain=CERTAIN_FAILURE;
	    break;
	  }
	} else if (sn == CGAL::POSITIVE) {
	  CGAL_assertion(csn != CGAL::NEGATIVE);
	  has_positive=true;
	}
	  
      } else if (csn == CGAL::POSITIVE) {
	has_positive=true;
      }
      //if (i==0) {
    }







    if ( has_zero) {
      if ((certain==CERTAIN_FAILURE 
	   && (lb != failures.first || starts_positive)) 
	  || rem_depth == 0
	  || lb == failures.first && ub == failures.second
	  /* || .9*(ub-lb) < (failures.second - failures.first)*/) {
	CGAL_UD_DEBUG(  "Not recursing with " 
			<< failures.first << " " << failures.second 
			<< " because " << has_zero << has_positive 
			<< certain << std::endl);
	ret= failures;
	return certain;
      } else {
	CGAL_UD_DEBUG( "Recursing with " 
		       << failures.first << " " << failures.second 
		       << " because " << has_zero << has_positive 
		       << certain << std::endl);
	double nwid = failures.second-failures.first;
	double wid = ub-lb;
	if (nwid > wid/2.0) NS*=2;
	if (nwid < wid/4.0) NS/=2;
	return isolate_failure(f, failures.first, failures.second,
			       has_positive,
			       rem_depth-1,
			       NS, ret);
      }
    } else {
      ++interval3_filtered_;
      //ret= std::pair<double,double>(1,-1);
      return NO_FAILURE;
    }
     
  }
 
  Isolate_result isolate_failure(INT *f, double lb, double ub,
				 bool starts_positive,
                                 std::pair<double,double> &ret) const {
    //CGAL::Protect_FPU_rounding<true> prot;
    INT v= evaluate_ipoly(f, lb);
    double c=lb;
    CGAL_UD_DEBUG("Const is " << v );
#ifndef NDEBUG
    double olb= lb;
    double oub=ub;
#endif
    if (v.inf() >= 0) {
      INT deriv= INT(0,1)*f[1]+INT(0,2)*f[2]+INT(0,3)*f[3]+INT(0,4)*f[4];
      CGAL_UD_DEBUG(" and derivitive is "<<  deriv 
		    << " making the step " << v/CGAL::abs(deriv) << std::endl;);
      
      INT slope;
      if (deriv.inf() >= 0) {
	++deriv_filtered_;
	return NO_FAILURE;
      }
      else slope= -deriv.inf();
      do {
        INT step= v/slope;
        lb=c;
        c= c+step.inf();
        if (c >= ub) {
          ++deriv_filtered_;
          return NO_FAILURE;
        }
        v= evaluate_ipoly(f,c);
      } while (v.inf() >0);
      ub= c;
    }
    CGAL_UD_DEBUG("lb changed from " << olb << " to " << lb << std::endl);
    CGAL_UD_DEBUG("ub changed from " << oub << " to " << lb << std::endl);

    return isolate_failure(Certificate_evaluator(f), lb, ub, starts_positive, 3, 20, ret);
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
      std::sort(changed_objects.begin(), changed_objects.end());
      for (typename Subscribers::iterator it= subscribers_.begin(); it != subscribers_.end(); ++it) {
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

#ifndef MOVE_ALL
  Active active_;
  std::vector<Coef_data> coef_cache_;
  double next_activation_;
  std::vector<Key> changed_objects_;
#endif

public:
  typename Kinetic_kernel::Positive_side_of_oriented_circle_2 soc_;
  Indirect_kernel ik_, fk_;
  

};

/*template <class V, class K >
  inline std::ostream &operator<<(std::ostream &out, const Active_objects_update_vector<V, K> &v) {
  return v.write(out);
  }


  template <class V, class K >
  inline std::istream &operator>>(std::istream &in, Active_objects_update_vector<V, K> &v) {
  return v.read(in);
  }*/

CGAL_END_NAMESPACE;
#endif
