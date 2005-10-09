#ifndef CGAL_CARTESIAN_INSTANTANEOUS_KERNEL_REP_2_H
#define  CGAL_CARTESIAN_INSTANTANEOUS_KERNEL_REP_2_H
#include <CGAL/KDS/basic.h>
#include <CGAL/Handle_for.h>
#include <map>
#include <iostream>
#include <CGAL/KDS/Instantaneous_adaptor.h>

#define MSA(Pred, pred, d) typedef Instantaneous_adaptor<typename Static_kernel::Pred##_##d, Handle, Point_2> Pred##_##d;\
Pred##_##d pred##_##d##_object() const {\
typename Static_kernel::Pred##_##d sp= static_kernel().pred##_##d##_object();\
return Pred##_##d(*this, sp);\
}

#define TSO(name) typedef typename Static_kernel::name name;

CGAL_KDS_BEGIN_NAMESPACE;


//! A Kernel that makes a snapshot of a Kinetic motion look like a normal Kernel
/*!
  This is only valid for two dimensional predicates and things
*/
template <class Moving_point_table_d, class Static_kernel_t> 
class Cartesian_instantaneous_kernel_rep_2 {
  typedef Cartesian_instantaneous_kernel_rep_2<Moving_point_table_d, Static_kernel_t> This;
public:
  typedef Moving_point_table_d Point_table;
  typedef Static_kernel_t Static_kernel;
  typedef typename Point_table::Object::Coordinate::NT Time;

  class Handle: public CGAL::Handle_for<This> {
  public:
    Handle(const Point_table *mot,
	   const Time &cur_time= Time(0),
	   const Static_kernel &sk= Static_kernel()):
      Handle_for<This>(This(mot, cur_time, sk)){
    }
    //! To get around the last of non-const geom_traits() function in triangulation
    void set_time(const Time &cur_time) const {
      if (Ptr()->time_ != cur_time){
	Ptr()->time_=cur_time;
	Ptr()->cache_.clear();
      }
    }
    
    const Time &time() const {
      return Ptr()->time_;
    }
    const Static_kernel &static_kernel() const {
      return Ptr()->skernel_;
    }
    typedef typename Moving_point_table_d::Key Point_2;

    typename Static_kernel::Point_2 
    to_static(const Point_2 &k) const {
      if (Ptr()->cache_.find(k) == Ptr()->cache_.end()){
	Ptr()->cache_[k]=typename Static_kernel::Point_2(Ptr()->mot_->at(k).x()(time()), 
							 Ptr()->mot_->at(k).y()(time()));
      }
      return Ptr()->cache_[k];
    } 
    MSA(Side_of_oriented_circle,side_of_oriented_circle, 2);
    //MSA(Side_of_oriented_sphere,side_of_oriented_sphere, 3);
    MSA(Orientation,orientation, 2);
    //MSA(Orientation,orientation, 3);
    MSA(Compare_x, compare_x, 2);
    MSA(Compare_y,compare_y, 2);
    MSA(Less_x, less_x, 2);
    MSA(Less_y, less_y, 2); 
    TSO(Segment_2);
    TSO(Triangle_2);

  };
    
  Cartesian_instantaneous_kernel_rep_2(const Point_table *mot,
				       const Time &cur_time,
				       const Static_kernel &sk):
    time_(cur_time),
    skernel_(sk),
    mot_(mot){
  }

 
protected:
  mutable Time time_;
  Static_kernel skernel_;
  const Point_table *mot_;
  mutable std::map<typename Point_table::Key,
		   typename Static_kernel::Point_2> cache_;
};

#undef MSA
#undef TSO

CGAL_KDS_END_NAMESPACE;

#endif
