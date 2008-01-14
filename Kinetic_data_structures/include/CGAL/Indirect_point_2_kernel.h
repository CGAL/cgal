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
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_INDIRECT_KERNEL_H
#define  CGAL_INDIRECT_KERNEL_H
#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Tools/Label.h>
#include <CGAL/Tools/Counter.h>
//#include <CGAL/Kinetic/Cartesian_static_converter.h>

#define CGAL_MSA(Pred, pred) typedef Adaptor<typename SK::Pred##_2> Pred##_2; \
  Pred##_2 pred##_2_object() const					\
  {									\
    typename SK::Pred##_2 sp= SK::pred##_2_object();			\
    return Pred##_2(rep_, sp);						\
  }

CGAL_BEGIN_NAMESPACE




template <class SK, class Key = CGAL::Label<int> >
class Indirect_point_2_kernel: public SK
{



  typedef Indirect_point_2_kernel<SK, Key> This;

  class Rep: public Kinetic::Ref_counted<Rep > {
  public:
    typedef typename SK::Point_2 Pt;
    typedef typename std::vector<Pt> Container;
    
    Rep(const Rep &o ): points_(o.points_){}
    Rep(){}

    void swap(Container &sc) {
      CGAL_precondition(sc.size() == points_.size());
      std::swap(sc, points_);
    }

    Key new_point(const Pt p) {
      points_.push_back();
      return Key(points_.size()-1);
    }
    
    template <class It>
    void set_points(It b, It e) {
      CGAL_precondition(std::distance(b,e) == points_.size());
      std::copy(b,e, points_.begin());
    }
    
    typedef CGAL::Counter<int, Key> Key_iterator;
    Key_iterator begin() const {
      return Key_iterator(0);
    }
    Key_iterator end() const {
      return Key_iterator(points_.size());
    }
    
    typedef std::pair<Key_iterator, Key_iterator> KP;
    
    template <class It>
    KP new_points(It b, It e) {
      int osz= points_.size();
      points_.resize(std::distance(b,e)+ osz);
      std::copy(b,e, points_.begin() + osz);
      return KP(Key_iterator(osz), Key_iterator(number_of_points()));
    }
  
    unsigned int number_of_points() const {
      return points_.size();
    }
    
    const Pt& point(Key k) const {
      CGAL_precondition(k.index() < points_.size());
      return points_[k.index()];
    }
    
    const Container& container() const {
      return points_;
    }
    
    void swap(Rep &o) {
      std::swap(points_, o.points_);
    }
  protected:
    Container points_;
  };




  template <class Predicate>
  class Adaptor {
  public:
    Adaptor(typename Rep::Handle r,
	    Predicate pred=Predicate()): ik_(r), pred_(pred) {
    }

    typedef typename Predicate::result_type result_type;
    typedef Key argument_type;
    typedef argument_type first_argument_type;
    typedef argument_type second_argument_type;
    typedef argument_type third_argument_type;
    typedef argument_type fourth_argument_type;
    typedef argument_type fifth_argument_type;
    typedef typename Arity_traits<Predicate>::Arity Arity;
    
    result_type operator()(const first_argument_type &arg0) const
    {
      return pred_(ik_->point(arg0));
    }
    
    result_type operator()(const first_argument_type &arg0,
			   const second_argument_type &arg1) const
    {
      return pred_(ik_->point(arg0), ik_->point(arg1));
    }
    
    result_type operator()(const first_argument_type &arg0,
			   const second_argument_type &arg1,
			   const third_argument_type &arg2) const
    {
      return pred_(ik_->point(arg0), ik_->point(arg1),
		   ik_->point(arg2));
    }
    
    result_type operator()(const first_argument_type &arg0,
			   const second_argument_type &arg1,
			   const third_argument_type &arg2,
			   const fourth_argument_type &arg3) const
    {
      return pred_(ik_->point(arg0), ik_->point(arg1),
		   ik_->point(arg2), ik_->point(arg3));
    }
    
    result_type operator()(const first_argument_type &arg0,
			   const second_argument_type &arg1,
			   const third_argument_type &arg2,
			   const fourth_argument_type &arg3,
			   const fifth_argument_type &arg4) const
    {
      return pred_(ik_->point(arg0), ik_->point(arg1),
		   ik_->point(arg2), ik_->point(arg3),
		   ik_->point(arg4));
    }
    
  protected:
    
    typename Rep::Handle ik_;
    Predicate pred_;
  };




public:

  Indirect_point_2_kernel(SK sk= SK()): SK(sk), rep_(new Rep()) {
  }

 
  typedef Key Point_2;
  typedef typename SK::Point_2 Geometric_point_2;
  typedef SK Direct_kernel;

  const Direct_kernel direct_kernel_object() const {
    return static_cast<Direct_kernel>(*this);
  }

  struct Current_coordinates {
    Current_coordinates(typename Rep::Handle rep): rep_(rep){}
    typedef Point_2 argument_type;
    typedef Geometric_point_2 result_type;
    const Geometric_point_2 & operator()(Point_2 k) const {
      return rep_->point(k);
    }
    typename Rep::Handle rep_;
  };

  Current_coordinates current_coordinates_object() const {
    return Current_coordinates(rep_);
  }
  

  Point_2 new_point_2(Geometric_point_2 pt) const {
    return rep_->new_point(pt);
  }

  template <class It>
  void set_point_2s(It b, It e) const {
    rep_->set_points(b,e);
  }

  typedef typename Rep::Key_iterator Key_iterator;
  typedef typename Rep::KP Key_range;
  
  Key_iterator point_2s_begin() const {
    return rep_->begin();
  }

  Key_iterator point_2s_end() const {
    return rep_->end();
  }

  template <class It>
  Key_range new_point_2s(It b, It e) const {
    return rep_->new_points(b,e);
  }

  Point_2 ith_point(unsigned int i) const {
    CGAL_precondition( i < rep_->number_of_points());
    return Point_2(i);
  }

  typedef typename Rep::Container Swapable_container;
  
  void swap(Swapable_container &sc){
    rep_->swap(sc);
  }

  const Swapable_container& container() const {
    return rep_->container();
  }

  unsigned int number_of_point_2s() const {
    return rep_->number_of_points();
  }

  void swap(This &o) {
    rep_->swap(*o.rep_);
  }

  // for compatibility with InstantaneousKernel
  template <class T>
  void set_time(T) const {
  }

  This clone() const {
    return This(new Rep(*rep_), *this);
  }

  CGAL_MSA(Side_of_oriented_circle,side_of_oriented_circle);
  CGAL_MSA(Orientation,orientation);
  CGAL_MSA(Compare_x, compare_x);
  CGAL_MSA(Compare_y,compare_y);
  CGAL_MSA(Less_x, less_x);
  CGAL_MSA(Less_y, less_y);
  CGAL_MSA(Compare_distance, compare_distance);


  //CGAL_MSA(Compare_distance, compare_distance, 3);

  //CGAL_MSA(Side_of_oriented_sphere,side_of_oriented_sphere, 3);
  //CGAL_MSA(Orientation,orientation, 3);
  //CGAL_MSA(Compare_x,compare_x, 3);
  //CGAL_MSA(Compare_y,compare_y, 3);
  //CGAL_MSA(Compare_z,compare_z, 3);
  //CGAL_MSA(Compare_xyz,compare_xyz, 3);
  //CGAL_MSA(Less_x, less_x, 3);
  //CGAL_MSA(Less_y, less_y, 3);
  //CGAL_MSA(Less_z, less_z, 3);
  //CGAL_MSA(Coplanar_orientation, coplanar_orientation, 3);
  //CGAL_MSA(Coplanar_side_of_bounded_circle, coplanar_side_of_bounded_circle, 3);

protected:
  Indirect_point_2_kernel(typename Rep::Handle rep, SK sk= SK()): SK(sk), rep_(rep) {
  }
  mutable typename Rep::Handle rep_;
};
#undef CGAL_MSA
CGAL_END_NAMESPACE
#endif
