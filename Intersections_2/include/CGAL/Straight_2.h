// Copyright (c) 2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_STRAIGHT_2_H
#define CGAL_STRAIGHT_2_H

#include <CGAL/Line_2_Line_2_intersection.h>
#include <CGAL/squared_distance_utils.h>
#include <CGAL/Kernel/global_functions_internal_2.h>

namespace CGAL {

namespace internal {

class Straight_2_base_ {
public:
    enum states {EMPTY, POINT, SEGMENT, RAY, LINE};
    enum bound_states {NO_UNBOUNDED=0, MIN_UNBOUNDED=1, MAX_UNBOUNDED=2,
                        BOTH_UNBOUNDED = 3, LINE_EMPTY = 4};
protected:
                            Straight_2_base_() ;
    int                     main_dir_;  // support_ is x or y directed (0/1).
    int                     dir_sign_;  // sign of main direction coord.
    unsigned int            bound_state_; // 0, 1, 2, 3, 4.
public:
    unsigned int            bound_state() const {return bound_state_;}
};

template <class K>
class Straight_2_: public Straight_2_base_ {
public:
                            Straight_2_() ;
                            Straight_2_(typename K::Point_2 const &point) ;
                            Straight_2_(typename K::Line_2 const &line) ;
                            Straight_2_(typename K::Ray_2 const &ray) ;
                            Straight_2_(typename K::Ray_2 const &ray,bool cooriented) ;
                            Straight_2_(typename K::Segment_2 const &seg) ;
    void                    cut_right_off(typename K::Line_2 const & cutter) ;
    int                     collinear_order(typename K::Point_2 const & p1,
                                            typename K::Point_2 const &p2) const ;
    void                    current(typename K::Line_2 &line) const;
    void                    current(typename K::Ray_2 &ray) const;
    void                    current(typename K::Segment_2 &seg) const;
    void                    current(typename K::Point_2 &point) const;
    states                  current_state() const;
    bool                    operator==(const Straight_2_<K>&) const;
    bool                    operator!=(const Straight_2_<K>&other) const
                { return !(*this == other);}
protected:
    typename K::Line_2          support_;   // The supporting line.
    typename K::Point_2         min_;
    typename K::Point_2         max_;
};



inline
Straight_2_base_::
Straight_2_base_()
{
    bound_state_ = LINE_EMPTY;
}

template <class K>
Straight_2_base_::states
Straight_2_<K>::
current_state() const
{
    switch (bound_state_) {
    case BOTH_UNBOUNDED:
        return LINE;
    case MIN_UNBOUNDED:
    case MAX_UNBOUNDED:
        return RAY;
    case NO_UNBOUNDED:
        return (collinear_order(min_, max_) == 0) ? POINT : SEGMENT;
    case LINE_EMPTY:
    default:
        return EMPTY;
    }
}

template <class K>
Straight_2_<K>::
Straight_2_()
{
    bound_state_ = LINE_EMPTY;
}

template <class K>
Straight_2_<K>::
Straight_2_(typename K::Line_2 const &line)
{
    support_ = line;
    typename K::Vector_2 const &dir = support_.direction().to_vector();
    main_dir_ = (CGAL_NTS abs(dir.x()) > CGAL_NTS abs(dir.y()) ) ? 0 : 1;
    dir_sign_ =
        CGAL_NTS sign(support_.direction().to_vector().cartesian(main_dir_));
    bound_state_ = BOTH_UNBOUNDED;
}

template <class K>
Straight_2_<K>::
Straight_2_(typename K::Point_2 const &point)
{
  typedef typename K::Direction_2 Direction_2;
  typedef typename K::Line_2 Line_2;
    support_ = Line_2(point, Direction_2(K::RT(1), K::RT(0)));
    main_dir_ = 0;
    dir_sign_ = 1;
    bound_state_ = NO_UNBOUNDED;
    min_ = point;
    max_ = point;
}

template <class K>
Straight_2_<K>::
Straight_2_(typename K::Ray_2 const &ray)
{
    support_ = ray.supporting_line();
    typename K::Vector_2 const &dir = ray.direction().to_vector();
    main_dir_ = (CGAL_NTS abs(dir.x()) > CGAL_NTS abs(dir.y()) ) ? 0 : 1;
    dir_sign_ =
        CGAL_NTS sign(support_.direction().to_vector().cartesian(main_dir_));
    bound_state_ = MAX_UNBOUNDED;
    min_ = ray.source();
}

template <class K>
Straight_2_<K>::
Straight_2_(typename K::Ray_2 const &ray_,bool cooriented)
{
        typename K::Ray_2 const &ray = cooriented ? ray_ : ray_.opposite();
        support_ = ray.supporting_line();
        /* the supporting line is cooriented with the input ray iff
        cooriented is true */
    typename K::Vector_2 const &dir = ray.direction().to_vector();
    main_dir_ = (CGAL_NTS abs(dir.x()) > CGAL_NTS abs(dir.y()) ) ? 0 : 1;
    dir_sign_ =
        CGAL_NTS sign(support_.direction().to_vector().cartesian(main_dir_));
    if (cooriented)
        {
            bound_state_ = MAX_UNBOUNDED;
                min_ = ray.source();
        }
        else
        {
                bound_state_ = MIN_UNBOUNDED;
                max_ = ray.source();
        }
}

template <class K>
Straight_2_<K>::
Straight_2_(typename K::Segment_2 const &seg)
{
    support_ = seg.supporting_line();
    typename K::Vector_2 const &dir = support_.direction().to_vector();
    main_dir_ = (CGAL_NTS abs(dir.x()) > CGAL_NTS abs(dir.y()) ) ? 0 : 1;
    dir_sign_ =
        CGAL_NTS sign(support_.direction().to_vector().cartesian(main_dir_));
    bound_state_ = NO_UNBOUNDED;
    min_ = seg.source();
    max_ = seg.target();
}


template <class K>
void
Straight_2_<K>::
current(typename K::Line_2 &line) const
{
    CGAL_kernel_assertion(bound_state_ == BOTH_UNBOUNDED);
    line = support_;
}

template <class K>
void
Straight_2_<K>::
current(typename K::Ray_2 &ray) const
{
  typedef typename K::Ray_2 Ray_2;
    CGAL_kernel_assertion(bound_state_ == MIN_UNBOUNDED
                        || bound_state_ == MAX_UNBOUNDED);
    if (bound_state_ == MIN_UNBOUNDED) {
        ray = Ray_2(max_, -support_.direction());
    } else {
        ray = Ray_2(min_, support_.direction());
    }
}

template <class K>
void
Straight_2_<K>::
current(typename K::Segment_2 &seg) const
{
  typedef typename K::Segment_2 Segment_2;
    CGAL_kernel_assertion(bound_state_ == NO_UNBOUNDED
                && collinear_order(min_, max_) != 0);
    seg = Segment_2(min_, max_);
}

template <class K>
void
Straight_2_<K>::
current(typename K::Point_2 &pt) const
{
    CGAL_kernel_assertion(bound_state_ == NO_UNBOUNDED
                && collinear_order(min_, max_) == 0);
    pt = min_;
}


template <class K>
bool Straight_2_<K>::operator==(const Straight_2_<K>& s) const
{
  typedef Straight_2_<K> Straight_2;
  //    enum bound_states {NO_UNBOUNDED=0, MIN_UNBOUNDED=1, MAX_UNBOUNDED=2,
  //                    BOTH_UNBOUNDED = 3, LINE_EMPTY = 4};
  if (bound_state_!=s.bound_state()) return false;
  if (bound_state_==Straight_2::LINE_EMPTY) return true; // empty
  if (support_!=s.support_)
    return false; // on different lines, even if both are points.
  switch (bound_state_)
    {
    case Straight_2::NO_UNBOUNDED:
      return min_==s.min_ && max_==s.max_;
    case Straight_2::MAX_UNBOUNDED:
      return min_==s.min_;
    case Straight_2::MIN_UNBOUNDED:
      return max_==s.max_;
    case Straight_2::BOTH_UNBOUNDED:
      return true;
    }
  return false;
}




template <class K>
int
sign_of_cross(typename K::Direction_2 const &dir1,
	      typename K::Direction_2 const &dir2,
	      const K&)
{
    return static_cast<int>(internal::orientation(dir1.to_vector(),
                                               dir2.to_vector(), K()));
}

template <class K>
void
Straight_2_<K>::
cut_right_off(typename K::Line_2 const & cutter)
// cut off any part of this straight that is to the right of the cutter.
{
    if (bound_state_ == LINE_EMPTY)
        return;
    Line_2_Line_2_pair<K> pair(&support_, &cutter);
    switch (pair.intersection_type()) {
    case Line_2_Line_2_pair<K>::NO_INTERSECTION:
        if (cutter.has_on_negative_side(support_.point()))
            bound_state_ = LINE_EMPTY;
        break;
    case Line_2_Line_2_pair<K>::LINE:
        break;
    case Line_2_Line_2_pair<K>::POINT:
        typename K::Point_2 ispoint = pair.intersection_point();
        bool new_point = false;
        switch (sign_of_cross(support_.direction(), cutter.direction(), K())) {
        case -1: // new minimum.
            if (bound_state_ & MIN_UNBOUNDED) {
                new_point = true;
                bound_state_ ^= MIN_UNBOUNDED;  // exclusive or removes flag.
            } else {
                if (collinear_order(ispoint, min_) == -1)
                    new_point = true;
            }
            if (new_point) {
                if (!(bound_state_ & MAX_UNBOUNDED)
                    && collinear_order(ispoint, max_) == -1)
                    bound_state_ = LINE_EMPTY;
                else
                    min_ = ispoint;
            }
            break;
        case 0: // should not happen
            CGAL_kernel_warning_msg(false, "Internal CGAL error.");
            break;
        case 1: // new maximum
            if (bound_state_ & MAX_UNBOUNDED) {
                new_point = true;
                bound_state_ ^= MAX_UNBOUNDED;  // exclusive or removes flag.
            } else {
                if (collinear_order(ispoint, max_) == 1)
                    new_point = true;
            }
            if (new_point) {
                if (!(bound_state_ & MIN_UNBOUNDED)
                    && collinear_order(ispoint, min_) == 1)
                    bound_state_ = LINE_EMPTY;
                else
                    max_ = ispoint;
            }
            break;
        }
        break;
    }
}

template <class K>
int
Straight_2_<K>::
collinear_order(typename K::Point_2 const &pt1, typename K::Point_2 const & pt2) const
// Compare two points on the support_ line.
// If the second point lies in the direction of the direction vector from
// the first point, the result is 1.
// Other results are -1 (other side) and 0 (coincidence).
{
  typename K::Construct_cartesian_const_iterator_2 construct_cccit;
  typename K::Cartesian_const_iterator_2 ccit1 = construct_cccit(pt1) + main_dir_;
  typename K::Cartesian_const_iterator_2 ccit2 = construct_cccit(pt2) + main_dir_;
    int diffsign;
    diffsign = CGAL_NTS sign(*ccit2 - *ccit1);
    if (diffsign == 0)
        return 0;
    return (diffsign == dir_sign_) ? 1 : -1;
}

} // namespace internal

} //namespace CGAL






#endif
