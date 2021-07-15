// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)


// custom point container
#ifndef CGAL_POINT_CONTAINER_H
#define CGAL_POINT_CONTAINER_H

#include <CGAL/license/Spatial_searching.h>


#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/internal/Get_dimension_tag.h>

#include <boost/optional.hpp>

namespace CGAL {

template <class Traits>
class Point_container {

private:
  typedef typename Traits::Point_d Point_d;
  typedef std::vector<const Point_d*> Point_vector;

public:
  typedef typename Traits::FT FT;

  typedef typename Point_vector::iterator iterator;
  typedef typename Point_vector::const_iterator const_iterator;
  typedef typename internal::Get_dimension_tag<Traits>::Dimension D;
private:
  Traits traits;
  // the iterator range of the Point_container
  boost::optional<iterator> m_b ;
  boost::optional<iterator> m_e ;

  int built_coord;    // a coordinate for which the pointer list is built
  Kd_tree_rectangle<FT,D> bbox;       // bounding box, i.e. rectangle of node
  Kd_tree_rectangle<FT,D> tbox;       // tight bounding box,
  // i.e. minimal enclosing bounding
  // box of points

  using Construct_cartesian_const_iterator_d = typename Traits::Construct_cartesian_const_iterator_d;
  using FTP = typename Construct_cartesian_const_iterator_d::result_type;

  std::vector< std::vector< std::pair<FTP, const Point_d*> > > m_references;
  std::vector<FTP> m_temporary;

  long m_start = -1;
  long m_end = -1;

  std::size_t m_dim = -1;
  std::size_t m_depth = -1;

public:

  void initialize_references(
    const iterator begin, const iterator end,
    const Construct_cartesian_const_iterator_d& construct_it,
    std::vector< std::pair<FTP, const Point_d*> >& reference) const {

    // std::cout << std::endl;
    std::size_t count = 0;
    for (auto it = begin; it != end; ++it) {
      const auto bit = construct_it(**it);
      reference[count] = std::make_pair(bit, *it);
      // std::cout << *reference[count] << std::endl;
      ++count;
    }
    // CGAL_assertion_msg(false, "TODO: FINISH INITIALIZE REFERENCE!");
  }

  FT compare_keys(
    const FTP a, const FTP b, const std::size_t p, const int dim) const {

    FT diff = FT(0);
    for (std::size_t i = 0; i < dim; ++i) {
      std::size_t r = i + p;
      r = (r < dim) ? r : r - dim;
      diff = a[r] - b[r];
      if (diff != FT(0)) break;
    }
    // CGAL_assertion_msg(false, "TODO: FINISH SUPER KEY COMPARE!");
    return diff;
  }

  void apply_merge_sort(
    std::vector< std::pair<FTP, const Point_d*> >& reference,
    std::vector<FTP>& temporary,
    const std::size_t low, const std::size_t high,
    const std::size_t p, const int dim) const {

    std::size_t i, j, k;
    if (high > low) {

      const std::size_t median = low + ( (high - low) >> 1 );
      apply_merge_sort(reference, temporary, low, median, p, dim);
      apply_merge_sort(reference, temporary, median + 1, high, p, dim);

      for (i = median + 1; i > low; --i) {
        temporary[i - 1] = reference[i - 1].first;
      }

      for (j = median; j < high; ++j) {
        temporary[median + (high - j)] = reference[j + 1].first;
      }

      for (k = low; k <= high; ++k) {
        reference[k].first = (
          compare_keys(temporary[i], temporary[j], p, dim) < FT(0) ) ?
          temporary[i++] :
          temporary[j--] ;
      }
    }

    // CGAL_assertion_msg(false, "TODO: FINISH MERGE SORT!");
  }

  std::size_t remove_duplicates(
    std::vector< std::pair<FTP, const Point_d*> >& reference,
    const std::size_t i, const int dim) const {

    std::size_t end = 0;
    for (std::size_t j = 1; j < reference.size(); ++j) {
      const FT compare = compare_keys(reference[j].first, reference[j - 1].first, i, dim);

      if (compare < FT(0)) {
        CGAL_assertion_msg(false, "ERROR: NEGATIVE SUPER KEY COMPARE RESULT!");
      } else if (compare > FT(0)) {
        reference[++end].first = reference[j].first;
      }
    }

    // CGAL_assertion_msg(false, "TODO: FINISH REMOVE DUPLICATES!");
    return end;
  }

  void print_references() const {

    for (const auto& reference : m_references) {
      std::cout << std::endl;
      for (const auto& item : reference) {
        std::cout << *(item.first) << std::endl;
      }
    }
  }

  struct Balanced_cmp {

    const int m_dim;
    const long m_median;
    const std::size_t m_axis;
    std::vector< std::pair<FTP, const Point_d*> >& m_reference;

    FT compare_keys( // TODO: Remove this duplicated method!
      const FTP a, const FTP b, const std::size_t p, const int dim) const {

      FT diff = FT(0);
      for (std::size_t i = 0; i < dim; ++i) {
        std::size_t r = i + p;
        r = (r < dim) ? r : r - dim;
        diff = a[r] - b[r];
        if (diff != FT(0)) break;
      }
      // CGAL_assertion_msg(false, "TODO: FINISH SUPER KEY COMPARE!");
      return diff;
    }

    Balanced_cmp(
      std::vector< std::pair<FTP, const Point_d*> >& reference,
      const std::size_t axis, const long median, const int dim) :
    m_dim(dim), m_median(median), m_axis(axis), m_reference(reference)
    { }

    bool operator()(const Point_d* pt) const {

      FTP a; bool duplicates_found = true;
      for (const auto& item : m_reference) {
        if (item.second == pt) {
          a = item.first;
          duplicates_found = false;
          break;
        }
      }
      if (duplicates_found) {
        CGAL_assertion_msg(false, "ERROR: DUPLICATES ARE FOUND!");
        return false;
      }

      const FT compare = compare_keys(
        a, m_reference[m_median].first, m_axis, m_dim);
      return compare < FT(0);
    }
  };

public:

  void set_data(
    const std::vector< std::vector< std::pair<FTP, const Point_d*> > >& references,
    const std::vector<FTP>& temporary) {

    m_references = references;
    m_temporary = temporary;
  }

  void set_data(
    const std::size_t dim, const std::size_t depth,
    const long start, const long end) {

    m_dim = dim;
    m_depth = depth;
    m_start = start;
    m_end = end;
  }

  void set_initial_depth(const std::size_t depth) {
    m_depth = depth;
  }

  void balanced_split(Point_container<Traits>& c) {

    // CGAL_assertion(dimension() == c.dimension());
    // CGAL_assertion(is_valid());

    // NEW CODE!
    CGAL_assertion(m_dim   != static_cast<std::size_t>(-1));
    CGAL_assertion(m_depth != static_cast<std::size_t>(-1));

    CGAL_assertion(m_start != -1);
    CGAL_assertion(m_end   != -1);

    const long median = m_start + ((m_end - m_start) / 2);
    const std::size_t axis = m_depth % m_dim;

    // std::cout << "data: "   << size()  << std::endl;
    // std::cout << "start: "  << m_start << std::endl;
    // std::cout << "end: "    << m_end   << std::endl;
    // std::cout << "median: " << median  << std::endl;

    CGAL_assertion(m_end > m_start + 2);
    for (std::size_t i = m_start; i <= m_end; ++i) {
      m_temporary[i] = m_references[0][i].first;
    }

    std::size_t lower, upper;
    for (std::size_t i = 1; i < m_dim; ++i) {
      lower = m_start-1;
      upper = median;
      for (std::size_t j = m_start; j <= m_end; ++j) {
        const FT compare = compare_keys(
          m_references[i][j].first, m_references[0][median].first, axis, m_dim);
        if (compare < FT(0)) {
          m_references[i-1][++lower].first = m_references[i][j].first;
        } else if (compare > FT(0)) {
          m_references[i-1][++upper].first = m_references[i][j].first;
        }
      }
    }

    for (std::size_t i = m_start; i <= m_end; ++i) {
      m_references[m_dim-1][i].first = m_temporary[i];
    }

    // std::cout << "DEPTH: " << m_depth << std::endl;
    // print_references();

    ///////////////////////

    // OLD CODE!

    // c.bbox = bbox;
    // const int split_coord = sep.cutting_dimension();
    // FT cutting_value = sep.cutting_value();

    // built_coord = split_coord;
    // c.built_coord = split_coord;

    // auto construct_it = traits.construct_cartesian_const_iterator_d_object();
    // Cmp<Traits> cmp(split_coord, cutting_value, construct_it);
    // iterator it = std::partition(begin(), end(), cmp);

    ///////////////////////

    // NEW VERSION!
    Balanced_cmp cmp(m_references[m_dim-1], axis, median, m_dim);
    iterator it = std::partition(begin(), end(), cmp);

    c.set_range(begin(), it);
    set_range(it, end());

    ///////////////////////

    // OLD VERSION!

    // std::vector<const Point_d*> data1, data2;
    // for (std::size_t i = m_start; i <= m_end; ++i) {
    //   if (i < median) {
    //     data1.push_back(m_references[m_dim-1][i].second);
    //   } else if (i > median) {
    //     data2.push_back(m_references[m_dim-1][i].second);
    //   }
    // }

    // c.set_range(data1.begin(), data1.end());
    // set_range(data2.begin(), data2.end());

    // std::cout << "data 1: " << std::endl;
    // for (const auto& d1 : data1) {
    //   std::cout << *d1 << std::endl;
    // }

    // std::cout << "data 2: " << std::endl;
    // for (const auto& d2 : data2) {
    //   std::cout << *d2 << std::endl;
    // }

    ///////////////////////

    // std::cout << "pre size c0: " << c.size() << std::endl;
    // std::cout << "pre size c1: " << size()   << std::endl;

    c.set_data(m_references, m_temporary);
    c.set_data(m_dim, m_depth + 1, m_start, lower);
    set_data(m_dim, m_depth + 1, median + 1, upper);

    // Adjust boxes.
    // bbox.set_lower_bound(split_coord, cutting_value);
    // tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(begin(), end(), construct_it);
    // c.bbox.set_upper_bound(split_coord, cutting_value);
    // c.tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(c.begin(), c.end(), construct_it);
    // CGAL_assertion(is_valid());
    // CGAL_assertion(c.is_valid());

    // CGAL_assertion_msg(false, "TODO: FINISH BALANCED SPLIT!");
  }

  inline const Kd_tree_rectangle<FT,D>&
  bounding_box() const
  {
    return bbox;
  }

  inline const Kd_tree_rectangle<FT,D>&
  tight_bounding_box() const
  {
    return tbox;
  }

  inline int
  dimension() const
  {
    return bbox.dimension();
  }

  inline int
  built_coordinate() const
  {
    return built_coord;
  }

  // coordinate of the maximal span
  inline int
  max_span_coord() const
  {
    return bbox.max_span_coord();
  }

  // coordinate of the maximal tight span
  inline int
  max_tight_span_coord() const
  {
    return tbox.max_span_coord();
  }

  inline FT
  max_span_lower() const
  {
    return bbox.min_coord(max_span_coord());
  }

  inline FT
  max_tight_span_lower() const
  {
    return tbox.min_coord(max_tight_span_coord());
  }

  inline FT
  max_span_upper() const
  {
    return bbox.max_coord(max_span_coord());
  }

  inline FT
  max_tight_span_upper() const
  {
    return tbox.max_coord(max_tight_span_coord());
  }

  inline FT
  max_spread() const
  {
    return  max_span_upper() -  max_span_lower();
  }

  inline FT
  max_tight_spread() const
  {
    return  max_tight_span_upper() -  max_tight_span_lower();
  }


  int
  max_tight_span_coord_balanced(FT Aspect_ratio) const
  {
    int cut_dim(-1);
    FT max_spread_points(FT(-1));
    FT max_length = max_spread();  // length of longest side of box
    int dim = dimension();
    for (int d=0; d<dim; d++) {
      FT length=bbox.max_coord(d)-bbox.min_coord(d);

      if (FT(2)*max_length/length <= Aspect_ratio) {
        FT spread=tbox.max_coord(d)-tbox.min_coord(d);

        if (spread > max_spread_points) {
          max_spread_points = spread;
          cut_dim = d;
        }
      }
    }
    // CGAL_assertion(cut_dim >= 0);
    return cut_dim;
  }

  FT
  max_span_upper_without_dim(int d) const
  {
    FT max_span(FT(0));
    int dim=dimension();
    for (int i=0; i<dim; i++) {
      FT span = bbox.max_coord(i)-bbox.min_coord(i);
      if (d != i && span > max_span) max_span=span;
    }
    return max_span;
  }

  FT
  balanced_fair(int d, FT Aspect_ratio)
  {
    FT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
    FT low_cut = bbox.min_coord(d) + small_piece; // lowest legal cut;
    FT high_cut = bbox.max_coord(d) - small_piece; //highest legal cut;
    // CGAL_assertion (high_cut >= low_cut);
    FT split_value = median(d);
    if (split_value < low_cut) split_value = low_cut;
    if (split_value > high_cut) split_value = high_cut;
    return split_value;
  }

  FT
  balanced_sliding_fair(int d, FT Aspect_ratio)
  {
    FT small_piece = max_span_upper_without_dim(d) / Aspect_ratio;
    FT low_cut = bbox.min_coord(d) + small_piece; // lowest legal cut;
    FT high_cut = bbox.max_coord(d) - small_piece; //highest legal cut;
    // CGAL_assertion (high_cut >= low_cut);
    FT split_value = median(d);
    FT max_span_lower = tbox.min_coord(d);
    FT max_span_upper = tbox.max_coord(d);
    if (split_value < low_cut) split_value= max_span_lower;
    if (split_value > high_cut) split_value = max_span_upper;
    return split_value;
  }

  //  points
  inline std::size_t
  size() const
  {
    return *m_e - *m_b;
  }

  inline const_iterator
  begin() const {
    return *m_b;
  }

  inline const_iterator
  end() const
  {
    return *m_e;
  }

  inline iterator
  begin()
  {
    return *m_b;
  }

  inline iterator
  end()
  {
    return *m_e;
  }

  inline bool
  empty() const
  {
    return !m_b || !m_e || (*m_b == *m_e ) ;
  }

  // building the container from a sequence of Point_d*
  Point_container(const int d, iterator begin, iterator end, const Traits& traits_) :
    traits(traits_), m_b(begin), m_e(end), bbox(d, begin, end, traits.construct_cartesian_const_iterator_d_object()), tbox(bbox)
  {
    built_coord = max_span_coord();
  }

  void
  set_range(iterator begin, iterator end)
  {
    m_b = begin;
    m_e = end;
  }


  // building an empty container
  Point_container(const int d,const Traits& traits_) :
    traits(traits_),bbox(d), tbox(d)
  {}

  template <class Traits2>
  struct Cmp {
    typedef typename Traits2::FT FT;
    typedef typename Traits2::Point_d Point_d;
    typedef std::vector<const Point_d*> Point_vector;

    int split_coord;
    FT value;
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;

    Cmp(int s, FT c,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it)
      : split_coord(s), value(c), construct_it(cst_it)
    {}

    bool
    operator()(const Point_d* pt) const
    {
      typename Traits2::Cartesian_const_iterator_d ptit;
      ptit = construct_it(*pt);
      return  *(ptit+split_coord) < value;
    }
  };


  template <class Traits2>
  struct Between {
    typedef typename Traits2::FT FT;
    typedef typename Traits2::Point_d Point_d;
    typedef std::vector<const Point_d*> Point_vector;

    int split_coord;
    FT low, high;
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;

    Between(int s, FT l, FT h,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it)
      : split_coord(s), low(l), high(h), construct_it(cst_it)
    {}

    bool
    operator()(const Point_d* pt) const
    {
      typename Traits2::Cartesian_const_iterator_d ptit;
      ptit = construct_it(*pt);
      if(! ( *(ptit+split_coord) <= high ) ){
        //        std::cerr << "Point " << *pt << " exceeds " << high << " in dimension " << split_coord << std::endl;
        return false;
      }
      if(! ( *(ptit+split_coord) >= low ) ){
        //std::cerr << "Point " << *pt << " below " << low << " in dimension " << split_coord << std::endl;
        return false;
      }
      return true;
    }
  };


  void recompute_tight_bounding_box()
  {
    tbox.template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(begin(), end(),traits.construct_cartesian_const_iterator_d_object());
  }


  bool
  is_valid() const
  {
    if(empty()) return true;
    bool b = true;
    for (int i = 0; i < dimension(); i++){
      CGAL_assertion( b = (b && (bbox.min_coord(i) <= tbox.min_coord(i))));
      CGAL_assertion( b = (b && (bbox.max_coord(i) >= tbox.max_coord(i))));

      typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
      Between<Traits> between(i,tbox.min_coord(i), tbox.max_coord(i), construct_it);
      for(const_iterator it = begin(); it != end(); it++){
        b = (b && between(*it));
      }
    }
    return b;
  }


  // note that splitting is restricted to the built coordinate
  template <class Separator>
  void split(Point_container<Traits>& c, Separator& sep,
             bool sliding=false)
  {
    CGAL_assertion(dimension()==c.dimension());
    CGAL_assertion(is_valid());
    c.bbox=bbox;

    const int split_coord = sep.cutting_dimension();
    FT cutting_value = sep.cutting_value();

    built_coord=split_coord;
    c.built_coord=split_coord;


    typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();

    Cmp<Traits> cmp(split_coord, cutting_value,construct_it);
    iterator it = std::partition(begin(), end(), cmp);
    // now [begin,it) are lower and [it,end) are upper
    if (sliding) { // avoid empty lists

      if (it == begin()) {
        iterator minelt = std::min_element(begin(),end(),comp_coord_val<Traits,int>(split_coord,construct_it));
        if(minelt != it){
          std::iter_swap(minelt,it);
        }
        cutting_value = *(construct_it(**it)+split_coord);
        sep.set_cutting_value(cutting_value);
        it++;
      }
      if (it == end()) {
        iterator maxelt = std::max_element(begin(),end(),comp_coord_val<Traits,int>(split_coord,construct_it));
        it--;
        if(maxelt != it){
          std::iter_swap(maxelt,it);
        }
        cutting_value = *(construct_it(**it)+split_coord);
        sep.set_cutting_value(cutting_value);
      }
    }

    c.set_range(begin(), it);
    set_range(it, end());
    // adjusting boxes
    bbox.set_lower_bound(split_coord, cutting_value);
    tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(begin(),end(),construct_it);
    c.bbox.set_upper_bound(split_coord, cutting_value);
    c.tbox. template update_from_point_pointers<typename Traits::Construct_cartesian_const_iterator_d>(c.begin(),c.end(),construct_it);
    CGAL_assertion(is_valid());
    CGAL_assertion(c.is_valid());
  }



  template <class Traits2, class Value>
  struct comp_coord_val {

  private:
    Value coord;
    const typename Traits2::Construct_cartesian_const_iterator_d& construct_it;

    typedef typename Traits2::Point_d Point_d;
  public:
    comp_coord_val (const Value& coordinate,const typename Traits2::Construct_cartesian_const_iterator_d& cst_it)
      : coord(coordinate), construct_it(cst_it)
    {}

    bool
    operator()(const Point_d *a, const Point_d *b) const
    {
      typename Traits2::Cartesian_const_iterator_d ait = construct_it(*a),
        bit = construct_it(*b);
      return *(ait+coord) < *(bit+coord);
    }
  };


  FT
  median(const int split_coord)
  {
    typename Traits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
    iterator mid = begin() + (end() - begin())/2;
    std::nth_element(begin(), mid, end(),comp_coord_val<Traits,int>(split_coord,construct_it));

    typename Traits::Cartesian_const_iterator_d mpit = construct_it((*(*mid)));
    FT val1 = *(mpit+split_coord);
    mid++;
    mpit = construct_it((*(*mid)));
    FT val2 = *(mpit+split_coord);
    return (val1+val2)/FT(2);
  }



private:
  explicit Point_container()
  {} // disable default constructor

};

  template <class Point>
  std::ostream&
  operator<< (std::ostream& s, Point_container<Point>& c)
  {
    s << "Points container of size " << c.size() << "\n cell:";
    s << c.bounding_box();
    s << "\n minimal box enclosing points:"; s << c.tight_bounding_box();
    return s;
  }

} // namespace CGAL

#endif // CGAL_POINT_CONTAINER_H
