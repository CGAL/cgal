#ifndef CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H
#define CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H

#include <CGAL/Iterator_range.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_segment_traits_2.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/join.hpp>

namespace CGAL
{

template <typename Kernel, typename PointIterator>
void construct_next_x_monotone_curve_in_range
(PointIterator& current, PointIterator end,
 typename Arr_segment_traits_2<Kernel>
 ::X_monotone_curve_2& curve,
 const Arr_segment_traits_2<Kernel>& traits)
{
  using Arr_traits = Arr_segment_traits_2<Kernel>;
  using Point_2 = typename std::iterator_traits<PointIterator>::value_type;
  typename Arr_traits::Construct_x_monotone_curve_2 curve_construct
    = traits.construct_x_monotone_curve_2_object();

  PointIterator ita = current;
  PointIterator itb = current;
  itb ++;

  if (itb == end)
  {
    current = PointIterator();
    return;
  }

  curve = curve_construct (*ita, *itb);
  current = itb;
}

template <typename Kernel, typename PointIterator>
void construct_next_x_monotone_curve_in_range
(PointIterator& current, PointIterator end,
 typename Arr_polyline_traits_2<Arr_segment_traits_2<Kernel> >
 ::X_monotone_curve_2& curve,
 const Arr_polyline_traits_2<Arr_segment_traits_2<Kernel> >& traits)
{
  using Arr_traits = Arr_polyline_traits_2<Arr_segment_traits_2<Kernel> >;
  using Point_2 = typename std::iterator_traits<PointIterator>::value_type;
  typename Arr_traits::Construct_x_monotone_curve_2 curve_construct
    = traits.construct_x_monotone_curve_2_object();

  // If end is reached, just put the index out of bound to match end() iterator
  if (current == end)
  {
    current = PointIterator();
    return;
  }

  PointIterator ita = current;
  PointIterator itb = current;
  itb ++;

  CGAL::Comparison_result comp = compare_xy (*ita, *itb);
  CGAL::Comparison_result comp_x = compare_x (*ita, *itb);

  for (; itb != end; ++ ita, ++ itb)
  {
    const Point_2& pa = *ita;
    const Point_2& pb = *itb;
    CGAL::Comparison_result new_comp = compare_xy (pa, pb);
    CGAL::Comparison_result new_comp_x = compare_x (pa, pb);
    if (new_comp != comp || new_comp_x != comp_x)
    {
      curve = curve_construct (current, itb);
      current = ita;
      return;
    }
  }

  curve = curve_construct (current, end);
  current = end;
}

template <typename PointIterator, typename ArrTraits>
class Points_to_x_monotone_curves_iterator
  : public boost::iterator_facade <Points_to_x_monotone_curves_iterator<PointIterator, ArrTraits>,
                                   typename ArrTraits::X_monotone_curve_2,
                                   std::input_iterator_tag>
{
  using X_monotone_curve = typename ArrTraits::X_monotone_curve_2;
  using Point_2 = typename std::iterator_traits<PointIterator>::value_type;

  using PointRange = CGAL::Iterator_range<PointIterator>;
  using Singleton = std::array<Point_2, 1>;
  using Range = boost::range::joined_range<const PointRange, const Singleton>;
  using Range_iterator = typename Range::const_iterator;

  PointRange m_input_range;
  Singleton m_first_element;
  Range m_range;

  Range_iterator m_current;
  X_monotone_curve m_curve;

public:

  Points_to_x_monotone_curves_iterator ()
    : m_input_range(PointIterator(), PointIterator())
    , m_range (m_input_range, m_first_element)
    , m_current()
  { }

  Points_to_x_monotone_curves_iterator (PointIterator begin, PointIterator end)
    : m_input_range(begin, end)
    , m_first_element ({*begin})
    , m_range (m_input_range, m_first_element)
    , m_current(m_range.begin())
  {
    construct_next_x_monotone_curve_in_range
      (m_current, m_range.end(), m_curve, ArrTraits());
  }


private:

  friend class boost::iterator_core_access;

  void increment()
  {
    construct_next_x_monotone_curve_in_range
      (m_current, m_range.end(), m_curve, ArrTraits());
  }

  bool equal (const Points_to_x_monotone_curves_iterator& other) const
  {
    return this->m_current == other.m_current;
  }

  X_monotone_curve& dereference() const
  {
    return const_cast<X_monotone_curve&>(m_curve);
  }
};

template <typename ArrTraits, typename PointRange>
Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator, ArrTraits>
points_to_x_monotone_curves_begin (const PointRange& points)
{
  return Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator,
                                              ArrTraits>(points.begin(), points.end());
}
template <typename ArrTraits, typename PointRange>
Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator, ArrTraits>
points_to_x_monotone_curves_end (const PointRange& points)
{
  return Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator,
                                              ArrTraits>();
}

}

#endif // CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H
