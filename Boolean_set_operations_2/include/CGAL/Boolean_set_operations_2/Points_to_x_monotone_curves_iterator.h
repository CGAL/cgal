#ifndef CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H
#define CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H

#include <CGAL/Iterator_range.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/join.hpp>

namespace CGAL
{

template <typename PointIterator, typename ArrPolylineTraits>
class Points_to_x_monotone_curves_iterator
  : public boost::iterator_facade <Points_to_x_monotone_curves_iterator<PointIterator, ArrPolylineTraits>,
                                   typename ArrPolylineTraits::X_monotone_curve_2,
                                   std::input_iterator_tag>
{
  using X_monotone_curve = typename ArrPolylineTraits::X_monotone_curve_2;
  using Point_2 = typename std::iterator_traits<PointIterator>::value_type;

  using PointRange = CGAL::Iterator_range<PointIterator>;
  using Singleton = std::array<Point_2, 1>;
  using Range = boost::range::joined_range<const PointRange, const Singleton>;
  using Range_iterator = typename Range::const_iterator;

  PointRange m_input_range;
  Singleton m_first_element;
  Range m_range;

  std::size_t m_index;
  X_monotone_curve m_curve;
  typename ArrPolylineTraits::Construct_x_monotone_curve_2 m_curve_construct;

public:

  Points_to_x_monotone_curves_iterator ()
    : m_input_range(PointIterator(), PointIterator())
    , m_range (m_input_range, m_first_element)
    , m_index(std::size_t(-1))
    , m_curve_construct(ArrPolylineTraits().construct_x_monotone_curve_2_object())
  { }

  Points_to_x_monotone_curves_iterator (PointIterator begin, PointIterator end)
    : m_input_range(begin, end)
    , m_first_element ({*begin})
    , m_range (m_input_range, m_first_element)
    , m_index(0)
    , m_curve_construct(ArrPolylineTraits().construct_x_monotone_curve_2_object())
  {
    construct_next_x_monotone_curve();
  }


private:

  friend class boost::iterator_core_access;

  Range_iterator iter (std::size_t index) const { return m_range.begin() + index; }
  const Point_2& value (std::size_t index) const { return *(iter(index)); }

  void increment()
  {
    construct_next_x_monotone_curve();
  }

  bool equal (const Points_to_x_monotone_curves_iterator& other) const
  {
    return this->m_index == other.m_index;
  }

  X_monotone_curve& dereference() const
  {
    return const_cast<X_monotone_curve&>(m_curve);
  }

  void construct_next_x_monotone_curve()
  {
    // If end is reached, just put the index out of bound to match end() iterator
    if (m_index == m_range.size())
    {
      m_index = std::size_t(-1);
      return;
    }

    CGAL::Comparison_result comp = compare_xy (value(m_index), value(m_index + 1));
    CGAL::Comparison_result comp_x = compare_x (value(m_index), value(m_index + 1));

    for (std::size_t i = m_index + 2; i < m_range.size(); ++ i)
    {
      const Point_2& pa = value(i-1);
      const Point_2& pb = value(i);
      CGAL::Comparison_result new_comp = compare_xy (pa, pb);
      CGAL::Comparison_result new_comp_x = compare_x (pa, pb);
      if (new_comp != comp || new_comp_x != comp_x)
      {
        m_curve = m_curve_construct (iter(m_index), iter(i));
        m_index = i - 1;
        return;
      }
    }

    m_curve = m_curve_construct (iter(m_index), m_range.end());
    m_index = m_range.size();
  }
};

template <typename ArrPolylineTraits, typename PointRange>
Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator, ArrPolylineTraits>
points_to_x_monotone_curves_begin (const PointRange& points)
{
  return Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator,
                                              ArrPolylineTraits>(points.begin(), points.end());
}
template <typename ArrPolylineTraits, typename PointRange>
Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator, ArrPolylineTraits>
points_to_x_monotone_curves_end (const PointRange& points)
{
  return Points_to_x_monotone_curves_iterator<typename PointRange::const_iterator,
                                              ArrPolylineTraits>();
}

}

#endif // CGAL_BSO_2_POINTS_TO_X_MONOTONE_CURVES_ITERATOR_H
