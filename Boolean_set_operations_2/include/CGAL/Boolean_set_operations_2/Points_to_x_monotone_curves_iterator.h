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

  Range_iterator m_current;
  X_monotone_curve m_curve;
  typename ArrPolylineTraits::Construct_x_monotone_curve_2 m_curve_construct;

public:

  Points_to_x_monotone_curves_iterator ()
    : m_input_range(PointIterator(), PointIterator())
    , m_range (m_input_range, m_first_element)
    , m_current()
    , m_curve_construct(ArrPolylineTraits().construct_x_monotone_curve_2_object())
  { }

  Points_to_x_monotone_curves_iterator (PointIterator begin, PointIterator end)
    : m_input_range(begin, end)
    , m_first_element ({*begin})
    , m_range (m_input_range, m_first_element)
    , m_current(m_range.begin())
    , m_curve_construct(ArrPolylineTraits().construct_x_monotone_curve_2_object())
  {
    construct_next_x_monotone_curve();
  }


private:

  friend class boost::iterator_core_access;

  void increment()
  {
    construct_next_x_monotone_curve();
  }

  bool equal (const Points_to_x_monotone_curves_iterator& other) const
  {
    return this->m_current == other.m_current;
  }

  X_monotone_curve& dereference() const
  {
    return const_cast<X_monotone_curve&>(m_curve);
  }

  void construct_next_x_monotone_curve()
  {
    // If end is reached, just put the index out of bound to match end() iterator
    if (m_current == m_range.end())
    {
      m_current = Range_iterator();
      return;
    }

    Range_iterator ita = m_current;
    Range_iterator itb = m_current;
    itb ++;

    CGAL::Comparison_result comp = compare_xy (*ita, *itb);
    CGAL::Comparison_result comp_x = compare_x (*ita, *itb);

    for (; itb != m_range.end(); ++ ita, ++ itb)
    {
      const Point_2& pa = *ita;
      const Point_2& pb = *itb;
      CGAL::Comparison_result new_comp = compare_xy (pa, pb);
      CGAL::Comparison_result new_comp_x = compare_x (pa, pb);
      if (new_comp != comp || new_comp_x != comp_x)
      {
        m_curve = m_curve_construct (m_current, itb);
        m_current = ita;
        return;
      }
    }

    m_curve = m_curve_construct (m_current, m_range.end());
    m_current = m_range.end();
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
