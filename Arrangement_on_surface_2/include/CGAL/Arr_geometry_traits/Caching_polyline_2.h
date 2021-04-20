#ifndef CGAL_CACHING_POLYLINE_2_H
#define CGAL_CACHING_POLYLINE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/iterator.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/iterator/zip_iterator.hpp>


namespace CGAL {

namespace internal {

template <typename Kernel_, typename Range_>
class Caching_polyline_2_iterator;

template <typename Kernel_, typename Range_>
class Caching_polyline_2
{
public:

  using Kernel = Kernel_;
  using Range = Range_;
  using Self = Caching_polyline_2<Kernel, Range>;

  using Point_2 = typename Kernel::Point_2;
  using Point_ptr = std::shared_ptr<Point_2>;
  using Line_2 = typename Kernel::Line_2;
  using Line_ptr = std::shared_ptr<Line_2>;
  using Line_cache = std::vector<Line_ptr>;
  using Extreme_point = std::pair<Point_ptr, Line_ptr>;

  using Size = std::size_t;
  using size_type = std::size_t;
  using Index = std::uint32_t;

  using iterator = Caching_polyline_2_iterator<Kernel, Range>;
  friend iterator;

  using Subcurve_type_2 = iterator;
  using Subcurve_iterator = Prevent_deref<iterator>;
  using Subcurve_const_iterator = Prevent_deref<iterator>;

protected:

  // Using -2 as null index as -1 can be the index of `m_first`
  inline Index null_idx() const { return Index(-2); }

  const Range* m_range;
  std::shared_ptr<Line_cache> m_line_cache;
  Extreme_point m_first;
  Extreme_point m_last;
  Index m_begin;
  Index m_end;
  bool m_reverse;
  bool m_is_directed_right;

public:

  Caching_polyline_2()
    : m_range(nullptr)
    , m_begin(null_idx())
    , m_end(null_idx())
    , m_reverse(false)
  { }

  // Create a caching polyline from a range
  Caching_polyline_2(const Kernel& kernel, const Range& range, bool force_closure = false)
    : m_range(&range)
    , m_reverse(false)
  {
    m_line_cache = std::make_shared<Line_cache>(range.size(), nullptr);
    m_begin = 0;
    m_end = range.size();
    if (force_closure)
      m_last = Extreme_point(std::make_shared<Point_2>(*range.begin()), nullptr);
    compute_direction(kernel);
  }

  // Create an isolated segment with no range support
  Caching_polyline_2(const Kernel& kernel, const Point_2& first, const Point_2& last)
    : m_first(std::make_shared<Point_2>(first), nullptr)
    , m_last(std::make_shared<Point_2>(last), nullptr)
    , m_begin(null_idx())
    , m_end(null_idx())
    , m_reverse(false)
  {
    compute_direction(kernel);
  }

  // Create a polyline from a subset of another polyline
  Caching_polyline_2(const Kernel& kernel, iterator begin, iterator end)
    : m_range(begin.support().m_range)
    , m_line_cache(begin.support().m_line_cache)
    , m_reverse(false)
  {
    const Self& support = begin.support();
    CGAL_assertion(&support == &end.support());

    if (std::distance(begin, end) < 2) // Polyline with less than 2 points is empty
      return;

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion(support.m_first.first != nullptr);
      m_first = support.m_first;
      m_begin = support.m_begin;
    }
    else
      m_begin = begin.base();

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion(support.m_last.first != nullptr);
      m_last = support.m_last;
      m_end = support.m_end;
    }
    else
      m_end = end.base();

    CGAL_assertion(number_of_subcurves() > 0);

    compute_direction(kernel);
  }

  // Create a polyline from a subset of another polyline and one or two new extreme points
  Caching_polyline_2(const Kernel& kernel, Extreme_point first, iterator begin, iterator end, Extreme_point last)
    : m_range(begin.support().m_range)
    , m_line_cache(begin.support().m_line_cache)
    , m_reverse(false)
  {
    const Self& support = begin.support();
    CGAL_assertion(&support == &end.support());

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion(first.first == nullptr);
      m_first = support.m_first;
      m_begin = support.m_begin;
    }
    else
    {
      m_first = first;
      m_begin = begin.base();
    }

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion(last.first == nullptr);
      m_last = support.m_last;
      m_end = support.m_end;
    }
    else
    {
      m_last = last;
      m_end = end.base();
    }

    compute_direction(kernel);
  }


  // Create the opposite polyline
  Caching_polyline_2 opposite() const
  {
    Caching_polyline_2 out(*this);
    out.m_reverse = !out.m_reverse;
    out.m_is_directed_right = !out.m_is_directed_right;
    return out;
  }

  // Create a new extreme point located on the `index`th segment of
  // the polyline (keep cached line)
  Extreme_point extreme_point(const Point_2& p, std::size_t index) const
  {
    Extreme_point out(std::make_shared<Point_2>(p), nullptr);
    if (index == 0)
      out.second = m_first.second;
    else if (index == number_of_subcurves() + 1)
      out.second = m_last.second;
    else
      out.second = (*m_line_cache)[m_begin + index - 1];
    return out;
  }

  void compute_direction(const Kernel& kernel)
  {
    m_is_directed_right = (kernel.compare_xy_2_object()(*points_begin(), *std::prev(points_end())) == SMALLER);
  }

  bool is_directed_right() const { return m_is_directed_right; }

  Bbox_2 bbox() const
  {
    return bbox_2(points_begin(), points_end());
  }

  iterator points_begin() const { return iterator(this, Tag_true()); }
  iterator points_end() const { return iterator(this, Tag_false()); }

  Subcurve_const_iterator subcurves_begin() const { return Subcurve_const_iterator(points_begin()); }
  Subcurve_const_iterator subcurves_end() const { return Subcurve_const_iterator(std::prev(points_end())); }

  size_type number_of_subcurves() const
  {
    return (m_end - m_begin) + (m_first.first ? 1 : 0) + (m_last.first ? 1 : 0) - 1;
  }

  inline Subcurve_type_2 operator[](const std::size_t i) const
  {
    return std::next(points_begin(), i);
  }

  inline void clear()
  {
    m_begin = null_idx();
    m_end = null_idx();
    m_first = nullptr;
    m_last = nullptr;
  }

  friend std::ostream& operator<<(std::ostream& os, const Self& p)
  {
    os << p.number_of_subcurves();
    for (iterator it = p.points_begin(); it != p.points_end(); ++ it)
      os << " " << *it;
    return os;
  }

private:

  const Point_2& point(const Index& idx) const
  {
    if (idx == m_begin - 1)
    {
      CGAL_assertion(m_first.first != nullptr);
      return *m_first.first;
    }
    // else
    if (idx == m_end)
    {
      CGAL_assertion(m_last.first != nullptr);
      return *m_last.first;
    }
    // else
    CGAL_assertion(m_begin <= idx && idx < m_end);
    return *std::next(m_range->begin(), idx);
  }

  Line_ptr& line_ptr(const Index& idx) const
  {
    if (idx == m_begin - 1)
      return const_cast<Line_ptr&>(m_first.second);
    // else
    if (idx == m_end)
      return const_cast<Line_ptr&>(m_last.second);
    // else
    CGAL_assertion(m_begin <= idx && idx < m_end);
    return const_cast<Line_ptr&>((*m_line_cache)[idx]);
  }

  const Line_2& line(const Kernel& kernel, const Index& index, const Point_2& a, const Point_2& b) const
  {
    Line_ptr& l = line_ptr(index);
    CGAL_BRANCH_PROFILER("Caching polyline line-cache access", br);
    if (!l)
    {
      CGAL_BRANCH_PROFILER_BRANCH(br);
      l = std::make_shared<Line_2>(kernel.construct_line_2_object()(a, b));
    }
    return *l;
  }

};

// This derived class is just there to make sure `Curve_2` and `X_monotone_curve_2` are different types
template <typename Kernel_, typename Range_>
class X_monotone_caching_polyline_2 : public Caching_polyline_2<Kernel_, Range_>
{
public:
  using Base = Caching_polyline_2<Kernel_, Range_>;
  using Self = X_monotone_caching_polyline_2<Kernel_, Range_>;
  using Point_2 = typename Base::Point_2;
  using iterator = typename Base::iterator;
  using Extreme_point = typename Base::Extreme_point;

  X_monotone_caching_polyline_2() : Base() { }
  X_monotone_caching_polyline_2(const Base& base) : Base(base) { }
  X_monotone_caching_polyline_2(const Kernel_& kernel, const Point_2& first, const Point_2& last) : Base(kernel, first, last) { }
  X_monotone_caching_polyline_2(const Kernel_& kernel, iterator begin, iterator end) : Base(kernel, begin, end) { }
  X_monotone_caching_polyline_2(const Kernel_& kernel, Extreme_point first, iterator begin, iterator end, Extreme_point last)
    : Base(kernel, first, begin, end, last) { }
};


template <typename Kernel_, typename Range_>
class Caching_polyline_2_iterator
  : public boost::iterator_facade<Caching_polyline_2_iterator<Kernel_, Range_>,
                                  typename Kernel_::Point_2,
                                  std::random_access_iterator_tag>
{
public:

  using Kernel = Kernel_;
  using Range = Range_;
  using Self = Caching_polyline_2_iterator<Kernel, Range>;
  using Polyline = Caching_polyline_2<Kernel, Range>;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;
  using Index = typename Polyline::Index;

protected:
  const Polyline* m_support;
  Index m_base;

public:

  Caching_polyline_2_iterator() = delete;

  Caching_polyline_2_iterator(const Polyline* support, const Tag_true&) // begin
    : m_support(support)
  {
    if (m_support->m_reverse)
    {
      if (m_support->m_last.first)
        m_base = m_support->m_end;
      else
        m_base = m_support->m_end - 1;
    }
    else
    {
      if (m_support->m_first.first)
        m_base = m_support->m_begin - 1;
      else
        m_base = m_support->m_begin;
    }
  }

  Caching_polyline_2_iterator(const Polyline* support, const Tag_false&) // end
    : m_support(support)
  {
    if (m_support->m_reverse)
    {
      if (m_support->m_first.first)
        m_base = m_support->m_begin - 2;
      else
        m_base = m_support->m_begin - 1;
    }
    else
    {
      if (m_support->m_last.first)
        m_base = m_support->m_end + 1;
      else
        m_base = m_support->m_end;
    }
  }

  // interoperability
  template <typename K, typename I>
  Caching_polyline_2_iterator(const Caching_polyline_2_iterator<K, I>& other)
    : m_support(other.m_support)
    , m_base(other.m_base)
  { }

  const Polyline& support() const { return *m_support; }
  const Index& base() const { return m_base; }

  // The iterator is also used as a wrapper for a segment (using it as
  // source and using the immediate following iterator as target)
  const Point_2& source() const { return const_dereference(); }
  const Point_2& target() const { return std::next(*this).const_dereference(); }

  bool is_directed_right() const
  {
    return m_support->is_directed_right();
  }

  const Point_2& left() const
  {
    if (is_directed_right())
      return source();
    // else
    return target();
  }

  const Point_2& right() const
  {
    if (is_directed_right())
      return target();
    // else
    return source();
  }

  const Line_2& line(const Kernel& kernel) const
  {
    return m_support->line(kernel, m_base, source(), target());
  }

private:
  friend class boost::iterator_core_access;

  // bool reverse = false; -> int r = 1 - 2 * int(reverse) = 1;
  // bool reverse = true; ->  int r = 1 - 2 * int(reverse) = -1;
  void increment()
  {
    move(1 - 2 * int(m_support->m_reverse));
  }

  void decrement()
  {
    move(-(1 - 2 * int(m_support->m_reverse)));
  }

  void advance(std::ptrdiff_t n)
  {
    move((1 - 2 * int(m_support->m_reverse)) * n);
  }

  void move(std::ptrdiff_t n)
  {
    CGAL_assertion(m_support != nullptr);
    m_base += n;
  }

  // interoperability
  template <typename K, typename I>
  std::ptrdiff_t distance_to(const Caching_polyline_2_iterator<K,I>& other) const
  {
    CGAL_assertion(m_support != nullptr);
    CGAL_assertion(other.m_support != nullptr);
    if (m_support->m_reverse)
      return std::ptrdiff_t(m_base) - std::ptrdiff_t(other.m_base);
    // else
    return std::ptrdiff_t(other.m_base) - std::ptrdiff_t(m_base);
  }

  // interoperability
  template <typename K, typename I>
  bool equal(const Caching_polyline_2_iterator<K, I>& other) const
  {
    CGAL_assertion(m_support != nullptr);
    CGAL_assertion(other.m_support != nullptr);
    return m_base == other.m_base;
  }

  Point_2& dereference() const
  {
    return const_cast<Point_2&>(const_dereference());
  }

  const Point_2& const_dereference() const
  {
    CGAL_assertion(m_support != nullptr);
    return m_support->point(m_base);
  }
};

} // namespace internal

} //namespace CGAL

#endif
