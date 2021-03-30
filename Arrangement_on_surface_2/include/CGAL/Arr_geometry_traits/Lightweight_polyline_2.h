#ifndef CGAL_LIGHTWEIGHT_POLYLINE_2_H
#define CGAL_LIGHTWEIGHT_POLYLINE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/iterator.h>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/iterator/zip_iterator.hpp>


namespace CGAL {

namespace internal {

template <typename Kernel_, typename Iterator>
class Lightweight_polyline_2_iterator;

template <typename Kernel_, typename Iterator_>
class Lightweight_polyline_2
{
public:

  using Kernel = Kernel_;
  using Base_iterator = Iterator_;
  using Self = Lightweight_polyline_2<Kernel, Base_iterator>;

  using Point_2 = typename Kernel::Point_2;
  using Point_ptr = std::shared_ptr<Point_2>;
  using Line_2 = typename Kernel::Line_2;
  using Line_ptr = std::shared_ptr<Line_2>;
  using Line_cache = std::vector<Line_ptr>;
  using Line_cache_iterator = typename Line_cache::const_iterator;
  using Extreme_point = std::pair<Point_ptr, Line_ptr>;

  using Size = std::size_t;
  using size_type = std::size_t;

  using Zip_iterator = boost::zip_iterator<std::pair<Base_iterator, Line_cache_iterator> >;
  using Internal_point = typename std::iterator_traits<Zip_iterator>::value_type;

  using iterator = Lightweight_polyline_2_iterator<Kernel_, Base_iterator>;
  friend iterator;

  using Subcurve_type_2 = iterator;
  using Subcurve_iterator = Prevent_deref<iterator>;
  using Subcurve_const_iterator = Prevent_deref<iterator>;

protected:

  Zip_iterator m_begin;
  Zip_iterator m_end;
  Extreme_point m_first;
  Extreme_point m_last;
  std::shared_ptr<Line_cache> m_line_cache;
  bool m_reverse;
  bool m_is_directed_right;

public:

  Lightweight_polyline_2() : m_reverse(false) { }

  Lightweight_polyline_2(const Point_2& first, const Point_2& last)
    : m_first(init(first))
    , m_last(init(last))
    , m_reverse(false)
  {
    compute_direction();
  }

  Lightweight_polyline_2 (Base_iterator begin, Base_iterator end, bool force_closure = false)
    : m_reverse(false)
  {
    m_line_cache = std::make_shared<Line_cache>(std::distance(begin, end), nullptr);
    m_begin = Zip_iterator(std::make_pair(begin, m_line_cache->begin()));
    m_end = Zip_iterator(std::make_pair(end, m_line_cache->end()));
    if (force_closure)
      m_last = init(*m_begin);
    compute_direction();
  }

  Lightweight_polyline_2 (iterator begin, iterator end)
    : m_line_cache (begin.support().m_line_cache), m_reverse(false)
  {
    const Self& support = begin.support();
    CGAL_assertion (&support == &end.support());

    if (std::distance(begin, end) < 2) // Polyline with less than 2 points is empty
      return;

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion (is_init(support.m_first));
      m_first = support.m_first;
      m_begin = support.m_begin;
    }
    else
      m_begin = begin.base();

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion (is_init(support.m_last));
      m_last = support.m_last;
      m_end = support.m_end;
    }
    else
      m_end = end.base();

    CGAL_assertion (number_of_subcurves() > 0);
    // This constructor should only create x-monotone polylines
    CGAL_assertion (is_x_monotone());

    compute_direction();
  }

  Lightweight_polyline_2 (Extreme_point first, iterator begin, iterator end, Extreme_point last)
    : m_line_cache (begin.support().m_line_cache), m_reverse(false)
  {
    const Self& support = begin.support();
    CGAL_assertion (&support == &end.support());

    if (begin.base() == support.m_begin - 1)
    {
      CGAL_assertion (!is_init(first));
      m_first = init(support.m_first);
      m_begin = support.m_begin;
    }
    else
    {
      m_first = first;
      m_begin = begin.base();
    }

    if (end.base() == support.m_end + 1)
    {
      CGAL_assertion (!is_init(last));
      m_last = init(support.m_last);
      m_end = support.m_end;
    }
    else
    {
      m_last = last;
      m_end = end.base();
    }

    // This constructor should only create x-monotone polylines
    CGAL_assertion (is_x_monotone());

    compute_direction();
  }

  Lightweight_polyline_2 opposite() const
  {
    Lightweight_polyline_2 out (*this);
    out.m_reverse = !out.m_reverse;
    out.m_is_directed_right = !out.m_is_directed_right;
    return out;
  }

  Extreme_point init (const Point_2& p, std::size_t index) const
  {
    Extreme_point out = this->init(p);
    if (index == 0)
      out.second = m_first.second;
    else if (index == number_of_subcurves() + 1)
      out.second = m_last.second;
    else
      out.second = (m_begin + index - 1)->second;
    return out;
  }

  void compute_direction()
  {
#ifdef CGAL_PROFILE
    std::size_t nb_curves = number_of_subcurves();
    if (nb_curves == 1)
    {
      if (is_init(m_first) && is_init(m_last))
      {
        CGAL_PROFILER ("Polyline with 1 isolated subcurve");
      }
      else
      {
        CGAL_PROFILER ("Polyline with 1 subcurve");
      }
    }
    else
    {
      CGAL_PROFILER ("Polyline with >=2 subcurves");
    }
#endif

    // In all this class we use boost::prior instead of std::prev for compatibility with zip iterators
    m_is_directed_right = (Kernel().compare_xy_2_object()(*points_begin(), *boost::prior(points_end())) == SMALLER);
  }

  bool is_directed_right() const { return m_is_directed_right; }

  bool is_x_monotone() const
  {
    auto compare_x_2 = Kernel().compare_x_2_object();
    iterator b = points_begin(), e = points_end() - 1;
    Comparison_result comp = Kernel().compare_x_2_object()(*b, *(b+1));
    for (iterator it = b + 1; it < e; ++ it)
      if (comp != Kernel().compare_x_2_object()(*it, *(it+1)))
        return false;
    return true;
  }

  Bbox_2 bbox() const
  {
    return bbox_2 (points_begin(), points_end());
  }

  iterator points_begin() const { return iterator (this, Tag_true()); }
  iterator points_end() const { return iterator (this, Tag_false()); }

  Subcurve_const_iterator subcurves_begin() const { return Subcurve_const_iterator(points_begin()); }
  Subcurve_const_iterator subcurves_end() const { return Subcurve_const_iterator(boost::prior(points_end())); }

  size_type number_of_subcurves() const
  {
    return std::distance (m_begin, m_end) + (is_init(m_first) ? 1 : 0) + (is_init(m_last) ? 1 : 0) - 1;
  }

  inline Subcurve_type_2 operator[](const std::size_t i) const
  {
    return std::next(points_begin(), i);
  }

  inline void clear()
  {
    m_begin = Base_iterator();
    m_end = Base_iterator();
    m_first = nullptr;
    m_last = nullptr;
  }

  friend std::ostream& operator<< (std::ostream& os, const Self& p)
  {
    os << p.number_of_subcurves();
    for (iterator it = p.points_begin(); it != p.points_end(); ++ it)
      os << " " << *it;
    return os;
  }

  Extreme_point uninit() const { return Extreme_point(nullptr, nullptr); }

private:

  bool is_init (const Extreme_point& p) const { return (p.first != nullptr); }
  Extreme_point init (const Point_2& p) const { return Extreme_point(std::make_shared<Point_2>(p), nullptr); }
  Extreme_point init (const Internal_point& p) const { return Extreme_point(std::make_shared<Point_2>(p.first), p.second); }
  Extreme_point init (const Extreme_point& p) const { return p; }

  const Point_2& point (const Internal_point& p) const { return p.first; }
  const Point_2& point (const Extreme_point& p) const { return *p.first; }

  const Line_2& line (std::shared_ptr<Line_2>& line, const Point_2& a, const Point_2& b) const
  {
    CGAL_BRANCH_PROFILER("Cache acces", br);
    if (!line)
    {
      CGAL_BRANCH_PROFILER_BRANCH(br);
      line = std::make_shared<Line_2>(a, b);
    }
    return *line;
  }

};

template <typename Kernel_, typename Iterator>
class Lightweight_polyline_2_iterator
  : public boost::iterator_facade<Lightweight_polyline_2_iterator<Kernel_, Iterator>,
                                  typename Kernel_::Point_2,
                                  typename std::iterator_traits
                                  <typename Lightweight_polyline_2<Kernel_, Iterator>::Zip_iterator>::iterator_category>
{
public:

  using Kernel = Kernel_;
  using Base_iterator = Iterator;
  using Self = Lightweight_polyline_2_iterator<Kernel, Base_iterator>;
  using Polyline = Lightweight_polyline_2<Kernel, Base_iterator>;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;
  using Zip_iterator = typename Polyline::Zip_iterator;

  friend Polyline;

protected:
  const Polyline* m_support;
  Zip_iterator m_base;

public:

  Lightweight_polyline_2_iterator () = delete;

  Lightweight_polyline_2_iterator (const Polyline* support, const Tag_true&) // begin
    : m_support (support)
  {
    if (m_support->m_reverse)
    {
      if (m_support->is_init(m_support->m_last))
        m_base = m_support->m_end;
      else
        m_base = m_support->m_end - 1;
    }
    else
    {
      if (m_support->is_init(m_support->m_first))
        m_base = m_support->m_begin - 1;
      else
        m_base = m_support->m_begin;
    }
  }

  Lightweight_polyline_2_iterator (const Polyline* support, const Tag_false&) // end
    : m_support (support)
  {
    if (m_support->m_reverse)
    {
      if (m_support->is_init(m_support->m_first))
        m_base = m_support->m_begin - 2;
      else
        m_base = m_support->m_begin - 1;
    }
    else
    {
      if (m_support->is_init(m_support->m_last))
        m_base = m_support->m_end + 1;
      else
        m_base = m_support->m_end;
    }
  }

  // interoperability
  template <typename K, typename I>
  Lightweight_polyline_2_iterator (const Lightweight_polyline_2_iterator<K, I>& other)
    : m_support (other.m_support)
    , m_base(other.m_base)
  { }

  const Polyline& support() const { return *m_support; }
  Zip_iterator base() const { return m_base; }

  // The iterator is also used as a wrapper for a segment (using it as
  // source and using the immediate following iterator as target)

  const Point_2& source() const { return const_dereference(); }
  const Point_2& target() const { return std::next(*this).const_dereference(); }

  bool is_vertical() const
  {
    const Point_2& ps = source();
    const Point_2& pt = target();

    return (Kernel().compare_x_2_object()(ps, pt) == EQUAL);
  }

  bool is_directed_right() const
  {
    CGAL_assertion ((Kernel().compare_xy_2_object()(source(), target()) == SMALLER) == m_support->is_directed_right());
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

  const Line_2& line() const
  {
    return m_support->line (cached_line(), source(), target());
  }

  std::shared_ptr<Line_2>& cached_line () const
  {
    CGAL_assertion (m_support != nullptr);
    if (m_base == m_support->m_begin - 1)
    {
      CGAL_assertion (m_support->is_init(m_support->m_first));
      return const_cast<std::shared_ptr<Line_2>&>(m_support->m_first.second);
    }
    else if (m_base == m_support->m_end)
    {
      CGAL_assertion (m_support->is_init(m_support->m_last));
      return const_cast<std::shared_ptr<Line_2>&>(m_support->m_last.second);
    }

    // else
    CGAL_assertion (std::distance (m_support->m_begin, m_base) >= 0);
    CGAL_assertion (std::distance (m_base, m_support->m_end) > 0);
    return const_cast<std::shared_ptr<Line_2>&>(m_base->second);
  }

private:
  friend class boost::iterator_core_access;

  void increment()
  {
    CGAL_assertion (m_support != nullptr);
    if (m_support->m_reverse)
      -- m_base;
    else
      ++ m_base;
  }

  void decrement()
  {
    CGAL_assertion (m_support != nullptr);
    if (m_support->m_reverse)
      ++ m_base;
    else
      -- m_base;
  }

  void advance (std::ptrdiff_t n)
  {
    CGAL_assertion (m_support != nullptr);
    if (m_support->m_reverse)
      m_base -= n;
    else
      m_base += n;
  }

  // interoperability
  template <typename K, typename I>
  std::ptrdiff_t distance_to (const Lightweight_polyline_2_iterator<K,I>& other) const
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (other.m_support != nullptr);
    if (m_support->m_reverse)
      return std::distance (other.m_base, m_base);
    // else
    return std::distance (m_base, other.m_base);
  }

  // interoperability
  template <typename K, typename I>
  bool equal (const Lightweight_polyline_2_iterator<K, I>& other) const
  {
    CGAL_assertion (m_support != nullptr);
    CGAL_assertion (other.m_support != nullptr);
    return m_base == other.m_base;
  }

  Point_2& dereference() const
  {
    return const_cast<Point_2&>(const_dereference());
  }

  const Point_2& const_dereference() const
  {
    CGAL_assertion (m_support != nullptr);
    if (m_base == m_support->m_begin - 1)
    {
      CGAL_assertion (m_support->is_init(m_support->m_first));
      return m_support->point(m_support->m_first);
    }
    else if (m_base == m_support->m_end)
    {
      CGAL_assertion (m_support->is_init(m_support->m_last));
      return m_support->point(m_support->m_last);
    }

    // else
    CGAL_assertion (std::distance (m_support->m_begin, m_base) >= 0);
    CGAL_assertion (std::distance (m_base, m_support->m_end) > 0);
    return m_support->point(*m_base);
  }
};

} // namespace internal

template <typename Curve>
class Indexed_sweep_curve_accessor;

template <typename Kernel, typename Iterator>
class Indexed_sweep_curve_accessor<internal::Lightweight_polyline_2<Kernel, Iterator> >
{
private:

  using Curve = internal::Lightweight_polyline_2<Kernel, Iterator>;
  Curve* m_first;
  std::size_t m_nb_vertices;

public:
  typedef Tag_true valid;

  template <typename It>
  Indexed_sweep_curve_accessor (It begin, It end)
    : m_nb_vertices (std::distance(begin, end) + 1)
    , m_first (&*begin)
  { }

  // Not used
  std::size_t nb_vertices() const { return m_nb_vertices; }
  std::size_t min_end_index (const Curve& c) const
  {
    if (c.subcurves_begin()->is_directed_right())
      return curve_index(c);
    // else
    return curve_index(c) + 1;
  }
  std::size_t max_end_index (const Curve& c) const
  {
    if (c.subcurves_begin()->is_directed_right())
      return curve_index(c) + 1;
    // else
    return curve_index(c);
  }
  const Curve& curve (const Curve& c) const { return c; }

  void before_init() const { }
  void after_init() const { }

private:

  std::size_t curve_index (const Curve& c) const
  {
    std::size_t out = std::size_t(&c - m_first);
//    std::cerr << out << " ";
    return out;
  }
};

} //namespace CGAL

#endif
