#ifndef CGAL_LIGHTWEIGHT_POLYLINE_2_H
#define CGAL_LIGHTWEIGHT_POLYLINE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/iterator.h>

#include <boost/iterator/iterator_facade.hpp>


#ifdef CGAL_USE_LIGHTWEIGHT_POLYLINES_WITH_CACHE
#define CGAL_LIGHTWEIGHT_POLYLINE_EMBED_CACHE
#endif

namespace CGAL {

namespace internal {

template <typename Kernel_, typename Iterator>
class Lightweight_polyline_2_iterator;

template <typename Kernel_, typename IsPointIterator>
class Lightweight_polyline_2_base;

template <typename Kernel_>
class Lightweight_polyline_2_base<Kernel_, Tag_true>
{
public:
  using Kernel = Kernel_;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;

  using Point_ref = const Point_2&;
  using Line_ref = Line_2;
  using Input_point = Point_2;
  using Hanging_point = std::shared_ptr<Point_2>;

  Hanging_point uninit() const { return nullptr; }
  bool is_init (const Hanging_point& p) const { return (p != nullptr); }
  Hanging_point init (const Input_point& p) const { return std::make_shared<Point_2>(p); }
  Hanging_point init (const Hanging_point& p) const { return p; }

  Point_ref point (const Input_point& p) const { return p; }
  Point_ref point (const Hanging_point& p) const { return *p; }

  template <typename PointA, typename PointB>
  Line_ref line (const std::nullptr_t&, const PointA& a, const PointB& b) const { return Line_2(point(a),point(b)); }
};

template <typename Kernel_>
class Lightweight_polyline_2_base<Kernel_, Tag_false>
{
public:
  using Kernel = Kernel_;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;

  using Point_ref = const Point_2&;
//  using Line_ref = const Line_2&;
  using Line_ref = Line_2;
  using Input_point = std::pair<const Point_2&, const std::shared_ptr<Line_2>& >;

  using Hanging_point = std::pair<std::shared_ptr<Point_2>, std::shared_ptr<Line_2> >;

  Hanging_point uninit() const { return Hanging_point(nullptr, nullptr); }
  bool is_init (const Hanging_point& p) const { return (p.first != nullptr); }
  Hanging_point init (const Point_2& p) const { return Hanging_point(std::make_shared<Point_2>(p), nullptr); }
  Hanging_point init (const Input_point& p) const { return Hanging_point(std::make_shared<Point_2>(p.first), p.second); }
  Hanging_point init (const Hanging_point& p) const { return p; }

  Point_ref point (const Point_2& p) const { return p; }
  Point_ref point (const Input_point& p) const { return p.first; }
  Point_ref point (const Hanging_point& p) const { return *p.first; }

  template <typename PointA, typename PointB>
  Line_ref line (std::shared_ptr<Line_2>& line, const PointA& a, const PointB& b) const
  {
    if (!line)
      line = std::make_shared<Line_2>(point(a), point(b));
    return *line;
  }
};


template <typename Kernel_, typename Iterator>
class Lightweight_polyline_2
  : public Lightweight_polyline_2_base
<Kernel_,
 Boolean_tag<std::is_same<typename std::iterator_traits<Iterator>::value_type,
                          typename Kernel_::Point_2>::value> >
{
public:

  using Kernel = Kernel_;
  using Base_iterator = Iterator;

  using Is_point_iterator = Boolean_tag<std::is_same<typename std::iterator_traits<Iterator>::value_type,
                                                     typename Kernel_::Point_2>::value>;
  using Base = Lightweight_polyline_2_base<Kernel, Is_point_iterator>;
  using Self = Lightweight_polyline_2<Kernel, Base_iterator>;

  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;

  using Size = std::size_t;
  using size_type = std::size_t;

  using iterator = Lightweight_polyline_2_iterator<Kernel_, Base_iterator>;
  friend iterator;

  using Input_point = typename Base::Input_point;
  using Hanging_point = typename Base::Hanging_point;

  using Subcurve_type_2 = iterator;
  using Subcurve_iterator = Prevent_deref<iterator>;
  using Subcurve_const_iterator = Prevent_deref<iterator>;

  using Base::is_init;
  using Base::init;

protected:

  Base_iterator m_begin;
  Base_iterator m_end;
  Hanging_point m_first;
  Hanging_point m_last;
  bool m_reverse;
  bool m_is_directed_right;

#ifdef CGAL_LIGHTWEIGHT_POLYLINE_EMBED_CACHE
  using Line_ptr = std::shared_ptr<Line_2>;
  std::shared_ptr<std::vector<Line_ptr> > m_cache;
#endif

public:

  Lightweight_polyline_2() : m_reverse(false) { }

  Lightweight_polyline_2(const Point_2& first, const Point_2& last)
    : m_first(init(first))
    , m_last(init(last))
    , m_reverse(false)
  {
    compute_direction();
  }

#ifdef CGAL_LIGHTWEIGHT_POLYLINE_EMBED_CACHE
  template <typename PointIterator>
  Lightweight_polyline_2 (PointIterator begin, PointIterator end, bool force_closure = false)
    : m_reverse(false)
  {
    m_cache = std::make_shared<std::vector<Line_ptr> >(std::distance(begin, end), nullptr);
    m_begin = Base_iterator(std::make_pair(begin, m_cache->begin()));
    m_end = Base_iterator(std::make_pair(end, m_cache->end()));
    if (force_closure)
      m_last = init(*m_begin);
    compute_direction();
  }
#endif

  Lightweight_polyline_2 (Base_iterator begin, Base_iterator end, bool force_closure = false)
    : m_begin(begin), m_end(end)
    , m_reverse(false)
  {
    CGAL_assertion (std::distance (begin, end) >= 2);
    if (force_closure)
      m_last = init(*m_begin);
    compute_direction();
  }

  Lightweight_polyline_2 (iterator begin, iterator end)
    : m_reverse(false)
#ifdef CGAL_LIGHTWEIGHT_POLYLINE_EMBED_CACHE
    , m_cache (begin.support().m_cache)
#endif
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

  Lightweight_polyline_2 (Hanging_point first, iterator begin, iterator end, Hanging_point last)
    : m_reverse(false)
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

  Hanging_point init (const Point_2& p, std::size_t index) const
  {
    Hanging_point out = this->init(p);

#ifdef CGAL_LIGHTWEIGHT_POLYLINE_EMBED_CACHE
    if (index == 0)
      out.second = m_first.second;
    else if (index == number_of_subcurves() + 1)
      out.second = m_last.second;
    else
      out.second = (m_begin + index - 1)->second;
#endif
    return out;
  }

  void compute_direction()
  {
#ifdef CGAL_PROFILE
    std::size_t nb_curves = number_of_subcurves();
    if (nb_curves == 1)
    {
      if (m_first && m_last)
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

};

template <typename Kernel_, typename Iterator>
class Lightweight_polyline_2_iterator
  : public boost::iterator_facade<Lightweight_polyline_2_iterator<Kernel_, Iterator>,
                                  typename Kernel_::Point_2,
                                  typename std::iterator_traits<Iterator>::iterator_category>
{
public:

  using Kernel = Kernel_;
  using Base_iterator = Iterator;
  using Self = Lightweight_polyline_2_iterator<Kernel, Base_iterator>;
  using Polyline = Lightweight_polyline_2<Kernel, Base_iterator>;
  using Line_is_cached = Boolean_tag<!Polyline::Is_point_iterator::value>;
  using Point_2 = typename Kernel::Point_2;
  using Line_2 = typename Kernel::Line_2;
  using Line_ref = typename Polyline::Line_ref;
  using Input_point = typename Polyline::Input_point;
  using Hanging_point = typename Polyline::Hanging_point;
  friend Polyline;

protected:
  const Polyline* m_support;
  Base_iterator m_base;

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
  Base_iterator base() const { return m_base; }

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

  Line_ref line() const
  {
    return m_support->line (cached_line(Line_is_cached()), source(), target());
  }

  std::shared_ptr<Line_2>& cached_line (const Tag_true&) const
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

  std::nullptr_t cached_line (const Tag_false&) const
  {
    return nullptr;
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
