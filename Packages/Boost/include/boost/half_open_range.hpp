// (C) Copyright Jeremy Siek and David Abrahams 2000-2001. Permission to copy,
// use, modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided "as is"
// without express or implied warranty, and with no claim as to its suitability
// for any purpose.
//
// Revision History:
// 11 Feb 2001  Use new iterator_adaptor interface, Fixes for Borland.
//              (Dave Abrahams)
// 04 Feb 2001  Support for user-defined iterator categories (Dave Abrahams)
// 30 Jan 2001  Initial Checkin (Dave Abrahams)

#ifndef BOOST_HALF_OPEN_RANGE_HPP_
# define BOOST_HALF_OPEN_RANGE_HPP_

# include <boost/counting_iterator.hpp>
# include <functional>
# include <cassert>
# include <boost/operators.hpp>
# include <string>
# include <stdexcept>
# include <iterator>

namespace boost {

namespace detail {

  // Template class choose_finish -- allows us to maintain the invariant that
  // start() <= finish() on half_open_range specializations that support random
  // access.
#ifdef __MWERKS__
  template <class T>
  const T& choose_finish(const T&, const T& finish, std::input_iterator_tag)
  {
      return finish;
  }

  template <class T>
  const T& choose_finish(const T&, const T& finish, std::output_iterator_tag)
  {
      return finish;
  }

  template <class T>
  const T& choose_finish(const T& start, const T& finish, std::random_access_iterator_tag)
  {
      return finish < start ? start : finish;
  }
#else
  template <bool is_random_access> struct finish_chooser;

  template <>
  struct finish_chooser<false>
  {
      template <class T>
      struct rebind
      {
          static T choose(const T&, const T& finish)
              { return finish; }
      };
  };

  template <>
  struct finish_chooser<true>
  {
      template <class T>
      struct rebind
      {
          static T choose(const T& start, const T& finish)
              { return finish < start ? start : finish; }
      };
  };

  template <class Category, class Incrementable>
  struct choose_finish
  {
      static const Incrementable choose(const Incrementable& start, const Incrementable& finish)
      {
          return finish_chooser<(
              ::boost::is_convertible<Category*,std::random_access_iterator_tag*>::value
              )>::template rebind<Incrementable>::choose(start, finish);
      }
  };
#endif
}

template <class Incrementable>
struct half_open_range
{
    typedef typename counting_iterator_generator<Incrementable>::type iterator;

 private: // utility type definitions
    // Using iter_t prevents compiler confusion with boost::iterator
    typedef typename counting_iterator_generator<Incrementable>::type iter_t;

    typedef std::less<Incrementable> less_value;
    typedef typename iter_t::iterator_category category;
    typedef half_open_range<Incrementable> self;

 public:
    typedef iter_t const_iterator;
    typedef typename iterator::value_type value_type;
    typedef typename iterator::difference_type difference_type;
    typedef typename iterator::reference reference;
    typedef typename iterator::reference const_reference;
    typedef typename iterator::pointer pointer;
    typedef typename iterator::pointer const_pointer;
    
    // It would be nice to select an unsigned type, but this is appropriate
    // since the library makes an attempt to select a difference_type which can
    // hold the difference between any two iterators.
    typedef typename iterator::difference_type size_type;

    half_open_range(Incrementable start, Incrementable finish)
        : m_start(start),
          m_finish(
#ifndef __MWERKS__
            detail::choose_finish<category,Incrementable>::choose(start, finish)
#else
            detail::choose_finish(start, finish, category())
#endif
              )
        {}

    // Implicit conversion from std::pair<Incrementable,Incrementable> allows us
    // to accept the results of std::equal_range(), for example.
    half_open_range(const std::pair<Incrementable,Incrementable>& x)
        : m_start(x.first),
          m_finish(
#ifndef __MWERKS__
              detail::choose_finish<category,Incrementable>::choose(x.first, x.second)
#else
            detail::choose_finish(x.first, x.second, category())
#endif
              )
        {}

    half_open_range& operator=(const self& x)
    {
        m_start = x.m_start;
        m_finish = x.m_finish;
        return *this;
    }
    
    half_open_range& operator=(const std::pair<Incrementable,Incrementable>& x)
    {
        m_start = x.first;
        m_finish =
#ifndef __MWERKS__
            detail::choose_finish<category,Incrementable>::choose(x.first, x.second);
#else
            detail::choose_finish(x.first, x.second, category();
#endif
    }
        
    iterator begin() const { return iterator(m_start); }
    iterator end() const { return iterator(m_finish); }
    
    Incrementable front() const { assert(!this->empty()); return m_start; }
    Incrementable back() const { assert(!this->empty()); return boost::prior(m_finish); }
                                     
    Incrementable start() const { return m_start; }
    Incrementable finish() const { return m_finish; }
    
    size_type size() const { return boost::detail::distance(begin(), end()); }

    bool empty() const
    {
        return m_finish == m_start;
    }

    void swap(half_open_range& x) {
        std::swap(m_start, x.m_start);
        std::swap(m_finish, x.m_finish);
    }
    
 public: // functions requiring random access elements
    
    // REQUIRES: x is reachable from this->front()
    bool contains(const value_type& x) const
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));
        return !less_value()(x, m_start) && less_value()(x, m_finish);
    }

    bool contains(const half_open_range& x) const
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));
        return x.empty() || !less_value()(x.m_start, m_start) && !less_value()(m_finish, x.m_finish);
    }

    bool intersects(const half_open_range& x) const
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));
        return less_value()(
            less_value()(this->m_start, x.m_start) ? x.m_start : this->m_start,
            less_value()(this->m_finish, x.m_finish) ? this->m_finish : x.m_finish);
    }

    half_open_range& operator&=(const half_open_range& x)
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));
        
        if (less_value()(this->m_start, x.m_start))
            this->m_start = x.m_start;
        
        if (less_value()(x.m_finish, this->m_finish))
            this->m_finish = x.m_finish;

        if (less_value()(this->m_finish, this->m_start))
            this->m_start = this->m_finish;

        return *this;
    }

    half_open_range& operator|=(const half_open_range& x)
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));

        if (!x.empty())
        {
            if (this->empty())
            {
                *this = x;
            }
            else
            {
                if (less_value()(x.m_start, this->m_start))
                    this->m_start = x.m_start;
        
                if (less_value()(this->m_finish, x.m_finish))
                    this->m_finish = x.m_finish;
            }
        }
        return *this;
    }

    // REQUIRES: x is reachable from this->front()
    const_iterator find(const value_type& x) const
    {
        BOOST_STATIC_ASSERT((boost::is_same<category, std::random_access_iterator_tag>::value));
        
        return const_iterator(this->contains(x) ? x : m_finish);
    }

    // REQUIRES: index >= 0 && index < size()
    value_type operator[](size_type index) const
    {
        assert(index >= 0 && index < size());
        return m_start + index;
    }
    
    value_type at(size_type index) const
    {
        if (index < 0 || index >= size())
            throw std::out_of_range(std::string("half_open_range"));
        return m_start + index;
    }
    
 private: // data members
    Incrementable m_start, m_finish;
};

template <class Incrementable>
half_open_range<Incrementable> operator|(
    half_open_range<Incrementable> x,
    const half_open_range<Incrementable>& y)
{
    return x |= y;
}

template <class Incrementable>
half_open_range<Incrementable> operator&(
    half_open_range<Incrementable> x,
    const half_open_range<Incrementable>& y)
{
    return x &= y;
}

template <class Incrementable>
inline bool operator==(
    const half_open_range<Incrementable>& x,
    const half_open_range<Incrementable>& y)
{
    const bool y_empty = y.empty();
    return x.empty() ? y_empty : !y_empty && x.start() == y.start() && x.finish() == y.finish();
}

template <class Incrementable>
inline bool operator!=(
    const half_open_range<Incrementable>& x,
    const half_open_range<Incrementable>& y)
{
    return !(x == y);
}

template <class Incrementable>
inline half_open_range<Incrementable>
make_half_open_range(Incrementable first, Incrementable last)
{
  return half_open_range<Incrementable>(first, last);
}

template <class Incrementable>
bool intersects(
    const half_open_range<Incrementable>& x,
    const half_open_range<Incrementable>& y)
{
    return x.intersects(y);
}
    
template <class Incrementable>
bool contains(
    const half_open_range<Incrementable>& x,
    const half_open_range<Incrementable>& y)
{
    return x.contains(y);
}
    
} // namespace boost

#ifndef BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION

namespace std {
template <class Incrementable> struct less<boost::half_open_range<Incrementable> >
        : binary_function<
            boost::half_open_range<Incrementable>,
            boost::half_open_range<Incrementable>,bool>
{
    bool operator()(
        const boost::half_open_range<Incrementable>& x,
        const boost::half_open_range<Incrementable>& y) const
    {
        less<Incrementable> cmp;
        return !y.empty() && (
            cmp(x.start(), y.start())
            || !cmp(y.start(), x.start())
               && cmp(x.finish(), y.finish()));
    }
};

template <class Incrementable> struct less_equal<boost::half_open_range<Incrementable> >
        : binary_function<
            boost::half_open_range<Incrementable>,
            boost::half_open_range<Incrementable>,bool>
{
    bool operator()(
        const boost::half_open_range<Incrementable>& x,
        const boost::half_open_range<Incrementable>& y) const
    {
        typedef boost::half_open_range<Incrementable> range;
        less<range> cmp;
        return !cmp(y,x);
    }
};
template <class Incrementable> struct greater<boost::half_open_range<Incrementable> >
        : binary_function<
            boost::half_open_range<Incrementable>,
            boost::half_open_range<Incrementable>,bool>
{
    bool operator()(
        const boost::half_open_range<Incrementable>& x,
        const boost::half_open_range<Incrementable>& y) const
    {
        typedef boost::half_open_range<Incrementable> range;
        less<range> cmp;
        return cmp(y,x);
    }
};

template <class Incrementable> struct greater_equal<boost::half_open_range<Incrementable> >
        : binary_function<
            boost::half_open_range<Incrementable>,
            boost::half_open_range<Incrementable>,bool>
{
    bool operator()(
        const boost::half_open_range<Incrementable>& x,
        const boost::half_open_range<Incrementable>& y) const
    {
        typedef boost::half_open_range<Incrementable> range;
        less<range> cmp;
        return !cmp(x,y);
    }
};
} // namespace std

#else

namespace boost {
// Can't partially specialize std::less et al, so we must provide the operators
template <class Incrementable>
bool operator<(const half_open_range<Incrementable>& x,
               const half_open_range<Incrementable>& y)
{
    return !y.empty() && (
        x.empty() || std::less<Incrementable>()(x.start(), y.start())
        || !std::less<Incrementable>()(y.start(), x.start())
                && std::less<Incrementable>()(x.finish(), y.finish()));
}

template <class Incrementable>
bool operator>(const half_open_range<Incrementable>& x,
               const half_open_range<Incrementable>& y)
{
    return y < x;
}

template <class Incrementable>
bool operator<=(const half_open_range<Incrementable>& x,
               const half_open_range<Incrementable>& y)
{
    return !(y < x);
}

template <class Incrementable>
bool operator>=(const half_open_range<Incrementable>& x,
               const half_open_range<Incrementable>& y)
{
    return !(x < y);
}
} // namespace boost

#endif
    

#endif // BOOST_HALF_OPEN_RANGE_HPP_
