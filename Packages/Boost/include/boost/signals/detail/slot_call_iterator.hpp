// Boost.Signals library

// Copyright Douglas Gregor 2001-2004. Use, modification and
// distribution is subject to the Boost Software License, Version
// 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

#ifndef BOOST_SIGNALS_SLOT_CALL_ITERATOR
#define BOOST_SIGNALS_SLOT_CALL_ITERATOR

#include <functional>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/signals/detail/config.hpp>
#include <boost/signals/connection.hpp>

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost {
  namespace BOOST_SIGNALS_NAMESPACE {
    namespace detail {
      // A cached return value from a slot
      template<typename T>
      struct cached_return_value {
        cached_return_value(const T& t) : value(t) {}

        T value;
      };

      // Generates a slot call iterator. Essentially, this is an iterator that:
      //   - skips over disconnected slots in the underlying list
      //   - calls the connected slots when dereferenced
      //   - caches the result of calling the slots
      template<typename Function, typename Iterator>
      class slot_call_iterator
        : public iterator_facade<slot_call_iterator<Function, Iterator>,
                                 typename Function::result_type,
                                 single_pass_traversal_tag,
                                 typename Function::result_type const&>
      {
        typedef iterator_facade<slot_call_iterator<Function, Iterator>,
                                typename Function::result_type,
                                single_pass_traversal_tag,
                                typename Function::result_type const&>
          inherited;

        typedef typename Function::result_type result_type;

        friend class iterator_core_access;

      public:
        slot_call_iterator() {}

        slot_call_iterator(Iterator iter_in, Iterator end_in, Function f)
          : iter(iter_in), end(end_in), f(f), cache()
        {
          iter = std::find_if(iter, end, std::not1(is_disconnected()));
        }

        typename inherited::reference
        dereference() const
        {
          if (!cache.get()) {
            cache.reset(new cached_return_value<result_type>(f(*iter)));
          }

          return cache->value;
        }

        void increment()
        {
          iter = std::find_if(++iter, end, std::not1(is_disconnected()));
          cache.reset();
        }

        bool equal(const slot_call_iterator& other) const
        {
          iter = std::find_if(iter, end, std::not1(is_disconnected()));
          other.iter = std::find_if(other.iter, other.end,
                                    std::not1(is_disconnected()));
          return iter == other.iter;
        }

      private:
        mutable Iterator iter;
        Iterator end;
        Function f;
        mutable shared_ptr< cached_return_value<result_type> > cache;
      };
    } // end namespace detail
  } // end namespace BOOST_SIGNALS_NAMESPACE
} // end namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif // BOOST_SIGNALS_SLOT_CALL_ITERATOR
