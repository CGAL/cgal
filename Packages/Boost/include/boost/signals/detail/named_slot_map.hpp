// Boost.Signals library

// Copyright Douglas Gregor 2001-2004. Use, modification and
// distribution is subject to the Boost Software License, Version
// 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

#ifndef BOOST_SIGNALS_NAMED_SLOT_MAP_HPP
#define BOOST_SIGNALS_NAMED_SLOT_MAP_HPP

#include <boost/signals/detail/config.hpp>
#include <boost/signals/detail/signals_common.hpp>
#include <boost/signals/connection.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/any.hpp>
#include <boost/utility.hpp>
#include <boost/function/function2.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <memory>
#include <utility>

namespace boost { namespace BOOST_SIGNALS_NAMESPACE {

enum connect_position { at_back, at_front };

namespace detail {

typedef function2<bool, any, any> compare_type;

// Used to delimit the front and back of the list for O(1) insertion.
struct front_type {};
struct back_type {};

// This function object bridges from a pair of any objects that hold
// values of type Key to the underlying function object that compares
// values of type Key.
template<typename Compare, typename Key>
class any_bridge_compare {
public:
  typedef bool result_type;
  typedef const any& first_argument_type;
  typedef const any& second_argument_type;

  any_bridge_compare(const Compare& c) : comp(c) {}

  bool operator()(const any& k1, const any& k2) const
  {
    if (k1.type() == typeid(front_type))
      return !(k2.type() == typeid(front_type));
    if (k1.type() == typeid(back_type))
      return false;
    if (k2.type() == typeid(front_type))
      return false;
    if (k2.type() == typeid(back_type))
      return true;

    // Neither is empty, so compare their values to order them
    // The strange */& is so that we will get a reference to the
    // value stored in the any object instead of a copy
    return comp(*any_cast<Key>(&k1), *any_cast<Key>(&k2));
  }

private:
  Compare comp;
};

class BOOST_SIGNALS_DECL named_slot_map_iterator :
  public iterator_facade<named_slot_map_iterator,
                         connection_slot_pair,
                         forward_traversal_tag>
{
  class impl;

  typedef iterator_facade<named_slot_map_iterator,
                          connection_slot_pair,
                          forward_traversal_tag> inherited;
public:
  named_slot_map_iterator();
  named_slot_map_iterator(const named_slot_map_iterator& other);
  ~named_slot_map_iterator();
  named_slot_map_iterator& operator=(const named_slot_map_iterator& other);

  connection_slot_pair& dereference() const;
  void increment();
  bool equal(const named_slot_map_iterator& other) const;

#if BOOST_WORKAROUND(BOOST_MSVC, <= 0x1701)
  void decrement();
  void advance(difference_type);
#endif

private:
  named_slot_map_iterator(std::auto_ptr<impl>);

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
  shared_ptr<impl> impl_;
#else
  scoped_ptr<impl> impl_;
#endif

  friend class named_slot_map;
};

class BOOST_SIGNALS_DECL named_slot_map
{
public:
  typedef named_slot_map_iterator iterator;

  named_slot_map(const compare_type& compare);
  ~named_slot_map();

  void clear();
  iterator begin();
  iterator end();
  iterator insert(const any& name, const connection& con, const any& slot,
                  connect_position at);
  void disconnect(const any& name);
  void erase(iterator pos);
  void remove_disconnected_slots();

private:
  class impl;

#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x564))
  shared_ptr<impl> impl_;
#else
  scoped_ptr<impl> impl_;
#endif
};

} } }

#endif // BOOST_SIGNALS_NAMED_SLOT_MAP_HPP
