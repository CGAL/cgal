// Boost.Signals library

// Copyright Doug Gregor 2001-2003. Use, modification and
// distribution is subject to the Boost Software License, Version
// 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

#ifndef BOOST_SIGNALS_SIGNAL_BASE_HEADER
#define BOOST_SIGNALS_SIGNAL_BASE_HEADER

#include <boost/signals/detail/config.hpp>
#include <boost/signals/detail/signals_common.hpp>
#include <boost/signals/connection.hpp>
#include <boost/signals/trackable.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/any.hpp>
#include <boost/utility.hpp>
#include <boost/function/function2.hpp>
#include <map>
#include <utility>
#include <vector>

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost {
  namespace BOOST_SIGNALS_NAMESPACE {
    namespace detail {
      // Forward declaration for the mapping from slot names to connections
      class named_slot_map;

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
          // if k1 is empty, then it precedes nothing
          if (k1.empty())
            return false;

          // if k2 is empty, then k1 must precede it
          if (k2.empty())
            return true;

          // Neither is empty, so compare their values to order them
          // The strange */& is so that we will get a reference to the
          // value stored in the any object instead of a copy
          return comp(*any_cast<Key>(&k1), *any_cast<Key>(&k2));
        }

      private:
        Compare comp;
      };

      // Must be constructed before calling the slots, because it safely
      // manages call depth
      class BOOST_SIGNALS_DECL call_notification {
      public:
        call_notification(const shared_ptr<signal_base_impl>&);
        ~call_notification();

        shared_ptr<signal_base_impl> impl;
      };

      // Implementation of base class for all signals. It handles the
      // management of the underlying slot lists.
      class BOOST_SIGNALS_DECL signal_base_impl {
      public:
        friend class call_notification;

        typedef function2<bool, any, any> compare_type;

        // Make sure that an exception does not cause the "clearing" flag to
        // remain set
        class temporarily_set_clearing {
        public:
          temporarily_set_clearing(signal_base_impl* b) : base(b)
          {
            base->flags.clearing = true;
          }

          ~temporarily_set_clearing()
          {
            base->flags.clearing = false;
          }

        private:
          signal_base_impl* base;
        };

        friend class temporarily_set_clearing;

        signal_base_impl(const compare_type&);
        ~signal_base_impl();

        // Disconnect all slots connected to this signal
        void disconnect_all_slots();

        // Are there any connected slots?
        bool empty() const;

        // The number of connected slots
        std::size_t num_slots() const;

        // Disconnect all slots in the given group
        void disconnect(const any&);

        // We're being notified that a slot has disconnected
        static void slot_disconnected(void* obj, void* data);

        connection connect_slot(const any& slot,
                                const any& name,
                                const std::vector<const trackable*>&);

      private:
        // Remove all of the slots that have been marked "disconnected"
        void remove_disconnected_slots() const;

      public:
        // Our call depth when invoking slots (> 1 when we have a loop)
        mutable int call_depth;

        struct {
          // True if some slots have disconnected, but we were not able to
          // remove them from the list of slots because there are valid
          // iterators into the slot list
          mutable bool delayed_disconnect:1;

          // True if we are disconnecting all disconnected slots
          bool clearing:1;
        } flags;

        // Slots
        typedef std::multimap<any, connection_slot_pair, compare_type>
          slot_container_type;
        typedef slot_container_type::iterator slot_iterator;
        typedef slot_container_type::value_type stored_slot_type;
        mutable slot_container_type slots_;
      };

      class BOOST_SIGNALS_DECL signal_base : public noncopyable {
      public:
        typedef signal_base_impl::compare_type compare_type;

        friend class call_notification;

        signal_base(const compare_type& comp);
        ~signal_base();

      public:
        // Disconnect all slots connected to this signal
        void disconnect_all_slots() { impl->disconnect_all_slots(); }

        // Are there any connected slots?
        bool empty() const { return impl->empty(); }

        // How many slots are connected?
        std::size_t num_slots() const { return impl->num_slots(); }

      protected:
        connection connect_slot(const any& slot,
                                const any& name,
                                const std::vector<const trackable*>& bound)
        {
          return impl->connect_slot(slot, name, bound);
        }

        typedef signal_base_impl::slot_iterator slot_iterator;
        typedef signal_base_impl::stored_slot_type stored_slot_type;

        shared_ptr<signal_base_impl> impl;
      };
    } // end namespace detail
  } // end namespace BOOST_SIGNALS_NAMESPACE
} // end namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif // BOOST_SIGNALS_SIGNAL_BASE_HEADER
