// Boost.Signals library

// Copyright Doug Gregor 2001-2003. Use, modification and
// distribution is subject to the Boost Software License, Version
// 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

#ifndef BOOST_SIGNALS_CONNECTION_HPP
#define BOOST_SIGNALS_CONNECTION_HPP

#include <boost/signals/detail/signals_common.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/operators.hpp>
#include <boost/any.hpp>
#include <list>
#include <cassert>
#include <utility>

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

namespace boost {
  namespace BOOST_SIGNALS_NAMESPACE {
    class trackable;

    namespace detail {
      // Represents an object that has been bound as part of a slot, and how
      // to notify that object of a disconnect
      struct bound_object {
        void* obj;
        void* data;
        void (*disconnect)(void*, void*);

        bool operator==(const bound_object& other) const
          { return obj == other.obj && data == other.data; }
        bool operator<(const bound_object& other) const
          { return obj < other.obj; }
      };

      // Describes the connection between a signal and the objects that are
      // bound for a specific slot. Enables notification of the signal and the
      // slots when a disconnect is requested.
      struct basic_connection {
        void* signal;
        void* signal_data;
        void (*signal_disconnect)(void*, void*);

        std::list<bound_object> bound_objects;
      };
    } // end namespace detail

    // The user may freely pass around the "connection" object and terminate
    // the connection at any time using disconnect().
    class BOOST_SIGNALS_DECL connection :
      private less_than_comparable1<connection>,
      private equality_comparable1<connection>
    {
    public:
      connection() : con(), controlling_connection(false) {}
      connection(const connection&);
      ~connection();

      // Disconnect the signal and slot, if they are connected
      void disconnect() const;

      // Returns true if the signal and slot are connected
      bool connected() const { return con.get() && con->signal_disconnect; }

      // Comparison of connections
      inline bool operator==(const connection& other) const;
      inline bool operator<(const connection& other) const;

      // Connection assignment
      connection& operator=(const connection& other) ;

      // Swap connections
      inline void swap(connection& other);

    public: // TBD: CHANGE THIS
      // Set whether this connection object is controlling or not
      void set_controlling() { controlling_connection = true; }

    private:
      friend class detail::signal_base_impl;
      friend class detail::slot_base;
      friend class trackable;

      // Reset this connection to refer to a different actual connection
      void reset(BOOST_SIGNALS_NAMESPACE::detail::basic_connection*);

      // Add a bound object to this connection (not for users)
      void add_bound_object(const BOOST_SIGNALS_NAMESPACE::detail::bound_object& b);

      friend class BOOST_SIGNALS_NAMESPACE::detail::bound_objects_visitor;

      // Pointer to the actual contents of the connection
      shared_ptr<BOOST_SIGNALS_NAMESPACE::detail::basic_connection> con;

      // True if the destruction of this connection object should disconnect
      bool controlling_connection;
    };

    // Similar to connection, but will disconnect the connection when it is
    // destroyed unless release() has been called.
    class BOOST_SIGNALS_DECL scoped_connection : public connection {
    public:
      scoped_connection() : connection(), released(false) {}
      scoped_connection(const connection&);
      scoped_connection(const scoped_connection&);
      ~scoped_connection();

      connection release();

      inline void swap(scoped_connection&);

      scoped_connection& operator=(const connection&);
      scoped_connection& operator=(const scoped_connection&);

    private:
      bool released;
    };

    inline connection::connection(const connection& other) :
      con(other.con), controlling_connection(other.controlling_connection)
    {
    }

    inline connection::~connection()
    {
      if (controlling_connection) {
        disconnect();
      }
    }

    inline void
    connection::reset(BOOST_SIGNALS_NAMESPACE::detail::basic_connection* new_con)
    {
      con.reset(new_con);
    }

    inline bool connection::operator==(const connection& other) const
    {
      return con.get() == other.con.get();
    }

    inline bool connection::operator<(const connection& other) const
    {
      return con.get() < other.con.get();
    }

    inline connection& connection::operator=(const connection& other)
    {
      connection(other).swap(*this);
      return *this;
    }

    inline void connection::swap(connection& other)
    {
      this->con.swap(other.con);
      std::swap(this->controlling_connection, other.controlling_connection);
    }

    inline void swap(connection& c1, connection& c2)
    {
      c1.swap(c2);
    }

    inline scoped_connection::scoped_connection(const connection& other) :
      connection(other),
      released(false)
    {
    }

    inline
    scoped_connection::scoped_connection(const scoped_connection& other) :
      connection(other),
      released(other.released)
    {
    }

    inline scoped_connection::~scoped_connection()
    {
      if (!released) {
        this->disconnect();
      }
    }

    inline connection scoped_connection::release()
    {
      released = true;
      return *this;
    }

    inline void scoped_connection::swap(scoped_connection& other)
    {
      this->connection::swap(other);
      bool other_released = other.released;
      other.released = this->released;
      this->released = other_released;
    }

    inline void swap(scoped_connection& c1, scoped_connection& c2)
    {
      c1.swap(c2);
    }

    inline scoped_connection&
    scoped_connection::operator=(const connection& other)
    {
      scoped_connection(other).swap(*this);
      return *this;
    }

    inline scoped_connection&
    scoped_connection::operator=(const scoped_connection& other)
    {
      scoped_connection(other).swap(*this);
      return *this;
    }

    namespace detail {
      struct connection_slot_pair {
        connection first;
        any second;

        connection_slot_pair() {}

        connection_slot_pair(const connection& c, const any& a)
          : first(c), second(a)
        {
        }

        // Dummys to allow explicit instantiation to work
        bool operator==(const connection_slot_pair&) const { return false; }
        bool operator<(const connection_slot_pair&) const { return false;}
      };

      // Determines if the underlying connection is disconnected
      struct is_disconnected {
        typedef std::pair<const any, connection_slot_pair> argument_type;
        typedef bool result_type;

        inline bool operator()(const argument_type& c) const
        {
          return !c.second.first.connected();
        }
      };

      // Autodisconnects the bound object when it is destroyed unless the
      // release method is invoked.
      class auto_disconnect_bound_object {
      public:
        auto_disconnect_bound_object(const bound_object& b) :
          binding(b), auto_disconnect(true)
        {
        }

        ~auto_disconnect_bound_object()
        {
          if (auto_disconnect)
            binding.disconnect(binding.obj, binding.data);
        }

        void release() { auto_disconnect = false; }

      private:
        bound_object binding;
        bool auto_disconnect;
      };
    } // end namespace detail
  } // end namespace BOOST_SIGNALS_NAMESPACE
} // end namespace boost

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif

#endif // BOOST_SIGNALS_CONNECTION_HPP
