// Boost.Signals library

// Copyright Doug Gregor 2001-2003. Use, modification and
// distribution is subject to the Boost Software License, Version
// 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

// This file intentionally does not have include guards, because it is meant
// to be included multiple times (one for each signalN class). The
// BOOST_SIGNALS_SIGNAL_TEMPLATE_HEADER_INCLUDED macro merely serves to
// suppress reinclusion of the files that this header depends on.

#ifndef BOOST_SIGNALS_SIGNAL_TEMPLATE_HEADER_INCLUDED
#define BOOST_SIGNALS_SIGNAL_TEMPLATE_HEADER_INCLUDED
#  include <boost/config.hpp>
#  include <boost/signals/connection.hpp>
#  include <boost/utility.hpp>
#  include <boost/ref.hpp>
#  include <boost/signals/slot.hpp>
#  include <boost/last_value.hpp>
#  include <boost/signals/detail/signal_base.hpp>
#  include <boost/signals/detail/slot_call_iterator.hpp>
#  include <cassert>
#  include <functional>
#  include <memory>
#  include <string>
#endif // !BOOST_SIGNALS_SIGNAL_TEMPLATE_HEADER_INCLUDED

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_PREFIX
#endif

// Include the appropriate functionN header
#define BOOST_SIGNAL_FUNCTION_N_HEADER BOOST_JOIN(<boost/function/function,BOOST_SIGNALS_NUM_ARGS.hpp>)
#include BOOST_SIGNAL_FUNCTION_N_HEADER

// Determine if a comma should follow a listing of the arguments/parameters
#if BOOST_SIGNALS_NUM_ARGS == 0
#  define BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
#else
#  define BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS ,
#endif // BOOST_SIGNALS_NUM_ARGS > 0

// Define class names used
#define BOOST_SIGNALS_SIGNAL BOOST_JOIN(signal,BOOST_SIGNALS_NUM_ARGS)
#define BOOST_SIGNALS_FUNCTION BOOST_JOIN(function,BOOST_SIGNALS_NUM_ARGS)
#define BOOST_SIGNALS_ARGS_STRUCT BOOST_JOIN(args,BOOST_SIGNALS_NUM_ARGS)
#define BOOST_SIGNALS_CALL_BOUND BOOST_JOIN(call_bound,BOOST_SIGNALS_NUM_ARGS)

// Define commonly-used instantiations
#define BOOST_SIGNALS_ARGS_STRUCT_INST \
  BOOST_SIGNALS_NAMESPACE::detail::BOOST_SIGNALS_ARGS_STRUCT<BOOST_SIGNALS_TEMPLATE_ARGS>

namespace boost {
  namespace BOOST_SIGNALS_NAMESPACE {
    namespace detail {
      // Holds the arguments for a bound slot call in a single place
      template<BOOST_SIGNALS_TEMPLATE_PARMS
               BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
               typename Dummy = int>
      struct BOOST_SIGNALS_ARGS_STRUCT {
        BOOST_SIGNALS_ARGS_STRUCT(BOOST_SIGNALS_COPY_PARMS)
          BOOST_SIGNALS_INIT_ARGS
        {
        }

        BOOST_SIGNALS_ARGS_AS_MEMBERS
      };

      // Function object that calls the function object given to it, passing
      // the bound arguments along to that underlying function object
      template<typename R>
      struct BOOST_SIGNALS_CALL_BOUND {
        template<BOOST_SIGNALS_TEMPLATE_PARMS
                 BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
                 typename F>
        struct caller {
          typedef BOOST_SIGNALS_ARGS_STRUCT<BOOST_SIGNALS_TEMPLATE_ARGS>*
            args_type;

          args_type args;

          typedef R result_type;

          caller() {}
          caller(args_type a) : args(a) {}

          template<typename Pair>
          R operator()(const Pair& slot) const
          {
            F* target = const_cast<F*>(any_cast<F>(&slot.second.second));
            return (*target)(BOOST_SIGNALS_BOUND_ARGS);
          }
        };
      };

      template<>
      struct BOOST_SIGNALS_CALL_BOUND<void> {
        template<BOOST_SIGNALS_TEMPLATE_PARMS
                 BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
                 typename F>
        struct caller {
          typedef BOOST_SIGNALS_ARGS_STRUCT<BOOST_SIGNALS_TEMPLATE_ARGS>*
            args_type;

          args_type args;

          typedef unusable result_type;

          caller(args_type a) : args(a) {}

          template<typename Pair>
          unusable operator()(const Pair& slot) const
          {
            F* target = const_cast<F*>(any_cast<F>(&slot.second.second));
            (*target)(BOOST_SIGNALS_BOUND_ARGS);
            return unusable();
          }
        };
      };
    } // namespace detail
  } // namespace BOOST_SIGNALS_NAMESPACE

  // The actual signalN class
  template<
    typename R,
    BOOST_SIGNALS_TEMPLATE_PARMS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    typename Combiner = last_value<R>,
    typename Group = int,
    typename GroupCompare = std::less<Group>,
    typename SlotFunction = BOOST_SIGNALS_FUNCTION<
                              R BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
                              BOOST_SIGNALS_TEMPLATE_ARGS>
  >
  class BOOST_SIGNALS_SIGNAL :
    public BOOST_SIGNALS_NAMESPACE::detail::signal_base, // management of slot list
    public BOOST_SIGNALS_NAMESPACE::trackable // signals are trackable
  {
  public:
    // The slot function type
    typedef SlotFunction slot_function_type;

    // Result type of a slot
    typedef typename BOOST_SIGNALS_NAMESPACE::detail::slot_result_type<R>::type
      slot_result_type;

    // Argument types
    BOOST_SIGNALS_ARG_TYPES

#if BOOST_SIGNALS_NUM_ARGS == 1
    typedef T1 argument_type;
#elif BOOST_SIGNALS_NUM_ARGS == 2
    typedef T1 first_argument_type;
    typedef T2 second_argument_type;
#endif

  private:
    // The real slot name comparison object type
    typedef BOOST_SIGNALS_NAMESPACE::detail::any_bridge_compare<GroupCompare, Group>
      real_group_compare_type;

    // The function object passed to the slot call iterator that will call
    // the underlying slot function with its arguments bound
    typedef BOOST_SIGNALS_NAMESPACE::detail::BOOST_SIGNALS_CALL_BOUND<R>
      outer_bound_slot_caller;
    typedef typename outer_bound_slot_caller::template
              caller<BOOST_SIGNALS_TEMPLATE_ARGS
                     BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
                     slot_function_type>
      call_bound_slot;

  public:
    // Combiner's result type
    typedef typename Combiner::result_type result_type;

    // Combiner type
    typedef Combiner combiner_type;

    // Slot type
    typedef slot<slot_function_type> slot_type;

    // Slot name type and comparison
    typedef Group group_type;
    typedef GroupCompare group_compare_type;

    typedef BOOST_SIGNALS_NAMESPACE::detail::slot_call_iterator<
              call_bound_slot, slot_iterator> slot_call_iterator;

    explicit
    BOOST_SIGNALS_SIGNAL(const Combiner& c = Combiner(),
                         const GroupCompare& comp = GroupCompare()) :
      BOOST_SIGNALS_NAMESPACE::detail::signal_base(real_group_compare_type(comp)),
      combiner(c)
    {
    }

    // Connect a slot to this signal
    BOOST_SIGNALS_NAMESPACE::connection connect(const slot_type&);
    BOOST_SIGNALS_NAMESPACE::connection connect(const group_type&, const slot_type&);

    // Disconnect a named slot
    void disconnect(const group_type& group)
    {
      impl->disconnect(group);
    }

    // Emit the signal
    result_type operator()(BOOST_SIGNALS_PARMS);
    result_type operator()(BOOST_SIGNALS_PARMS) const;

  private:
    Combiner combiner;
  };

  template<
    typename R,
    BOOST_SIGNALS_TEMPLATE_PARMS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    typename Combiner,
    typename Group,
    typename GroupCompare,
    typename SlotFunction
  >
  BOOST_SIGNALS_NAMESPACE::connection
  BOOST_SIGNALS_SIGNAL<
    R, BOOST_SIGNALS_TEMPLATE_ARGS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    Combiner, Group, GroupCompare, SlotFunction
  >::connect(const slot_type& in_slot)
  {
    // If the slot has been disconnected, just return a disconnected
    // connection
    if (!in_slot.is_active()) {
      return BOOST_SIGNALS_NAMESPACE::connection();
    }

    return impl->connect_slot(in_slot.get_slot_function(),
                              any(),
                              in_slot.get_bound_objects());
  }

  template<
    typename R,
    BOOST_SIGNALS_TEMPLATE_PARMS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    typename Combiner,
    typename Group,
    typename GroupCompare,
    typename SlotFunction
  >
  BOOST_SIGNALS_NAMESPACE::connection
  BOOST_SIGNALS_SIGNAL<
    R, BOOST_SIGNALS_TEMPLATE_ARGS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    Combiner, Group, GroupCompare, SlotFunction
  >::connect(const group_type& group,
             const slot_type& in_slot)
  {
    return impl->connect_slot(in_slot.get_slot_function(),
                              group,
                              in_slot.get_bound_objects());
  }

  template<
    typename R,
    BOOST_SIGNALS_TEMPLATE_PARMS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    typename Combiner,
    typename Group,
    typename GroupCompare,
    typename SlotFunction
  >
  typename BOOST_SIGNALS_SIGNAL<
             R, BOOST_SIGNALS_TEMPLATE_ARGS
             BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
             Combiner, Group, GroupCompare, SlotFunction>::result_type
  BOOST_SIGNALS_SIGNAL<
    R, BOOST_SIGNALS_TEMPLATE_ARGS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    Combiner, Group, GroupCompare, SlotFunction
  >::operator()(BOOST_SIGNALS_PARMS)
  {
    // Notify the slot handling code that we are making a call
    BOOST_SIGNALS_NAMESPACE::detail::call_notification notification(this->impl);

    // Construct a function object that will call the underlying slots
    // with the given arguments.
#if BOOST_SIGNALS_NUM_ARGS == 0
    BOOST_SIGNALS_ARGS_STRUCT_INST args;
#else
    BOOST_SIGNALS_ARGS_STRUCT_INST args(BOOST_SIGNALS_ARGS);
#endif // BOOST_SIGNALS_NUM_ARGS > 0
    call_bound_slot f(&args);

    // Let the combiner call the slots via a pair of input iterators
    return combiner(slot_call_iterator(notification.impl->slots_.begin(),
                                       impl->slots_.end(), f),
                    slot_call_iterator(notification.impl->slots_.end(),
                                       impl->slots_.end(), f));
  }

  template<
    typename R,
    BOOST_SIGNALS_TEMPLATE_PARMS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    typename Combiner,
    typename Group,
    typename GroupCompare,
    typename SlotFunction
  >
  typename BOOST_SIGNALS_SIGNAL<
             R, BOOST_SIGNALS_TEMPLATE_ARGS
             BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
             Combiner, Group, GroupCompare, SlotFunction>::result_type
  BOOST_SIGNALS_SIGNAL<
    R, BOOST_SIGNALS_TEMPLATE_ARGS
    BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS
    Combiner, Group, GroupCompare, SlotFunction
  >::operator()(BOOST_SIGNALS_PARMS) const
  {
    // Notify the slot handling code that we are making a call
    BOOST_SIGNALS_NAMESPACE::detail::call_notification notification(this->impl);

    // Construct a function object that will call the underlying slots
    // with the given arguments.
#if BOOST_SIGNALS_NUM_ARGS == 0
    BOOST_SIGNALS_ARGS_STRUCT_INST args;
#else
    BOOST_SIGNALS_ARGS_STRUCT_INST args(BOOST_SIGNALS_ARGS);
#endif // BOOST_SIGNALS_NUM_ARGS > 0

    call_bound_slot f(&args);

    // Let the combiner call the slots via a pair of input iterators
    return combiner(slot_call_iterator(notification.impl->slots_.begin(),
                                       impl->slots_.end(), f),
                    slot_call_iterator(notification.impl->slots_.end(),
                                       impl->slots_.end(), f));
  }
} // namespace boost

#undef BOOST_SIGNAL_FUNCTION_N_HEADER
#undef BOOST_SIGNALS_ARGS_STRUCT_INST
#undef BOOST_SIGNALS_CALL_BOUND
#undef BOOST_SIGNALS_ARGS_STRUCT
#undef BOOST_SIGNALS_FUNCTION
#undef BOOST_SIGNALS_SIGNAL
#undef BOOST_SIGNALS_COMMA_IF_NONZERO_ARGS

#ifdef BOOST_HAS_ABI_HEADERS
#  include BOOST_ABI_SUFFIX
#endif
