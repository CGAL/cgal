// Boost.Function library

//  Copyright Doug Gregor 2001-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

#ifndef BOOST_FUNCTION_BASE_HEADER
#define BOOST_FUNCTION_BASE_HEADER

#include <stdexcept>
#include <string>
#include <memory>
#include <new>
#include <boost/config.hpp>
#include <boost/assert.hpp>
#include <boost/type_traits/arithmetic_traits.hpp>
#include <boost/type_traits/composite_traits.hpp>
#include <boost/type_traits/is_stateless.hpp>
#include <boost/ref.hpp>
#include <boost/pending/ct_if.hpp>
#include <boost/detail/workaround.hpp>

#ifndef BOOST_NO_SFINAE
#  include "boost/utility/enable_if.hpp"
#endif

#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300 || defined(__ICL) && __ICL <= 600 || defined(__MWERKS__) && __MWERKS__ < 0x2406 && !defined(BOOST_STRICT_CONFIG)
#  define BOOST_FUNCTION_TARGET_FIX(x) x
#else
#  define BOOST_FUNCTION_TARGET_FIX(x)
#endif // not MSVC

#if defined(__sgi) && defined(_COMPILER_VERSION) && _COMPILER_VERSION <= 730 && !defined(BOOST_STRICT_CONFIG)
// Work around a compiler bug.
// boost::python::objects::function has to be seen by the compiler before the
// boost::function class template.
namespace boost { namespace python { namespace objects {
  class function;
}}}
#endif

#if defined (BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)                    \
 || defined(BOOST_BCB_PARTIAL_SPECIALIZATION_BUG)                         \
 || !(BOOST_STRICT_CONFIG || !defined(__SUNPRO_CC) || __SUNPRO_CC > 0x540)
#  define BOOST_FUNCTION_NO_FUNCTION_TYPE_SYNTAX
#endif

#if !defined(BOOST_FUNCTION_NO_FUNCTION_TYPE_SYNTAX)
namespace boost {

#if defined(__sgi) && defined(_COMPILER_VERSION) && _COMPILER_VERSION <= 730 && !defined(BOOST_STRICT_CONFIG)
// The library shipping with MIPSpro 7.3.1.3m has a broken allocator<void>
class function_base;

template<typename Signature,
         typename Allocator = std::allocator<function_base> >
class function;
#else
template<typename Signature, typename Allocator = std::allocator<void> >
class function;
#endif

template<typename Signature, typename Allocator>
inline void swap(function<Signature, Allocator>& f1,
                 function<Signature, Allocator>& f2)
{
  f1.swap(f2);
}

} // end namespace boost
#endif // have partial specialization

namespace boost {
  namespace detail {
    namespace function {
      /**
       * A union of a function pointer and a void pointer. This is necessary
       * because 5.2.10/6 allows reinterpret_cast<> to safely cast between
       * function pointer types and 5.2.9/10 allows static_cast<> to safely
       * cast between a void pointer and an object pointer. But it is not legal
       * to cast between a function pointer and a void* (in either direction),
       * so function requires a union of the two. */
      union any_pointer
      {
        void* obj_ptr;
        const void* const_obj_ptr;
        void (*func_ptr)();
        char data[1];
      };

      inline any_pointer make_any_pointer(void* o)
      {
        any_pointer p;
        p.obj_ptr = o;
        return p;
      }

      inline any_pointer make_any_pointer(const void* o)
      {
        any_pointer p;
        p.const_obj_ptr = o;
        return p;
      }

      inline any_pointer make_any_pointer(void (*f)())
      {
        any_pointer p;
        p.func_ptr = f;
        return p;
      }

      /**
       * The unusable class is a placeholder for unused function arguments
       * It is also completely unusable except that it constructable from
       * anything. This helps compilers without partial specialization to
       * handle Boost.Function objects returning void.
       */
      struct unusable
      {
        unusable() {}
        template<typename T> unusable(const T&) {}
      };

      /* Determine the return type. This supports compilers that do not support
       * void returns or partial specialization by silently changing the return
       * type to "unusable".
       */
      template<typename T> struct function_return_type { typedef T type; };

      template<>
      struct function_return_type<void>
      {
        typedef unusable type;
      };

      // The operation type to perform on the given functor/function pointer
      enum functor_manager_operation_type {
        clone_functor_tag,
        destroy_functor_tag
      };

      // Tags used to decide between different types of functions
      struct function_ptr_tag {};
      struct function_obj_tag {};
      struct member_ptr_tag {};
      struct function_obj_ref_tag {};
      struct stateless_function_obj_tag {};

      template<typename F>
      class get_function_tag
      {
        typedef typename ct_if<(is_pointer<F>::value),
                            function_ptr_tag,
                            function_obj_tag>::type ptr_or_obj_tag;

        typedef typename ct_if<(is_member_pointer<F>::value),
                            member_ptr_tag,
                            ptr_or_obj_tag>::type ptr_or_obj_or_mem_tag;

        typedef typename ct_if<(is_reference_wrapper<F>::value),
                             function_obj_ref_tag,
                             ptr_or_obj_or_mem_tag>::type or_ref_tag;

      public:
        typedef typename ct_if<(is_stateless<F>::value),
                            stateless_function_obj_tag,
                            or_ref_tag>::type type;
      };

      // The trivial manager does nothing but return the same pointer (if we
      // are cloning) or return the null pointer (if we are deleting).
      inline any_pointer trivial_manager(any_pointer f,
                                         functor_manager_operation_type op)
      {
        if (op == clone_functor_tag)
          return f;
        else
          return make_any_pointer(reinterpret_cast<void*>(0));
      }

      /**
       * The functor_manager class contains a static function "manage" which
       * can clone or destroy the given function/function object pointer.
       */
      template<typename Functor, typename Allocator>
      struct functor_manager
      {
      private:
        typedef Functor functor_type;

        // For function pointers, the manager is trivial
        static inline any_pointer
        manager(any_pointer function_ptr,
                functor_manager_operation_type op,
                function_ptr_tag)
        {
          if (op == clone_functor_tag)
            return function_ptr;
          else
            return make_any_pointer(static_cast<void (*)()>(0));
        }

        // For function object pointers, we clone the pointer to each
        // function has its own version.
        static inline any_pointer
        manager(any_pointer function_obj_ptr,
                functor_manager_operation_type op,
                function_obj_tag)
        {
#ifndef BOOST_NO_STD_ALLOCATOR
        typedef typename Allocator::template rebind<functor_type>::other
          allocator_type;
        typedef typename allocator_type::pointer pointer_type;
#else
        typedef functor_type* pointer_type;
#endif // BOOST_NO_STD_ALLOCATOR

#  ifndef BOOST_NO_STD_ALLOCATOR
          allocator_type allocator;
#  endif // BOOST_NO_STD_ALLOCATOR

          if (op == clone_functor_tag) {
            functor_type* f =
              static_cast<functor_type*>(function_obj_ptr.obj_ptr);

            // Clone the functor
#  ifndef BOOST_NO_STD_ALLOCATOR
            pointer_type copy = allocator.allocate(1);
            allocator.construct(copy, *f);

            // Get back to the original pointer type
            functor_type* new_f = static_cast<functor_type*>(copy);
#  else
            functor_type* new_f = new functor_type(*f);
#  endif // BOOST_NO_STD_ALLOCATOR
            return make_any_pointer(static_cast<void*>(new_f));
          }
          else {
            /* Cast from the void pointer to the functor pointer type */
            functor_type* f =
              reinterpret_cast<functor_type*>(function_obj_ptr.obj_ptr);

#  ifndef BOOST_NO_STD_ALLOCATOR
            /* Cast from the functor pointer type to the allocator's pointer
               type */
            pointer_type victim = static_cast<pointer_type>(f);

            // Destroy and deallocate the functor
            allocator.destroy(victim);
            allocator.deallocate(victim, 1);
#  else
            delete f;
#  endif // BOOST_NO_STD_ALLOCATOR

            return make_any_pointer(static_cast<void*>(0));
          }
        }
      public:
        /* Dispatch to an appropriate manager based on whether we have a
           function pointer or a function object pointer. */
        static any_pointer
        manage(any_pointer functor_ptr, functor_manager_operation_type op)
        {
          typedef typename get_function_tag<functor_type>::type tag_type;
          return manager(functor_ptr, op, tag_type());
        }
      };

      // A type that is only used for comparisons against zero
      struct useless_clear_type {};
    } // end namespace function
  } // end namespace detail

/**
 * The function_base class contains the basic elements needed for the
 * function1, function2, function3, etc. classes. It is common to all
 * functions (and as such can be used to tell if we have one of the
 * functionN objects).
 */
class function_base
{
public:
  function_base() : manager(0)
  {
    functor.obj_ptr = 0;
  }

  // Is this function empty?
  bool empty() const { return !manager; }

public: // should be protected, but GCC 2.95.3 will fail to allow access
  detail::function::any_pointer (*manager)(
    detail::function::any_pointer,
    detail::function::functor_manager_operation_type);
  detail::function::any_pointer functor;
};

/**
 * The bad_function_call exception class is thrown when a boost::function
 * object is invoked
 */
class bad_function_call : public std::runtime_error
{
public:
  bad_function_call() : std::runtime_error("call to empty boost::function") {}
};

/* Poison comparison between Boost.Function objects (because it is
 * meaningless). The comparisons would otherwise be allowed because of the
 * conversion required to allow syntax such as:
 *   boost::function<int, int> f;
 *   if (f) { f(5); }
 */
void operator==(const function_base&, const function_base&);
void operator!=(const function_base&, const function_base&);

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1310)
inline bool operator==(const function_base& f,
                       detail::function::useless_clear_type*)
{
  return f.empty();
}

inline bool operator!=(const function_base& f,
                       detail::function::useless_clear_type*)
{
  return !f.empty();
}

inline bool operator==(detail::function::useless_clear_type*,
                       const function_base& f)
{
  return f.empty();
}

inline bool operator!=(detail::function::useless_clear_type*,
                       const function_base& f)
{
  return !f.empty();
}
#endif

namespace detail {
  namespace function {
    inline bool has_empty_target(const function_base* f)
    {
      return f->empty();
    }

    inline bool has_empty_target(...)
    {
      return false;
    }
  } // end namespace function
} // end namespace detail
} // end namespace boost

#endif // BOOST_FUNCTION_BASE_HEADER
