// (C) Copyright David Abrahams 2002.
// (C) Copyright Jeremy Siek    2002.
// (C) Copyright Thomas Witt    2002.
// Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef BOOST_ITERATOR_FACADE_23022003THW_HPP
#define BOOST_ITERATOR_FACADE_23022003THW_HPP

#include <boost/static_assert.hpp>

#include <boost/iterator.hpp>
#include <boost/iterator/interoperable.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include <boost/iterator/detail/facade_iterator_category.hpp>
#include <boost/iterator/detail/enable_if.hpp>

#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/type_traits/add_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_pod.hpp>

#include <boost/mpl/apply_if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/always.hpp>
#include <boost/mpl/apply.hpp>

#if BOOST_WORKAROUND(BOOST_MSVC, == 1200)
# include <boost/mpl/if.hpp>
#endif

#include <boost/iterator/detail/config_def.hpp> // this goes last

namespace boost
{
  // This forward declaration is required for the friend declaration
  // in iterator_core_access
  template <class I, class V, class TC, class R, class D> class iterator_facade;

  namespace detail
  {
    // A binary metafunction class that always returns bool.  VC6
    // ICEs on mpl::always<bool>, probably because of the default
    // parameters.
    struct always_bool2
    {
        template <class T, class U>
        struct apply
        {
            typedef bool type;
        };
    };
    
    //
    // enable if for use in operator implementation.
    //
    template <
        class Facade1
      , class Facade2
      , class Return
    >
    struct enable_if_interoperable
      : ::boost::iterators::enable_if<
           mpl::or_<
               is_convertible<Facade1, Facade2>
             , is_convertible<Facade2, Facade1>
           >
         , Return
        >
    {
    };

    //
    // Generates associated types for an iterator_facade with the
    // given parameters.
    //
    template <
        class ValueParam
      , class CategoryOrTraversal
      , class Reference 
      , class Difference
    >
    struct iterator_facade_types
    {
        typedef typename facade_iterator_category<
            CategoryOrTraversal, ValueParam, Reference
        >::type iterator_category;
        
        typedef typename remove_const<ValueParam>::type value_type;
        
        typedef typename mpl::apply_if<
            detail::iterator_writability_disabled<ValueParam,Reference>
          , add_pointer<typename add_const<value_type>::type>
          , add_pointer<value_type>
        >::type pointer;
      
# if defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)                          \
    && (BOOST_WORKAROUND(_STLPORT_VERSION, BOOST_TESTED_AT(0x452))              \
        || BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, BOOST_TESTED_AT(310)))     \
    || BOOST_WORKAROUND(BOOST_RWSTD_VER, BOOST_TESTED_AT(0x20101))              \
    || BOOST_WORKAROUND(BOOST_DINKUMWARE_STDLIB, <= 310)

        // To interoperate with some broken library/compiler
        // combinations, user-defined iterators must be derived from
        // std::iterator.  It is possible to implement a standard
        // library for broken compilers without this limitation.
#  define BOOST_ITERATOR_FACADE_NEEDS_ITERATOR_BASE 1

        typedef
           iterator<iterator_category, value_type, Difference, pointer, Reference>
        base;
# endif
    };


    // operator->() needs special support for input iterators to strictly meet the
    // standard's requirements. If *i is not a reference type, we must still
    // produce a (constant) lvalue to which a pointer can be formed. We do that by
    // returning an instantiation of this special proxy class template.

    template <class T>
    struct operator_arrow_proxy
    {
        operator_arrow_proxy(T const* px) : m_value(*px) {}
        const T* operator->() const { return &m_value; }
        // This function is needed for MWCW and BCC, which won't call operator->
        // again automatically per 13.3.1.2 para 8
        operator const T*() const { return &m_value; }
        T m_value;
    };

    // A metafunction that gets the result type for operator->.  Also
    // has a static function make() which builds the result from a
    // Reference
    template <class Value, class Reference, class Pointer>
    struct operator_arrow_result
    {
        // CWPro8.3 won't accept "operator_arrow_result::type", and we
        // need that type below, so metafunction forwarding would be a
        // losing proposition here.
        typedef typename mpl::if_<
            is_reference<Reference>
          , Pointer
          , operator_arrow_proxy<Value>
        >::type type;

        static type make(Reference x)
        {
            return type(&x);
        }
    };

# if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
    // Deal with ETI
    template<>
    struct operator_arrow_result<int, int, int>
    {
        typedef int type;
    };
# endif

    //
    // Iterator is actually an iterator_facade, so we do not have to
    // go through iterator_traits to access the traits.
    //
    template <class Iterator>
    class operator_brackets_proxy
    {
        typedef typename Iterator::reference  reference;
        typedef typename Iterator::value_type value_type;

     public:
        operator_brackets_proxy(Iterator const& iter)
          : m_iter(iter)
        {}

        operator reference() const
        {
            return *m_iter;
        }

        operator_brackets_proxy& operator=(value_type const& val)
        {
            *m_iter = val;
            return *this;
        }

     private:
        Iterator m_iter;
    };

    template <class Value, class Reference>
    struct use_operator_brackets_proxy
      : mpl::and_<
            // Really we want an is_copy_constructible trait here,
            // but is_POD will have to suffice in the meantime.
            boost::is_POD<Value>
          , iterator_writability_disabled<Value,Reference>
        >
    {};
        
    template <class Iterator, class Value, class Reference>
    struct operator_brackets_result
    {
        typedef typename mpl::if_<
            use_operator_brackets_proxy<Value,Reference>
          , Value 
          , operator_brackets_proxy<Iterator>
        >::type type;
    };

    template <class Iterator>
    operator_brackets_proxy<Iterator> make_operator_brackets_result(Iterator const& iter, mpl::false_)
    {
        return operator_brackets_proxy<Iterator>(iter);
    }

    template <class Iterator>
    typename Iterator::value_type make_operator_brackets_result(Iterator const& iter, mpl::true_)
    {
      return *iter;
    }

    struct choose_difference_type
    {
        template <class I1, class I2>
        struct apply
          :
# ifdef BOOST_NO_ONE_WAY_ITERATOR_INTEROP
          iterator_difference<I1>
# elif BOOST_WORKAROUND(BOOST_MSVC, == 1200)
          mpl::if_<
              is_convertible<I2,I1>
            , typename I1::difference_type
            , typename I2::difference_type
          >
# else 
          mpl::apply_if<
              is_convertible<I2,I1>
            , iterator_difference<I1>
            , iterator_difference<I2>
          >
# endif 
        {};

    };
  } // namespace detail


  // Macros which describe the declarations of binary operators
# ifdef BOOST_NO_STRICT_ITERATOR_INTEROPERABILITY
#  define BOOST_ITERATOR_FACADE_INTEROP_HEAD(prefix, op, result_type)   \
    template <                                                          \
        class Derived1, class V1, class TC1, class R1, class D1         \
      , class Derived2, class V2, class TC2, class R2, class D2         \
    >                                                                   \
    prefix typename mpl::apply2<result_type,Derived1,Derived2>::type    \
    operator op(                                                        \
        iterator_facade<Derived1, V1, TC1, R1, D1> const& lhs           \
      , iterator_facade<Derived2, V2, TC2, R2, D2> const& rhs)
# else 
#  define BOOST_ITERATOR_FACADE_INTEROP_HEAD(prefix, op, result_type)   \
    template <                                                          \
        class Derived1, class V1, class TC1, class R1, class D1         \
      , class Derived2, class V2, class TC2, class R2, class D2         \
    >                                                                   \
    prefix typename detail::enable_if_interoperable<                    \
        Derived1, Derived2                                              \
      , typename mpl::apply2<result_type,Derived1,Derived2>::type       \
    >::type                                                             \
    operator op(                                                        \
        iterator_facade<Derived1, V1, TC1, R1, D1> const& lhs           \
      , iterator_facade<Derived2, V2, TC2, R2, D2> const& rhs)
# endif 

#  define BOOST_ITERATOR_FACADE_PLUS_HEAD(prefix,args)              \
    template <class Derived, class V, class TC, class R, class D>   \
    prefix Derived operator+ args

  //
  // Helper class for granting access to the iterator core interface.
  //
  // The simple core interface is used by iterator_facade. The core
  // interface of a user/library defined iterator type should not be made public
  // so that it does not clutter the public interface. Instead iterator_core_access
  // should be made friend so that iterator_facade can access the core
  // interface through iterator_core_access.
  //
  class iterator_core_access
  {
# if defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS)                  \
    || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
      // Tasteless as this may seem, making all members public allows member templates
      // to work in the absence of member template friends.
   public:
# else
      
      template <class I, class V, class TC, class R, class D> friend class iterator_facade;

#  define BOOST_ITERATOR_FACADE_RELATION(op)                                \
      BOOST_ITERATOR_FACADE_INTEROP_HEAD(friend,op, detail::always_bool2);

      BOOST_ITERATOR_FACADE_RELATION(==)
      BOOST_ITERATOR_FACADE_RELATION(!=)

      BOOST_ITERATOR_FACADE_RELATION(<)
      BOOST_ITERATOR_FACADE_RELATION(>)
      BOOST_ITERATOR_FACADE_RELATION(<=)
      BOOST_ITERATOR_FACADE_RELATION(>=)
#  undef BOOST_ITERATOR_FACADE_RELATION

      BOOST_ITERATOR_FACADE_INTEROP_HEAD(
          friend, -, detail::choose_difference_type)
      ;

      BOOST_ITERATOR_FACADE_PLUS_HEAD(
          friend                                
          , (iterator_facade<Derived, V, TC, R, D> const&
           , typename Derived::difference_type)
      )
      ;

      BOOST_ITERATOR_FACADE_PLUS_HEAD(
          friend
        , (typename Derived::difference_type
           , iterator_facade<Derived, V, TC, R, D> const&)
      )
      ;

# endif

      template <class Facade>
      static typename Facade::reference dereference(Facade const& f)
      {
          return f.dereference();
      }

      template <class Facade>
      static void increment(Facade& f)
      {
          f.increment();
      }

      template <class Facade>
      static void decrement(Facade& f)
      {
          f.decrement();
      }

      template <class Facade1, class Facade2>
      static bool equal(Facade1 const& f1, Facade2 const& f2, mpl::true_)
      {
          return f1.equal(f2);
      }

      template <class Facade1, class Facade2>
      static bool equal(Facade1 const& f1, Facade2 const& f2, mpl::false_)
      {
          return f2.equal(f1);
      }

      template <class Facade>
      static void advance(Facade& f, typename Facade::difference_type n)
      {
          f.advance(n);
      }

      template <class Facade1, class Facade2>
      static typename Facade1::difference_type distance_from(
          Facade1 const& f1, Facade2 const& f2, mpl::true_)
      {
          return -f1.distance_to(f2);
      }

      template <class Facade1, class Facade2>
      static typename Facade2::difference_type distance_from(
          Facade1 const& f1, Facade2 const& f2, mpl::false_)
      {
          return f2.distance_to(f1);
      }

   private:
      // objects of this class are useless
      iterator_core_access(); //undefined
  };

  //
  // iterator_facade - use as a public base class for defining new
  // standard-conforming iterators.
  //
  template <
      class Derived             // The derived iterator type being constructed
    , class Value
    , class CategoryOrTraversal
    , class Reference   = Value&
    , class Difference  = std::ptrdiff_t
  >
  class iterator_facade
# ifdef BOOST_ITERATOR_FACADE_NEEDS_ITERATOR_BASE
    : public detail::iterator_facade_types<
         Value, CategoryOrTraversal, Reference, Difference
      >::base
#  undef BOOST_ITERATOR_FACADE_NEEDS_ITERATOR_BASE
# endif
  {
   private:
      //
      // Curiously Recurring Template interface.
      //
      typedef Derived derived_t;

      Derived& derived()
      {
          return static_cast<Derived&>(*this);
      }

      Derived const& derived() const
      {
          return static_cast<Derived const&>(*this);
      }

      typedef detail::iterator_facade_types<
         Value, CategoryOrTraversal, Reference, Difference
      > associated_types;
      
   public:

      typedef typename associated_types::value_type value_type;
      typedef Reference reference;
      typedef Difference difference_type;
      typedef typename associated_types::pointer pointer;
      typedef typename associated_types::iterator_category iterator_category;

      reference operator*() const
      {
          return iterator_core_access::dereference(this->derived());
      }

      typename detail::operator_arrow_result<
          value_type
        , reference
        , pointer
      >::type
      operator->() const
      {
          return detail::operator_arrow_result<
              value_type
            , reference
            , pointer
          >::make(*this->derived());
      }
        
      typename detail::operator_brackets_result<Derived,Value,Reference>::type
      operator[](difference_type n) const
      {
          typedef detail::use_operator_brackets_proxy<Value,Reference> use_proxy;
          
          return detail::make_operator_brackets_result<Derived>(
              this->derived() + n
            , use_proxy()
          );
      }

      Derived& operator++()
      {
          iterator_core_access::increment(this->derived());
          return this->derived();
      }

      Derived operator++(int)
      {
          Derived tmp(this->derived());
          ++*this;
          return tmp;
      }

      Derived& operator--()
      {
          iterator_core_access::decrement(this->derived());
          return this->derived();
      }

      Derived operator--(int)
      {
          Derived tmp(this->derived());
          --*this;
          return tmp;
      }

      Derived& operator+=(difference_type n)
      {
          iterator_core_access::advance(this->derived(), n);
          return this->derived();
      }

      Derived& operator-=(difference_type n)
      {
          iterator_core_access::advance(this->derived(), -n);
          return this->derived();
      }

      Derived operator-(difference_type x) const
      {
          Derived result(this->derived());
          return result -= x;
      }

# if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
      // There appears to be a bug which trashes the data of classes
      // derived from iterator_facade when they are assigned unless we
      // define this assignment operator.  This bug is only revealed
      // (so far) in STLPort debug mode, but it's clearly a codegen
      // problem so we apply the workaround for all MSVC6.
      iterator_facade& operator=(iterator_facade const&)
      {
          return *this;
      }
# endif
  };

  //
  // Operator implementation. The library supplied operators
  // enables the user to provide fully interoperable constant/mutable
  // iterator types. I.e. the library provides all operators
  // for all mutable/constant iterator combinations.
  //
  // Note though that this kind of interoperability for constant/mutable
  // iterators is not required by the standard for container iterators.
  // All the standard asks for is a conversion mutable -> constant.
  // Most standard library implementations nowadays provide fully interoperable
  // iterator implementations, but there are still heavily used implementations
  // that do not provide them. (Actually it's even worse, they do not provide
  // them for only a few iterators.)
  //
  // ?? Maybe a BOOST_ITERATOR_NO_FULL_INTEROPERABILITY macro should
  //    enable the user to turn off mixed type operators
  //
  // The library takes care to provide only the right operator overloads.
  // I.e.
  //
  // bool operator==(Iterator,      Iterator);
  // bool operator==(ConstIterator, Iterator);
  // bool operator==(Iterator,      ConstIterator);
  // bool operator==(ConstIterator, ConstIterator);
  //
  //   ...
  //
  // In order to do so it uses c++ idioms that are not yet widely supported
  // by current compiler releases. The library is designed to degrade gracefully
  // in the face of compiler deficiencies. In general compiler
  // deficiencies result in less strict error checking and more obscure
  // error messages, functionality is not affected.
  //
  // For full operation compiler support for "Substitution Failure Is Not An Error"
  // (aka. enable_if) and boost::is_convertible is required.
  //
  // The following problems occur if support is lacking.
  //
  // Pseudo code
  //
  // ---------------
  // AdaptorA<Iterator1> a1;
  // AdaptorA<Iterator2> a2;
  //
  // // This will result in a no such overload error in full operation
  // // If enable_if or is_convertible is not supported
  // // The instantiation will fail with an error hopefully indicating that
  // // there is no operator== for Iterator1, Iterator2
  // // The same will happen if no enable_if is used to remove
  // // false overloads from the templated conversion constructor
  // // of AdaptorA.
  //
  // a1 == a2;
  // ----------------
  //
  // AdaptorA<Iterator> a;
  // AdaptorB<Iterator> b;
  //
  // // This will result in a no such overload error in full operation
  // // If enable_if is not supported the static assert used
  // // in the operator implementation will fail.
  // // This will accidently work if is_convertible is not supported.
  //
  // a == b;
  // ----------------
  //

# ifdef BOOST_NO_ONE_WAY_ITERATOR_INTEROP
#  define BOOST_ITERATOR_CONVERTIBLE(a,b) mpl::true_()
# else
#  define BOOST_ITERATOR_CONVERTIBLE(a,b) is_convertible<a,b>()
# endif

# define BOOST_ITERATOR_FACADE_INTEROP(op, result_type, return_prefix, base_op) \
  BOOST_ITERATOR_FACADE_INTEROP_HEAD(inline, op, result_type)                   \
  {                                                                             \
      /* For those compilers that do not support enable_if */                   \
      BOOST_STATIC_ASSERT((                                                     \
          is_interoperable< Derived1, Derived2 >::value                         \
      ));                                                                       \
      return_prefix iterator_core_access::base_op(                              \
          static_cast<Derived1 const&>(lhs)                                     \
        , static_cast<Derived2 const&>(rhs)                                     \
        , BOOST_ITERATOR_CONVERTIBLE(Derived2,Derived1)                         \
      );                                                                        \
  }

# define BOOST_ITERATOR_FACADE_RELATION(op, return_prefix, base_op) \
  BOOST_ITERATOR_FACADE_INTEROP(                                    \
      op                                                            \
    , detail::always_bool2                                          \
    , return_prefix                                                 \
    , base_op                                                       \
  )

  BOOST_ITERATOR_FACADE_RELATION(==, return, equal)
  BOOST_ITERATOR_FACADE_RELATION(!=, return !, equal)

  BOOST_ITERATOR_FACADE_RELATION(<, return 0 >, distance_from)
  BOOST_ITERATOR_FACADE_RELATION(>, return 0 <, distance_from)
  BOOST_ITERATOR_FACADE_RELATION(<=, return 0 >=, distance_from)
  BOOST_ITERATOR_FACADE_RELATION(>=, return 0 <=, distance_from)
# undef BOOST_ITERATOR_FACADE_RELATION

  // operator- requires an additional part in the static assertion
  BOOST_ITERATOR_FACADE_INTEROP(
      -
    , detail::choose_difference_type
    , return
    , distance_from
  )
# undef BOOST_ITERATOR_FACADE_INTEROP
# undef BOOST_ITERATOR_FACADE_INTEROP_HEAD

# define BOOST_ITERATOR_FACADE_PLUS(args)           \
  BOOST_ITERATOR_FACADE_PLUS_HEAD(inline, args)     \
  {                                                 \
      Derived tmp(static_cast<Derived const&>(i));  \
      return tmp += n;                              \
  }

BOOST_ITERATOR_FACADE_PLUS((
  iterator_facade<Derived, V, TC, R, D> const& i
  , typename Derived::difference_type n
))

BOOST_ITERATOR_FACADE_PLUS((
    typename Derived::difference_type n
    , iterator_facade<Derived, V, TC, R, D> const& i
))
# undef BOOST_ITERATOR_FACADE_PLUS
# undef BOOST_ITERATOR_FACADE_PLUS_HEAD

} // namespace boost

#include <boost/iterator/detail/config_undef.hpp>

#endif // BOOST_ITERATOR_FACADE_23022003THW_HPP
