// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CLASS_DWA200216_HPP
# define CLASS_DWA200216_HPP

# include <boost/python/detail/prefix.hpp>

# include <boost/python/class_fwd.hpp>
# include <boost/python/object/class.hpp>

# include <boost/python/bases.hpp>
# include <boost/python/object.hpp>
# include <boost/python/type_id.hpp>
# include <boost/python/data_members.hpp>
# include <boost/python/make_function.hpp>
# include <boost/python/signature.hpp>
# include <boost/python/init.hpp>
# include <boost/python/args_fwd.hpp>

# include <boost/python/object/select_holder.hpp>
# include <boost/python/object/class_wrapper.hpp>
# include <boost/python/object/make_instance.hpp>
# include <boost/python/object/pickle_support.hpp>
# include <boost/python/object/add_to_namespace.hpp>
# include <boost/python/object/class_converters.hpp>

# include <boost/python/detail/overloads_fwd.hpp>
# include <boost/python/detail/operator_id.hpp>
# include <boost/python/detail/def_helper.hpp>
# include <boost/python/detail/force_instantiate.hpp>

# include <boost/type_traits/is_same.hpp>
# include <boost/type_traits/is_convertible.hpp>
# include <boost/type_traits/is_member_function_pointer.hpp>
# include <boost/type_traits/is_polymorphic.hpp>

# include <boost/mpl/size.hpp>
# include <boost/mpl/for_each.hpp>
# include <boost/mpl/bool.hpp>
# include <boost/mpl/not.hpp>
# include <boost/mpl/or.hpp>
# include <boost/mpl/vector/vector10.hpp>

# include <boost/utility.hpp>
# include <boost/detail/workaround.hpp>

# if BOOST_WORKAROUND(__MWERKS__, <= 0x3004)                        \
    /* pro9 reintroduced the bug */                                 \
    || (BOOST_WORKAROUND(__MWERKS__, > 0x3100)                      \
        && BOOST_WORKAROUND(__MWERKS__, BOOST_TESTED_AT(0x3201)))   \
    || BOOST_WORKAROUND(__GNUC__, < 3)

#  define BOOST_PYTHON_NO_MEMBER_POINTER_ORDERING 1

# endif

# ifdef BOOST_PYTHON_NO_MEMBER_POINTER_ORDERING
#  include <boost/mpl/and.hpp>
#  include <boost/type_traits/is_member_pointer.hpp>
# endif

namespace boost { namespace python {

template <class DerivedVisitor> class def_visitor;

enum no_init_t { no_init };

namespace detail
{
  // This function object is used with mpl::for_each to write the id
  // of the type a pointer to which is passed as its 2nd compile-time
  // argument. into the iterator pointed to by its runtime argument
  struct write_type_id
  {
      write_type_id(type_info**p) : p(p) {}

      // Here's the runtime behavior
      template <class T>
      void operator()(T*) const
      {
          *(*p)++ = type_id<T>();
      }

      type_info** p;
  };

  template <class T, class Prev = detail::not_specified>
  struct select_held_type;

  template <class T1, class T2, class T3>
  struct has_noncopyable;

  // Register to_python converters for a class T.  The first argument
  // will be mpl::true_ unless noncopyable was specified as a
  // class_<...> template parameter. The 2nd argument is a pointer to
  // the type of holder that must be created. The 3rd argument is a
  // reference to the Python type object to be created.
  template <class T, class SelectHolder>
  inline void register_class_to_python(mpl::true_, SelectHolder, T* = 0)
  {
      typedef typename SelectHolder::type holder;
      force_instantiate(objects::class_cref_wrapper<T, objects::make_instance<T,holder> >());
      SelectHolder::register_();
  }

  template <class T, class SelectHolder>
  inline void register_class_to_python(mpl::false_, SelectHolder, T* = 0)
  {
      SelectHolder::register_();
  }

  //
  // register_wrapper_class -- register the relationship between a
  // virtual function callback wrapper class and the class being
  // wrapped.
  //
  template <class T>
  inline void register_wrapper_class_impl(T*, T*, int) {}
  
  template <class Wrapper, class T>
  inline void register_wrapper_class_impl(Wrapper*, T*, ...)
  {
      objects::register_class_from_python<Wrapper, mpl::vector1<T> >();
      objects::copy_class_object(type_id<T>(), type_id<Wrapper>());
  }
  
  template <class Held, class T>
  inline void register_wrapper_class(Held* = 0, T* = 0)
  {
      register_wrapper_class_impl((Held*)0, (T*)0,  0);
  }
  
  template <class T>
  struct is_data_member_pointer
      : mpl::and_<
            is_member_pointer<T>
          , mpl::not_<is_member_function_pointer<T> >
        >
  {};
  
# ifdef BOOST_PYTHON_NO_MEMBER_POINTER_ORDERING
#  define BOOST_PYTHON_DATA_MEMBER_HELPER(D) , detail::is_data_member_pointer<D>()
#  define BOOST_PYTHON_YES_DATA_MEMBER , mpl::true_
#  define BOOST_PYTHON_NO_DATA_MEMBER , mpl::false_
# elif defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)
#  define BOOST_PYTHON_DATA_MEMBER_HELPER(D) , 0
#  define BOOST_PYTHON_YES_DATA_MEMBER , int
#  define BOOST_PYTHON_NO_DATA_MEMBER , ...
# else 
#  define BOOST_PYTHON_DATA_MEMBER_HELPER(D)
#  define BOOST_PYTHON_YES_DATA_MEMBER
#  define BOOST_PYTHON_NO_DATA_MEMBER
# endif
  
  namespace error
  {
    //
    // A meta-assertion mechanism which prints nice error messages and
    // backtraces on lots of compilers. Usage:
    //
    //      assertion<C>::failed
    //
    // where C is an MPL metafunction class
    //
    
    template <class C> struct assertion_failed { };
    template <class C> struct assertion_ok { typedef C failed; };

    template <class C>
    struct assertion
        : mpl::if_<C, assertion_ok<C>, assertion_failed<C> >::type
    {};

    //
    // Checks for validity of arguments used to define virtual
    // functions with default implementations.
    //
    
    template <class Default>
    void not_a_derived_class_member(Default) {}
    
    template <class T, class Fn>
    struct virtual_function_default
    {
        template <class Default>
        static void
        must_be_derived_class_member(Default const&)
        {
            typedef typename assertion<mpl::not_<is_same<Default,Fn> > >::failed test0;
# if !BOOST_WORKAROUND(__MWERKS__, <= 0x2407)
            typedef typename assertion<is_polymorphic<T> >::failed test1;
# endif 
            typedef typename assertion<is_member_function_pointer<Fn> >::failed test2;
            not_a_derived_class_member<Default>(Fn());
        }
    };
  }
}

// This is the primary mechanism through which users will expose
// C++ classes to Python.
template <
    class T // class being wrapped
    , class X1 // = detail::not_specified
    , class X2 // = detail::not_specified
    , class X3 // = detail::not_specified
    >
class class_ : public objects::class_base
{
 public: // types
    typedef objects::class_base base;
    typedef T wrapped_type;
    
    typedef class_<T,X1,X2,X3> self;
    BOOST_STATIC_CONSTANT(bool, is_copyable = (!detail::has_noncopyable<X1,X2,X3>::value));

    // held_type - either T, a class derived from T or a smart pointer
    // to a (class derived from) T.    
    typedef typename detail::select_held_type<
        X1, typename detail::select_held_type<
        X2, typename detail::select_held_type<
        X3
    >::type>::type>::type held_type;

    typedef objects::select_holder<T,held_type> select_holder;
    
 private: // types

    typedef typename detail::select_bases<X1
            , typename detail::select_bases<X2
              , typename boost::python::detail::select_bases<X3>::type
              >::type
            >::type bases;


    // A helper class which will contain an array of id objects to be
    // passed to the base class constructor
    struct id_vector
    {
        id_vector()
        {
            // Stick the derived class id into the first element of the array
            ids[0] = type_id<T>();

            // Write the rest of the elements into succeeding positions.
            type_info* p = ids + 1;
            mpl::for_each(detail::write_type_id(&p), (bases*)0, (add_pointer<mpl::_>*)0);
        }

        BOOST_STATIC_CONSTANT(
            std::size_t, size = mpl::size<bases>::value + 1);
        type_info ids[size];
    };
    friend struct id_vector;

 public: // constructors
    
    // Construct with the class name, with or without docstring, and default __init__() function
    class_(char const* name, char const* doc = 0);

    // Construct with class name, no docstring, and an uncallable __init__ function
    class_(char const* name, no_init_t);

    // Construct with class name, docstring, and an uncallable __init__ function
    class_(char const* name, char const* doc, no_init_t);

    // Construct with class name and init<> function
    template <class DerivedT>
    inline class_(char const* name, init_base<DerivedT> const& i)
        : base(name, id_vector::size, id_vector().ids)
    {
        this->register_holder();
        this->def(i);
    }

    // Construct with class name, docstring and init<> function
    template <class DerivedT>
    inline class_(char const* name, char const* doc, init_base<DerivedT> const& i)
        : base(name, id_vector::size, id_vector().ids, doc)
    {
        this->register_holder();
        this->def(i);
    }

 public: // member functions
    
    // Generic visitation
    template <class Derived>
    self& def(def_visitor<Derived> const& visitor)
    {
        visitor.visit(*this);
        return *this;
    }

    // Wrap a member function or a non-member function which can take
    // a T, T cv&, or T cv* as its first parameter, a callable
    // python object, or a generic visitor.
    template <class F>
    self& def(char const* name, F f)
    {
        this->def_impl(name, f, detail::def_helper<char const*>(0), &f);
        return *this;
    }

    template <class A1, class A2>
    self& def(char const* name, A1 a1, A2 const& a2)
    {
        this->def_maybe_overloads(name, a1, a2, &a2);
        return *this;
    }

    template <class Fn, class A1, class A2>
    self& def(char const* name, Fn fn, A1 const& a1, A2 const& a2)
    {
        //  The arguments are definitely:
        //      def(name, function, policy, doc_string)
        //      def(name, function, doc_string, policy)

        this->def_impl(
            name, fn
            , detail::def_helper<A1,A2>(a1,a2)
            , &fn);

        return *this;
    }

    template <class Fn, class A1, class A2, class A3>
    self& def(char const* name, Fn fn, A1 const& a1, A2 const& a2, A3 const& a3)
    {
        this->def_impl(
            name, fn
            , detail::def_helper<A1,A2,A3>(a1,a2,a3)
            , &fn);

        return *this;
    }

    //
    // Data member access
    //
    template <class D>
    self& def_readonly(char const* name, D const& d)
    {
        return this->def_readonly_impl(name, d BOOST_PYTHON_DATA_MEMBER_HELPER(D));
    }

    template <class D>
    self& def_readwrite(char const* name, D const& d)
    {
        return this->def_readwrite_impl(name, d BOOST_PYTHON_DATA_MEMBER_HELPER(D));
    }
    
    template <class D>
    self& def_readonly(char const* name, D& d)
    {
        return this->def_readonly_impl(name, d BOOST_PYTHON_DATA_MEMBER_HELPER(D));
    }

    template <class D>
    self& def_readwrite(char const* name, D& d)
    {
        return this->def_readwrite_impl(name, d BOOST_PYTHON_DATA_MEMBER_HELPER(D));
    }

    // Property creation
    template <class Get>
    self& add_property(char const* name, Get fget)
    {
        base::add_property(name, this->make_getter(fget));
        return *this;
    }

    template <class Get, class Set>
    self& add_property(char const* name, Get fget, Set fset)
    {
        base::add_property(name, this->make_getter(fget), this->make_setter(fset));
        return *this;
    }
        
    template <class Get>
    self& add_static_property(char const* name, Get fget)
    {
        base::add_static_property(name, object(fget));
        return *this;
    }

    template <class Get, class Set>
    self& add_static_property(char const* name, Get fget, Set fset)
    {
        base::add_static_property(name, object(fget), object(fset));
        return *this;
    }
        
    template <class U>
    self& setattr(char const* name, U const& x)
    {
        this->base::setattr(name, object(x));
        return *this;
    }

    // Pickle support
    template <typename PickleSuiteType>
    self& def_pickle(PickleSuiteType const& x)
    {
      error_messages::must_be_derived_from_pickle_suite(x);
      detail::pickle_suite_finalize<PickleSuiteType>::register_(
        *this,
        &PickleSuiteType::getinitargs,
        &PickleSuiteType::getstate,
        &PickleSuiteType::setstate,
        PickleSuiteType::getstate_manages_dict());
      return *this;
    }

    self& staticmethod(char const* name)
    {
        this->make_method_static(name);
        return *this;
    }
 private: // helper functions

    // Builds a method for this class around the given [member]
    // function pointer or object, appropriately adjusting the type of
    // the first signature argument so that if f is a member of a
    // (possibly not wrapped) base class of T, an lvalue argument of
    // type T will be required.
    //
    // @group PropertyHelpers {
    template <class F>
    object make_getter(F f)
    {
        typedef typename api::is_object_operators<F>::type is_obj_or_proxy;
        
        return this->make_fn_impl(
            f, is_obj_or_proxy(), (char*)0, detail::is_data_member_pointer<F>()
        );
    }
    
    template <class F>
    object make_setter(F f)
    {
        typedef typename api::is_object_operators<F>::type is_obj_or_proxy;
        
        return this->make_fn_impl(
            f, is_obj_or_proxy(), (int*)0, detail::is_data_member_pointer<F>()
        );
    }
    
    template <class F>
    object make_fn_impl(F const& f, mpl::false_, void*, mpl::false_)
    {
        return python::make_function(f, default_call_policies(), detail::get_signature(f, (T*)0));
    }

    template <class D, class B>
    object make_fn_impl(D B::*pm_, mpl::false_, char*, mpl::true_)
    {
        D T::*pm = pm_;
        return python::make_getter(pm);
    }

    template <class D, class B>
    object make_fn_impl(D B::*pm_, mpl::false_, int*, mpl::true_)
    {
        D T::*pm = pm_;
        return python::make_setter(pm);
    }

    template <class F>
    object make_fn_impl(F const& x, mpl::true_, void*, mpl::false_)
    {
        return x;
    }
    // }
    
    template <class D, class B>
    self& def_readonly_impl(
        char const* name, D B::*pm_ BOOST_PYTHON_YES_DATA_MEMBER)
    {
        return this->add_property(name, pm_);
    }

    template <class D, class B>
    self& def_readwrite_impl(
        char const* name, D B::*pm_ BOOST_PYTHON_YES_DATA_MEMBER)
    {
        return this->add_property(name, pm_, pm_);
    }

    template <class D>
    self& def_readonly_impl(
        char const* name, D& d BOOST_PYTHON_NO_DATA_MEMBER)
    {
        return this->add_static_property(name, python::make_getter(d));
    }

    template <class D>
    self& def_readwrite_impl(
        char const* name, D& d  BOOST_PYTHON_NO_DATA_MEMBER)
    {
        return this->add_static_property(name, python::make_getter(d), python::make_setter(d));
    }

    inline void register_() const;
    inline void register_holder();
    
    //
    // These two overloads discriminate between def() as applied to a
    // generic visitor and everything else.
    //
    // @group def_impl {
    template <class Helper, class LeafVisitor, class Visitor>
    inline void def_impl(
        char const* name
      , LeafVisitor
      , Helper const& helper
      , def_visitor<Visitor> const* v
    )
    {
        v->visit(*this, name,  helper);
    }

    template <class Fn, class Helper>
    inline void def_impl(
        char const* name
      , Fn fn
      , Helper const& helper
      , ...
    )
    {
        objects::add_to_namespace(
            *this
          , name
          , make_function(
                fn
              , helper.policies()
              , helper.keywords()
              , detail::get_signature(fn, (T*)0)
            )
          , helper.doc()
        );

        this->def_default(name, fn, helper, mpl::bool_<Helper::has_default_implementation>());
    }
    // }

    //
    // These two overloads handle the definition of default
    // implementation overloads for virtual functions. The second one
    // handles the case where no default implementation was specified.
    //
    // @group def_default {
    template <class Fn, class Helper>
    inline void def_default(
        char const* name
        , Fn fn
        , Helper const& helper
        , mpl::bool_<true>)
    {
        detail::error::virtual_function_default<T,Fn>::must_be_derived_class_member(
            helper.default_implementation());
            
        objects::add_to_namespace(
            *this, name,
            make_function(
                helper.default_implementation(), helper.policies(), helper.keywords())
            );
    }
    
    template <class Fn, class Helper>
    inline void def_default(char const*, Fn, Helper const&, mpl::bool_<false>)
    { }
    // }
    
    //
    // These two overloads discriminate between def() as applied to
    // regular functions and def() as applied to the result of
    // BOOST_PYTHON_FUNCTION_OVERLOADS(). The final argument is used to
    // discriminate.
    //
    // @group def_maybe_overloads {
    template <class OverloadsT, class SigT>
    void def_maybe_overloads(
        char const* name
        , SigT sig
        , OverloadsT const& overloads
        , detail::overloads_base const*)

    {
        //  convert sig to a type_list (see detail::get_signature in signature.hpp)
        //  before calling detail::define_with_defaults.
        detail::define_with_defaults(
            name, overloads, *this, detail::get_signature(sig));
    }

    template <class Fn, class A1>
    void def_maybe_overloads(
        char const* name
        , Fn fn
        , A1 const& a1
        , ...)
    {
        this->def_impl(
            name
          , fn
          , detail::def_helper<A1>(a1)
          , &fn
        );

    }
    // }
};


//
// implementations
//

// register converters
template <class T, class X1, class X2, class X3>
inline void class_<T,X1,X2,X3>::register_() const
{
    objects::register_class_from_python<T,bases>();

    typedef BOOST_DEDUCED_TYPENAME select_holder::held_type held_t;
    detail::register_wrapper_class<held_t,T>();
    
    detail::register_class_to_python<T>(
        mpl::bool_<is_copyable>()
      , select_holder()
    );
}

template <class T, class X1, class X2, class X3>
inline void class_<T,X1,X2,X3>::register_holder()
{
    this->register_();
    typedef typename select_holder::type holder;
    this->set_instance_size(
        objects::additional_instance_size<holder>::value
    );
}

template <class T, class X1, class X2, class X3>
inline class_<T,X1,X2,X3>::class_(char const* name, char const* doc)
    : base(name, id_vector::size, id_vector().ids, doc)
{
    this->register_holder();
    select_holder::assert_default_constructible();
    this->def(init<>());
}

template <class T, class X1, class X2, class X3>
inline class_<T,X1,X2,X3>::class_(char const* name, no_init_t)
    : base(name, id_vector::size, id_vector().ids)
{
    this->register_();
    this->def_no_init();
}

template <class T, class X1, class X2, class X3>
inline class_<T,X1,X2,X3>::class_(char const* name, char const* doc, no_init_t)
    : base(name, id_vector::size, id_vector().ids, doc)
{
    this->register_();
    this->def_no_init();
}

namespace detail
{
  template <class T1, class T2, class T3>
  struct has_noncopyable
      : mpl::or_<
          is_same<T1,noncopyable>
        , is_same<T2,noncopyable>
        , is_same<T3,noncopyable>
      >
  {};


  template <class T, class Prev>
  struct select_held_type
    : mpl::if_<
          mpl::or_<
              specifies_bases<T>
            , is_same<T,noncopyable>
              >
        , Prev
        , T
      >
  {
  };
}

}} // namespace boost::python

# undef BOOST_PYTHON_DATA_MEMBER_HELPER
# undef BOOST_PYTHON_YES_DATA_MEMBER
# undef BOOST_PYTHON_NO_DATA_MEMBER
# undef BOOST_PYTHON_NO_MEMBER_POINTER_ORDERING

#endif // CLASS_DWA200216_HPP
