// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef SELECT_HOLDER_DWA2002322_HPP
# define SELECT_HOLDER_DWA2002322_HPP

# include <boost/python/has_back_reference.hpp>
# include <boost/python/detail/not_specified.hpp>
# include <boost/python/pointee.hpp>
# include <boost/python/object/value_holder.hpp>
# include <boost/python/object/pointer_holder.hpp>
# include <boost/python/object/class_wrapper.hpp>
# include <boost/python/object/make_ptr_instance.hpp>
# include <boost/python/object/instance.hpp>
# include <boost/python/detail/force_instantiate.hpp>

# include <boost/type.hpp>

# include <boost/mpl/bool.hpp>
# include <boost/mpl/if.hpp>
# include <boost/mpl/or.hpp>
# include <boost/mpl/not.hpp>

# include <boost/type_traits/same_traits.hpp>
# include <boost/type_traits/is_base_and_derived.hpp>
# include <boost/type_traits/alignment_traits.hpp>

# include <cstddef>

namespace boost { namespace python { namespace objects {

namespace detail
{

  // check_default_constructible --
  //
  // Used to give a clean error message when the user doesn't specify
  // any __init__ functions, but when the class being wrapped doesn't
  // have an appropriate default constructor (or if
  // has_back_reference<T> is true, a constructor taking PyObject*).
  
  // A helpful compile-time assertion which gives a reasonable error
  // message if T can't be default-constructed.
  template <class T>
  static int specify_init_arguments_or_no_init_for_class_(T const&);

  // U is expected to take an initial hidden PyObject* in its
  // constructor. Normally this means U is a virtual function
  // dispatcher subclass for T.
  template <class T, class U>
  void check_default_constructible(T*, U*, mpl::true_)
  {
      python::detail::force_instantiate(
          sizeof(specify_init_arguments_or_no_init_for_class_<T>(U((::PyObject*)0)))
          );
  }

  // Handles the "normal" case where T is held directly and
  // has_back_reference<T> is not specialized.
  template <class T>
  void check_default_constructible(T*, T*, mpl::false_)
  {
      python::detail::force_instantiate(
          sizeof(specify_init_arguments_or_no_init_for_class_<T>(T()))
          );
  }

  //
  // select_value_holder/select_pointer_holder --
  //
  //   An instantiation of one of these data-free class templates is
  //   returned by select_holder::execute(), below.  Each provides the
  //   following public interface:
  //
  //    static void assert_default_constructible() -- called when no
  //    init<...> arguments are specified in class_<T, ...>'s
  //    constructor; causes a compile-time error when T has no
  //    corresponding default constructor.
  //
  //    typedef ... type -- the class derived from instance_holder
  //    which will manage a Held object in Python class instances
  //
  //    static type* get() { return 0; } -- just a way to access the
  //    computed type at runtime.
  //
  //    static void register_() -- forces registration of any
  //    to_python converters corresponding to Held.
  
  template <class T, class Held>
  struct select_value_holder
  {
   private:
      typedef mpl::or_<
          mpl::not_<is_same<T,Held> >
        , has_back_reference<T>
      > use_back_ref;
  
   public:
      static void assert_default_constructible()
      {
          detail::check_default_constructible((T*)0,(Held*)0, use_back_ref());
      }
  
      typedef typename mpl::if_<
          use_back_ref
        , value_holder_back_reference<T,Held>
        , value_holder<T>
      >::type type;

      typedef Held held_type;
      
      static inline void register_() {}

      static type* get() { return 0; }
  };

  template <class T,class Ptr>
  struct select_pointer_holder
  {
   private:
      typedef typename python::pointee<Ptr>::type wrapper;
      typedef mpl::or_<
          mpl::not_<is_same<T,wrapper> >
        , has_back_reference<T>
      > use_back_ref;

   public:
      static void assert_default_constructible()
      {
          detail::check_default_constructible((T*)0,(wrapper*)0, use_back_ref());
      }
  
      typedef typename mpl::if_<
            use_back_ref
          , pointer_holder_back_reference<Ptr,T>
          , pointer_holder<Ptr,T>
      >::type type;
      
      typedef typename pointee<Ptr>::type held_type;
      
      static inline void register_()
      {
          select_pointer_holder::register_(use_back_ref());
      }

      static type* get() { return 0; }
      
   private:
      static inline void register_(mpl::true_)
      {
      }

      struct construct_from_pointer
      {
          static type* execute(PyObject*, Ptr x)
          {
              return new type(x);
          }
      };
      
      static inline void register_(mpl::false_)
      {
          python::detail::force_instantiate(
              objects::class_value_wrapper<Ptr, make_ptr_instance<T,type> >());
      }
  };
}

//    select_holder<T,Held>::execute((Held*)0)
//
// Returns an instantiation of
// detail::select_value_holder or detail::select_pointer_holder, as
// appropriate for class_<T,Held>
template <class T, class Held>
struct select_holder
  : mpl::if_<
        is_same<Held, python::detail::not_specified>
      , detail::select_value_holder<T,T>
      , typename mpl::if_<
            mpl::or_<
                is_same<T,Held>
              , is_base_and_derived<T, Held>
            >
          , detail::select_value_holder<T,Held>
          , detail::select_pointer_holder<T, Held>
        >::type
    >::type
{
};

}}} // namespace boost::python::objects

#endif // SELECT_HOLDER_DWA2002322_HPP
