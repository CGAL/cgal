// Boost.Function library

//  Copyright Douglas Gregor 2001-2003. Use, modification and
//  distribution is subject to the Boost Software License, Version
//  1.0. (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

// For more information, see http://www.boost.org

// Note: this header is a header template and must NOT have multiple-inclusion
// protection.
#include <boost/function/detail/prologue.hpp>

#define BOOST_FUNCTION_TEMPLATE_PARMS BOOST_PP_ENUM_PARAMS(BOOST_FUNCTION_NUM_ARGS, typename T)

#define BOOST_FUNCTION_TEMPLATE_ARGS BOOST_PP_ENUM_PARAMS(BOOST_FUNCTION_NUM_ARGS, T)

#define BOOST_FUNCTION_PARM(J,I,D) BOOST_PP_CAT(T,I) BOOST_PP_CAT(a,I)

#define BOOST_FUNCTION_PARMS BOOST_PP_ENUM(BOOST_FUNCTION_NUM_ARGS,BOOST_FUNCTION_PARM,BOOST_PP_EMPTY)

#define BOOST_FUNCTION_ARGS BOOST_PP_ENUM_PARAMS(BOOST_FUNCTION_NUM_ARGS, a)

#define BOOST_FUNCTION_ARG_TYPE(J,I,D) \
  typedef BOOST_PP_CAT(T,I) BOOST_PP_CAT(arg, BOOST_PP_CAT(BOOST_PP_INC(I),_type));

#define BOOST_FUNCTION_ARG_TYPES BOOST_PP_REPEAT(BOOST_FUNCTION_NUM_ARGS,BOOST_FUNCTION_ARG_TYPE,BOOST_PP_EMPTY)

// Type of the default allocator
#ifndef BOOST_NO_STD_ALLOCATOR
#  define BOOST_FUNCTION_DEFAULT_ALLOCATOR std::allocator<function_base>
#else
#  define BOOST_FUNCTION_DEFAULT_ALLOCATOR int
#endif // BOOST_NO_STD_ALLOCATOR

// Comma if nonzero number of arguments
#if BOOST_FUNCTION_NUM_ARGS == 0
#  define BOOST_FUNCTION_COMMA
#else
#  define BOOST_FUNCTION_COMMA ,
#endif // BOOST_FUNCTION_NUM_ARGS > 0

// Class names used in this version of the code
#define BOOST_FUNCTION_FUNCTION BOOST_JOIN(function,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_FUNCTION_INVOKER \
  BOOST_JOIN(function_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_VOID_FUNCTION_INVOKER \
  BOOST_JOIN(void_function_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_VOID_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(void_function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_STATELESS_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(stateless_function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_STATELESS_VOID_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(stateless_void_function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_GET_FUNCTION_INVOKER \
  BOOST_JOIN(get_function_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_GET_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(get_function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)
#define BOOST_FUNCTION_GET_STATELESS_FUNCTION_OBJ_INVOKER \
  BOOST_JOIN(get_stateless_function_obj_invoker,BOOST_FUNCTION_NUM_ARGS)

#ifndef BOOST_NO_VOID_RETURNS
#  define BOOST_FUNCTION_VOID_RETURN_TYPE void
#  define BOOST_FUNCTION_RETURN(X) X
#else
#  define BOOST_FUNCTION_VOID_RETURN_TYPE ::boost::detail::function::unusable
#  define BOOST_FUNCTION_RETURN(X) X; return BOOST_FUNCTION_VOID_RETURN_TYPE ()
#endif

namespace boost {
  namespace detail {
    namespace function {
      template<
        typename FunctionPtr,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
        >
      struct BOOST_FUNCTION_FUNCTION_INVOKER
      {
        static R invoke(any_pointer function_ptr BOOST_FUNCTION_COMMA
                        BOOST_FUNCTION_PARMS)
        {
          FunctionPtr f = reinterpret_cast<FunctionPtr>(function_ptr.func_ptr);
          return f(BOOST_FUNCTION_ARGS);
        }
      };

      template<
        typename FunctionPtr,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
        >
      struct BOOST_FUNCTION_VOID_FUNCTION_INVOKER
      {
        static BOOST_FUNCTION_VOID_RETURN_TYPE
        invoke(any_pointer function_ptr BOOST_FUNCTION_COMMA
               BOOST_FUNCTION_PARMS)

        {
          FunctionPtr f = reinterpret_cast<FunctionPtr>(function_ptr.func_ptr);
          BOOST_FUNCTION_RETURN(f(BOOST_FUNCTION_ARGS));
        }
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
      >
      struct BOOST_FUNCTION_FUNCTION_OBJ_INVOKER
      {
        static R invoke(any_pointer function_obj_ptr BOOST_FUNCTION_COMMA
                        BOOST_FUNCTION_PARMS)

        {
          FunctionObj* f = (FunctionObj*)(function_obj_ptr.obj_ptr);
          return (*f)(BOOST_FUNCTION_ARGS);
        }
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
      >
      struct BOOST_FUNCTION_VOID_FUNCTION_OBJ_INVOKER
      {
        static BOOST_FUNCTION_VOID_RETURN_TYPE
        invoke(any_pointer function_obj_ptr BOOST_FUNCTION_COMMA
               BOOST_FUNCTION_PARMS)

        {
          FunctionObj* f = (FunctionObj*)(function_obj_ptr.obj_ptr);
          BOOST_FUNCTION_RETURN((*f)(BOOST_FUNCTION_ARGS));
        }
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
      >
      struct BOOST_FUNCTION_STATELESS_FUNCTION_OBJ_INVOKER
      {
        static R invoke(any_pointer BOOST_FUNCTION_COMMA BOOST_FUNCTION_PARMS)
        {
          FunctionObj f = FunctionObj();
          return f(BOOST_FUNCTION_ARGS);
        }
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
      >
      struct BOOST_FUNCTION_STATELESS_VOID_FUNCTION_OBJ_INVOKER
      {
        static BOOST_FUNCTION_VOID_RETURN_TYPE
        invoke(any_pointer BOOST_FUNCTION_COMMA BOOST_FUNCTION_PARMS)

        {
          FunctionObj f = FunctionObj();
          BOOST_FUNCTION_RETURN(f(BOOST_FUNCTION_ARGS));
        }
      };

      template<
        typename FunctionPtr,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
      >
      struct BOOST_FUNCTION_GET_FUNCTION_INVOKER
      {
        typedef typename ct_if<(is_void<R>::value),
                            BOOST_FUNCTION_VOID_FUNCTION_INVOKER<
                            FunctionPtr,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >,
                          BOOST_FUNCTION_FUNCTION_INVOKER<
                            FunctionPtr,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >
                       >::type type;
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
       >
      struct BOOST_FUNCTION_GET_FUNCTION_OBJ_INVOKER
      {
        typedef typename ct_if<(is_void<R>::value),
                            BOOST_FUNCTION_VOID_FUNCTION_OBJ_INVOKER<
                            FunctionObj,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >,
                          BOOST_FUNCTION_FUNCTION_OBJ_INVOKER<
                            FunctionObj,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >
                       >::type type;
      };

      template<
        typename FunctionObj,
        typename R BOOST_FUNCTION_COMMA
        BOOST_FUNCTION_TEMPLATE_PARMS
       >
      struct BOOST_FUNCTION_GET_STATELESS_FUNCTION_OBJ_INVOKER
      {
        typedef typename ct_if<(is_void<R>::value),
                            BOOST_FUNCTION_STATELESS_VOID_FUNCTION_OBJ_INVOKER<
                            FunctionObj,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >,
                          BOOST_FUNCTION_STATELESS_FUNCTION_OBJ_INVOKER<
                            FunctionObj,
                            R BOOST_FUNCTION_COMMA
                            BOOST_FUNCTION_TEMPLATE_ARGS
                          >
                       >::type type;
      };

    } // end namespace function
  } // end namespace detail

  template<
    typename R BOOST_FUNCTION_COMMA
    BOOST_FUNCTION_TEMPLATE_PARMS,
    typename Allocator = BOOST_FUNCTION_DEFAULT_ALLOCATOR
  >
  class BOOST_FUNCTION_FUNCTION : public function_base
  {
  public:
#ifndef BOOST_NO_VOID_RETURNS
    typedef R         result_type;
#else
    typedef  typename detail::function::function_return_type<R>::type
      result_type;
#endif // BOOST_NO_VOID_RETURNS

  private:
    struct clear_type {};

  public:
    BOOST_STATIC_CONSTANT(int, args = BOOST_FUNCTION_NUM_ARGS);

    // add signature for boost::lambda
    template<typename Args>
    struct sig
    {
      typedef result_type type;
    };

#if BOOST_FUNCTION_NUM_ARGS == 1
    typedef T0 argument_type;
#elif BOOST_FUNCTION_NUM_ARGS == 2
    typedef T0 first_argument_type;
    typedef T1 second_argument_type;
#endif

    BOOST_STATIC_CONSTANT(int, arity = BOOST_FUNCTION_NUM_ARGS);
    BOOST_FUNCTION_ARG_TYPES

    typedef Allocator allocator_type;
    typedef BOOST_FUNCTION_FUNCTION self_type;

    BOOST_FUNCTION_FUNCTION() : function_base()
                              , invoker(0) {}

    // MSVC chokes if the following two constructors are collapsed into
    // one with a default parameter.
    template<typename Functor>
    BOOST_FUNCTION_FUNCTION(Functor BOOST_FUNCTION_TARGET_FIX(const &) f
#ifndef BOOST_NO_SFINAE
                            ,typename enable_if_c<
                            (::boost::type_traits::ice_not<
                             (is_integral<Functor>::value)>::value),
                                        int>::type = 0
#endif // BOOST_NO_SFINAE
                            ) :
      function_base(),
      invoker(0)
    {
      this->assign_to(f);
    }

#ifndef BOOST_NO_SFINAE
    BOOST_FUNCTION_FUNCTION(clear_type*) : function_base(), invoker(0) {}
#else
    BOOST_FUNCTION_FUNCTION(int zero) : function_base(), invoker(0)
    {
      BOOST_ASSERT(zero == 0);
    }
#endif

    BOOST_FUNCTION_FUNCTION(const BOOST_FUNCTION_FUNCTION& f) :
      function_base(),
      invoker(0)
    {
      this->assign_to_own(f);
    }

    ~BOOST_FUNCTION_FUNCTION() { clear(); }

#if BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
    // MSVC 6.0 and prior require all definitions to be inline, but
    // these definitions can become very costly.
    result_type operator()(BOOST_FUNCTION_PARMS) const
    {
      if (this->empty())
        boost::throw_exception(bad_function_call());

      return invoker(this->functor BOOST_FUNCTION_COMMA BOOST_FUNCTION_ARGS);
    }
#else
    result_type operator()(BOOST_FUNCTION_PARMS) const;
#endif

    // The distinction between when to use BOOST_FUNCTION_FUNCTION and
    // when to use self_type is obnoxious. MSVC cannot handle self_type as
    // the return type of these assignment operators, but Borland C++ cannot
    // handle BOOST_FUNCTION_FUNCTION as the type of the temporary to
    // construct.
    template<typename Functor>
#ifndef BOOST_NO_SFINAE
    typename enable_if_c<
               (::boost::type_traits::ice_not<
                 (is_integral<Functor>::value)>::value),
               BOOST_FUNCTION_FUNCTION&>::type
#else
    BOOST_FUNCTION_FUNCTION&
#endif
    operator=(Functor BOOST_FUNCTION_TARGET_FIX(const &) f)
    {
      self_type(f).swap(*this);
      return *this;
    }

#ifndef BOOST_NO_SFINAE
    BOOST_FUNCTION_FUNCTION& operator=(clear_type*)
    {
      this->clear();
      return *this;
    }
#else
    BOOST_FUNCTION_FUNCTION& operator=(int zero)
    {
      BOOST_ASSERT(zero == 0);
      this->clear();
      return *this;
    }
#endif

    // Assignment from another BOOST_FUNCTION_FUNCTION
    BOOST_FUNCTION_FUNCTION& operator=(const BOOST_FUNCTION_FUNCTION& f)
    {
      if (&f == this)
        return *this;

      self_type(f).swap(*this);
      return *this;
    }

    void swap(BOOST_FUNCTION_FUNCTION& other)
    {
      if (&other == this)
        return;

      std::swap(this->manager, other.manager);
      std::swap(this->functor, other.functor);
      std::swap(invoker, other.invoker);
    }

    // Clear out a target, if there is one
    void clear()
    {
      if (this->manager) {
        function_base::functor =
          this->manager(this->functor, detail::function::destroy_functor_tag);
      }

      this->manager = 0;
      invoker = 0;
    }

#if (defined __SUNPRO_CC) && (__SUNPRO_CC <= 0x530) && !(defined BOOST_NO_COMPILER_CONFIG)
    // Sun C++ 5.3 can't handle the safe_bool idiom, so don't use it
    operator bool () const { return !this->empty(); }
#else
  private:
    struct dummy {
      void nonnull() {};
    };

    typedef void (dummy::*safe_bool)();

  public:
    operator safe_bool () const
      { return (this->empty())? 0 : &dummy::nonnull; }

    bool operator!() const
      { return this->empty(); }
#endif

  private:
    void assign_to_own(const BOOST_FUNCTION_FUNCTION& f)
    {
      if (!f.empty()) {
        invoker = f.invoker;
        this->manager = f.manager;
        this->functor =
          f.manager(f.functor, detail::function::clone_functor_tag);
      }
    }

    template<typename Functor>
    void assign_to(Functor f)
    {
      typedef typename detail::function::get_function_tag<Functor>::type tag;
      this->assign_to(f, tag());
    }

    template<typename FunctionPtr>
    void assign_to(FunctionPtr f, detail::function::function_ptr_tag)
    {
      clear();

      if (f) {
        typedef typename detail::function::BOOST_FUNCTION_GET_FUNCTION_INVOKER<
                           FunctionPtr,
                           R BOOST_FUNCTION_COMMA
                           BOOST_FUNCTION_TEMPLATE_ARGS
                         >::type
          actual_invoker_type;

        invoker = &actual_invoker_type::invoke;
        this->manager =
          &detail::function::functor_manager<FunctionPtr, Allocator>::manage;
        this->functor =
          this->manager(detail::function::make_any_pointer(
                            // should be a reinterpret cast, but some compilers
                            // insist on giving cv-qualifiers to free functions
                            (void (*)())(f)
                          ),
                          detail::function::clone_functor_tag);
      }
    }

#if BOOST_FUNCTION_NUM_ARGS > 0
    template<typename MemberPtr>
    void assign_to(MemberPtr f, detail::function::member_ptr_tag)
    {
      this->assign_to(mem_fn(f));
    }
#endif // BOOST_FUNCTION_NUM_ARGS > 0

    template<typename FunctionObj>
    void assign_to(FunctionObj f, detail::function::function_obj_tag)
    {
      if (!detail::function::has_empty_target(boost::addressof(f))) {
        typedef
          typename detail::function::BOOST_FUNCTION_GET_FUNCTION_OBJ_INVOKER<
                                       FunctionObj,
                                       R BOOST_FUNCTION_COMMA
                                       BOOST_FUNCTION_TEMPLATE_ARGS
                                     >::type
          actual_invoker_type;

        invoker = &actual_invoker_type::invoke;
        this->manager = &detail::function::functor_manager<
                                    FunctionObj, Allocator>::manage;
#ifndef BOOST_NO_STD_ALLOCATOR
        typedef typename Allocator::template rebind<FunctionObj>::other
          rebound_allocator_type;
        typedef typename rebound_allocator_type::pointer pointer_type;
        rebound_allocator_type allocator;
        pointer_type copy = allocator.allocate(1);
        allocator.construct(copy, f);

        // Get back to the original pointer type
        FunctionObj* new_f = static_cast<FunctionObj*>(copy);
#else
        FunctionObj* new_f = new FunctionObj(f);
#endif // BOOST_NO_STD_ALLOCATOR
        this->functor =
          detail::function::make_any_pointer(static_cast<void*>(new_f));
      }
    }

    template<typename FunctionObj>
    void assign_to(const reference_wrapper<FunctionObj>& f,
                   detail::function::function_obj_ref_tag)
    {
      if (!detail::function::has_empty_target(f.get_pointer())) {
        typedef
          typename detail::function::BOOST_FUNCTION_GET_FUNCTION_OBJ_INVOKER<
                                       FunctionObj,
                                       R BOOST_FUNCTION_COMMA
                                       BOOST_FUNCTION_TEMPLATE_ARGS
                                     >::type
          actual_invoker_type;

        invoker = &actual_invoker_type::invoke;
        this->manager = &detail::function::trivial_manager<FunctionObj>::get;
        this->functor =
          this->manager(
            detail::function::make_any_pointer(
              const_cast<FunctionObj*>(f.get_pointer())),
            detail::function::clone_functor_tag);
      }
    }

    template<typename FunctionObj>
    void assign_to(FunctionObj, detail::function::stateless_function_obj_tag)
    {
      typedef
          typename detail::function::
                     BOOST_FUNCTION_GET_STATELESS_FUNCTION_OBJ_INVOKER<
                       FunctionObj,
                       R BOOST_FUNCTION_COMMA
                       BOOST_FUNCTION_TEMPLATE_ARGS
                     >::type
          actual_invoker_type;
      invoker = &actual_invoker_type::invoke;
      this->manager = &detail::function::trivial_manager<FunctionObj>::get;
      this->functor = detail::function::make_any_pointer(this);
    }

    typedef result_type (*invoker_type)(detail::function::any_pointer
                                        BOOST_FUNCTION_COMMA
                                        BOOST_FUNCTION_TEMPLATE_ARGS);

    invoker_type invoker;
  };

  template<typename R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_PARMS ,
           typename Allocator>
  inline void swap(BOOST_FUNCTION_FUNCTION<
                     R BOOST_FUNCTION_COMMA
                     BOOST_FUNCTION_TEMPLATE_ARGS ,
                     Allocator
                   >& f1,
                   BOOST_FUNCTION_FUNCTION<
                     R BOOST_FUNCTION_COMMA
                     BOOST_FUNCTION_TEMPLATE_ARGS,
                     Allocator
                   >& f2)
  {
    f1.swap(f2);
  }

#if !BOOST_WORKAROUND(BOOST_MSVC, <= 1200)
  template<typename R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_PARMS,
           typename Allocator>
  typename BOOST_FUNCTION_FUNCTION<
      R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_ARGS, 
      Allocator>::result_type
   BOOST_FUNCTION_FUNCTION<R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_ARGS,

                           Allocator>
  ::operator()(BOOST_FUNCTION_PARMS) const
  {
    if (this->empty())
      boost::throw_exception(bad_function_call());
    
    return invoker(this->functor BOOST_FUNCTION_COMMA BOOST_FUNCTION_ARGS);
  }
#endif

// Poison comparisons between boost::function objects of the same type.
template<typename R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_PARMS ,
         typename Allocator>
  void operator==(const BOOST_FUNCTION_FUNCTION<
                          R BOOST_FUNCTION_COMMA
                          BOOST_FUNCTION_TEMPLATE_ARGS ,
                          Allocator>&,
                  const BOOST_FUNCTION_FUNCTION<
                          R BOOST_FUNCTION_COMMA
                          BOOST_FUNCTION_TEMPLATE_ARGS ,
                  Allocator>&);
template<typename R BOOST_FUNCTION_COMMA BOOST_FUNCTION_TEMPLATE_PARMS ,
         typename Allocator>
  void operator!=(const BOOST_FUNCTION_FUNCTION<
                          R BOOST_FUNCTION_COMMA
                          BOOST_FUNCTION_TEMPLATE_ARGS ,
                          Allocator>&,
                  const BOOST_FUNCTION_FUNCTION<
                          R BOOST_FUNCTION_COMMA
                          BOOST_FUNCTION_TEMPLATE_ARGS ,
                  Allocator>&);

#if !defined(BOOST_FUNCTION_NO_FUNCTION_TYPE_SYNTAX)

#if BOOST_FUNCTION_NUM_ARGS == 0
#define BOOST_FUNCTION_PARTIAL_SPEC R (void)
#else
#define BOOST_FUNCTION_PARTIAL_SPEC R (BOOST_PP_ENUM_PARAMS(BOOST_FUNCTION_NUM_ARGS,T))
#endif

template<typename R BOOST_FUNCTION_COMMA
         BOOST_FUNCTION_TEMPLATE_PARMS,
         typename Allocator>
class function<BOOST_FUNCTION_PARTIAL_SPEC, Allocator>
  : public BOOST_FUNCTION_FUNCTION<R, BOOST_FUNCTION_TEMPLATE_ARGS
                                   BOOST_FUNCTION_COMMA Allocator>
{
  typedef BOOST_FUNCTION_FUNCTION<R, BOOST_FUNCTION_TEMPLATE_ARGS
                                  BOOST_FUNCTION_COMMA Allocator> base_type;
  typedef function self_type;

  struct clear_type {};

public:
  typedef typename base_type::allocator_type allocator_type;

  function() : base_type() {}

  template<typename Functor>
  function(Functor f
#ifndef BOOST_NO_SFINAE
           ,typename enable_if_c<
                            (::boost::type_traits::ice_not<
                          (is_integral<Functor>::value)>::value),
                       int>::type = 0
#endif
           ) :
    base_type(f)
  {
  }

#ifndef BOOST_NO_SFINAE
  function(clear_type*) : base_type() {}
#endif

  function(const self_type& f) : base_type(static_cast<const base_type&>(f)){}

  function(const base_type& f) : base_type(static_cast<const base_type&>(f)){}

  self_type& operator=(const self_type& f)
  {
    self_type(f).swap(*this);
    return *this;
  }

  template<typename Functor>
#ifndef BOOST_NO_SFINAE
  typename enable_if_c<
                            (::boost::type_traits::ice_not<
                         (is_integral<Functor>::value)>::value),
                      self_type&>::type
#else
  self_type&
#endif
  operator=(Functor f)
  {
    self_type(f).swap(*this);
    return *this;
  }

#ifndef BOOST_NO_SFINAE
  self_type& operator=(clear_type*)
  {
    this->clear();
    return *this;
  }
#endif

  self_type& operator=(const base_type& f)
  {
    self_type(f).swap(*this);
    return *this;
  }
};

#undef BOOST_FUNCTION_PARTIAL_SPEC
#endif // have partial specialization

} // end namespace boost

// Cleanup after ourselves...
#undef BOOST_FUNCTION_DEFAULT_ALLOCATOR
#undef BOOST_FUNCTION_COMMA
#undef BOOST_FUNCTION_FUNCTION
#undef BOOST_FUNCTION_FUNCTION_INVOKER
#undef BOOST_FUNCTION_VOID_FUNCTION_INVOKER
#undef BOOST_FUNCTION_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_VOID_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_STATELESS_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_STATELESS_VOID_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_GET_FUNCTION_INVOKER
#undef BOOST_FUNCTION_GET_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_GET_STATELESS_FUNCTION_OBJ_INVOKER
#undef BOOST_FUNCTION_GET_MEM_FUNCTION_INVOKER
#undef BOOST_FUNCTION_TEMPLATE_PARMS
#undef BOOST_FUNCTION_TEMPLATE_ARGS
#undef BOOST_FUNCTION_PARMS
#undef BOOST_FUNCTION_PARM
#undef BOOST_FUNCTION_ARGS
#undef BOOST_FUNCTION_ARG_TYPE
#undef BOOST_FUNCTION_ARG_TYPES
#undef BOOST_FUNCTION_VOID_RETURN_TYPE
#undef BOOST_FUNCTION_RETURN
