// Copyright David Abrahams 2002. Permission to copy, use,
// modify, sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.
#ifndef CLASS_CONVERTERS_DWA2002119_HPP
# define CLASS_CONVERTERS_DWA2002119_HPP

# include <boost/python/converter/registry.hpp>
# include <boost/python/converter/shared_ptr_from_python.hpp>

# include <boost/python/object/inheritance.hpp>

# include <boost/python/detail/force_instantiate.hpp>

# include <boost/type_traits/add_pointer.hpp>
# include <boost/type_traits/is_polymorphic.hpp>

# include <boost/mpl/for_each.hpp>

# include <boost/detail/workaround.hpp>

namespace boost { namespace python { namespace objects { 

//////////////////////////////////////////////////////////////////////
//
// register_base_of<T> -
//      A BinaryMetaFunction object which registers a single base
//      class of T, and the corresponding cast(s)
//


// register_downcast/do_nothing -
//      Helpers for register_base_of<> which take care of registering
//      down-casts
template <class Base, class Derived>
struct register_downcast
{
    static void execute()
    {
        register_conversion<Base, Derived>(true);
    }
};

struct do_nothing
{
    static void execute() { }
};

// Here's where the real work gets done:
template <class Derived>
struct register_base_of
{
    // Here's the runtime part:
    template <class Base>
    void operator()(Base*) const
    {
        // Register the Base class
        register_dynamic_id<Base>();
        // Register the up-cast
        register_conversion<Derived,Base>(false);

        // Register the down-cast, if appropriate.
        mpl::if_<
# if BOOST_WORKAROUND(__MWERKS__, <= 0x2407)
            mpl::true_
# else
            is_polymorphic<Base>
# endif 
          , register_downcast<Base,Derived>
          , do_nothing
        >::type::execute();
    }
};

// Brings into existence all converters associated with a class. Bases
// is expected to be an mpl sequence of base types.
template <class Derived, class Bases>
inline void register_class_from_python(Derived* = 0, Bases* = 0)
{
    // Static object constructor performs registration
    static converter::shared_ptr_from_python<Derived> shared_ptr_registration;

    // register all up/downcasts here
    register_dynamic_id<Derived>();

    // register each base in the sequence
    mpl::for_each(register_base_of<Derived>(), (Bases*)0, (add_pointer<mpl::_>*)0);
}

}}} // namespace boost::python::object

#endif // CLASS_CONVERTERS_DWA2002119_HPP
