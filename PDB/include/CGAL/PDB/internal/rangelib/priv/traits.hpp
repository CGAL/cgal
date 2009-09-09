// Boost. Iterable Range Library (rangelib)
//
// Copyright 2003-2004 John Torjo (john@torjo.com) and Matthew Wilson (matthew@synesis.com.au)
//
// Permission to copy, use, sell and distribute this software is granted
// provided this copyright notice appears in all copies.
// Permission to modify the code and to distribute modified code is granted
// provided this copyright notice appears in all copies, and a notice
// that the code was modified is included with the copyright notice.
//
// This software is provided "as is" without express or implied warranty,
// and with no claim as to its suitability for any purpose.
 
// See http://www.boost.org for updates, documentation, and revision history.


#ifndef BOOST_TRAITS_HPP_INCLUDED
#define BOOST_TRAITS_HPP_INCLUDED

#include <boost/mpl/has_xxx.hpp>
#include <boost/mpl/if.hpp>

#include <iterator>

#include <boost/type_traits/ice.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/is_pointer.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <stddef.h> // for size_t

#include <CGAL/PDB/internal/rangelib/priv/defs.hpp>

// rangelib = Range Library
namespace CGAL { namespace PDB { namespace internal { namespace rangelib { 

namespace detail {
    // creates classes has_iterator_category & has_value_type
    BOOST_MPL_HAS_XXX_TRAIT_DEF( iterator_category)
    BOOST_MPL_HAS_XXX_TRAIT_DEF( value_type)
    BOOST_MPL_HAS_XXX_TRAIT_DEF( reference)
    BOOST_MPL_HAS_XXX_TRAIT_DEF( is_range_type)

    BOOST_MPL_HAS_XXX_TRAIT_DEF( key_type)
    
    // allow extracting information from a given type
    struct get_container_details_const {
        template< class container_type> struct underlying_type {
            typedef typename container_type::const_iterator iterator;
          // moved const
            typedef const typename container_type::value_type reference;
            typedef const typename container_type::value_type value_type;
            typedef value_type* pointer;

//            typedef iterator const_iterator;
//            typedef reference const_reference;
//            typedef pointer const_pointer;
            //typedef typename container_type::const_pointer pointer;
        };
    };

    struct get_container_details {
        template< class container_type> struct underlying_type {
            typedef typename container_type::iterator iterator;
            typedef typename container_type::value_type reference;
            typedef typename container_type::value_type value_type;
            //typedef typename container_type::pointer pointer;
            typedef value_type* pointer;

//            typedef typename container_type::const_iterator const_iterator;
//            typedef typename container_type::const_reference const_reference;
//            typedef typename container_type::const_pointer const_pointer;
        };
    };

    struct get_ptr_details {
        template< class ptr_type> struct underlying_type {
            typedef ptr_type iterator;
#ifdef __BORLANDC__
            typedef /* typename */ ::boost::remove_pointer<ptr_type>::type value_type;
#else /* ? __BORLANDC__ */
            typedef typename ::boost::remove_pointer<ptr_type>::type value_type;
#endif /* __BORLANDC__ */
            typedef value_type & reference;
            typedef value_type * pointer;

//            typedef const iterator const_iterator;
//            typedef const reference const_reference;
//            typedef const pointer const_pointer;
        };
    };
    struct get_iterator_details {
            struct has_ref_ptr {
                template< class i_type> struct find {
                    typedef typename i_type::reference reference;
                    typedef typename i_type::pointer pointer;
                };
            };        
            struct has_no_ref_ptr {
                template< class i_type> struct find {
                    typedef typename i_type::value_type value_type;
                    typedef value_type& reference;
                    typedef value_type* pointer;
                };
            };


        template< class i_type> struct underlying_type {
            typedef typename i_type::value_type value_type;
#if !defined(CGAL_PDB_BOOST_RTL_WORKAROUND_VC6)
            // FIXME if value_type is const, append const to the reference and pointer
            typedef typename i_type::reference reference;
            typedef typename i_type::pointer pointer;
#else
            // VC6 - default STL implementation does not have these defined
            //
            // however, note the our filter iterators have these defined;
            // if so, use them
            BOOST_STATIC_CONSTANT( bool, has_ref = ::CGAL::PDB::internal::rangelib::detail::has_reference< i_type>::value);
            typedef typename ::boost::mpl::if_c< has_ref, 
                ::CGAL::PDB::internal::rangelib::detail::get_iterator_details::has_ref_ptr, ::CGAL::PDB::internal::rangelib::detail::get_iterator_details::has_no_ref_ptr>::type refptr_find_type;

            typedef typename refptr_find_type::template find<i_type>::reference reference;
            typedef typename refptr_find_type::template find<i_type>::pointer pointer;
#endif

//            typedef iterator const_iterator;
//            typedef const reference const_reference;
//            typedef const pointer const_pointer;
            typedef i_type iterator;
        };
    };


    template< class iterator_type>
    class iterator_traits {
    private:
        BOOST_STATIC_CONSTANT( bool,  is_ptr = ::boost::is_pointer< iterator_type>::value);

        // * wrapper, so that we can derive from - needed by has_iterator_category/ has_value_type!
        template< class any_type> struct wrapper {};

        // * find if it's an iterator or a pointer
#ifdef __BORLANDC__
        typedef ::boost::mpl::if_c< is_ptr, 
#else /* ? __BORLANDC__ */
        typedef typename ::boost::mpl::if_c< is_ptr, 
#endif /* __BORLANDC__ */
            ::CGAL::PDB::internal::rangelib::detail::get_ptr_details, ::CGAL::PDB::internal::rangelib::detail::get_iterator_details >::type ask_type;
    public:
        typedef typename ask_type::template underlying_type< iterator_type>::iterator iterator;
        typedef typename ask_type::template underlying_type< iterator_type>::value_type value_type;
        typedef typename ask_type::template underlying_type< iterator_type>::reference reference;
        typedef typename ask_type::template underlying_type< iterator_type>::pointer pointer;

        typedef iterator const_iterator;
        typedef reference const_reference;
    };


    template< class container_type >
    class container_traits {
    private:    
        BOOST_STATIC_CONSTANT( bool, is_c = ::boost::is_const< container_type>::value);
        typedef typename ::boost::mpl::if_c< is_c,
            ::CGAL::PDB::internal::rangelib::detail::get_container_details_const, 
            ::CGAL::PDB::internal::rangelib::detail::get_container_details>::type chosen_type;
    public:
        typedef typename chosen_type::template underlying_type< container_type>::iterator iterator;
        typedef typename chosen_type::template underlying_type< container_type>::value_type value_type;
        typedef typename chosen_type::template underlying_type< container_type>::reference reference;
        typedef typename chosen_type::template underlying_type< container_type>::pointer pointer;
    };

}


// forward declaration
template< class container_type> struct crange;


// finds the real range type, given a range/ container or const container
template< class r>
struct range_finder {

    typedef typename ::boost::remove_const<r>::type nonconst_r;
    
    BOOST_STATIC_CONSTANT( bool, is_range = ::CGAL::PDB::internal::rangelib::detail::has_is_range_type<nonconst_r>::value);
    typedef typename ::boost::mpl::if_c< is_range, nonconst_r, ::CGAL::PDB::internal::rangelib::crange<r> >::type range_type;

    typedef typename range_type::iterator iterator_type;
    typedef typename ::std::iterator_traits<iterator_type>::difference_type difference_type;
};



namespace detail {
    template< class function>
    struct first_argument_finder {
        // FIXME - in the future, allow function to be a C++ function,
        //         which I will automatically wrap into a binary_function (partial specialization)
        typedef typename function::first_argument_type ref_type;
        typedef typename ::boost::remove_reference<ref_type>::type argument_type;
    };
}

}}}}


#endif // CGAL_PDB_BOOST_RTL_HPP_INCLUDED
