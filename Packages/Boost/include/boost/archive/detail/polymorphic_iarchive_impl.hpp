#ifndef BOOST_ARCHIVE_DETAIL_POLYMORPHIC_IARCHIVE_IMPL_HPP
#define BOOST_ARCHIVE_DETAIL_POLYMORPHIC_IARCHIVE_IMPL_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// polymorphic_iarchive_impl.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cstddef>
#include <string>
#include <ostream>
#include <boost/noncopyable.hpp>
#include <boost/cstdint.hpp>

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/archive/polymorphic_iarchive.hpp>

namespace boost { 
namespace archive {
namespace detail{

class basic_iserializer;
class basic_pointer_iserializer;

template<class ArchiveImplementation>
class polymorphic_iarchive_impl : 
    public polymorphic_iarchive,
    // note: gcc dynamic cross cast fails if the the derivation below is 
    // not public.  I think this is a mistake.
    public /*protected*/ ArchiveImplementation,
    private boost::noncopyable
{
private:
    // these are used by the serialization library.
    virtual void load_object(
        void *t, 
        const basic_iserializer & bis
    ){
        ArchiveImplementation::load_object(t, bis);
    }
    virtual const basic_pointer_iserializer * load_pointer(
        void * & t, 
        const basic_pointer_iserializer * bpis_ptr,
        const basic_pointer_iserializer * (*finder)(
            const boost::serialization::extended_type_info & type
        )
    ){
        return ArchiveImplementation::load_pointer(t, bpis_ptr, finder);
    }
    virtual void delete_created_pointers(){
        ArchiveImplementation::delete_created_pointers();
    }
    virtual unsigned int library_version() const{
        return ArchiveImplementation::library_version();
    }
    // primitive types the only ones permitted by polymorphic archives
    virtual void load(bool & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(char & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(signed char & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(unsigned char & t){
        ArchiveImplementation::load(t);
    }
    #ifndef BOOST_NO_CWCHAR
    #ifndef BOOST_NO_INTRINSIC_WCHAR_T
    virtual void load(wchar_t & t){
        ArchiveImplementation::load(t);
    }
    #endif
    #endif
    virtual void load(short & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(unsigned short & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(int & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(unsigned int & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(long & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(unsigned long & t){
        ArchiveImplementation::load(t);
    }
    #if !defined(BOOST_NO_INTRINSIC_INT64_T)
    virtual void load(boost::int64_t & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(boost::uint64_t & t){
        ArchiveImplementation::load(t);
    }
    #endif
    virtual void load(float & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(double & t){
        ArchiveImplementation::load(t);
    }
    virtual void load(std::string & t){
        ArchiveImplementation::load(t);
    }
    #ifndef BOOST_NO_STD_WSTRING
    virtual void load(std::wstring & t){
        ArchiveImplementation::load(t);
    }
    #endif
    virtual void load_binary(void * t, std::size_t size){
        ArchiveImplementation::load(t);
    }

    // used for xml and other tagged formats default does nothing
    virtual void load_start(const char * name){
        ArchiveImplementation::load_start(name);
    }
    virtual void load_end(const char * name){
        ArchiveImplementation::load_end(name);
    }

    virtual void register_basic_serializer(const detail::basic_iserializer & bis){
        ArchiveImplementation::register_basic_serializer(bis);
    }
public:
    // define operators for non-const arguments.  Don't depend one the const
    // ones below because the compiler MAY make a temporary copy to
    // create the const parameter (Though I havn't seen this happen). 
    // the >> operator
    template<class T>
    polymorphic_iarchive & operator>>(T & t){
        // if this assertion trips. It means we're trying to load a
        // const object with a compiler that doesn't have correct
        // funtion template ordering.  On other compilers, this is
        // handled below.
        BOOST_STATIC_ASSERT(! boost::is_const<T>::value);
        return polymorphic_iarchive::operator>>(t);
    }

    // the & operator 
    template<class T>
    polymorphic_iarchive & operator&(T & t){
        return polymorphic_iarchive::operator&(t);
    }

    // define the following pair in order to permit passing of const and non_const
    // temporary objects. These are needed to properly implement serialization
    // wrappers.

    #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
    // the >> operator
    template<class T>
    polymorphic_iarchive & operator>>(const T & t){
        // this should only be used for wrappers.  Check that here
        return polymorphic_iarchive::operator>>(t);
    }
    // the & operator 
    template<class T>
    polymorphic_iarchive & operator&(const T & t){
        return polymorphic_iarchive::operator&(t);
    }
    #endif
    // all current archives take a stream as constructor argument
    template <class _Elem, class _Tr>
    polymorphic_iarchive_impl(
        std::basic_istream<_Elem, _Tr> & is, 
        unsigned int flags = 0
    ) :
        ArchiveImplementation(is, flags | boost::archive::no_header)
    {
        // postpone archive initialization until build is complete
        if(0 == flags & no_header)
            ArchiveImplementation::init();
    }
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_DETAIL_POLYMORPHIC_IARCHIVE_IMPL_HPP
