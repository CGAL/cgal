#ifndef BOOST_ARCHIVE_DETAIL_POLYMORPHIC_OARCHIVE_IMPL_HPP
#define BOOST_ARCHIVE_DETAIL_POLYMORPHIC_OARCHIVE_IMPL_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// polymorphic_oarchive_impl.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <string>
#include <ostream>
#include <boost/noncopyable.hpp>
#include <boost/cstdint.hpp>
#include <cstddef> // size_t

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#include <boost/archive/polymorphic_oarchive.hpp>

namespace boost { 
namespace archive {
namespace detail{

class basic_oserializer;
class basic_pointer_oserializer;

template<class ArchiveImplementation>
class polymorphic_oarchive_impl : 
    public polymorphic_oarchive,
    // note: gcc dynamic cross cast fails if the the derivation below is 
    // not public.  I think this is a mistake.
    public /*protected*/ ArchiveImplementation,
    private boost::noncopyable
{
private:
    // these are used by the serialization library.
    virtual void save_object(
        const void *x, 
        const detail::basic_oserializer & bos
    ){
        ArchiveImplementation::save_object(x, bos);
    }
    virtual void save_pointer(
        const void * t, 
        const detail::basic_pointer_oserializer * bpos_ptr
    ){
        ArchiveImplementation::save_pointer(t, bpos_ptr);
    }
    virtual void save_null_pointer(){
        ArchiveImplementation::save_null_pointer();
    }
    // primitive types the only ones permitted by polymorphic archives
    virtual void save(const bool t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const char t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const signed char t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const unsigned char t){
        ArchiveImplementation::save(t);
    }
    #ifndef BOOST_NO_CWCHAR
    #ifndef BOOST_NO_INTRINSIC_WCHAR_T
    virtual void save(const wchar_t t){
        ArchiveImplementation::save(t);
    }
    #endif
    #endif
    virtual void save(const short t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const unsigned short t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const int t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const unsigned int t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const long t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const unsigned long t){
        ArchiveImplementation::save(t);
    }
    #if !defined(BOOST_NO_INTRINSIC_INT64_T)
    virtual void save(const boost::int64_t t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const boost::uint64_t t){
        ArchiveImplementation::save(t);
    }
    #endif
    virtual void save(const float t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const double t){
        ArchiveImplementation::save(t);
    }
    virtual void save(const std::string & t){
        ArchiveImplementation::save(t);
    }
    #ifndef BOOST_NO_STD_WSTRING
    virtual void save(const std::wstring & t){
        ArchiveImplementation::save(t);
    }
    #endif
    virtual unsigned int library_version() const{
        return ArchiveImplementation::library_version();
    }
    virtual void save_binary(const void * t, std::size_t size){
        ArchiveImplementation::save(t);
    }

    // used for xml and other tagged formats default does nothing
    virtual void save_start(const char * name){
        ArchiveImplementation::save_start(name);
    }
    virtual void save_end(const char * name){
        ArchiveImplementation::save_end(name);
    }
    virtual void end_preamble(){
        ArchiveImplementation::end_preamble();
    }
    virtual void register_basic_serializer(const detail::basic_oserializer & bos){
        ArchiveImplementation::register_basic_serializer(bos);
    }
public:
    // to avoie ambiguities when using this class directly, trap an pass one
    // to the implemenation these operations.
    // note: we presume that older compilers will never create a const
    // argument from a non-const by copyiing
    template<class T>
    polymorphic_oarchive & operator<<(const T & t){
        return polymorphic_oarchive::operator<<(t);
    }

    // the & operator 
    template<class T>
    polymorphic_oarchive & operator&(const T & t){
        return polymorphic_oarchive::operator&(t);
    }

    // define operators for non-const arguments.  Don't depend one the const
    // ones below because the compiler MAY make a temporary copy to
    // create the const parameter (Though I havn't seen this happen). 
    #ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
        // the << operator
        template<class T>
        polymorphic_oarchive & operator<<(T & t){
            // if trap here, we're saving a tracted non-const
            // value - this could be a stack variable with the same
            // address for multiple items. This would be the source of very 
            // subtle errors and should be double checked
            // BOOST_STATIC_WARNING(
            //     serialization::tracking_level == serialization::track_never
            // );
            return polymorphic_oarchive::operator<<(t);
        }
        // the & operator 
        template<class T>
        polymorphic_oarchive & operator&(T & t){
            return polymorphic_oarchive::operator&(t);
        }
    #endif

    // all current archives take a stream as constructor argument
    template <class _Elem, class _Tr>
    polymorphic_oarchive_impl(
        std::basic_ostream<_Elem, _Tr> & os, 
        unsigned int flags = 0
    ) :
        ArchiveImplementation(os, flags | no_header)
    {
        // postpone archive initialization until build is complete
        if(0 == (flags & no_header))
            ArchiveImplementation::init();
    }
};

} // namespace detail
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_DETAIL_POLYMORPHIC_OARCHIVE_IMPL_HPP
