/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// xml_woarchive_impl.ipp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com .
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/config.hpp>
#ifndef BOOST_NO_STD_WSTREAMBUF

#include <ostream>
#include <string>
#include <algorithm>
#include <locale>

#include <cstring> // strlen
#include <boost/config.hpp> // msvc 6.0 needs this to suppress warnings
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::strlen; 
} // namespace std
#endif

#include <boost/throw_exception.hpp>
#include <boost/utf8_codecvt_facet.hpp>

#include <cstring>
#include <cstdlib> // mbtowc

#include <boost/config.hpp> // for BOOST_DEDUCED_TYPENAME
#if defined(BOOST_NO_STDC_NAMESPACE) && ! defined(__LIBCOMO__)
namespace std{ 
    using ::strlen; 
    using ::mbtowc; 
} //std
#endif

#include <boost/pfto.hpp>

#include <boost/archive/iterators/xml_escape.hpp>
#include <boost/archive/iterators/wchar_from_mb.hpp>
#include <boost/archive/iterators/ostream_iterator.hpp>
#include <boost/archive/iterators/dataflow_exception.hpp>

#include <boost/archive/add_facet.hpp>

namespace boost {
namespace archive {

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// implemenations of functions specific to wide char archives

std::wostream & operator<<(std::wostream &os, const char *t){
    for(;;){
        wchar_t wc;
        int result = std::mbtowc(&wc, t, 10 /* max number */);
        if(0 < result)
            os.put(wc);
        else
        if(0 == result)
            break;
        else
            boost::throw_exception(
                iterators::dataflow_exception(
                    iterators::dataflow_exception::invalid_conversion
                )
            );
    }
    return os;
}

std::wostream & operator<<(std::wostream &os, const char t){
    wchar_t wc;
    std::mbtowc(&wc, &t, 1);
    os.put(wc);
    return os;
}

// copy chars to output escaping to xml and widening characters as we go
template<class InputIterator>
void save_iterator(std::wostream &os, InputIterator begin, InputIterator end){
    typedef iterators::wchar_from_mb<
        iterators::xml_escape<InputIterator>
    > xmbtows;
    std::copy(
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(begin)),
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(end)),
        boost::archive::iterators::ostream_iterator<wchar_t>(os)
    );
}

template<class Archive>
void xml_woarchive_impl<Archive>::save(const std::string & s){
    // note: we don't use s.begin() and s.end() because dinkumware
    // doesn't have string::value_type defined. So use a wrapper
    // around these values to implement the definitions.
    const char * begin = s.data();
    const char * end = begin + s.size();
    save_iterator(os, begin, end);
}

#ifndef BOOST_NO_STD_WSTRING
template<class Archive>
void xml_woarchive_impl<Archive>::save(const std::wstring & ws){
    typedef iterators::xml_escape<std::wstring::const_iterator> xmbtows;
    std::copy(
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(ws.data())),
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(ws.data() + ws.size())),
        boost::archive::iterators::ostream_iterator<wchar_t>(os)
    );
}
#endif //BOOST_NO_STD_WSTRING

template<class Archive>
void xml_woarchive_impl<Archive>::save(const char * s){
   save_iterator(os, s, s + std::strlen(s));
}

#ifndef BOOST_NO_INTRINSIC_WCHAR_T
template<class Archive>
void xml_woarchive_impl<Archive>::save(const wchar_t * ws){
    os << ws;
    typedef iterators::xml_escape<const wchar_t *> xmbtows;
    std::copy(
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(ws)),
        xmbtows(BOOST_MAKE_PFTO_WRAPPER(ws + std::wcslen(ws))),
        boost::archive::iterators::ostream_iterator<wchar_t>(os)
    );
}
#endif

template<class Archive>
xml_woarchive_impl<Archive>::xml_woarchive_impl(
    std::wostream & os_,
    unsigned int flags
) :
    basic_text_oprimitive<std::wostream>(
        os_,
        true // don't change the codecvt - use the one below
    ),
    basic_xml_oarchive<Archive>(flags)
{
    // Standard behavior is that imbue can be called
    // a) before output is invoked or
    // b) after flush has been called.  This prevents one-to-many
    // transforms (such as one to many transforms from getting
    // mixed up.  Unfortunately, STLPort doesn't respect b) above
    // so the restoration of the original archive locale done by
    // the locale_saver doesn't get processed,
    // before the current one is destroyed.
    // so the codecvt doesn't get replaced with the orginal
    // so closing the stream invokes codecvt::do_unshift
    // so it crashes because the corresponding locale that contained
    // the codecvt isn't around any more.
    // we can hack around this by using a static codecvt that never
    // gets destroyed.
    if(0 == (flags & no_codecvt)){
        utf8_codecvt_facet_wchar_t *pfacet;
        #if defined(__SGI_STL_PORT) || defined(_STLPORT_VERSION)
            static utf8_codecvt_facet_wchar_t facet(static_cast<size_t>(1));
            pfacet = & facet;
        #else
            pfacet = new utf8_codecvt_facet_wchar_t;
        #endif
        archive_locale.reset(add_facet(std::locale::classic(), pfacet));
        os.imbue(* archive_locale);
    }
    if(0 == (flags & no_header))
        this->init();
}

} // namespace archive
} // namespace boost

#endif //BOOST_NO_STD_WSTREAMBUF
