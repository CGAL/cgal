/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// text_oarchive_impl.ipp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com .
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <string>
#include <boost/config.hpp>
#include <locale>
#include <cstddef> // size_t

#include <boost/config.hpp>
#if defined(BOOST_NO_STDC_NAMESPACE)
namespace std{ 
    using ::size_t; 
} // namespace std
#endif

#ifndef BOOST_NO_CWCHAR
#include <cwchar>
#ifdef BOOST_NO_STDC_NAMESPACE
namespace std{ using ::wcslen; }
#endif
#endif

#include <boost/archive/codecvt_null.hpp>
#include <boost/archive/add_facet.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace boost { 
namespace archive {

//////////////////////////////////////////////////////////////////////
// implementation of basic_text_oprimitive overrides for the combination
// of template parameters used to create a text_oprimitive

template<class Archive>
void text_oarchive_impl<Archive>::save(const char * s)
{
    unsigned len = std::ostream::traits_type::length(s);
    *this->This() << len;
    this->This()->newtoken();
    os << s;
}

template<class Archive>
void text_oarchive_impl<Archive>::save(const std::string &s)
{
    unsigned size = s.size();
    *this->This() << size;
    this->This()->newtoken();
    os << s;
}

#ifndef BOOST_NO_CWCHAR
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
template<class Archive>
void text_oarchive_impl<Archive>::save(const wchar_t * ws)
{
    std::size_t l = std::wcslen(ws);
    * this->This() << l;
    this->This()->newtoken();
    os.write((const char *)ws, l * sizeof(wchar_t)/sizeof(char));
}
#endif

#ifndef BOOST_NO_STD_WSTRING
template<class Archive>
void text_oarchive_impl<Archive>::save(const std::wstring &ws)
{
    std::size_t l = ws.size();
    * this->This() << l;
    this->This()->newtoken();
    os.write((const char *)(ws.data()), l * sizeof(wchar_t)/sizeof(char));
}
#endif
#endif // BOOST_NO_CWCHAR

} // namespace archive
} // namespace boost

