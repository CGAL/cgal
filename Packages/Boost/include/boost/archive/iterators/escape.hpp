#ifndef BOOST_ARCHIVE_ITERATORS_ESCAPE_HPP
#define BOOST_ARCHIVE_ITERATORS_ESCAPE_HPP

// MS compatible compilers support #pragma once
#if defined(_MSC_VER) && (_MSC_VER >= 1020)
# pragma once
#endif

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// escape.hpp

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com . 
// Use, modification and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <cassert>

#include <boost/config.hpp> // for BOOST_DEDUCED_TYPENAME
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/iterator_traits.hpp>

namespace boost { 
namespace archive {
namespace iterators {

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// insert escapes into text

template<class Derived, class Base>
class escape : 
    public boost::iterator_adaptor<
        Derived, 
        Base, 
        BOOST_DEDUCED_TYPENAME boost::iterator_value<Base>::type,
        single_pass_traversal_tag,
        BOOST_DEDUCED_TYPENAME boost::iterator_value<Base>::type
    >{
    typedef BOOST_DEDUCED_TYPENAME boost::iterator_value<Base>::type base_value_type;
    typedef BOOST_DEDUCED_TYPENAME boost::iterator_reference<Base>::type reference_type;
    friend class boost::iterator_core_access;

    typedef BOOST_DEDUCED_TYPENAME boost::iterator_adaptor<
        Derived, 
        Base, 
        base_value_type,
        single_pass_traversal_tag,
        base_value_type
    > super_t;

    typedef escape<Derived, Base> this_t;

    bool equal(const this_t & rhs) const {
        return 
            NULL == m_bnext
            && NULL == m_bend
            && this->base_reference() == rhs.base_reference()
        ;
    }

    //Access the value referred to 
    reference_type dereference() const {
        return m_current_value;
    }

   void increment(){
        if(++m_bnext < m_bend){
            m_current_value = *m_bnext;
            return;
        }
        ++(this->base_reference());
        m_bnext = NULL;
        m_bend = NULL;
        m_current_value = (static_cast<Derived *>(this))->fill(m_bnext, m_bend);
    }

    // buffer to handle pending characters
    const base_value_type *m_bnext;
    const base_value_type *m_bend;
    BOOST_DEDUCED_TYPENAME boost::iterator_value<Base>::type m_current_value;
    bool m_full;
public:
    escape(Base base) : 
        super_t(base),
        m_bnext(NULL),
        m_bend(NULL)
    {
        m_current_value = static_cast<Derived *>(this)->fill(m_bnext, m_bend);
    }
};

} // namespace iterators
} // namespace archive
} // namespace boost

#endif // BOOST_ARCHIVE_ITERATORS_ESCAPE_HPP
