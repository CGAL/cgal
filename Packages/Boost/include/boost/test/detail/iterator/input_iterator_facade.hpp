//  (C) Copyright Gennadiy Rozental 2004.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.
//
//  File        : $RCSfile$
//
//  Version     : $Revision$
//
//  Description : 
// ***************************************************************************

#ifndef BOOST_INPUT_ITERATOR_FACADE_HPP_071894GER
#define BOOST_INPUT_ITERATOR_FACADE_HPP_071894GER

// Boost
#include <boost/iterator/iterator_facade.hpp>

namespace boost {

namespace unit_test {

// ************************************************************************** //
// **************          input_iterator_core_access          ************** //
// ************************************************************************** //

class input_iterator_core_access
{
#if defined(BOOST_NO_MEMBER_TEMPLATE_FRIENDS) || BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x551))
public:
#else
    template <class I, class V, class R, class TC> friend class input_iterator_facade;
#endif

    template <class Facade>
    static bool get( Facade& f )
    {
        return f.get();
    }

private:
    // objects of this class are useless
    input_iterator_core_access(); //undefined
};

// ************************************************************************** //
// **************            input_iterator_facade             ************** //
// ************************************************************************** //

template<typename Derived,
         typename ValueType,
         typename Reference = ValueType const&,
         typename Traversal = single_pass_traversal_tag>
class input_iterator_facade : public iterator_facade<Derived,ValueType,Traversal,Reference>
{
public:
    // Constructor
    input_iterator_facade() : m_valid( false ), m_value() {}

protected: // provide access to the Derived
    void                init()
    {
        m_valid = true;
        increment();
    }

    // Data members
    bool                m_valid;
    ValueType           m_value;

private:
    friend class boost::iterator_core_access;

    // iterator facade interface implementation
    void                increment()
    {
        // we make post-end incrementation indefinetly safe 
        if( m_valid )
            m_valid = input_iterator_core_access::get( *static_cast<Derived*>(this) );
    }
    Reference dereference() const
    {
        return m_value;
    }

    // iterator facade interface implementation
    bool                equal( input_iterator_facade const& rhs ) const
    {
        // two invalid iterator equals, inequal otherwise
        return !m_valid && !rhs.m_valid;
    }
};

} // namespace unit_test

} // namespace boost

// ***************************************************************************
//  Revision History :
//  
//  $Log$
//  Revision 1.1  2004/11/20 10:52:23  spion
//  Initial revision
//
//  Revision 1.6  2004/10/05 02:17:14  rogeeff
//  como fix
//
//  Revision 1.5  2004/07/19 12:29:57  rogeeff
//  guard rename
//  mingw port
//
//  Revision 1.4  2004/06/07 07:33:50  rogeeff
//  detail namespace renamed
//
//  Revision 1.3  2004/06/05 11:03:12  rogeeff
//  input_iterator_adaptor simplified
//  token_iterator added
//
//  Revision 1.2  2004/05/25 10:29:09  rogeeff
//  use standard getline
//  eliminate initialize
//  proper handle \n in wide case
//
//  Revision 1.1  2004/05/21 06:30:10  rogeeff
//  ifstream_line_iterator added
//
// ***************************************************************************

#endif // BOOST_INPUT_ITERATOR_FACADE_HPP_071894GER

