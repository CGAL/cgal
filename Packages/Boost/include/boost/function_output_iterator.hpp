// (C) Copyright Jeremy Siek 2001. Permission to copy, use, modify,
// sell and distribute this software is granted provided this
// copyright notice appears in all copies. This software is provided
// "as is" without express or implied warranty, and with no claim as
// to its suitability for any purpose.

// Revision History:

// 27 Feb 2001   Jeremy Siek
//      Initial checkin.

#ifndef BOOST_FUNCTION_OUTPUT_ITERATOR_HPP
#define BOOST_FUNCTION_OUTPUT_ITERATOR_HPP

#include <iterator>

namespace boost {

  template <class UnaryFunction>
  class function_output_iterator {
    typedef function_output_iterator self;
  public:
    typedef std::output_iterator_tag iterator_category;
    typedef void                value_type;
    typedef void                difference_type;
    typedef void                pointer;
    typedef void                reference;

    explicit function_output_iterator(const UnaryFunction& f = UnaryFunction())
      : m_f(f) {}

    struct output_proxy {
      output_proxy(UnaryFunction& f) : m_f(f) { }
      template <class T> output_proxy& operator=(const T& value) {
        m_f(value); 
        return *this; 
      }
      UnaryFunction& m_f;
    };
    output_proxy operator*() { return output_proxy(m_f); }
    self& operator++() { return *this; } 
    self& operator++(int) { return *this; }
  private:
    UnaryFunction m_f;
  };

  template <class UnaryFunction>
  inline function_output_iterator<UnaryFunction>
  make_function_output_iterator(const UnaryFunction& f = UnaryFunction()) {
    return function_output_iterator<UnaryFunction>(f);
  }

} // namespace boost

#endif // BOOST_FUNCTION_OUTPUT_ITERATOR_HPP
