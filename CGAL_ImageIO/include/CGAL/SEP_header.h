// Copyright (c) 2016, 2017 GeometryFactory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_SEP_HEADER_HPP
#define CGAL_SEP_HEADER_HPP

#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>

#include <boost/tuple/tuple.hpp>
#include <boost/fusion/adapted/boost_tuple.hpp>
#include <boost/array.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#ifdef CGAL_SEP_READER_DEBUG
#  include <boost/tuple/tuple_io.hpp>
// The debug the Boost.Spirit.Qi parser, this is needed:
namespace std {
template <typename A, typename B>
inline std::ostream& operator<<(std::ostream& out, const std::pair<A, B>& id) {
  return out << id.first << " " << id.second;
}
} //end namespace std
#endif // CGAL_SEP_READER_DEBUG

namespace CGAL {

struct SEP_header_aux
{
  typedef std::map<std::string, int> Int_dict;
  typedef std::map<std::string, double> Double_dict;
  typedef std::map<std::string, std::string> String_dict;
  Int_dict int_dict;
  Double_dict double_dict;
  String_dict string_dict;
  
public:
  template <typename Tuple_string_variant>
  SEP_header_aux& operator<<(const Tuple_string_variant& tuple)
  {
    using boost::get;
    visitor vis(this, get<0>(tuple));
    boost::apply_visitor(vis, get<1>(tuple));
    return *this;
  }

private:
  struct visitor : public boost::static_visitor<> {
    SEP_header_aux* self;
    std::string key;
    visitor(SEP_header_aux* header, std::string key)
      : self(header)
      , key(key)
    {}

    template <typename T>
    void operator()(const T& t) {
      // std::cerr << "My assignement ("
      //           << typeid(t).name() << "): "
      //           << key << "=" << t << std::endl;
      self->add(key, t);
    }
  };

  void add(const std::string& key, int v) {
    int_dict[key] = v;
    double_dict[key] = v;
  }

  void add(const std::string& key, const double& v) {
    double_dict[key] = v;
  }

  void add(const std::string& key, const std::string& v) {
    string_dict[key] = v;
  }
}; // end struct SEP_header_aux

} // end namespace CGAL

BOOST_FUSION_ADAPT_STRUCT(
    CGAL::SEP_header_aux,
    ( CGAL::SEP_header_aux::Int_dict,     int_dict)
    ( CGAL::SEP_header_aux::Double_dict,  double_dict)
    ( CGAL::SEP_header_aux::String_dict,  string_dict)
)

namespace CGAL {

class SEP_header {
  
  boost::array<std::size_t, 3> _n;
  boost::array<double, 3> _d;
  boost::array<double, 3> _o;

  SEP_header_aux::String_dict _string_dict;

  int _dim;

public:
  /// constructor
  SEP_header(std::string fileName) : _dim(-1) {
    std::ifstream input(fileName.c_str());
    if(!input) {
      std::cerr << "Error: cannot open the header file \"" 
                << fileName << "\"!\n";
      return;
    }
    SEP_header_aux header;
    if(parseHeader(input, header))
    {
      get_n(1) = header.int_dict["n1"];
      get_n(2) = header.int_dict["n2"];
      get_n(3) = header.int_dict["n3"];
      _dim = 0;
      for(int i = 1; i <= 3; ++i) {
        if(get_n(i) == 0) get_n(i) = 1;
        if(get_n(i) != 1) ++_dim;
      }
      get_d(1) = header.double_dict["d1"];
      get_d(2) = header.double_dict["d2"];
      get_d(3) = header.double_dict["d3"];
      get_o(1) = header.double_dict["o1"];
      get_o(2) = header.double_dict["o2"];
      get_o(3) = header.double_dict["o3"];
      _string_dict = header.string_dict;
    }
  }

  /// copy-constructor
  SEP_header(const SEP_header& other)
    : _n(other._n)
    , _d(other._d)
    , _o(other._o)
    , _dim(other._dim)
  {}

  /// const getters
  /// @{
  double o(int i) const { return _o[i-1]; }
  double d(int i) const { return _d[i-1]; }
  std::size_t n(int i) const { return _n[i-1]; }
  int dimension() const { return _dim;}
  std::string string_field(std::string key) const {
    SEP_header_aux::String_dict::const_iterator it = _string_dict.find(key);
    if(it == _string_dict.end()) return std::string();
    else return it->second;
  }
  /// @}
  
  /// non-const getters
  /// @{
  double& get_o(int i) { return _o[i-1]; }
  double& get_d(int i) { return _d[i-1]; }
  std::size_t& get_n(int i) { return _n[i-1]; }
  int& get_dim() { return _dim; }
  /// @}

  std::ostream& display_information(std::string fileName,
                                   std::ostream& out) const {
    out << "Header file: " << fileName << std::endl;
    out << "Parameters:\n";
    out << "dimension = " << dimension() << "\n";
    for(int i = 1; i <= 3 ; ++i) {
      out << "d" << i << " = " << d(i) << "\n";
      out << "n" << i << " = " << n(i) << "\n";
      out << "o" << i << " = " << o(i) << "\n";
    }
    out << "in = " << string_field("in") << std::endl;
    return out;
  }

private:
  template <typename Iterator>
  struct sep_header_grammar : boost::spirit::qi::grammar<Iterator, SEP_header_aux()>
  {
    sep_header_grammar()
      : sep_header_grammar::base_type(header,
                                      "sep_header_grammar")
    {
      namespace qi = boost::spirit::qi;
      namespace ascii = boost::spirit::ascii;
      using qi::lexeme;
      using ascii::char_;
      using qi::int_;
      using qi::real_parser;
      using boost::phoenix::at_c;
      using qi::_1;
      using qi::_2;
      using qi::_val;

      qi::real_parser<double,
                      qi::strict_real_policies<double> > double_with_dot;

      skip_comment = ('#' >> *(qi::char_ - qi::eol) >> qi::eol);
      skip_comment.name("my skip comment");

      skip = skip_comment | qi::space;
      skip.name("my skipper");

      identifier = +((qi::alnum|char_('_'))[_val += _1]);
      identifier.name("my identifier");

      string = +(char_ - qi::space)[_val += _1];
      string.name("my string");

      quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
      quoted_string.name("my quoted string");

      variant_of_value %= ( double_with_dot
                            |
                            int_
                            |
                            (quoted_string | string) );
      variant_of_value.name("my value");

      assignment %= (identifier >> '=' >> variant_of_value);
      assignment.name("my assignment");

      header = (*skip) >> +( assignment[_val << _1] >> (*skip));
      header.name("SEP header");

      using boost::phoenix::construct;
      using boost::phoenix::val;
      qi::on_error<qi::fail>
        (
         header,
         std::cout
         << val("Error! Expecting ")
         << qi::_4                               // what failed?
         << val(" here: \"")
         << construct<std::string>(qi::_3, qi::_2)   // iterators to error-pos, end
         << val("\"")
         << std::endl
         );
#ifdef CGAL_SEP_READER_DEBUG
      qi::debug(skip);
      qi::debug(skip_comment);
      qi::debug(identifier);
      qi::debug(header);
      qi::debug(variant_of_value);
      qi::debug(assignment);
      qi::debug(string);
      qi::debug(quoted_string);
#endif // CGAL_SEP_READER_DEBUG
    } // end constructor of sep_header_grammar

    typedef boost::variant<double,
                           int,
                           std::string> value;
    typedef boost::tuple<std::string, value> entry_type;

    boost::spirit::qi::rule<Iterator>                    skip_comment;
    boost::spirit::qi::rule<Iterator>                    skip;
    boost::spirit::qi::rule<Iterator, std::string()>     identifier;
    boost::spirit::qi::rule<Iterator, std::string()>     quoted_string;
    boost::spirit::qi::rule<Iterator, std::string()>     string;
    boost::spirit::qi::rule<Iterator, value()>           variant_of_value;
    boost::spirit::qi::rule<Iterator, entry_type() >     assignment;
    boost::spirit::qi::rule<Iterator, SEP_header_aux() > header;
  }; // end class template sep_header_grammar<Iterator>

  bool parseHeader(std::ifstream& input, SEP_header_aux& header) {
    std::string file_content;

    input.seekg(0, std::ios::end);
    file_content.reserve(input.tellg());
    input.seekg(0, std::ios::beg);

    file_content.assign(std::istreambuf_iterator<char>(input),
                        std::istreambuf_iterator<char>());
    if(!input) return false;

    typedef std::string::const_iterator iterator;

    iterator begin(file_content.begin()), end(file_content.end());

    sep_header_grammar<iterator> grammar;

    bool b = boost::spirit::qi::parse(begin, end, grammar, header);
    return b && (begin == end);
  }
}; // end class SEP_header

} // end namespace CGAL

#endif // CGAL_SEP_HEADER_HPP
