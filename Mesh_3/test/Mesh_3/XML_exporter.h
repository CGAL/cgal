// Copyright (c) 2012 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
//
// Author(s)     : Clement Jamin
//

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>

template<typename value_type = std::string>
class Simple_XML_exporter
{
public:
  typedef value_type                        Value_type;
  typedef std::vector<value_type>           Element;
  typedef std::map<std::string, value_type> Element_with_map;
  typedef std::vector<Element>              List_of_elements;

  Simple_XML_exporter(
    const std::string &list_name,
    const std::string &element_name,
    const std::vector<std::string> &subelement_names,
    bool add_timestamp = true)
  : m_list_name(list_name),
    m_element_name(element_name),
    m_subelement_names(subelement_names),
    m_add_timestamp(add_timestamp)
  {}

  bool add_element(const Element &element)
  {
    if (element.size() == m_subelement_names.size())
    {
      m_list_of_elements.push_back(element);
      return true;
    }
    else
    {
      std::cerr << "ERROR: element.size() == m_subelement_names.size()" << std::endl;
      return false;
    }
  }

  bool add_element(Element_with_map &element)
  {
    Element elt;

    std::vector<std::string>::const_iterator
      it_subelement_name = m_subelement_names.begin();
    std::vector<std::string>::const_iterator
      it_subelement_name_end = m_subelement_names.end();
    for ( ; it_subelement_name != it_subelement_name_end ; ++it_subelement_name)
    {
      elt.push_back(element[*it_subelement_name]);
    }

    return add_element(elt);
  }

  bool export_to_xml(const std::string &filename) const
  {
    std::ofstream xmlfile;
    xmlfile.open(filename.c_str());
    xmlfile << "<?xml version='1.0'?>" << std::endl;
    xmlfile << "<" << m_list_name << ">" << std::endl;

    typename List_of_elements::const_iterator it_element = m_list_of_elements.begin();
    typename List_of_elements::const_iterator it_element_end = m_list_of_elements.end();
    for (int id = 1 ; it_element != it_element_end ; ++it_element, ++id)
    {
      xmlfile << "  <" << m_element_name << ">" << std::endl;
      std::vector<std::string>::const_iterator
        it_subelement_name = m_subelement_names.begin();
      std::vector<std::string>::const_iterator
        it_subelement_name_end = m_subelement_names.end();

      if (m_add_timestamp)
        xmlfile << "    <id> " << time(NULL) << " </id>" << std::endl;

      for (int i = 0 ;
           it_subelement_name != it_subelement_name_end ;
           ++it_subelement_name, ++i)
      {
        xmlfile
          << "    <" << *it_subelement_name << "> "
          << (*it_element)[i]
          << " </" << *it_subelement_name << ">" << std::endl;
      }
      xmlfile << "  </" << m_element_name << ">" << std::endl;
    }

    xmlfile << "</" << m_list_name << ">" << std::endl;
    xmlfile.close();
    return 0;

  }

protected:
  std::string                       m_list_name;
  std::string                       m_element_name;
  std::vector<std::string>          m_subelement_names;
  List_of_elements                  m_list_of_elements;
  bool                              m_add_timestamp;
};





template<typename value_type = std::string>
class Streaming_XML_exporter
{
public:
  typedef value_type                        Value_type;
  typedef std::vector<value_type>           Element;
  typedef std::map<std::string, value_type> Element_with_map;
  typedef std::vector<Element>              List_of_elements;

  Streaming_XML_exporter(
    const std::string &filename,
    const std::string &list_name,
    const std::string &element_name,
    const std::vector<std::string> &subelement_names,
    bool add_timestamp = true)
  : m_list_name(list_name),
    m_element_name(element_name),
    m_subelement_names(subelement_names),
    m_add_timestamp(add_timestamp)
  {
    m_xml_fstream.open(filename.c_str());
    if (m_xml_fstream.good())
    {
      m_xml_fstream << "<?xml version='1.0'?>" << std::endl;
      m_xml_fstream << "<" << m_list_name << ">" << std::endl;
    }
    else
    {
      std::cerr << "Could not open file '" << filename << "'." << std::endl;
    }
  }

  virtual ~Streaming_XML_exporter()
  {
    close_file();
  }

  void close_file()
  {
    m_xml_fstream.close();
  }

  bool add_element(const Element &element)
  {
    if (element.size() == m_subelement_names.size())
    {
      m_xml_fstream << "  <" << m_element_name << ">" << std::endl;
      std::vector<std::string>::const_iterator
        it_subelement_name = m_subelement_names.begin();
      std::vector<std::string>::const_iterator
        it_subelement_name_end = m_subelement_names.end();

      if (m_add_timestamp)
      {
        m_xml_fstream << "    <id> " << time(NULL) << " </id>" << std::endl;
      }

      for (int i = 0 ;
           it_subelement_name != it_subelement_name_end ;
           ++it_subelement_name, ++i)
      {
        m_xml_fstream
          << "    <" << *it_subelement_name << "> "
          << element[i]
          << " </" << *it_subelement_name << ">" << std::endl;
      }
      m_xml_fstream << "  </" << m_element_name << ">" << std::endl;

      // Save current pointer position
      std::ofstream::streampos pos = m_xml_fstream.tellp();
      // Close the XML file (temporarily) so that the XML file is always correct
      m_xml_fstream << "</" << m_list_name << ">" << std::endl;
      // Restore the pointer position so that the next "add_element" will overwrite
      // the end of the file
      m_xml_fstream.seekp(pos);

      m_xml_fstream.flush();
      return true;
    }
    else
    {
      std::cerr << "ERROR: element.size() == m_subelement_names.size()" << std::endl;
      return false;
    }
  }

  bool add_element(Element_with_map &element)
  {
    Element elt;

    std::vector<std::string>::const_iterator
      it_subelement_name = m_subelement_names.begin();
    std::vector<std::string>::const_iterator
      it_subelement_name_end = m_subelement_names.end();
    for ( ; it_subelement_name != it_subelement_name_end ; ++it_subelement_name)
    {
      elt.push_back(element[*it_subelement_name]);
    }

    return add_element(elt);
  }

protected:
  std::ofstream                     m_xml_fstream;
  std::string                       m_list_name;
  std::string                       m_element_name;
  std::vector<std::string>          m_subelement_names;
  bool                              m_add_timestamp;
};
