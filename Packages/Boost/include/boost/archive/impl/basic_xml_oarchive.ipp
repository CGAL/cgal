/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// basic_xml_oarchive.ipp:

// (C) Copyright 2002 Robert Ramey - http://www.rrsd.com .
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.

#include <boost/archive/basic_xml_archive.hpp>
#include <boost/archive/basic_xml_oarchive.hpp>

namespace boost {
namespace archive {

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8
// implemenations of functions common to both types of xml output

template<class Archive>
void basic_xml_oarchive<Archive>::write_attribute(
    const char *attribute_name, 
    int t,
    const char *conjunction
){
    this->This()->put(' ');
    this->This()->put(attribute_name);
    this->This()->put(conjunction);
    this->This()->save(t);
    this->This()->put('"');
}
    
template<class Archive>
void basic_xml_oarchive<Archive>::write_attribute(
    const char *attribute_name, 
    const char *key
){
    this->This()->put(' ');
    this->This()->put(attribute_name);
    this->This()->put("=\"");
    this->This()->put(key);
    this->This()->put('"');
}
    
template<class Archive>
void basic_xml_oarchive<Archive>::indent(){
    int i;
    for(i = depth; i-- > 0;)
        this->This()->put('\t');
}

template<class Archive>
void basic_xml_oarchive<Archive>::save_start(const char *name)
{
    if(NULL == name)
        return;
    end_preamble();
    if(depth > 0){
        this->This()->put('\n');
        indent();
    }
    ++depth;
    this->This()->put('<');
    this->This()->save(name);
    pending_preamble = true;
    indent_next = false;
}

template<class Archive>
void basic_xml_oarchive<Archive>::save_end(const char *name)
{
    if(NULL == name)
        return;
    end_preamble();
    --depth;
    if(indent_next){
        this->This()->put('\n');
        indent();
    }
    indent_next = true;
    this->This()->put("</");
    this->This()->save(name);
    this->This()->put('>');
    if(0 == depth)
        this->This()->put('\n');
}

template<class Archive>
void basic_xml_oarchive<Archive>::end_preamble(){
    if(pending_preamble){
        this->This()->put('>');
        pending_preamble = false;
    }
}

template<class Archive>
void basic_xml_oarchive<Archive>::save_override(const object_id_type & t, int)
{
    int i = t.t; // extra .t is for borland
    write_attribute(OBJECT_ID, i, "=\"_"); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(
    const object_reference_type & t,
    int
){
    int i = t.t; // extra .t is for borland
    write_attribute(OBJECT_REFERENCE, i, "=\"_"); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(const version_type & t, int)
{
    int i = t.t; // extra .t is for borland
    write_attribute(VERSION, i); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(const class_id_type & t, int)
{
    write_attribute(CLASS_ID, t); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(
    const class_id_reference_type & t,
    int
){
    write_attribute(CLASS_ID_REFERENCE, t); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(
    const class_id_optional_type & t,
    int
){
    write_attribute(CLASS_ID, t); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(const class_name_type & t, int)
{
    const char * key = t;
    if(NULL == key)
        return;
    write_attribute(CLASS_NAME, key); 
}
template<class Archive>
void basic_xml_oarchive<Archive>::save_override(const tracking_type & t, int)
{
    write_attribute(TRACKING, t.t); // extra .t is for borland 
}

template<class Archive>
void basic_xml_oarchive<Archive>::init(){
    // xml header
    header = true;
    this->This()->put("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");
    this->This()->put("<!DOCTYPE boost_serialization>\n");
    // xml document wrapper - outer root
    this->This()->put("<boost_serialization");
    write_attribute("signature", ARCHIVE_SIGNATURE); 
    write_attribute("version", static_cast<unsigned int>(ARCHIVE_VERSION)); 
    this->This()->put(">\n");
}

template<class Archive>
basic_xml_oarchive<Archive>::basic_xml_oarchive(unsigned int flags) :
    depth(0),
    indent_next(false),
    pending_preamble(false),
    header(false)
{
}

template<class Archive>
basic_xml_oarchive<Archive>::~basic_xml_oarchive(){
    if(!header)
        return;
    this->This()->put("</boost_serialization>\n");
}

} // namespace archive
} // namespace boost
