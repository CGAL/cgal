#include <regex_dictionary.h>
#include <regex++.h>
#include <iostream>

using namespace std;

typedef hash_map< string, Regular_expression > Regex_map;

static Regex_map regex_map;

void
regex_register( const string& name, const string &regex ) {
  regex_map[name].init( regex.c_str() );
  //std::cerr << "!! Warning: registered name=[" << name << "] regex=[" << regex << "]" << std::endl;
}

bool
regex_does_match( const string& name, const string& text ) {
  bool result = regex_map[name].match( text );
  //std::cerr << "!! checking " << name << " against " << text << " yields: " << result << std::endl;
  //std::cerr << "match1 : " << regex_get_submatch(1) << std::endl;
  //std::cerr << "match2 : " << regex_get_submatch(2) << std::endl;
  //std::cerr << "match3 : " << regex_get_submatch(3) << std::endl;
  return result;
}

const string&
regex_get_submatch( unsigned int num ) {
  return Regular_expression::get_match( num );
}
