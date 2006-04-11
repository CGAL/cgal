#ifndef REGEXPP_H
#define REGEXPP_H

#include <sys/types.h>
#include <regex.h>
#include <mstring.h>
#include <vector>

struct Regular_expression {
  struct Exception {};
  Regular_expression() : regmatch(NULL) {}
  Regular_expression( const string& rx );
  Regular_expression( const char*   rx );
  ~Regular_expression();

  void init( const char *rx );

  bool match( const string& s );
  bool match( const char*   s );

  static const string& get_match( unsigned int num );
protected:

  static vector <string>  submatch;
  regex_t                 regex;
  regmatch_t             *regmatch;
};

#endif








