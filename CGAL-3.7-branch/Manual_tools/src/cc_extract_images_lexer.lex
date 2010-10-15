%{
#include <iostream>
#include "mstring.h"

extern bool cc_option;

static string remove_quotes( string s ) {
  s = s.substr( 1, s.size() - 2 );
  if( s.size() >= 2 && s[0] == '.' && s[1] == '/' )
    s = s.substr( 2, s.size() -2 );
  return s;
}

static void report_filename( string s ) {
  s = remove_quotes( s );
  const bool cc_file = s.find( "cc_" ) == std::size_t(0);
  if( cc_file == cc_option )
    std::cout << s << std::endl;
}

%}

%option caseless
%option noyywrap
%option nodefault

%x FilenameMode

ws               [ \t\n\r]
ext              gif|png|jpg|jpeg
letter           [^<>\"]
filename         \"{letter}+"."{ext}\"
unquotedfilename {letter}+"."{ext}

/* src{ws}*={ws}*\"{ws}*       BEGIN( FilenameMode ); */

%%
src{ws}*={ws}*{unquotedfilename} {
                 std::cerr << "WARNING: unquoted images not allowed in XHTML"
                           << std::endl; }
src{ws}*={ws}*                  { BEGIN(FilenameMode); break; }
href{ws}*={ws}*                  { BEGIN(FilenameMode); break; }
<FilenameMode>{filename}        { report_filename(yytext);
                                  BEGIN(INITIAL); }
<*>{ws}                         /* ignore whitespace */
<*>.                            { BEGIN(INITIAL); }

%%

