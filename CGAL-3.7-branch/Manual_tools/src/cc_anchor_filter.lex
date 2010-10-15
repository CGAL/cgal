%{
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

extern "C" {
#include <errno.h>
#include <dirent.h>
}


#ifdef __GNUC__

#if ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

#include <ext/hash_map>
#include <ext/hash_set>

using __gnu_cxx::hash;
using __gnu_cxx::hash_map;

namespace __gnu_cxx {

#else // ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

#include <hash_map>
#include <hash_set>

namespace std {

#endif // ((__GNUC__ > 3) || ((__GNUC__ == 3) && (__GNUC_MINOR__ >= 1)))

template <>
struct hash<std::string>
{
  size_t operator()(const std::string& str) const
  {
    unsigned long h = 0;
    const char* s = str.data();
    for (size_t len = str.length(); len > 0; --len, ++s)
      h = 5*h + (unsigned long)(*s);
    return size_t(h);
  }
};
};

#endif

using namespace std;

/* Hack, to get rid of the yywrap. */
#define YY_SKIP_YYWRAP 1
#define yywrap() 1

class Dictionary {
  typedef hash_map< string, string > Map;
  Map map;
public:
  void add( const string& key, const string& value ) {
    if( !is_defined( key ) )
      map[key] = value;
  }

  bool is_defined( const string& key ) const {
   return map.find( key ) != map.end();
  }

  string operator[]( const string& key ) const {
   Map::const_iterator it = map.find( key );
   if( it != map.end() ) {
     return it->second;
   } else {
     cerr << " internal warning: unknown dictionary key [" << key << "]" << endl;
     return string();
   }
  }
};


/* program return code: 0 = o.k., 1 = error, 2 = warning */
int error_code = 0;

/* Is used to count nesting levels during PARAMMODE */
int nesting;

/* Is used to count nesting levels during NOLINKMODE */
int linknesting;


Dictionary
  dict_labels,
  dict_labels_text,
  dict_bib,
  dict_anchormode_bib,
  dict_cc,
  dict_internal;

const string reference_icon = "&#x261E;"; // Unicode for pointing finger.

string   output_path, reftext;
ofstream output_file;

void compress_spaces_in_string( string& s ) {
    for ( size_t i = 0; i < s.size(); ++i) {
    if ( isspace( s[i])) {
        size_t k = 1;
        while ( i + k < s.size() && isspace( s[ i + k]))
        ++k;
        s.replace( i, k, " ");
    }
    }
}


bool
match_cc_idfier( string s, ostream *out = NULL ) {
  string::size_type begin = 0;
  string::size_type len   = s.length();

  do {
    const string id = s.substr( begin, len );
    //cout << " !! Warning: checking [" << id << "]" << endl;

    string::size_type lt = id.find( "&lt;" );
    string::size_type gt = id.rfind( "&gt;" );
    bool template_param_matched = false;

    if( lt != string::npos && gt != string::npos ) {
      string tparam = s.substr( begin + lt + 4, gt - lt - 4 );
      //std::cout << tparam << std::endl;
      if( match_cc_idfier( tparam ) ) {
        template_param_matched = true;
      }
    }

    if( !template_param_matched && dict_cc.is_defined( id ) ) {
      // only match the whole id if _no_ substring inside < > matched
      if( out )
        *out << dict_cc[ id ];
      else
        return true;
    }
    else {
      string::size_type t = id.find_last_of( "&:, " );
      if( t != string::npos ) {
        if( t >= string::size_type(1) && id[t-1] == ':' )
          --t;
        len = t;
        continue;
      } else
        if( out ) *out << id;
    }

    begin += len;
    while( begin < s.length() ) {
      string x = s.substr(begin, 4);
      if( x == "&lt;" || x == "&gt;" ) {
        if( out ) *out << x;
        begin += 4;
     } else if( x[0] == ':' && x[1] == ':' ) {
        if( out ) *out << "::";
        begin += 2;
      } else if( x[0] == ',' ) {
        if( out ) *out << ',';
        begin += 1;
      } else if( x[0] == ' ' ) {
        if( out ) *out << ' ';
        begin += 1;
      } else
        break;
    }
    len   = s.length() - begin;
  } while( len > 0 );
  return false;
}

#define ECHO output_file << yytext

%}

/* Possibly nested blocks where cross-linking is suppressed. Use     */
/* carefully, since it blocks all kinds of replacements, such as     */
/* links to chapters, citations, references etc.                     */
%x NOLINKMODE

/* Avoid substitutions while within an HTML Tag. */
%x TAGMODE

%x ANCHORMODE

/* Avoid substitutions while within a Heading. */
%x HEADINGMODE

/* Avoid substitutions while within the Header. */
%x HEADERMODE

/* aggressive crosslinking: every CCidfier1 is crosslinked, when listed in dict */
%x CROSSLINKMODE

%x BIBMODE
%x BIBMODE2


letter          [a-zA-Z]
noletter        [^a-zA-Z]
digit           [0-9]
CCletter        [a-zA-Z_]
LT              "&lt;"
GT              "&gt;"
CCidfier        ({CCletter}({CCletter}|{digit})*)
CCidfier1       {CCidfier}({CCidfier}|[ ]?({LT}[ ]?)+{CCidfier}|([ ]?{GT})*[ ]?"::"[ ]?{CCidfier}|[ ]?[,][ ]?{CCidfier})*([ ]?{GT})*
intsign         signed[ ]+|unsigned[ ]+
intsize         long[ ]+long[ ]+|long[ ]+|short[ ]+
inttype         ({intsign}?{intsize}?"int")|({intsign}?"char")
floattype       float|double|long[ ]+double
CCidfier2       {CCidfier1}|{inttype}|{floattype}
ws              [ \t\n\r]*
par             {ws}("<"[pP]">"{ws})*
glue            {ws}("<"[pP]">"{ws})*"<!GLUE>"{ws}("<"[pP]">"{ws})*
head            "<"[hH][eE][aA][dD]
endhead         "</"[hH][eE][aA][dD]">"
nolinkbegin     "<!NoLinkBegin>"
nolinkend       "<!NoLinkEnd>"
nolinksync      "<!NoLinkSync>"
cccbegin        "[cccbegin]"
cccend          "[cccend]"

%%


{nolinkbegin}                  { BEGIN( NOLINKMODE); linknesting = 1; }
{nolinkend}                    { fputs( "\nERROR: Found \\lcNoLinkEnd without a matching \\lcNoLinkBegin.\n", stderr); error_code = 1; }
{nolinksync}                   { ;}
<NOLINKMODE>{nolinkbegin}      { ++linknesting; }
<NOLINKMODE>{nolinkend}        { --linknesting;
                                 if ( ! linknesting) BEGIN( INITIAL);
                               }
<NOLINKMODE>{nolinksync}       { fputs( "\nERROR: \\lcNoLinkSync used in \\lcNoLinkBegin/End block, which indicates missing or runaway \\lcNoLinkEnd.\n", stderr);
                                 linknesting = 0;
                                 BEGIN( INITIAL);
                                 error_code = 1;
                               }

<NOLINKMODE>.                  { ECHO;}

{cccbegin}                     { BEGIN( CROSSLINKMODE ); }
{cccend}                       { cerr << endl << "ERROR: found [cccend] withouth matching [cccbegin]" << endl; error_code = 1; }

{head}                         { BEGIN(HEADERMODE); ECHO; }
<HEADERMODE>"<I>"|"</I>"       { ; }
<HEADERMODE>([^[]|[\n])*{endhead} { BEGIN(INITIAL); ECHO; /* Avoid substitutions in the head */ }

"<A HREF="["][^"]*"Biblio_"[^"]+["]">" { ECHO; BEGIN(BIBMODE2); }
<BIBMODE2>[^<]+"</A>" {
  yytext [ yyleng - 4 ] = '\0'; // cut off trailing </A>
  if( dict_bib.is_defined( yytext ) )
    output_file  << dict_bib[ yytext ] << "</A>";
  else
    output_file  << yytext << "</A>";
  BEGIN(INITIAL);
}

"<A NAME="["]"Biblio_"[^"]+    { ECHO; BEGIN( BIBMODE ); } // we silently assume, that anchor == text
<BIBMODE>["]"></A><B>["[^]]+"]</B>" {
  yytext [ yyleng - 5 ] = '\0'; // cut off trailing </B>
  const char *my_yytext = yytext + 10;
  if( dict_bib.is_defined( my_yytext ) )
    output_file  << "\"></A><B>[" << dict_bib[ my_yytext ] << "]</B>";
  else {
    cerr << "warning: undefined bib-entry \"" << my_yytext << "\"" << endl;
    output_file  << "\"></A><B>[" << my_yytext << "]</B>"; // at least print something
  }
  BEGIN( INITIAL );
}

"<"[hH][1-6]                   { BEGIN( HEADINGMODE); ECHO;}
<HEADINGMODE>"</"[hH][1-6]">"  { BEGIN( INITIAL); ECHO;}
<NOLINKMODE,HEADERMODE,HEADINGMODE>{cccbegin}        { ; } // consume [cccbegin] markers in headers/headingss/nolinkmode
<NOLINKMODE,HEADERMODE,HEADINGMODE>{cccend}          { ; }
<HEADINGMODE>.                 { ECHO; }

&[^;]+[;]                      { ECHO; }

"<"[^ \t\n\r]                  { BEGIN(ANCHORMODE); ECHO; }
<ANCHORMODE>[^>]*[>]           { BEGIN(INITIAL); ECHO; }

[\[]reftext[:][^\]]+[\]] {
  yytext [ yyleng - 1 ] = '\0'; // cut off trailing ]
  reftext = string(yytext + 9); // skip \reftext:
}

[\[]ref[:][^\]]+[\]] {
   yytext [ yyleng - 1 ] = '\0'; // cut off trailing ]
   const char *my_yytext = yytext + 5;
   if( dict_labels.is_defined( my_yytext ) ) {
     //cerr << " !! label [" << my_yytext << "] is defined as [" << dict_labels[ my_yytext ] << " [" << dict_labels_text[ my_yytext ] << "]" << endl;
     output_file  << "<A HREF=\"" << dict_labels[ my_yytext ] << "\">";
     if( reftext != string() ) {
       if( reftext == "\\icon" )
         output_file << reference_icon;
       else
         output_file << reftext;
     }
     else if( dict_labels_text.is_defined( my_yytext ) )
       output_file << dict_labels_text[ my_yytext ];
     else
       output_file << reference_icon;
     output_file << "</A>";
   } else {
     cerr << " !! Warning: undefined label \"" << my_yytext << "\"" << endl;
     if( reftext != string() )
       output_file << reftext;
   }
   reftext = string();
}

[\[]internal[:][^\]]+[\]] {
   yytext [ yyleng - 1 ] = '\0';    // cut off trailing "]
   string my_yytext( yytext + 10 ); // skip [internal:"

   if( dict_internal.is_defined( my_yytext ) )
     output_file  << dict_internal[ my_yytext ];
}

<CROSSLINKMODE>{CCidfier2} {
  //output_file << "matched [" << yytext << "] as ccidfieR" << std::endl;
  string my_yytext( yytext );
  // so that "long  long       int" is properly crosslinked
  compress_spaces_in_string( my_yytext );
  //std::cout << "[" << my_yytext << "]" << std::endl;

  match_cc_idfier( my_yytext, &output_file );
}
<CROSSLINKMODE>[<][^>]+[>] { ECHO; } // do not crosslink inside html tags
<CROSSLINKMODE>.        { ECHO; }
<CROSSLINKMODE>{cccend} { BEGIN( INITIAL ); }

%%

void
process( const string &filename ) {
  // open input file
  FILE *in = fopen( filename.c_str(), "r" );
  if( !in ) {
    cerr << " !! Error: could not open file " << filename << " for reading." << endl;
    return;
  }

  // open output file
  output_file.open( (output_path + filename).c_str() );
  if( !output_file ) {
    cerr << " !! Error: could not open file " << (output_path + filename) << " for writing." << endl;
    return;
  }

  // create new flex buffer
  yyin = in;
  yy_switch_to_buffer( yy_create_buffer( yyin, YY_BUF_SIZE));


  // actual lex
  yylex();

  // cleanup
  YY_BUFFER_STATE buf = YY_CURRENT_BUFFER;
  yy_delete_buffer( buf);
  fclose( in );
  output_file.close();
}

int main( int argc, char** argv ) {
  if( argc < 3 ) {
    cerr << "usage: " << argv[0] << " ruleset outputpath" << endl;
    exit(1);
  }

  ifstream in( argv[1] );
  if( !in ) {
    cerr << "cannot find ruleset \"" << argv[1] << "\"" << endl;
  }
  //cerr << "loading ruleset : " << argv[1] << endl;
  string type, key, value, line;
  while( in ) {
    in >> type;
    std::getline( in, line );

    std::string::size_type start = line.find_first_not_of( " " );
    std::string::size_type end   = line.find( '\t', start + 1 );
    if( end == std::string::npos ) {
      cerr << " !! rule malformatted: " << type << line << endl;
      continue;
    }
    key   = line.substr( start, end -1 );
    value = line.substr( end + 1, line.find_last_not_of( " " ) - end );

    //cerr << "type: " << type << " key: [" << key << "] value: [" << value << "]" << endl;
    if( type == "l" ) {
      dict_labels.add(key,value);
    } else if( type == "lt" ) {
      dict_labels_text.add(key,value);
    } else if ( type == "b" ) {
      dict_bib.add(key,value);
    } else if ( type == "c" ) {
      dict_cc.add(key,value);
      //cerr << "c key: [" << key << "] value: [" << value << "]" << endl;
    } else if ( type == "i" ) {
      //cerr << "i key: [" << key << "] value: [" << value << "]" << endl;
      dict_internal.add(key,value);
    }
  }

  output_path = argv[2];

  if( output_path[ output_path.length() - 1 ] != '/' )
    output_path += '/';
  DIR *pdir;
  struct dirent *pent;

  pdir = opendir( "." );
  if( !pdir ) {
   exit(1);
  }
  errno = 0;
  // foreach file in *.html { process( file ) }
  while( ( pent = readdir(pdir) ) ) {
    string filename = pent->d_name;

    // filename.endsWith( ".html" )
    if( filename.length() > 5 && filename.compare( filename.length() - 5 , 5, ".html" ) == 0 )
      process( filename );
  }
  closedir( pdir );

  return error_code;
}
