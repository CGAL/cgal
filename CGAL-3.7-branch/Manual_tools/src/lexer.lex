/**************************************************************************
 
  html_lex.yy
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : lexical scanner for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : flex, bison, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

%{
#include <lexer.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string_conversion.h>
#include <internal_macros.h>
#include <error.h>
#include <config.h>
#include <syntax.tab.h>

#include <input.h>
#include <lex_include_impl.h>
#include <macro_dictionary.h>

int    old_state = 0;    // Old parser state before param. parsing.

// Used to communicate with the parser for DEFWITHARGS's tokens.
int number_of_args = 0;

string yy_string;             // Used like yytext, if we need a temp. string.

// Used to parse {} nested expressions as (La)TeX macro parameters
// using ParameterMode.
string current_macro;         // Current macro to be expanded
int    parameter_nesting = 0;
int    parameter_count  = 0;  // used to collect params during macro expansion.
                              // Changing from 1 to 0 stops parameter parsing.
string parameters[36];        // Parameter list incl. optional parameters
string parameter_format;      // String describing the parameter format
                              // m = mandatory, o = optional parameter
int    parameter_index  = 0;  // Index in parameter list. Incr. from 0.
int    parameter_options = 0; // counts number of optional parameters.
int    parameter_endopt = 0;  // number of optional parameters at end.

bool   parameter_option = false; // indicates whether we parse an opt.param.

bool is_parameter_parsing_done( bool optional, bool more);
bool expand_macro();

// Used to pass the parameter from IfParameterParsing to
// is_parameter_parsing_done.
string current_parameter;


// String to collect the result of a CCParameterMode parameter.
string cc_string;
// String to store the filename of ccReferenceFile intermediately
// The filename is already processed to not contain illegal characters anymore.
string cc_filename;

/* Set this flag to 1 to switch back to old_state. */
int set_old_state = 0;

/* remember the necessary stop character for \verb"..." */
/* or stop environment for \begin{verbatim}             */
char   stop_character;
string stop_envir;

void count_newlines( const char* s);

#define skipspaces() if ( skipspaces_eof()) yyterminate()
bool   skipspaces_eof();

#define skipseparator() if ( skipseparator_eof()) yyterminate()
bool   skipseparator_eof();

#define skiplimitedspaces() if ( skiplimitedspaces_eof()) yyterminate()
bool   skiplimitedspaces_eof();

int    skiplimitedspaces_param();

#define skiplimitedsepspaces() if ( skiplimitedsepspaces_eof()) yyterminate()
bool   skiplimitedsepspaces_eof();

#define skipnonspaces() if ( skipnonspaces_eof()) yyterminate()
bool   skipnonspaces_eof();

size_t removespaces( char* s);
void   inc_line();

#define IfParameterParsing( option)                                          \
	bool more_param = false;                                             \
	if ( (option && (parameter_count == 0) && (parameter_endopt > 1))    \
	  || (!option && (parameter_count == 1) && (parameter_endopt > 0))) {\
	    skiplimitedsepspaces();                                          \
	    int c = yyinput();                                               \
	    if ( c == EOF)                                                   \
		yyterminate();                                               \
	    if ( c == '[') {                                                 \
		more_param = true;                                           \
	    } else                                                           \
		unput( c);                                                   \
	}                                                                    \
	if ( is_parameter_parsing_done( (option), more_param))

#define BeginParameterStart()        \
    if ( YY_START == AllttMode)      \
        BEGIN( AllttParameterStart); \
    else                             \
        BEGIN( ParameterStart)

#define BeginParameterMode()         \
    if ( YY_START == AllttMode)      \
        BEGIN( AllttParameterMode);  \
    else                             \
        BEGIN( ParameterMode)

// Some inline function to process yytext.

inline char* next_char( char* s, char c) {
    while ( *s != c)
        ++s;
    return s;
}

inline char* next_alpha( char* s) {
    while ( ! isalpha(*s))
        ++s;
    return s;
}

inline char* next_non_alpha( char* s) {
    while ( isalpha(*s))
        ++s;
    return s;
}

inline char* next_digit( char* s) {
    while ( ! isdigit(*s))
        ++s;
    return s;
}

#define CC_Special(x) count_newlines(yytext); cc_string+=x; break

#define YY_BREAK  /* a do nothing */

/* --------------------------------------------------------------------
    Parsing Modes:

    --  INITIAL:          parses (La)TeX

    --  AllttMode:        Parses alltt environments. Copies all characters
                          literally except \macros and {}.

    --  DelimiterMode     parses LaTeX \verb"..." statements until a
                          stop character occurs. Another example is \path.

    --  EndTokenMode      parses LaTeX verbose environments until an
                          \end{envir} occurs. Another example is cprog.

    --  IncludeMode       parses lciInclude filename,

    --  ParameterStart:   starts a (La)TeX macro parameter, i.e. 
                          comments, braces, brackets, and escaped symbols,
			  are recognized. White spaces are removed. 
                          A brace starts a normal parameter, a bracket 
			  starts an optional parameter, both change then
			  to ParameterMode. All other inputs (macro or
			  single symbols) are returned as single token
			  parameter.

    --  ParameterMode:    parses (La)TeX macro parameters body, i.e. 
                          comments, braces and escaped braces are
                          recognized. All other characters are
                          literally copied and at the end returned as
                          PARAMETER, the surrounding {} stripped off.
			  In case of [] a PARAMETER_OPTION is returned.

    --  AllttParameterStart: same as ParameterStart, but in AllttMode

    --  AllttParameterMode:  same as ParameterMode, but in AllttMode

    --  CCParameterMode:  Parses C++ parameters, similar to ParameterMode.
                          Counts only parantheses. Does not treat escaped
                          parantheses or TeX comments (%). Expands the
                          few macros known for CC style (e.g. \tt).

-------------------------------------------------------------------- */

%}

%x AllttMode
%x IncludeMode
%x ParameterStart
%x ParameterMode
%x AllttParameterStart
%x AllttParameterMode
%x CCParameterMode
%x DelimiterMode
%x EndTokenMode

sep             [\001]
sepsp           [\t \001]
seps            {sepsp}*
sepnls          [\t \n\r\001]*
letter          [a-zA-Z]
digit           [0-9]
CCletter        [a-zA-Z_]
idfier          {letter}+
envir           [a-zA-Z0-9*]+
envirmacro      (("\\begin@")|("\\end")){envir}
texmacro        ([\\]((.)|([\\]"*")|({idfier}"*"?)))|{envirmacro}
deftexmacro     {texmacro}(("@"[mo]+)?)
texmacroskip    [\\]{idfier}"*"?
filename        [^ \t\n\\\{\}\[\]()\001]+
space           [\t ]
w               {space}*
ws              {space}+
number          {digit}+

%%
 /* Mode switching can be triggered from the parser */
 /* ----------------------------------------------- */
	if (set_old_state) {
	    BEGIN( old_state);
	    set_old_state = 0;
	}

 /* Separator, Newlines, paragraphs, comments, and EOF        */
 /* --------------------------------------------------------- */

<INITIAL,AllttMode>{sep}  { break; /* ignore separator */ }

 /* Count line numbers in all modes for better error messages */
{w}"\n"{w}("\n"{w})+   {  /* create a TeX paragraph */
                    count_newlines( yytext);
		    yylval.text = "\n<P>\n\n";
		    return STRING;
}
<INITIAL,AllttMode>[\n]	{
		    inc_line();
		    yylval.character = yytext[0];
		    return CHAR;
		}
<INITIAL,AllttMode>[\\]"\n"      {
		    inc_line();
		    yylval.character = ' ';
		    return CHAR;
}
"%".*[\n]{w}  { /* Match one line TeX comments remove spaces in next line */
		    inc_line();
		    break;
	        }
"%".*[\n]{w}("\n"{w})+  { /* Match comments with an empty line -> par */
                    count_newlines( yytext);
		    yylval.text = "\n<P>\n\n";
		    return STRING;
	        }
"%".*  {            /* Match TeX comments at the last line in a file (EOF) */ 
                    break;}

<<EOF>>         {   if ( YY_START != INITIAL)
                        printErrorMessage( ParsingStateError);
                    yyterminate();
}

 /* Handle include files      */
 /* ------------------------- */

<INITIAL,AllttMode>[\\](lciInclude){seps}[\{]{seps}   {  
		    old_state = YY_START;
		    BEGIN ( IncludeMode);
		    break;
}
<IncludeMode>{filename}          {
                    /* remove remaining characters before the '}' */
		    int c = yyinput();
		    while( c != EOF && c != '}') {
			if ( c == '\n')
			    inc_line();
			c = yyinput();
		    }
		    if ( c == EOF) {
			printErrorMessage( EOFInIncludeFilenameError);
			yyterminate();
		    }
		    BEGIN( old_state);
                    include_stack.push_file( yytext );
                    break;
}


 /* TeX Macro Definitions                         */
 /* LaTeX newcommand is as internal macro defined */
 /* --------------------------------------------- */

<INITIAL,AllttMode>[\\]lciBeginGroup    { return '{'; }
<INITIAL,AllttMode>[\\]lciEndGroup      { return '}'; }

<INITIAL,AllttMode>[\\]def{seps}{deftexmacro}{seps}([#][0-9])+ {
                    number_of_args = atoi( yytext + yyleng - 1);
		    old_state = YY_START;
                    BeginParameterStart();
                    char* s = yytext + 1;
		    s = next_char(s, '\\');
		    char* p = s;
		    ++s;
		    s = next_non_alpha( s);
		    *s = '\0';
	            yylval.text = newstr( p);
                    return DEFWITHARGS;
}
<INITIAL,AllttMode>[\\]gdef{seps}{deftexmacro}{seps}([#][0-9])+ {
                    number_of_args = atoi( yytext + yyleng - 1);
		    old_state = YY_START;
                    BeginParameterStart();
                    char* s = yytext + 1;
		    s = next_char(s, '\\');
		    char* p = s;
		    ++s;
		    s = next_non_alpha( s);
		    *s = '\0';
	            yylval.text = newstr( p);
                    return GDEFWITHARGS;
}
<INITIAL,AllttMode>[\\]g?def{seps}{deftexmacro}{seps}[#][^\{]* {
		    old_state = YY_START;
                    BeginParameterStart();
		    int i = 1;
		    while ( yytext[i] != '\\')
			++i;
	            yylval.text = newstr( yytext+i);
                    return DEFWITHUNKNOWNARGS;
}

 /* Macro Parameter Parsing   */
 /* ------------------------- */

<ParameterStart>[\\]"\n" { 
		    inc_line();
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
		    break;
}
<ParameterStart>"\n"     |
<ParameterStart>"%".*[\n]{w}  { /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    inc_line();
		    break;
}
<ParameterStart>{sepsp}+ { break;} /* ignore seps + sp before parameters */
<ParameterStart>"%".*    { break;} /* Match one line TeX comments */
                                   /* at the last line in file */

<ParameterStart>"\{"     {
		    current_parameter = string();
                    ++parameter_nesting;
                    BeginParameterMode();
		    break;
}
<ParameterStart>"\["     {
		    current_parameter = string();
                    ++parameter_nesting;
		    parameter_option = true;
                    BeginParameterMode();
		    break;
}
<ParameterStart>{texmacroskip}  {
		    current_parameter = yytext;
		    int c = skiplimitedspaces_param();
		    if ( c == EOF)
			yyterminate();
		    if ( isalpha(c) || c == '*' || c == '#')
			current_parameter += SEPARATOR;
		    IfParameterParsing( false)
		        return PARAMETER;
		    break;
}
<ParameterStart>{deftexmacro}  |
<ParameterStart>.              {
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
		    break;
}


<ParameterMode>[\\]"\n" |
<ParameterMode>"\n"     { /* keep newlines */
		    current_parameter += yytext;
                    inc_line();
		    break;
}
<ParameterMode>"%".*[\n]{w}("\n"{w})+  { 
                    /* Match comments with an empty line -> par */
		    current_parameter += "\n\n";
                    count_newlines( yytext);
		    break;
}
<ParameterMode>"%".*[\n]{w}  { /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    current_parameter += SEPARATOR;
                    inc_line();
		    break;
}
<ParameterMode>"%".*    { /* Match one line TeX comments */
                          /* at the last line in file */
		    current_parameter += SEPARATOR;
		    break;
}
<ParameterMode>"\{"     {
		    current_parameter += yytext;
                    ++parameter_nesting;
		    break;
}
<ParameterMode>"\}"     {
                    if ( --parameter_nesting == 0) {
			IfParameterParsing( false) {
			    if ( parameter_option) {
				printErrorMessage( ParameterOptionError);
				parameter_option = false;
				return PARAMETER_OPTION;
			    }
			    return PARAMETER;
			}
		    } else
			current_parameter += yytext;
		    break;
}
<ParameterMode>"]"      {
                    if ( parameter_nesting == 1 && parameter_option) {
			parameter_nesting--;
			parameter_option = false;
			IfParameterParsing( true)
			    return PARAMETER_OPTION;
		    } else
			current_parameter += yytext;
		    break;
}
<ParameterMode>{texmacroskip}  {
		    current_parameter += yytext;
		    int c = skiplimitedspaces_param();
		    if ( c == EOF)
			yyterminate();
		    if ( isalpha(c) || c == '*' || c == '#')
			current_parameter += SEPARATOR;
		    break;
}
<ParameterMode>{deftexmacro}   |
<ParameterMode>[^\\\{\}%\]\n]+ |
<ParameterMode>.               {
		    current_parameter += yytext;
		    break;
}

 /* Macro Parameter Parsing in AllttMode  */
 /* ------------------------------------- */

<AllttParameterStart>[\\]"\n" | 
<AllttParameterStart>"\n"     {
		    inc_line();
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
		    break;
}
<AllttParameterStart>{sep}+   { break;} /* ignore seps before parameters */
<AllttParameterStart>"\{"     {
		    current_parameter = string();
                    ++parameter_nesting;
                    BEGIN( AllttParameterMode);
		    break;
}
<AllttParameterStart>"\["     {
		    current_parameter = string();
                    ++parameter_nesting;
		    parameter_option = true;
                    BEGIN( AllttParameterMode);
		    break;
}
<AllttParameterStart>{deftexmacro}  |
<AllttParameterStart>.              {
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
		    break;
}


<AllttParameterMode>[\\]"\n" |
<AllttParameterMode>"\n"     { /* keep newlines */
		    current_parameter += yytext;
                    inc_line();
		    break;
}
<AllttParameterMode>"\{"     {
		    current_parameter += yytext;
                    ++parameter_nesting;
		    break;
}
<AllttParameterMode>"\}"     {
                    if ( --parameter_nesting == 0) {
			IfParameterParsing( false) {
			    if ( parameter_option) {
				printErrorMessage( ParameterOptionError);
				parameter_option = false;
				return PARAMETER_OPTION;
			    }
			    return PARAMETER;
			}
		    } else
			current_parameter += yytext;
		    break;
}
<AllttParameterMode>"]"      {
                    if ( parameter_nesting == 1 && parameter_option) {
			parameter_nesting--;
			parameter_option = false;
			IfParameterParsing( true)
			    return PARAMETER_OPTION;
		    } else
			current_parameter += yytext;
		    break;
}
<AllttParameterMode>{deftexmacro}   |
<AllttParameterMode>[^\\\{\}%\]\n]+ |
<AllttParameterMode>.               {
		    current_parameter += yytext;
		    break;
}

 /* C++ Parameter Parsing     */
 /* ------------------------- */
<CCParameterMode>[\\]"\n" |
<CCParameterMode>"\n"     {
                    inc_line();
		    cc_string += ' ';
		    break;
}
<CCParameterMode>"\{"     {
                    if ( parameter_nesting++ > 0)
                        cc_string += yytext;
		    break;
}
<CCParameterMode>"\}"     {
                    if ( --parameter_nesting == 0) {
			BEGIN( old_state);
			include_stack.push_string( in_string->name(),
						   current_macro,
						   in_string->line());
			current_macro = string();
		    } else
			cc_string += yytext;
		    break;
}

 /* Special treatment of several macros in C++ text */
<CCParameterMode>[\\][^a-zA-Z] { /* capture all quoted special symbols */
                    cc_string += yytext[1];
		    break;
}
<CCParameterMode>[\\]"tt"{sepnls} { CC_Special("|T|");}
<CCParameterMode>[\\]"bf"{sepnls} { CC_Special("|B|");}
<CCParameterMode>[\\]"em"{sepnls} { CC_Special("|I|");}
<CCParameterMode>[\\]"it"{sepnls} { CC_Special("|I|");}
<CCParameterMode>[\\]"sl"{sepnls} { CC_Special("|I|");}
<CCParameterMode>[\\]"ccFont"{sepnls}       { CC_Special("|I|"); }
<CCParameterMode>[\\]"l"?"dots"{sepnls}     { CC_Special("..."); }

<CCParameterMode>{deftexmacro}    { 
                    printErrorMessage( MacroInCModeError);
                    cc_string += yytext; 
		    break;
}

<CCParameterMode>{sep}              { break; /* ignore separator */ }
<CCParameterMode>[^\n\{\}\\\001]+   |
<CCParameterMode>.                  {
                    cc_string += yytext;
		    break;
}


<INITIAL,AllttMode>[\\]lciParseCC{seps}[\{][^\}]*[\}]   {
                    count_newlines(yytext);
		    yytext[ yyleng - 1] = '\0';
		    char* s = yytext + 1; 
		    s = next_char(s, '{');
		    current_macro  = s + 1;
                    skiplimitedsepspaces();
		    int c = yyinput();
		    if ( c != '{') {
			printErrorMessage( ParseCCError);
			current_macro = string();
		    } else {
			old_state = YY_START;
			BEGIN( CCParameterMode);
			parameter_nesting = 1;
			cc_string = string();
		    }
		    break;
}

 /* Special Parser Modes      */
 /* ------------------------- */
[\\]lciBeginAlltt { 
                    skiplimitedspaces();
                    BEGIN( AllttMode); 
		    break;
                }
<AllttMode>[\\]lciEndAlltt { 
                    BEGIN( INITIAL);
                    skiplimitedspaces();
		    break;
                }


<INITIAL,AllttMode>[\\]lciParseUntilDelimiter{seps}[\{]{seps}{texmacro}{seps}[\}]   {
		    old_state = YY_START;
		    BEGIN( DelimiterMode);
                    char* s = yytext + 1;
		    s = next_char(s, '\\');
		    char* p = s;
		    ++s;
		    s = next_non_alpha( s);
		    *s = '\0';
		    current_macro  = p;
		    skipseparator();
		    while ( (stop_character = yyinput()) == '\n')
			inc_line();
		    break;
                }
<DelimiterMode>"\n"	{
                    inc_line();
		    yymore();
		    break;
		}
<DelimiterMode>.	{
		    if ( yytext[yyleng-1] == stop_character) {
		        BEGIN( old_state);
			yytext[yyleng-1] = '\0';
			Macro_item item = fetchMacro( current_macro);
			parameters[0] = yytext;
			include_stack.push_string( current_macro + " in "
						   +item.filename, 
						   expandMacro(current_macro,
							       item, 
							       parameters, 
							       1, 0),
						   item.line);
			parameters[0] = string();
			current_macro = string();
                    } else {
		        yymore();
		    }
		    break;
		}

<INITIAL,AllttMode>[\\]lciParseFile{seps}[\{]{seps}{texmacro}{seps}[\}]{seps}[\{][^\}]+[\}]   {
                    count_newlines(yytext);
		    yytext[ yyleng - 1] = '\0';
		    char* s = yytext + 1;
		    s = next_char(s, '\\');
		    char* p = s;
		    ++s;
		    s = next_non_alpha(s);
		    *s = '\0';
		    while ( *s++ != '{');
		    string filename(s);
		    crop_string(filename);
		    remove_separator(filename);
		    parameters[0] = string();
		    append_file_to_string( filename, parameters[0]);
		    Macro_item item = fetchMacro( p);
		    include_stack.push_string( string(p)+" in "
					       +item.filename, 
					       expandMacro(p,
							   item, 
							   parameters,
							   1, 0),
					       item.line);
		    parameters[0] = string();
		    break;
                }

<INITIAL,AllttMode>[\\]lciParseUntilEndToken{seps}[\{]{seps}{texmacro}{seps}[\}]{seps}[\{]{seps}{envir}{seps}[\}]   {
		    old_state = YY_START;
		    BEGIN( EndTokenMode);
		    char* s = yytext + 1;
		    s = next_char(s, '\\');
		    char* p = s;
		    ++s;
		    s = next_non_alpha(s);
		    *s = '\0';
		    current_macro  = p;
		    s = next_alpha(s);
		    p = s;
		    s = next_non_alpha(s);
		    *s = '\0';
		    stop_envir = p;
		    break;
                }
<EndTokenMode>"\n"	{
                    parameters[0] += yytext[0];
		    inc_line();
		    break;
		}
<EndTokenMode>{sep}     { break; /* ignore separator */ }
<EndTokenMode>.	{   parameters[0] += yytext[0];	break; }

<EndTokenMode>[\\]end{seps}[\{]{seps}{envir}{seps}[\}]	{
                    parameters[1] = yytext;
		    char* s = yytext + 5;
		    s = next_alpha(s);
		    char* p = s;
		    s = next_non_alpha(s);
		    *s = '\0';
		    if ( stop_envir == p) {
		        BEGIN( old_state);
			include_stack.push_string( in_string->name(), 
						   string("\\end{") + p + '}',
						   in_string->line());
			Macro_item item = fetchMacro( current_macro);
			include_stack.push_string( current_macro + " in "
						   +item.filename, 
						   expandMacro(current_macro,
							       item, 
							       parameters, 
							       1, 0),
						   item.line);
			current_macro = string();
			stop_envir    = string();
			parameters[0] = string();
		    } else {
			parameters[0] += parameters[1];
		    }
		    parameters[1] = string();
		    break;
		}

<INITIAL,AllttMode>[\\]lcAsciiToHtml  {
		    old_state = YY_START;
                    BeginParameterStart();
		    return ASCIITOHTML;
                }

<INITIAL,AllttMode>[\\]lciRawOutput  {
		    old_state = YY_START;
                    BeginParameterStart();
		    return RAWOUTPUT;
                }
<INITIAL,AllttMode>[\\]lciRawOutputN{seps}[\{]{seps}{number}{seps}[\}]  {
                    char* s = yytext + 1;
		    while ( *s && ! isdigit(*s))
			++s;
	            int n = atoi( s);
                    if ( n > 0) {
			s = new char[n+1];
			yylval.text = s;
			while (n--)
			    if ( (*(s++) = yyinput()) == '\n')
				inc_line();
			*s = '\0';
		        return RAWOUTPUTN;
                    }
		    break;
                }
<INITIAL,AllttMode>[\\]lciAsciiOutputN{seps}[\{]{seps}{number}{seps}[\}]  {
                    char* s = yytext + 1;
		    while ( *s && ! isdigit(*s))
			++s;
	            int n = atoi( s);
                    if ( n > 0) {
			s = new char[n+1];
			char* p = s;
			while (n--)
			    if ( (*(s++) = yyinput()) == '\n')
				inc_line();
			*s = '\0';
			yylval.text = convert_ascii_to_html( p);
			delete[] p;
		        return RAWOUTPUTN;
                    }
		    break;
                }
<INITIAL,AllttMode>[\\]lciRawSkipN{seps}[\{]{seps}{number}{seps}[\}]  {
                    char* s = yytext + 1;
		    while ( *s && ! isdigit(*s))
			++s;
	            int n = atoi( s);
                    if ( n > 0) {
			while (n--)
			    if ( yyinput() == '\n')
				inc_line();
                    }
		    break;
                }


 /* Index      */
 /* ---------------------- */
<INITIAL>[\\]lciOpenFileforIndex  {
    OpenFileforIndex();      
    break;
}

<INITIAL>[\\]lciCloseFileforIndex  {
    CloseFileforIndex();      
    break;
}

<INITIAL>[\\]lciIndex  {
    handleIndex();      
    break;
}

<INITIAL>[\\]lciIndexTraitsClass  {
    handleIndexTraitsClass();      
    break;
}

<INITIAL>[\\]lciIndexRefName  {
    handleIndexRefName();      
    break;
}

 /* Special multicharacter sequences */
 /* -------------------------------- */


"---"           {
		    yylval.text = " - ";
		    return STRING;
                }

"--"            {
		    yylval.text = "-";
		    return STRING;
                }


 /* Skip unwished tokens */
 /* -------------------- */
<INITIAL,AllttMode>[\\]lciSkipWhiteSpace {
	            skipspaces();
		    break;
                }
<INITIAL,AllttMode>[\\]lciSkipNonWhiteSpace {
	            skipspaces();
		    skipnonspaces();
		    break;
                }

 /* Single letter quotes that cannot be put in the style file */
 /* --------------------------------------------------------- */
<INITIAL,AllttMode>[\\][\{\}%]   {
		    yylval.character = yytext[1];
		    return CHAR;
                }
<INITIAL,AllttMode>[\\]{space}   {
		    yylval.character = ' ';
		    return CHAR;
		}


 /* TeX macro expansion                    */
 /* -------------------------------------- */
<INITIAL,AllttMode>{texmacro}       {
                   yyleng = removespaces( yytext);
                   if ( ! expand_macro()) {
                      return STRING;
                   }
                   break;
                }


 /* Specials TeX sequence macro expansion:    */
 /* see also below for active char expansion  */
 /* ----------------------------------------- */
<INITIAL>"$$"       {
                    yyleng = removespaces( yytext);
		    if ( ! expand_macro())
			return STRING;
		    break;
                }


 /* The rest: spaces and single characters */
 /* -------------------------------------- */
[\\]?{ws}	{
		    if ( *yytext == '\\')
		        yylval.text = yytext + 1;
		    else
		        yylval.text = yytext;
	  	    return STRING;
		}
<AllttMode>{ws}	{
	            yylval.text = yytext;
	  	    return STRING;
		}

  /* Speed up for the usual case of plain text, numbers and spaces */

<INITIAL,AllttMode>[a-zA-Z0-9 \t]+ {
	            yylval.text = yytext;
	  	    return STRING;
}

<INITIAL,AllttMode>[\{\}]    |
.	                     {
                    if ( is_active_char( yytext[0])) {
			if ( ! expand_macro())
			    return STRING;
		    }
		    else {
			yylval.character = yytext[0];
			if ( is_html_multi_character( yylval.character)) {
			    yylval.text = html_multi_character(
				yylval.character);
			    return STRING;
			}
			return CHAR;
		    }
		    break;
		}

<AllttMode>.	{
		    yylval.character = yytext[0];
		    if ( is_html_multi_character( yylval.character)) {
			yylval.text = html_multi_character( yylval.character);
			return STRING;
		    }
		    return CHAR;
		}

%%

/* returns true if EOF has been detected */
bool skipspaces_eof() {
    int c = yyinput();
    while( c != EOF && (isspace(c) || c == SEPARATOR)) {
	if ( c == '\n')
	    inc_line();
        c = yyinput();
    }
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return true;    
    } 
    yyunput( c, yytext );
    return false;
}

/* returns true if EOF has been detected */
bool skipseparator_eof() {
    int c = yyinput();
    while( c != EOF && c == SEPARATOR) {
        c = yyinput();
    }
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return true;    
    } 
    yyunput( c, yytext );
    return false;
}

/* returns true if EOF has been detected */
bool skiplimitedspaces_eof() {
    if ( YY_START == AllttMode)
	return false;
    int nl_count = 0;
    int c = yyinput();
    while( c != EOF && isspace(c) && c != SEPARATOR) {
	if ( c == '\n') {
	    if ( nl_count > 0) {
		yyunput( c, yytext );
		yyunput( c, yytext );
		return false;
	    } else
		nl_count++;
	}
        c = yyinput();
    }
    if ( nl_count > 0)
	inc_line();
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return true;    
    } 
    yyunput( c, yytext );
    return false;
}

/* returns last character read. If its EOF, terminate, */
/* else the char has been put back to the input stream. */
int skiplimitedspaces_param() {
    int nl_count = 0;
    int c = yyinput();
    while( c != EOF && isspace(c) && c != SEPARATOR) {
	if ( c == '\n') {
	    if ( nl_count > 0) {
		yyunput( c, yytext );
		yyunput( c, yytext );
		return '\n';
	    } else
		nl_count++;
	}
        c = yyinput();
    }
    if ( nl_count > 0)
	inc_line();
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return c;
    } 
    yyunput( c, yytext );
    return c;
}

/* returns true if EOF has been detected */
bool skiplimitedsepspaces_eof() {
    if ( YY_START == AllttMode)
	return skipseparator_eof();
    int nl_count = 0;
    int c = yyinput();
    while( c != EOF && (isspace(c) || c == SEPARATOR)) {
	if ( c == '\n') {
	    if ( nl_count > 0) {
		yyunput( c, yytext );
		yyunput( c, yytext );
		return false;
	    } else
		nl_count++;
	}
        c = yyinput();
    }
    if ( nl_count > 0)
	inc_line();
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return true;    
    } 
    yyunput( c, yytext );
    return false;
}

/* returns true if EOF has been detected */
bool skipnonspaces_eof() {
    int c = yyinput();
    while( c != EOF && ( ! isspace(c) || c == SEPARATOR)) {
        c = yyinput();
    }
    if (c == EOF) {
	printErrorMessage( EOFInMacroExpansionError);
	return true;    
    } 
    yyunput( c, yytext );
    return false;
}

size_t removespaces( char* s) {
    char* begin = s;
    char* p = s;
    while ( *s) {
	if ( isspace( *s) || *s == SEPARATOR)
	    ++s;
	else
	    *p++ = *s++;
    }
    *p = '\0';
    return p - begin;
}


void inc_line() {
    in_string->line() ++;
    if ( line_switch)
	cerr << in_string->name() << " line " << in_string->line() << endl;
}

void count_newlines( const char* s) {
    while ( *s) {
	if ( *s == '\n')
	    inc_line();
	++s;
    }
}

bool is_parameter_parsing_done( bool option, bool more_param) {
    BEGIN( old_state);
    if ( parameter_count == 0 && parameter_endopt == 0) { 
	// Evaluation strategy for scanner defined commands
	yylval.text = newstr( current_parameter.c_str());
	current_parameter = string();
	BeginParameterStart();
	return true;
    }
    // Evaluation strategy for user defined commands
    parameter_format += ( option ? 'o' : 'm');
    parameters[ parameter_index + parameter_options] = current_parameter;
    current_parameter = string();
    if ( option) {
	parameter_options++;
	if ( parameter_endopt > 0 && parameter_count == 0)
	    parameter_endopt--;
    } else {
	parameter_index++;
	parameter_count--;
    }
    if ((parameter_count == 0) && (parameter_endopt > 0)) {
	if ( more_param) {
	    parameter_nesting = 1;
	    parameter_option = true;
	    BeginParameterMode();
	    return false;
	}
	parameter_endopt = 0;
    }
    if ( parameter_count == 0 && parameter_endopt == 0) {
	string macro = current_macro;
	if ( parameter_options > 0) {
	    macro += '@';
	    macro += parameter_format;
	}
	Macro_item item = fetchMacro( macro);
	include_stack.push_string( current_macro + " in "+item.filename, 
				   expandMacro(macro, item, parameters, 
					       parameter_index, 
					       parameter_options),
				   item.line);
	current_macro = string();
	for ( int i = 0; i < parameter_index + parameter_options; ++i) {
	    parameters[i] = string();
	}
	parameter_index   = 0;
	parameter_options = 0;
	parameter_format  = string();
    } else
	BeginParameterStart();
    return false;
}

// returns true if all went well. Return-value == false means return STRING.
bool expand_macro() {
    if ( ! definedMacro( yytext)) {
        printErrorMessage( MacroUndefinedError, yytext);
	yylval.text = yytext;
	return false;    
    }
    Macro_item item = fetchMacro( yytext);
    if ( item.n_param > 0 || item.n_opt_at_end > 0) {
	parameter_count     = item.n_param;
	parameter_index     = 0;
	parameter_options   = 0;
	parameter_option    = false;
	old_state = YY_START;
	if ( item.n_param == 0) {
	    string s( yytext);
	    if ( yytext[yyleng-1] == '*' || isalpha( yytext[yyleng-1]))
		skiplimitedsepspaces_eof();
	    int c = yyinput();
	    if ( c == '[') {
		parameter_nesting = 1;
		parameter_option  = true;
		parameter_endopt  = item.n_opt_at_end;
		current_macro     = s;
		BeginParameterMode();
	    } else {
		if (c != EOF)
		    yyunput( c, yytext );
		include_stack.push_string( s + " in " + item.filename,
					   item.body, item.line);
	    }
	} else {
	    parameter_endopt    = item.n_opt_at_end;
	    current_macro       = yytext;
	    BeginParameterStart();
	}
    } else {
	string s( yytext);
	if ( yytext[yyleng-1] == '*' || isalpha( yytext[yyleng-1]))
	    skiplimitedspaces_eof();
	if ( item.fct) {
	    if ( macro_exp_switch) {
		string repl = item.fct( s, parameters, 0, 0);
		cerr << '`' << s << "' Internally replaced with `" 
		     << repl << "'" <<endl;
		include_stack.push_string( s + " in " + item.filename,
					   repl, item.line);

	    } else
		include_stack.push_string( s + " in " + item.filename,
					   item.fct( s, parameters, 0, 0),
					   item.line);
	    return true;
	}
	if ( macro_exp_switch)
	    cerr << '`' << s << "' Replaced with `" << item.body << "'" <<endl;
	include_stack.push_string( s + " in " + item.filename,
				   item.body, item.line);
    }
    return true;
}

extern "C" int yywrap() {
	include_stack.pop();
	if ( include_stack.empty())
	    return 1;
	return 0;
}


void printScannerState( ostream& out) {
    switch ( YY_START) {
    case INITIAL: 
	out << "INITIAL"; 
	break;
    case AllttMode: 
	out << "AllttMode"; 
	break;
    case IncludeMode: 
	out << "IncludeMode"; 
	break;
    case ParameterStart: 
	out << "ParameterStart"; 
	break;
    case ParameterMode: 
	out << "ParameterMode"; 
	break;
    case CCParameterMode: 
	out << "CCParameterMode"; 
	break;
    case DelimiterMode: 
	out << "DelimiterMode"; 
	break;
    case EndTokenMode: 
	out << "EndTokenMode"; 
	break;
    default:
	out << "<unknow mode>";
    }
}
