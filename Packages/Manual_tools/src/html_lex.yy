/**************************************************************************
 
  html_lex.yy
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : lexical scanner for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : flex, bison, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <html_lex.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <database.h>
#include <html_config.h>
#include <string_conversion.h>
#include <internal_macros.h>
#include <html_error.h>
#include <html_syntax.tab.h>

#include <lex_include.h>
#include <lex_include_impl.h>
#include <macro_dictionary.h>

int    old_state = 0;    // Old parser state before param. parsing.

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
// String to store the filename of ccHtmlClassFile intermediately
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

/* --------------------------------------------------------------------
    Parsing Modes:

    --  INITIAL:          parses (La)TeX

    --  AllttMode:        Parses alltt environments. Copies all characters
                          literally except \macros and {}.

    --  DelimiterMode     parses LaTeX \verb"..." statements until a
                          stop character occurs. Another example is \path.

    --  EndTokenMode      parses LaTeX verbose environments until an
                          \end{envir} occurs. Another example is cprog.

    --  IncludeMode       parses input/include filename,

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
filename        [^ \t\n\\\{\}\[\]()\001]+
space           [\t ]
w               {space}*
ws              {space}+
number          {digit}+
rmblockintro    ([\{][\\](rm))|([\\]((text)|(math))rm[\{])
ttblockintro    ([\{][\\](tt))|([\\]((text)|(math))tt[\{])
emblockintro    ([\{][\\](em))|([\\]emph[\{])
itblockintro    ([\{][\\]((it)|(sl)))|([\\]((text)|(math))((it)|(sl))[\{])
scblockintro    ([\{][\\](sc))|([\\]textsc[\{])
sfblockintro    ([\{][\\](sf))|([\\]((text)|(math))sf[\{])
bfblockintro    ([\{][\\]((bf)|(mathbold)))|([\\]((text)|(math))bf[\{])

%%
 /* Mode switching can be triggered from the parser */
 /* ----------------------------------------------- */
	if (set_old_state) {
	    BEGIN( old_state);
	    set_old_state = 0;
	}

 /* Separator, Newlines, paragraphs, comments, and EOF        */
 /* --------------------------------------------------------- */

<INITIAL,AllttMode>{sep}  { /* ignore separator */ }

 /* Count line numbers in all modes for better error messages */
{w}"\n"{w}("\n"{w})+   {  /* create a TeX paragraph */
                    count_newlines( yytext);
		    yylval.string.text = "\n<P>\n\n";
		    yylval.string.len  = -1;
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
	        }
"%".*[\n]{w}("\n"{w})+  { /* Match comments with an empty line -> par */
                    count_newlines( yytext);
		    yylval.string.text = "\n<P>\n\n";
		    yylval.string.len  = -1;
		    return STRING;
	        }
"%".*  { /* Match one line TeX comments at the last line in a file (EOF) */ }

<<EOF>>         {   if ( YY_START != INITIAL)
                        printErrorMessage( ParsingStateError);
                    yyterminate();
}

 /* Handle include files      */
 /* ------------------------- */

<INITIAL,AllttMode>[\\]((include)|(input)){seps}[\{]{seps}   {  
		    old_state = YY_START;
		    BEGIN ( IncludeMode); 
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
		    include_stack.push_tex_file( yytext);
}

 /* TeX Macro Definitions                         */
 /* LaTeX newcommand is as internal macro defined */
 /* --------------------------------------------- */

<INITIAL,AllttMode>[\\]lciBeginGroup    { return '{'; }
<INITIAL,AllttMode>[\\]lciEndGroup      { return '}'; }

<INITIAL,AllttMode>[\\]def{seps}{deftexmacro}{seps}([#][0-9])+ {
		    old_state = YY_START;
                    BEGIN( ParameterStart);
		    int i = 1;
		    while ( yytext[i] != '\\')
			++i;
	            yylval.string.text = newstr( yytext+i);
		    yylval.string.len  = yyleng-i;
		    i = 1;
		    while ( isalpha( yylval.string.text[i]))
			++i;
		    ((char *)(yylval.string.text))[i] = '\0';
                    return DEFWITHARGS;
}
<INITIAL,AllttMode>[\\]gdef{seps}{deftexmacro}{seps}([#][0-9])+ {
		    old_state = YY_START;
                    BEGIN( ParameterStart);
		    int i = 1;
		    while ( yytext[i] != '\\')
			++i;
	            yylval.string.text = newstr( yytext+i);
		    yylval.string.len  = yyleng-i;
		    i = 1;
		    while ( isalpha( yylval.string.text[i]))
			++i;
		    ((char *)(yylval.string.text))[i] = '\0';
                    return GDEFWITHARGS;
}
<INITIAL,AllttMode>[\\]g?def{seps}{deftexmacro}{seps}[#][^\{]* {
		    old_state = YY_START;
                    BEGIN( ParameterStart);
		    int i = 1;
		    while ( yytext[i] != '\\')
			++i;
	            yylval.string.text = newstr( yytext+i);
		    yylval.string.len  = yyleng-i;
                    return DEFWITHUNKNOWNARGS;
}

 /* Macro Parameter Parsing   */
 /* ------------------------- */

<ParameterStart>[\\]"\n" { 
		    inc_line();
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
}
<ParameterStart>"\n"     |
<ParameterStart>"%".*[\n]{w}  { /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    inc_line();
}
<ParameterStart>{sepsp}+ | /* ignore separator and spaces before parameters */
<ParameterStart>"%".*    { /* Match one line TeX comments */
                           /* at the last line in file */
}
<ParameterStart>"\{"     {
		    current_parameter = string();
                    ++parameter_nesting;
                    BEGIN( ParameterMode);
}
<ParameterStart>"\["     {
		    current_parameter = string();
                    ++parameter_nesting;
		    parameter_option = true;
                    BEGIN( ParameterMode);
}
<ParameterStart>{deftexmacro}  |
<ParameterStart>.              {
		    current_parameter = yytext;
		    IfParameterParsing( false)
		        return PARAMETER;
}


<ParameterMode>[\\]"\n" { 
		    current_parameter += yytext;
                    inc_line();
}
<ParameterMode>"%".*[\n]{w}("\n"{w})+  { 
                    /* Match comments with an empty line -> par */
		    current_parameter += "\n\n";
                    count_newlines( yytext);
}
<ParameterMode>"%".*[\n]{w}  { /* Match one line TeX comments */
		    /* remove spaces in next line  */
                    inc_line();
}
<ParameterMode>"%".*    { /* Match one line TeX comments */
                          /* at the last line in file */
}
<ParameterMode>"\n"     { /* keep newlines */
		    current_parameter += yytext;
                    inc_line();
}
<ParameterMode>"\{"     {
		    current_parameter += yytext;
                    ++parameter_nesting;
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
}
<ParameterMode>"]"      {
                    if ( parameter_nesting == 1 && parameter_option) {
			parameter_nesting--;
			parameter_option = false;
			IfParameterParsing( true)
			    return PARAMETER_OPTION;
		    } else
			current_parameter += yytext;
}
<ParameterMode>{deftexmacro}   |
<ParameterMode>[^\\\{\}%\]\n]+ |
<ParameterMode>.               {
		    current_parameter += yytext;
}

 /* C++ Parameter Parsing     */
 /* ------------------------- */
<CCParameterMode>[\\]"\n" |
<CCParameterMode>"\n"     {
                    inc_line();
		    cc_string += '\n';
}
<CCParameterMode>"\{"     {
                    if ( parameter_nesting++ > 0)
                        cc_string += yytext;
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
}

 /* Special treatment of several macros in C++ text */
<CCParameterMode>[\\][^a-zA-Z] { /* capture all quoted special symbols */
                    cc_string += yytext[1];
}
<CCParameterMode>[\\]"tt"{sepnls} { count_newlines(yytext); cc_string+="|T|";}
<CCParameterMode>[\\]"bf"{sepnls} { count_newlines(yytext); cc_string+="|B|";}
<CCParameterMode>[\\]"em"{sepnls} { count_newlines(yytext); cc_string+="|I|";}
<CCParameterMode>[\\]"it"{sepnls} { count_newlines(yytext); cc_string+="|I|";}
<CCParameterMode>[\\]"sl"{sepnls} { count_newlines(yytext); cc_string+="|I|";}
<CCParameterMode>[\\]"ccFont"{sepnls}       { count_newlines(yytext); 
                                              cc_string += "|I|"; }
<CCParameterMode>[\\]"l"?"dots"{sepnls}     { count_newlines(yytext);
                                              cc_string += "..."; }

<CCParameterMode>{deftexmacro}    { 
                    printErrorMessage( MacroInCModeError);
                    cc_string += yytext; 
}

<CCParameterMode>{sep}              { /* ignore separator */ }
<CCParameterMode>[^\n\{\}\\\001]+   |
<CCParameterMode>.                  {
                    cc_string += yytext;
}


<INITIAL,AllttMode>[\\]lciParseCC[\{][^\}]*[\}]   {
                    count_newlines(yytext);
		    yytext[ yyleng - 1] = '\0';
		    current_macro  = yytext + 12;
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
}

 /* Special Parser Modes      */
 /* ------------------------- */
[\\]lciBeginAlltt { 
                    skiplimitedspaces();
                    BEGIN( AllttMode); 
                }
<AllttMode>[\\]lciEndAlltt { 
                    skiplimitedspaces();
                    BEGIN( INITIAL);
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
                }
<DelimiterMode>"\n"	{
                    inc_line();
		    yymore();
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
                }
<EndTokenMode>"\n"	{
                    parameters[0] += yytext[0];
		    inc_line();
		}
<EndTokenMode>{sep}     { /* ignore separator */ }
<EndTokenMode>.	{   parameters[0] += yytext[0];	}

<EndTokenMode>[\\]end[\{]{seps}{envir}{seps}[\}]	{
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
		}


<INITIAL,AllttMode>[\\]lcAsciiToHtml  {
		    old_state = YY_START;
                    BEGIN( ParameterStart);
		    return ASCIITOHTML;
                }

<INITIAL,AllttMode>[\\]lciRawOutput  {
		    old_state = YY_START;
                    BEGIN( ParameterStart);
		    return RAWOUTPUT;
                }
<INITIAL,AllttMode>[\\]lciRawOutputN{seps}[\{]{seps}{number}{seps}[\}]  {
                    char* s = yytext + 1;
		    while ( *s && ! isdigit(*s))
			++s;
	            int n = atoi( s);
                    if ( n > 0) {
			s = new char[n+1];
			yylval.string.text = s;
			yylval.string.len  = n;
			while (n--)
			    if ( (*(s++) = yyinput()) == '\n')
				inc_line();
			*s = '\0';
		        return RAWOUTPUTN;
                    }
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
			yylval.string.text = convert_ascii_to_html( p);
			delete[] p;
		        return RAWOUTPUTN;
                    }
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
                }


 /* Chapter and labels triggering new file and linking */
 /* -------------------------------------------------- */
<INITIAL,AllttMode>[\\]lciChapter[*]?{seps}[\{]{seps}  {
		    return CHAPTER;
}
<INITIAL,AllttMode>[\\]lciSection[*]?{seps}("["[^\]]*"]")?[\{]{seps}  {
                    count_newlines(yytext);
		    return SECTION;
}
<INITIAL,AllttMode>[\\]label{seps}[\{][^\}]+[\}]  {
                    count_newlines(yytext);
	            char* s = yytext + 6;
	            while ( *s++ != '{')
                        ;
	            while ( isspace(*s) || *s == SEPARATOR)
			++s;
		    yylval.string.text = s;
		    s = yytext + yyleng - 2;
	            while( isspace( *s) || *s == SEPARATOR)
                        --s;
		    *s = '\0';
		    yylval.string.len  = strlen( yylval.string.text);
                    if ( yylval.string.len > 0)
		        return LABEL;
}


 /* Different keywords from the manual style triggering C++ formatting */
 /* ------------------------------------------------------------------ */
<INITIAL,AllttMode>[\\]begin{seps}[\{]{seps}lciClass{seps}[\}]   {
		    return BEGINCLASS;
		}
<INITIAL,AllttMode>[\\]end{seps}[\{]{seps}lciClass{seps}[\}]   {
		    return ENDCLASS;
		}

 /* Flexibility for HTML class files. */
 /* -------------------------------------------------------------- */
<INITIAL,AllttMode>[\\]lciStoreHtmlFileName   {
                    cc_filename = cc_string;
}
<INITIAL,AllttMode>[\\]lciHtmlFileNameBegin   {
		    return HTMLBEGINCLASSFILE;
}

 /* Grouping symbols */
 /* ---------------- */

<INITIAL,AllttMode>{ttblockintro}  {  /* TeX styles like {\tt ... */
                    skiplimitedspaces();
		    return TTBLOCKINTRO;
		}
<INITIAL,AllttMode>{emblockintro}  { skiplimitedspaces();return EMBLOCKINTRO;}
<INITIAL,AllttMode>{itblockintro}  { skiplimitedspaces();return ITBLOCKINTRO;}
<INITIAL,AllttMode>{scblockintro}  { skiplimitedspaces();return SCBLOCKINTRO;}
<INITIAL,AllttMode>{bfblockintro}  { skiplimitedspaces();return BFBLOCKINTRO;}
<INITIAL,AllttMode>{rmblockintro}  { skiplimitedspaces();return RMBLOCKINTRO;}
<INITIAL,AllttMode>{sfblockintro}  { skiplimitedspaces();return SFBLOCKINTRO;}


"---" {
		        yylval.string.text = " - ";
		        yylval.string.len  = 3;
			return STRING;
                }

"--" {
		        yylval.string.text = "-";
		        yylval.string.len  = 1;
			return STRING;
                }


 /* Skip unwished tokens */
 /* -------------------- */
<INITIAL,AllttMode>[\\]lciSkipWhiteSpace {
	                skipspaces();
                }
<INITIAL,AllttMode>[\\]lciSkipNonWhiteSpace {
	                skipspaces();
	                skipnonspaces();
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
		    if ( ! expand_macro())
			return STRING;
                }


 /* The rest: spaces and single characters */
 /* -------------------------------------- */
[\\]?{ws}	{
		    if ( *yytext == '\\') {
		        yylval.string.text = yytext + 1;
		        yylval.string.len  = yyleng - 1;
		    } else {
		        yylval.string.text = yytext;
		        yylval.string.len  = yyleng;
		    }
	  	    return STRING;
		}
<AllttMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
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
			    yylval.string.text = html_multi_character(
						     yylval.character);
			    yylval.string.len  = strlen( yylval.string.text);
			    return STRING;
			}
			return CHAR;
		    }
		}

<AllttMode>.	{
	            yylval.character = yytext[0];
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
    unput( c);
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
    unput( c);
    return false;
}

// returns true if fine, false if a paragraph has been encountered.
bool skiplimitedspaces_eof() {
    int nl_count = 0;
    int c = yyinput();
    while( c != EOF && isspace(c) && c != SEPARATOR) {
	if ( c == '\n') {
	    if ( nl_count > 0) {
		unput( c);
		unput( c);
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
    unput( c);
    return false;
}

// returns true if fine, false if a paragraph has been encountered.
bool skiplimitedsepspaces_eof() {
    int nl_count = 0;
    int c = yyinput();
    while( c != EOF && (isspace(c) || c == SEPARATOR)) {
	if ( c == '\n') {
	    if ( nl_count > 0) {
		unput( c);
		unput( c);
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
    unput( c);
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
    unput( c);
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
    if ( parameter_count == 0 && parameter_endopt == 0) { 
	// Evaluation strategy for scanner defined commands
	yylval.string.text = newstr( current_parameter.c_str());
	yylval.string.len  = current_parameter.size();
	current_parameter = string();
	BEGIN( ParameterStart);
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
	    BEGIN( ParameterMode);
	    return false;
	}
	parameter_endopt = 0;
    }
    if ( parameter_count == 0 && parameter_endopt == 0) {
	BEGIN( old_state);
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
	BEGIN( ParameterStart);
    return false;
}

// returns true if all went well. Return-value == false means return STRING.
bool expand_macro() {
    if ( ! definedMacro( yytext)) {
	if ( ! quiet_switch) {
	    cerr << endl << "Unknown macro " << yytext
		 << " in `" << in_string->name()
		 << " in line " << in_string->line() << "'.";
	    if ( stack_trace_switch)
		printErrorMessage( MacroUndefinedError);
	}
	yylval.string.text = yytext;
	yylval.string.len  = -1;
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
		BEGIN( ParameterMode);
	    } else {
		if (c != EOF)
		    unput( c);
		include_stack.push_string( s + " in " + item.filename,
					   item.body, item.line);
	    }
	} else {
	    parameter_endopt    = item.n_opt_at_end;
	    current_macro       = yytext;
	    BEGIN( ParameterStart);
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
