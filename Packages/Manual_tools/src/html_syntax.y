/**************************************************************************
 
  html_syntax.y
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : grammatical parser for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <html_syntax.h>
#include <html_lex.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <string_conversion.h>
#include <html_error.h>
#include <lex_include.h>
#include <macro_dictionary.h>
#include <cpp_formatting.h>
#include <internal_macros.h>

#include <database.h>
#include <html_config.h>


/* Own prototypes */
/* ============== */
int yyerror( char *s);

%}

%union {
    const char* text;        /* a chunk of zero terminated text */
    char        character;   /* a character */
    int         number;      /* an integer  */
    TextToken*  pTextToken;
    Text*       pText;
}

/* Elementary data types */
/* --------------------- */
%token <text>      STRING
%token <character> CHAR

/* Newcommand and parameter parsing */
/* -------------------------------- */
%token <text>      DEFWITHARGS
%token <text>      GDEFWITHARGS
%token <text>      DEFWITHUNKNOWNARGS
%token <text>      PARAMETER
%token <text>      PARAMETER_OPTION

/* Special Scanning Modes */
/* ---------------------- */
%token             ASCIITOHTML
%token             RAWOUTPUT
%token <text>      RAWOUTPUTN

/* Special TeX Parsing    */
/* ---------------------- */
%token             CHAPTER

/* C++ Parsing            */
/* ---------------------- */
%token             BEGINCLASS
%token             ENDCLASS

/* File Management        */
/* ---------------------- */
%token             HTMLBEGINCLASSFILE


/* Type declarations for production rules       */
/* -------------------------------------------- */
%type  <pTextToken> string_token 
%type  <pText>      comment_group  comment_sequence  comment_token



%%
/* Grammar: Top Level */
/* ================== */

input:            /* empty */
                | input stmt
;

stmt:             error               {}
                | direct_string       {}
                | CHAPTER '{' comment_sequence '}'  {
		                          handleChapter( * $3); 
                                          delete $3;
                                      }
		| BEGINCLASS          {   handleClassEnvironment(); }
		| ENDCLASS            {   handleClassEnd(); }
		| HTMLBEGINCLASSFILE comment_group 
		                      {  handleHtmlClassFile( cc_filename,*$2);
		                          delete $2;
                                      }
		| '{' input '}'
;


/* Auxiliary Rules                 */
/* =============================== */

direct_string:    STRING        { handleString( $1); }
                | CHAR          { handleChar( $1); }
                | ASCIITOHTML PARAMETER  {
		                  char* s = convert_ascii_to_html($2);
				  handleString( s);
				  delete[] $2;
				  delete[] s;
				  set_old_state = 1;
				}
                | RAWOUTPUT PARAMETER  {
				  handleString( $2);
				  delete[] $2;
				  set_old_state = 1;
				}
                | RAWOUTPUTN   {
				  handleString( $1);
				  delete[] $1;
				  set_old_state = 1;
				}
                | DEFWITHARGS  PARAMETER {
				  if ( number_of_args < 1 ||
				       number_of_args > 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertMacro( $1,in_string->name(),
						   in_string->line(), 
						   $2, number_of_args);
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
                | GDEFWITHARGS  PARAMETER {
				  if ( number_of_args < 1 ||
				       number_of_args > 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertGlobalMacro( $1,
							 in_string->name(),
							 in_string->line(), 
							 $2,
							 number_of_args);
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
                | DEFWITHUNKNOWNARGS  PARAMETER {
				  if ( ! quiet_switch) {
				      cerr << endl 
					   << "Unknown macro definition: " 
					   << $1 << " in `" 
					   << in_string->name() << " in line "
					   << in_string->line() << "'.";
				      printErrorMessage( MacroDefUnknownError);
				  }
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
;


string_token:     STRING        { $$ = new TextToken( $1); }
                | CHAR          { $$ = new TextToken(); $$->add( $1); }
                | ASCIITOHTML PARAMETER  {
		                  char* s = convert_ascii_to_html($2);
				  $$ = new TextToken( s);
				  delete[] $2;
				  delete[] s;
				  set_old_state = 1;
				}
                | RAWOUTPUT PARAMETER  {
				  $$ = new TextToken( $2);
				  delete[] $2;
				  set_old_state = 1;
				}
                | RAWOUTPUTN   {
				  $$ = new TextToken( $1);
				  delete[] $1;
				  set_old_state = 1;
				}
                | DEFWITHARGS  PARAMETER {
				  $$ = new TextToken();
				  if ( number_of_args < 1 ||
				       number_of_args> 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertMacro( $1,in_string->name(),
						   in_string->line(), 
						   $2, number_of_args);
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
                | GDEFWITHARGS  PARAMETER {
				  $$ = new TextToken();
				  if ( number_of_args < 1 ||
				       number_of_args > 9)
				      printErrorMessage( NParamRangeError);
				  else
				      insertGlobalMacro( $1,
							 in_string->name(),
							 in_string->line(), 
							 $2,
							 number_of_args);
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
                | DEFWITHUNKNOWNARGS  PARAMETER {
				  $$ = new TextToken();
				  if ( ! quiet_switch) {
				      cerr << endl 
					   << "Unknown macro definition: " 
					   << $1 << " in `" 
					   << in_string->name() << " in line "
					   << in_string->line() << "'.";
				      printErrorMessage( MacroDefUnknownError);
				  }
				  delete[] $1;
				  delete[] $2;
				  set_old_state = 1;
				}
;


/* A sequence of words forming a comment */
/* ===================================== */
comment_group:      '{' comment_sequence '}'  { $$ = $2; }
;

comment_sequence:        /* empty */  { $$ = new Text( managed); }
                       | comment_sequence comment_token {
		                  $$ = $1;
				  $$->append( * $2);
				  delete $2;
		                }
;

comment_token:      string_token  { $$ = new Text( * $1, managed); }
                  | comment_group { $$ = $1; }
;


/* End of Grammar */
/* ============== */
%%

/* EOF */
