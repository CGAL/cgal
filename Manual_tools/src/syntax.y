/**************************************************************************
 
  html_syntax.y
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : grammatical parser for TeX and C++ code mixed files.
              Taylored for HTML manual generation.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1996 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

%{
#include <syntax.h>
#include <lexer.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <string_conversion.h>
#include <error.h>
#include <config.h>
#include <input.h>
#include <macro_dictionary.h>
#include <cpp_formatting.h>
#include <internal_macros.h>



/* Own prototypes */
/* ============== */
int yyerror( char *s);

%}

%union {
    const char*   text;        /* a chunk of zero terminated text */
    char          character;   /* a character */
    int           number;      /* an integer  */
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


%%
/* Grammar: Top Level */
/* ================== */

input:    /* empty */
        | input stmt
;

stmt:     error         {}
        | '{' input '}' {}
        | STRING        { handleString( $1); }
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
                              cerr << endl << ErrorText
                                   << ": Unknown macro definition: " << $1 
                                   << " in `" << in_string->name() 
                                   << " in line " << in_string->line() << "'." 
                                   << ResetColor;
                              printErrorMessage( MacroDefUnknownError);
                          }
                          delete[] $1;
                          delete[] $2;
                          set_old_state = 1;
                        }
;


/* End of Grammar */
/* ============== */
%%

/* EOF */
