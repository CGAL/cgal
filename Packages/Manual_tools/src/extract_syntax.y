/**************************************************************************
 
  extract_syntax.y
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : grammatical parser for TeX and C++ code mixed files.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Declarations from lex.yy */
/* ======================== */
extern int set_CCMode;
extern int set_NestingMode;
extern int set_INITIAL;
extern int line_number;

extern const char* file_name;
extern       char* creationvariable;

extern char *yytext;

int yylex( void);
void init_scanner( FILE* in);


/* Declarations for the parser */
/* =========================== */
/* This variable flags for bison that we are in the CCMode */
int CCMode = 0;

/* count the brace nesting level for debug purposes */
int nesting = 0;

/* Own prototypes */
int yyerror( char *s);

/* Datastructures for the parser */
/* ============================= */
#include <buffer.h>
#include <config.h>
%}

%union {
    struct {
        char*      text;    /* a (meaningless) chunk of zero terminated text */
        int        len;     /* its length */
    }            string;
    char         character;   /* a (meaningless) character */
    Buffer*      pBuffer;     /* Buffer for collected strings and characters */
    Buffer_list* pText;
}

/* Elementary data types */
/* --------------------- */
%token <string>    STRING
%token <string>    SPACE
%token <character> CHAR
%token             NEWLINE

/* Keywords to trigger on */
/* ---------------------- */
%token             BEGINCLASS
%token             ENDCLASS
%token             BEGINREFCLASS
%token             ENDREFCLASS
%token             BEGINREFPAGE
%token             ENDREFPAGE
%token <string>    CREATIONVARIABLE
%token             CONSTRUCTOR
%token             METHOD
%token             FUNCTION
%token             FUNCTIONTEMPLATE
%token             VARIABLE
%token             TYPEDEF
%token             NESTEDTYPE
%token             ENUM
%token             STRUCT
%token             GLOBALFUNCTION
%token             GLOBALFUNCTIONTEMPLATE
%token             GLOBALVARIABLE
%token             GLOBALTYPEDEF
%token             GLOBALENUM
%token             GLOBALSTRUCT
%token             DECLARATION

/* Special action keywords */
/* ----------------------- */
%token             CCSTYLE
%token             CCSECTION
%token             CCSUBSECTION
%token             INCLUDE
%token             HEADING
%token             COMMENTHEADING
%token             GOBBLETHREEPARAMS
%token             GOBBLETWOPARAMS
%token             GOBBLEONEPARAM
%token             GOBBLEAFTERONEPARAM
%token             BLOCKINTRO
%token             NEWCOMMAND

/* handle LALR(1) restriction */
/* -------------------------- */
%token             LALRRESTRICTION

%type  <pBuffer>   string
%type  <pBuffer>   whitespace
%type  <pBuffer>   cc_stmts
%type  <pBuffer>   cc_stmt
%type  <pBuffer>   cc_stmts_skip_space
%type  <pBuffer>   comment_token
%type  <pBuffer>   declaration
%type  <pBuffer>   classname
%type  <pBuffer>   template_params

%type  <pText>     whitespaces
%type  <pText>     optional_whitespaces
%type  <pText>     comment_group
%type  <pText>     comment_sequence
%type  <pText>     compound_comment
%type  <pText>     full_comment_sequence  
%type  <pText>     non_empty_comment_sequence
%type  <pText>     nested_token_sequence
%type  <pText>     nested_token

%%
/* Grammar: Top Level */
/* ================== */

input:            /* empty */
                | input stmt
;

stmt:             string              { delete $1;}
                | SPACE               {}
                | NEWLINE
		| BEGINCLASS
		  classname           {   handleClass( $2->begin());
		                          delete $2;}
		    decl_sequence 
		  ENDCLASS
                                      {
					  handleClassEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| BEGINREFCLASS
		  classname           {   handleClass( $2->begin());
		                          handleRefPage( $2->begin());
		                          delete $2;}
		  decl_sequence
		  ENDREFCLASS
                                      {
                                          handleRefPageEnd();
					  handleClassEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| BEGINREFCLASS
		  '[' nested_token_sequence ']'
		  classname           {   handleClass( $5->begin());
		                          handleRefPage( $5->begin());
		                          delete_list( $3);
		                          delete $5;}
		  decl_sequence
		  ENDREFCLASS
                                      {
                                          handleRefPageEnd();
					  handleClassEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| BEGINREFPAGE
		  classname           {   handleRefPage( $2->begin());
		                          delete $2; }
		| BEGINREFPAGE
		  '[' nested_token_sequence ']'
		  classname           {   handleRefPage( $5->begin());
		                          delete_list( $3);
		                          delete $5; }
                | ENDREFPAGE          {   handleRefPageEnd(); }
		| CREATIONVARIABLE    {}
		| CCSTYLE  '{'  nested_token_sequence '}'  { set_INITIAL = 1; 
                                                             delete_list($3);
                                                           }
                | HEADING '{' comment_sequence '}'         { delete_list($3); }
                | INCLUDE '{' comment_sequence '}'         { delete_list($3); }
                | COMMENTHEADING '{' comment_sequence '}'  { delete_list($3); }
                | gobble_parameters
                | GOBBLEAFTERONEPARAM reduced_group reduced_group
                | global_tagged_declarator
		| group
;

group:            '{'           { nesting++; }
                  input
		  '}'           { nesting--; }
                | BLOCKINTRO    { nesting++; }
                  input
		  '}'           { nesting--; }
;

reduced_group:    '{'           { nesting++; }
                  reduced_sequence
		  '}'           { nesting--; }
                | BLOCKINTRO    { nesting++; }
                  reduced_sequence
		  '}'           { nesting--; }
;

reduced_sequence:   /* empty */
		 |  reduced_sequence reduced_statement
;

reduced_statement:
		  string  { delete $1;}
                | SPACE   {}
                | NEWLINE
		| CREATIONVARIABLE    {}
		| CCSTYLE  '{' nested_token_sequence '}'   { set_INITIAL = 1;
                                                             delete_list($3);
                                                           }
                | HEADING '{' comment_sequence '}'         { delete_list($3); }
                | INCLUDE '{' comment_sequence '}'         { delete_list($3); }
                | COMMENTHEADING '{' comment_sequence '}'  { delete_list($3); }
                | gobble_parameters
                | GOBBLEAFTERONEPARAM reduced_group reduced_group
		| reduced_group
;

/* Auxiliary Rules                 */
/* =============================== */
string:           STRING        {
                                  $$ = new Buffer;
				  $$->add( $1.text, $1.len);
                                }
                | CHAR          {
                                  $$ = new Buffer;
				  $$->add( $1);
                                }
                | string STRING {
		                  $1->add( $2.text, $2.len);
				  $$ = $1;
		                }
                | string CHAR   {
		                  $1->add( $2);
				  $$ = $1;
		                }
;

optional_whitespaces:  /* empty */    { $$ = new Buffer_list(); }
                  | whitespaces       { $$ = $1; }
;

whitespaces:        whitespace  { $$ = new Buffer_list( 1, $1); }
                  | GOBBLEAFTERONEPARAM comment_group comment_group {
                                  $$ = $2;
				  delete_list( $3);
		                }
                  | whitespaces whitespace {
                                  $$ = $1;
                                  $$->push_back( $2);
                                }
                  | whitespaces GOBBLEAFTERONEPARAM comment_group 
                    comment_group {
                                  $$ = $3;
				  delete_list( $1);
				  delete_list( $4);
		                }
;

whitespace:         SPACE       { $$ = new Buffer( $1.text, $1.len); }
                  | NEWLINE     { $$ = new Buffer( "\n", 1); }
		  | gobble_parameters
                                { $$ = new Buffer( " ", 1); }
;



/* Class Declaration with Comments */
/* =============================== */
decl_sequence:    comment_sequence  {
				  handleMainComment( * $1);
				  delete_list($1);
		                }
		| decl_sequence
		  tagged_declarator 
		  comment_sequence {
				  handleMainComment( * $3);
				  delete_list($3);
		                }
;

tagged_declarator:
		  global_tagged_declarator
		| CONSTRUCTOR   declaration   comment_group {
		                  handleFunctionDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
		| METHOD        declaration   comment_group {
		                  handleMethodDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
;

global_tagged_declarator:
		  error               {}
                | FUNCTION      declaration   comment_group {
		                  handleFunctionDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
		| FUNCTIONTEMPLATE 
                      template_params
                      optional_whitespaces
		      declaration 
		      comment_group {
		                  handleFunctionTemplateDeclaration(
					  $2->begin(),
					  $4->begin());
				  delete $2;
				  delete_list( $3);
				  delete $4;
				  handleComment( * $5);
				  delete_list( $5);
		                }
 		| VARIABLE      declaration   comment_group {
		                  handleDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
 		| TYPEDEF       declaration   comment_group {
		                  handleDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
 		| NESTEDTYPE    declaration   comment_group {
		                  handleNestedType( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
 		| ENUM          declaration   comment_group {
		                  handleDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
 		| STRUCT        declaration   comment_group {
		                  handleDeclaration( $2->begin());
				  delete $2;
				  handleComment( * $3);
				  delete_list( $3);
		                }
		| GLOBALFUNCTION      declaration   {
		                  handleFunctionDeclaration( $2->begin());
				  delete $2;
		                }
		| GLOBALFUNCTIONTEMPLATE 
		      template_params 
		      optional_whitespaces 
		      declaration  {
		                  handleFunctionTemplateDeclaration( 
					  $2->begin(),
					  $4->begin());
				  delete $2;
				  delete_list( $3);
				  delete $4;
		                }
 		| GLOBALVARIABLE      declaration   {
		                  handleDeclaration( $2->begin());
				  delete $2;
		                }
 		| GLOBALTYPEDEF       declaration   {
		                  handleDeclaration( $2->begin());
				  delete $2;
		                }
 		| GLOBALENUM          declaration   {
		                  handleDeclaration( $2->begin());
				  delete $2;
		                }
 		| GLOBALSTRUCT        declaration   {
		                  handleDeclaration( $2->begin());
				  delete $2;
		                }
 		| DECLARATION   declaration {
		                  handleDeclaration( $2->begin());
				  delete $2;
		                }
;

/* A sequence of words forming a comment */
/* ===================================== */
comment_group:      optional_whitespaces '{' comment_sequence '}'  { 
                                  $$ = $3; 
				  delete_list( $1);
                                }
                  | optional_whitespaces BLOCKINTRO comment_sequence '}'  { 
                                  $$ = $3; 
				  delete_list( $1);
                                }
;

comment_sequence:   optional_whitespaces { $$ = new Buffer_list(); 
                                           delete_list( $1);
                                         }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces { $$ = $2; 
	                                   delete_list( $1);
	                                   delete_list( $3);
                                         }
;

full_comment_sequence:   /* empty */  { $$ = new Buffer_list(); }
                  | whitespaces       { $$ = $1; }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces  { 
		                  $$ = $1;
				  $$->insert($$->end(), $2->begin(),$2->end());
				  $$->insert($$->end(), $3->begin(),$3->end());
				  delete $2;
				  delete $3;
		                }
;

non_empty_comment_sequence:
		    comment_token     { $$ = new Buffer_list( 1, $1); }
                  | compound_comment  { $$ = $1; }
                  | non_empty_comment_sequence optional_whitespaces 
                    comment_token {
		                  $$ = $1;
				  $$->insert($$->end(), $2->begin(),$2->end());
				  delete $2;
				  $$->push_back( $3);
		                }
                  | non_empty_comment_sequence 
		    optional_whitespaces
		    compound_comment {
		                  $$ = $1;
				  $$->insert($$->end(), $2->begin(),$2->end());
				  delete $2;
				  $$->insert($$->end(), $3->begin(),$3->end());
				  delete $3;
		                }
;

comment_token:      string   { $$ = $1; }
                  | comment_token string {
		                  $$ = $1;
				  $$->add( * $2);
				  delete $2;
		                }
;

compound_comment:   '{' full_comment_sequence '}' {
				  $$ = $2;
		                  $$->push_front( new Buffer( "{", 1));
		                  $$->push_back(  new Buffer( "}", 1));
		                }
                  | BLOCKINTRO full_comment_sequence '}' {
				  $$ = $2;
		                }
                  | CCSTYLE '{' nested_token_sequence '}'  { 
				  $$ = $3;
		                  set_INITIAL = 1;
                                  if ( $$->front()->is_space())  // should not
		                      $$->push_front( new Buffer( "`", 1));
                                  else
		                      $$->front()->prepend( "`");
                                  if ( $$->back()->is_space())
		                      $$->push_back( new Buffer( "'", 1));
                                  else
		                      $$->back()->add( "'");
		  }
                  | CCSECTION '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->push_front( new Buffer( " ", 1));
		                  $$->push_front( new Buffer( "SECTION:"));
		                  $$->push_front( new Buffer( "\n", 1));
		                  $$->push_front( new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( 
				      "===================================="
				      "===================================="));
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                }
                  | CCSUBSECTION '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->push_front( new Buffer( " ", 1));
		                  $$->push_front( new Buffer( "Subsection:"));
		                  $$->push_front( new Buffer( "\n", 1));
		                  $$->push_front( new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( 
				      "------------------------------------"
				      "------------------------------------"));
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                }
                  | HEADING '{' comment_sequence '}'  {
				  $$ = $3;
		                }
                  | INCLUDE '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->push_front( new Buffer( "<"));
		                  $$->push_front( new Buffer( " ", 1));
		                  $$->push_front( new Buffer( "#include"));
		                  $$->push_back(  new Buffer(  ">"));
		                }
                  | COMMENTHEADING '{' comment_sequence '}'  {
				  $$ = $3;
		                }
		  | CREATIONVARIABLE   {
                                  $$ = new Buffer_list();
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer(
					"New creation variable is:"));
		                  $$->push_back(  new Buffer( " ", 1));
                                  Buffer* t = new Buffer( $1.text, $1.len);
		                  t->prepend( "`");
		                  t->add(     "'");
		                  $$->push_back(  t);
		                  $$->push_back(  new Buffer( "\n", 1));
		                  $$->push_back(  new Buffer( "\n", 1));
		                }
;

/* Parsing of a C++ expression/statement with nested expressions */
/* ============================================================= */
nested_token_sequence:
		    /* empty */ { $$ = new Buffer_list(); }
		  | nested_token_sequence nested_token
		                {
				  $1->insert($1->end(), $2->begin(),$2->end());
				  delete $2;
				  $$ = $1;
				}
;

nested_token:       string      { $$ = new Buffer_list( 1, $1); }
                  | SPACE       {
                                  $$ = new Buffer_list( 1, 
					     new Buffer( $1.text, $1.len));
		  }
                  | NEWLINE     {
                                  $$ = new Buffer_list( 1, 
							new Buffer( "\n", 1));
                                }
		  | '{' nested_token_sequence '}' {
		                  $2->push_front( new Buffer( "{", 1));
		                  $2->push_back(  new Buffer( "}", 1));
				  $$ = $2;
		                }
		  | BLOCKINTRO nested_token_sequence '}' {
				  $$ = $2;
		                }
		  | '[' nested_token_sequence ']' {
		                  $2->push_front( new Buffer( "[", 1));
		                  $2->push_back(  new Buffer( "]", 1));
				  $$ = $2;
		                }
		  | '(' nested_token_sequence ')' {
		                  $2->push_front( new Buffer( "(", 1));
		                  $2->push_back(  new Buffer( ")", 1));
				  $$ = $2;
		                }
                  | CCSTYLE '{' nested_token_sequence '}'  { 
				  $$ = $3;
		                  set_INITIAL = 1;
                                  if ( $$->front()->is_space())  // should not
		                      $$->push_front( new Buffer( "`", 1));
                                  else
		                      $$->front()->prepend( "`");
                                  if ( $$->back()->is_space())
		                      $$->push_back( new Buffer( "'", 1));
                                  else
		                      $$->back()->add( "'");
                                }
;

/* Parsing of a C++ Declaration (function, method ..., not class) */
/* ============================================================== */
declaration:      '{'           { nesting++; 
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  nesting--;
		                  set_INITIAL = 1;
				  CCMode = 0;
				  $$ = $3;
		                }
;

classname:        '{'           { nesting++; 
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  nesting--;
		                  set_INITIAL = 1;
				  CCMode = 0;
				  $$ = $3;
		                }
;

template_params: '{'           { nesting++; 
                                  CCMode = 1; 
                                }
                  cc_stmts_skip_space
		  '}'           { 
		                  nesting--;
		                  set_INITIAL = 1;
				  CCMode = 0;
				  $$ = $3;
		                }
;

cc_stmts:         /* empty */
                                { $$ = new Buffer;}
		  | cc_stmts cc_stmt {
		                  $$ = $1;
				  $$->add( * $2);
				  delete $2;
		                }
;

cc_stmt:          string        { $$ = $1;
                                }
                | SPACE         { $$ = new Buffer;
                                  $$->add( ' ');
		                }
                | NEWLINE       { $$ = new Buffer;
                                  $$->add( ' ');
		                }
		| '{'           { nesting++; }
                  cc_stmts
		  '}'           { nesting--; 
		                  $$ = $3;
				  $$->prepend( '{');
				  $$->add( '}');
		                }
;

cc_stmts_skip_space: 
		  /* empty */
                                { $$ = new Buffer;}
		| string
		  cc_stmts      { $$ = $1;
				  $$->add( * $2);
				  delete $2;
		                }
                | SPACE         cc_stmts { $$ = $2;}
                | NEWLINE       cc_stmts { $$ = $2;}
		| '{'           { nesting++;}
                  cc_stmts
		  '}'           { nesting--;}
                  cc_stmts      { 
		                  $$ = $3;
				  $$->prepend( '{');
				  $$->add( '}');
				  $$->add( * $6);
				  delete $6;
		                }
;

/* Parse over non useful parameters */
/* -------------------------------- */
gobble_parameters:
                  GOBBLETHREEPARAMS reduced_group reduced_group reduced_group
                | GOBBLETWOPARAMS   reduced_group reduced_group
                | GOBBLEONEPARAM    reduced_group
                | GOBBLEONEPARAM    CHAR
                | NEWCOMMAND        reduced_group reduced_group



/* End if Grammar */
/* ============== */
%%

int yyerror( char *s) {
    fprintf( stderr,
	     "error 1 in line %d in %s: in %s-code: nesting %d: %s.\n", 
	     line_number,
	     file_name,
	     (CCMode ? "CC" : "TeX"),
	     nesting,
	     s);
    return 0;
}

// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in parser.y

const char* errorMessage( ErrorNumber n) {
    switch ( n) {
    case NoError:
	return "NO ERROR";
    case ParseError:
	return "parse error";
    case VariableUsedError:
	return "The creationvariable was used but not defined";
    case ClassnameUsedError:
	return "The classname was used but not defined";
    case RefNameUsedError:
	return "The ref-name was used but not defined";
    case TemplateParamExpectedError:
        return "A template parameter is missing";
    case MalformedFunctionDeclaration:
        return "The function declaration is malformed";
    case MalformedTemplateParamError:
        return "The template parameter is malformed (<> nesting ..)";
    case SemicolonMissingError:
        return "The declaration does not end in a semicolon";
    }
    return "UNKNOWN ERROR MESSAGE NUMBER";
}

void  printErrorMessage( ErrorNumber n){
    cerr << "error " << int(n) << " in line " << line_number << " in `" 
         << file_name << "': " << errorMessage( n) << "." << endl;
}
















