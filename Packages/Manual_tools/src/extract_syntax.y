/**************************************************************************
 
  extract_syntax.y
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : grammatical parser for TeX and C++ code mixed files.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
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
#include <database.h>
#include <config.h>
%}

%union {
    struct {
        char*      text;    /* a (meaningless) chunk of zero terminated text */
        int        len;     /* its length */
    }          string;
    char       character;   /* a (meaningless) character */
    Buffer*    pBuffer;     /* Buffer for collected strings and characters */
    TextToken* pTextToken;
    Text*      pText;
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
%token             BEGINCLASSTEMPLATE
%token             ENDCLASSTEMPLATE
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
%token             BEGINTEXONLY
%token             ENDTEXONLY
%token             BEGINHTMLONLY
%token             ENDHTMLONLY

/* Special action keywords */
/* ----------------------- */
%token             CCSTYLE
%token             CCSECTION
%token             CCSUBSECTION
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

%type  <pBuffer>    string
%type  <pBuffer>    declaration      classname          template_params
%type  <pBuffer>    cc_stmts         cc_stmt            cc_stmts_skip_space
%type  <pText>      comment_group    comment_sequence  
%type  <pText>      nested_token_sequence               nested_token
%type  <pText>      compound_comment                    full_comment_sequence  
%type  <pText>      non_empty_comment_sequence

%type  <pText>      whitespaces      optional_whitespaces

%type  <pTextToken> comment_token    non_empty_token    whitespace



%%
/* Grammar: Top Level */
/* ================== */

input:            /* empty */
                | input stmt
;

stmt:             string  { delete $1;}
                | SPACE   {}
                | NEWLINE
		| BEGINCLASS
		  classname           {   handleClass( $2->string());
		                          delete $2;}
		    decl_sequence 
		  ENDCLASS
                                      {
					  handleClassEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| BEGINCLASSTEMPLATE
		  classname           {   handleClassTemplate( $2->string());
		                          delete $2;}
		  decl_sequence
		  ENDCLASSTEMPLATE
                                      {
					  handleClassTemplateEnd();
                                          free( creationvariable);
                                          creationvariable = NULL;
				      }
		| CREATIONVARIABLE    {}
		| CCSTYLE  '{'  nested_token_sequence '}'  { set_INITIAL = 1; 
                                                             delete $3;
                                                           }
                | HEADING '{' comment_sequence '}'         { delete $3; }
                | COMMENTHEADING '{' comment_sequence '}'  { delete $3; }
		| BEGINTEXONLY  nested_token_sequence ENDTEXONLY  {
                                                             delete $2;
                                                           }
		| BEGINHTMLONLY nested_token_sequence ENDHTMLONLY  { 
                                                             set_INITIAL = 1;
                                                             delete $2;
                                                           }
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
                                                             delete $3;
                                                           }
                | HEADING '{' comment_sequence '}'         { delete $3; }
                | COMMENTHEADING '{' comment_sequence '}'  { delete $3; }
		| BEGINTEXONLY  nested_token_sequence ENDTEXONLY   {
                                                             delete $2;
                                                           }
		| BEGINHTMLONLY nested_token_sequence ENDHTMLONLY  { 
                                                             set_INITIAL = 1;
                                                             delete $2;
                                                           }
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

non_empty_token:    string      { $$ = new TextToken( 
					       $1->string(), $1->length());
                                  delete $1;
                                }
;

optional_whitespaces:  /* empty */    { $$ = new Text( managed); }
                  | whitespaces       { $$ = $1; }
;

whitespaces:        whitespace  { $$ = new Text( * $1, managed); }
                  | GOBBLEAFTERONEPARAM comment_group comment_group {
                                  $$ = $2;
				  delete $3;
		                }
                  | whitespaces whitespace {
                                  $$ = $1;
                                  $$->append( * $2);
                                }
                  | whitespaces GOBBLEAFTERONEPARAM comment_group comment_group {
                                  $$ = $3;
				  delete $1;
				  delete $4;
		                }
;

whitespace:         SPACE       { $$ = new TextToken( $1.text, $1.len, true); }
                  | NEWLINE     { $$ = new TextToken( "\n", 1, true); }
		  | gobble_parameters
                                { $$ = new TextToken( " ", 1, true); }
;



/* Class Declaration with Comments */
/* =============================== */
decl_sequence:    comment_sequence  {
				  handleMainComment( * $1);
				  delete $1;
		                }
		| decl_sequence
		  tagged_declarator 
		  comment_sequence {
				  handleMainComment( * $3);
				  delete $3;
		                }
;

tagged_declarator:
		  global_tagged_declarator
		| CONSTRUCTOR   declaration   comment_group {
		                  handleFunctionDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
		| METHOD        declaration   comment_group {
		                  handleMethodDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
;

global_tagged_declarator:
		  FUNCTION      declaration   comment_group {
		                  handleFunctionDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
		| FUNCTIONTEMPLATE 
                      template_params
                      optional_whitespaces
		      declaration 
		      comment_group {
		                  handleFunctionTemplateDeclaration(
					  $2->string(),
					  $4->string());
				  delete $2;
				  delete $3;
				  delete $4;
				  handleComment( * $5);
				  delete $5;
		                }
 		| VARIABLE      declaration   comment_group {
		                  handleDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
 		| TYPEDEF       declaration   comment_group {
		                  handleDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
 		| NESTEDTYPE    declaration   comment_group {
		                  handleNestedType( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
 		| ENUM          declaration   comment_group {
		                  handleDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
 		| STRUCT        declaration   comment_group {
		                  handleDeclaration( $2->string());
				  delete $2;
				  handleComment( * $3);
				  delete $3;
		                }
		| GLOBALFUNCTION      declaration   {
		                  handleFunctionDeclaration( $2->string());
				  delete $2;
		                }
		| GLOBALFUNCTIONTEMPLATE 
		      template_params 
		      optional_whitespaces 
		      declaration  {
		                  handleFunctionTemplateDeclaration( 
					  $2->string(),
					  $4->string());
				  delete $2;
				  delete $3;
				  delete $4;
		                }
 		| GLOBALVARIABLE      declaration   {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
 		| GLOBALTYPEDEF       declaration   {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
 		| GLOBALENUM          declaration   {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
 		| GLOBALSTRUCT        declaration   {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
 		| DECLARATION   declaration {
		                  handleDeclaration( $2->string());
				  delete $2;
		                }
;

/* A sequence of words forming a comment */
/* ===================================== */
comment_group:      optional_whitespaces '{' comment_sequence '}'  { 
                                  $$ = $3; 
				  delete $1;
                                }
                  | optional_whitespaces BLOCKINTRO comment_sequence '}'  { 
                                  $$ = $3; 
				  delete $1;
                                }
;

comment_sequence:   optional_whitespaces { $$ = new Text(managed); 
                                           delete $1;
                                         }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces { $$ = $2; 
	                                   delete $1;
	                                   delete $3;
                                         }
;

full_comment_sequence:   /* empty */  { $$ = new Text(managed); }
                  | whitespaces       { $$ = $1; }
                  | optional_whitespaces
                    non_empty_comment_sequence
                    optional_whitespaces  { 
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
;

non_empty_comment_sequence:
		    comment_token     { $$ = new Text( * $1, managed); }
                  | compound_comment  { $$ = $1; }
                  | non_empty_comment_sequence optional_whitespaces comment_token {
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
                  | non_empty_comment_sequence 
		    optional_whitespaces
		    compound_comment {
		                  $$ = $1;
				  $$->append( * $2);
				  $$->append( * $3);
		                }
;

comment_token:      non_empty_token   { $$ = $1; }
                  | comment_token non_empty_token {
		                  $$ = $1;
				  $$->add( * $2);
				  delete $2;
		                }
;

compound_comment:   '{' full_comment_sequence '}' {
				  $$ = $2;
		                  $$->cons(   *new TextToken( "{", 1));
		                  $$->append( *new TextToken( "}", 1));
		                }
                  | BLOCKINTRO full_comment_sequence '}' {
				  $$ = $2;
		                }
                  | CCSTYLE '{' nested_token_sequence '}'  { 
				  $$ = $3;
		                  set_INITIAL = 1;
                                  if ( $$->head().isSpace)  // should not
		                      $$->cons(   *new TextToken( "`", 1));
                                  else
		                      $$->head().prepend( "`");
                                  InListFIter< TextToken> ix( * $$);
				  ForAll( ix) {
				      if ( ix.isLast())
					  if ( ix->isSpace)
					      $$->append( *new TextToken(
                                                                  "'", 1));
					  else
					      ix->add( "'");
				  }
                                }
		  | BEGINTEXONLY  nested_token_sequence ENDTEXONLY   {
                                  $$ = $2;
                                }
		  | BEGINHTMLONLY nested_token_sequence ENDHTMLONLY  { 
                                  set_INITIAL = 1;
				  delete $2;
                                }
                  | CCSECTION '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->cons(   *new TextToken( " ", 1, true));
		                  $$->cons(   *new TextToken( "SECTION:"));
		                  $$->cons(   *new TextToken( "\n", 1, true));
		                  $$->cons(   *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( 
				      "===================================="
				      "===================================="));
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                }
                  | CCSUBSECTION '{' comment_sequence '}'  {
				  $$ = $3;
		                  $$->cons(   *new TextToken( " ", 1, true));
		                  $$->cons(   *new TextToken( "Subsection:"));
		                  $$->cons(   *new TextToken( "\n", 1, true));
		                  $$->cons(   *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( 
				      "------------------------------------"
				      "------------------------------------"));
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                }
                  | HEADING '{' comment_sequence '}'  {
				  $$ = $3;
		                }
                  | COMMENTHEADING '{' comment_sequence '}'  {
				  $$ = $3;
		                }
		  | CREATIONVARIABLE   {
                                  $$ = new Text( managed);
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken(
					"New creation variable is:"));
		                  $$->append( *new TextToken( " ", 1, true));
                                  TextToken* t = new TextToken($1.text,$1.len);
		                  t->prepend( "`");
		                  t->add(     "'");
		                  $$->append( *t);
		                  $$->append( *new TextToken( "\n", 1, true));
		                  $$->append( *new TextToken( "\n", 1, true));
		                }
;

/* Parsing of a C++ expression/statement with nested expressions */
/* ============================================================= */
nested_token_sequence:
		    /* empty */ {
		                  $$ = new Text(managed);
		                }
		  | nested_token_sequence nested_token
		                {
				  $1->append( * $2);
				  $$ = $1;
				}
;

nested_token:       string      {
                                  $$ = new Text(*new TextToken( 
						    $1->string(),
						    $1->length()),
						managed);
				  delete $1;
                                }
                  | SPACE       {
                                  $$ = new Text(*new TextToken( 
						    $1.text,
						    $1.len,
						    true),
						managed);
		  }
                  | NEWLINE     {
                                  $$ = new Text(*new TextToken( "\n", 1, true),
						managed);
                                }
		  | '{' nested_token_sequence '}' {
		                  $2->cons(   *new TextToken( "{", 1));
		                  $2->append( *new TextToken( "}", 1));
				  $$ = $2;
		                }
		  | BLOCKINTRO nested_token_sequence '}' {
				  $$ = $2;
		                }
		  | '[' nested_token_sequence ']' {
		                  $2->cons(   *new TextToken( "[", 1));
		                  $2->append( *new TextToken( "]", 1));
				  $$ = $2;
		                }
		  | '(' nested_token_sequence ')' {
		                  $2->cons(   *new TextToken( "(", 1));
		                  $2->append( *new TextToken( ")", 1));
				  $$ = $2;
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
				  $$->add( $2);
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
				  $$->add( $2);
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
				  $$->add( $6);
				  delete $6;
		                }
;

/* Parse over non useful parameters */
/* -------------------------------- */
gobble_parameters:
                  GOBBLETHREEPARAMS reduced_group reduced_group reduced_group
                | GOBBLETWOPARAMS   reduced_group reduced_group
                | GOBBLEONEPARAM    reduced_group
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
    cerr << "error " << n << " in line " << line_number << " in `" 
         << file_name << "': " << errorMessage( n) << "." << endl;
}
















