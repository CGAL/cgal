/**************************************************************************
 
  extract_lex.yy
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : lexical scanner for TeX and C++ code mixed files.
  System    : flex, bison, C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

%{
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <buffer.h>
#include <config.h>
#include <extract_syntax.tab.h>

/* Set this flag to 1 to switch immediately to CCMode. */
int set_CCMode      = 0;
/* Set this flag to 1 to switch immediately to NestingMode. */
int set_NestingMode = 0;
/* Set this flag to 1 to switch back to INITIAL. */
int set_INITIAL     = 0;

/* Tag to mark the unchecked keyword */
int   unchecked_tag = 0;

/* Count the linenumber for better errormessages. */
int   line_number  = 1;

/* store the creationvariable */
char* creationvariable = NULL;

/* remember the necessary stop character for \verb"..." */
char stop_character;

/* The classname and other state variables */
char* global_classname       = 0;
char* global_template_params = 0;
char* global_ref_name        = 0;

/* Hack, to get rid of the yywrap. */
#define YY_SKIP_YYWRAP 1
#define yywrap() 1

void skipspaces( void);

#define YY_BREAK  /* a do nothing */

%}

/* The normal scanning mode parses TeX conventions.      */
/* In CCMode, it parses C++ conventions.                 */
/* The NestingMode parses only (){}[] nested expressions */
/* The VerbMode parses LaTeX \verb"..." statements as    */
/*     a sequence of characters                          */
%x CCMode
%x NestingMode
%x VerbMode
%x VerbatimMode
%x HtmlBlockMode

letter          [a-zA-Z]
noletter        [^a-zA-Z]
digit           [0-9]
CCletter        [a-zA-Z_]
idfier          {letter}+
CCidfier        ({CCletter}({CCletter}|{digit})*)
space           [\t ]
w               {space}*
ws              {space}+
escchar         [\\]
sign            [+-]
exp             [eE]
number          {digit}+
signNumber      ({sign}?{number})
floatNumber     ({signNumber}\.|{signNumber}\.{number})
expNumber       ({floatNumber}|{signNumber}){exp}{signNumber}
No              ({signNumber}|{floatNumber}|{expNumber})
operator        [^a-zA-Z_0-9\n\r\t \\]
blockintro      [\{][\\]((tt)|(em)|(it)|(sc)|(sl))

%%
 /* Mode switching can be triggered from the parser */
 /* ----------------------------------------------- */
	if (set_CCMode) {
	    BEGIN( CCMode);
	    set_CCMode = 0;
	} 
	if (set_NestingMode) {
	    BEGIN( NestingMode);
	    set_NestingMode = 0;
	} 
	if (set_INITIAL) {
	    BEGIN( INITIAL);
	    set_INITIAL = 0;
	}

 /* Count line numbers in all modes for better error messages */
 /* --------------------------------------------------------- */
<INITIAL,CCMode,NestingMode>[\n]	{
		    line_number++;
		    return NEWLINE;
		}

 /* Rules for TeX conventions */
 /* ------------------------- */
[\\]"%"         {   /* Avoid the quoted comment symbol */
		    yylval.character   = '%';
		    return CHAR;
		}
"%".*[\n]{w}    {   /* Match one line TeX comments */
		    /* remove spaces in next line  */
		    unput( '\n');
                    break;
		}
"%".*           {   /* Match one line TeX comments at the last line in file */
                    break;
		}
[\\]path[^a-zA-Z]   |   /* path */
[\\]verb[^a-zA-Z]   {   /* match LaTeX \verb"..." constructs */
		    BEGIN( VerbMode);
		    stop_character = yytext[ yyleng-1];
		    yylval.string.text = " ";
		    yylval.string.len  = 0;
	  	    return SPACE;
                }
<VerbMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return SPACE;
		}
<VerbMode>.	{
		    if ( yytext[0] == stop_character) {
		        BEGIN( INITIAL);
	                yylval.string.text = " ";
		        yylval.string.len  = 0;
	  	        return SPACE;
                    }
	            yylval.character = yytext[0];
	  	    return CHAR;
		}

 /* Trigger on different keywords */
 /* ----------------------------- */
[\\]begin{w}[\{]ccClass(Template)?[\}]{w}   {
		    BEGIN( CCMode);
		    return BEGINCLASS;
		}
[\\]end{w}[\{]ccClass(Template)?[\}]   {
		    return ENDCLASS;
		}
[\\]begin{w}[\{]ccRef((Class)|(Concept))[\}]{w}   {
		    BEGIN( CCMode);
		    return BEGINREFCLASS;
		}
[\\]end{w}[\{]ccRef((Class)|(Concept))[\}]   {
		    return ENDREFCLASS;
		}
[\\]begin{w}[\{]ccRef{idfier}[\}]{w}   {
		    BEGIN( CCMode);
		    return BEGINREFPAGE;
		}
[\\]end{w}[\{]ccRef{idfier}[\}]   {
		    return ENDREFPAGE;
		}
[\\]ccCreationVariable{w}[\{]{w}.*{w}[\}]   {
		    char *s = yytext + yyleng - 2;
		    while (( *s == ' ') || ( *s == '\t'))
		        s--;
		    char *r = yytext;
		    while ( *r != '{')
		        r++;
		    r++;
		    while (( *r == ' ') || ( *r == '\t'))
		        r++;
		    r[-1] = '`';
		    s[1]  = '\'';
		    s[2]  = 0;
		    free( creationvariable);
		    creationvariable = strdup( r - 1);
		    s[1] = 0;
	            yylval.string.text = r;
		    yylval.string.len  = s - r + 1;
		    return CREATIONVARIABLE;
		}		    
[\\]ccConstructor{w} {   /* constructor declaration: change to CCMode */
		    BEGIN( CCMode);
		    return CONSTRUCTOR;
		}
[\\]ccMemberFunction{w}   {   /* method declaration: change to CCMode */
		    BEGIN( CCMode);
		    return METHOD;
		}
[\\]ccMethod{w}   {   /* method declaration: change to CCMode */
		    BEGIN( CCMode);
		    return METHOD;
		}
[\\]ccFunction{w} {   /* function declaration: change to CCMode */
		    BEGIN( CCMode);
		    return FUNCTION;
		}
[\\]ccFunctionTemplate{w} {   /* function template declaration: 
			       change to CCMode */
		    BEGIN( CCMode);
		    return FUNCTIONTEMPLATE;
		}
[\\]ccVariable{w} {   /* variable declaration: change to CCMode */
		    BEGIN( CCMode);
		    return VARIABLE;
		}
[\\]ccTypedef{w} {   /* typedef declaration: change to CCMode */
		    BEGIN( CCMode);
		    return TYPEDEF;
		}
[\\]ccNestedType{w} {   /* nested type declaration: change to CCMode */
		    BEGIN( CCMode);
		    return NESTEDTYPE;
		}
[\\]ccEnum{w} {   /* enum declaration: change to CCMode */
		    BEGIN( CCMode);
		    return ENUM;
		}
[\\]ccStruct{w} { /* struct declaration: change to CCMode */
		    BEGIN( CCMode);
		    return STRUCT;
		}
[\\]ccGlobalFunction{w} {   /* function declaration: change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALFUNCTION;
		}
[\\]ccGlobalFunctionTemplate{w} {   /* function template declaration: 
			       change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALFUNCTIONTEMPLATE;
		}
[\\]ccGlobalVariable{w} {   /* variable declaration: change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALVARIABLE;
		}
[\\]ccGlobalTypedef{w} {   /* typedef declaration: change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALTYPEDEF;
		}
[\\]ccGlobalEnum{w} {   /* enum declaration: change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALENUM;
		}
[\\]ccGlobalStruct{w} {   /* struct declaration: change to CCMode */
		    BEGIN( CCMode);
		    return GLOBALSTRUCT;
		}
[\\]ccDeclaration{w} {   /* general declaration: change to CCMode */
		    BEGIN( CCMode);
		    return DECLARATION;
		}
[\\]ccHidden{w}   {   /* treat like a space */
		    yylval.string.text = "       ";
		    yylval.string.len  = 7;
	  	    return SPACE;
		}
[\\]ccUnchecked{w}  {  /* trigger a global boolean and treat it like a space */
		    unchecked_tag = 1;
		    yylval.string.text = "          ";
		    yylval.string.len  = 10;
	  	    return SPACE;
		}


 /* Specialized keywords, partly TeX, partly from the manual style */
 /* -------------------------------------------------------------- */
[\\]((ccc)|(ccStyle)){w}  { /* CCstyle formatting: change to NestingMode */
		    BEGIN( NestingMode);
		    return CCSTYLE;
		}
[\\]cc(Pure)?Var{w}      {
		    if ( creationvariable) {
	                yylval.string.text = creationvariable;
		        yylval.string.len  = strlen( creationvariable);
		    } else {
		        printErrorMessage( VariableUsedError);
	                yylval.string.text = "Unknown creationvariable";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
[\\]cc(Pure)?Class(Template)?Name{w}      {
		    if ( global_classname) {
	                yylval.string.text = global_classname;
		        yylval.string.len  = strlen( global_classname);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "[Unknown classname]";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
[\\]cc(Pure)?TemplateParameters{w}      {
		    if ( global_template_params) {
	                yylval.string.text = global_template_params;
		        yylval.string.len  = strlen( global_template_params);
		    } else {
		        printErrorMessage( ClassnameUsedError);
	                yylval.string.text = "Unknown template_params";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
[\\]cc(Pure)?RefName{w}      {
		    if ( global_ref_name) {
	                yylval.string.text = global_ref_name;
		        yylval.string.len  = strlen( global_ref_name);
		    } else {
		        printErrorMessage( RefNameUsedError);
	                yylval.string.text = "[Unknown RefName]";
		        yylval.string.len  = strlen( yylval.string.text);
		    }
		    return STRING;
                 }
[\\]ccSection{w} {  return CCSECTION; }
[\\]ccSubsection{w}  {  return CCSUBSECTION; }
[\\]ldots{w}     {
	            yylval.string.text = "...";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]RCSdef(Date)?{w}    {
		    return GOBBLETWOPARAMS;
                 }
[\\]CC{w}        {
	            yylval.string.text = "C++";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]gcc{w}        {
	            yylval.string.text = "g++";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]nat{w}       {
	            yylval.string.text = "|N";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]real{w}      {
	            yylval.string.text = "|R";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]N{w}         {
	            yylval.string.text = "|N";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]R{w}         {
	            yylval.string.text = "|R";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]Q{w}         {
	            yylval.string.text = "|Q";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]Z{w}         {
	            yylval.string.text = "Z";
		    yylval.string.len  = 1;
		    return STRING;
                 }
[\\]E{w}         {
	            yylval.string.text = "|E";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]stl{w}         {
	            yylval.string.text = "STL";
		    yylval.string.len  = 3;
		    return STRING;
                 }
[\\]leda{w}         {
	            yylval.string.text = "LEDA";
		    yylval.string.len  = 4;
		    return STRING;
                 }
[\\]cgal{w}      {
	            yylval.string.text = "CGAL";
		    yylval.string.len  = 4;
		    return STRING;
                 }
[\\]protocgal{w} {
	            yylval.string.text = "C++GAL";
		    yylval.string.len  = 6;
		    return STRING;
                 }
[\\]plageo{w}    {
	            yylval.string.text = "PLAGEO";
		    yylval.string.len  = 6;
		    return STRING;
                 }
[\\]ccDefinition{w} {
	            yylval.string.text = "DEFINITION";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccIsModel{w} {
	            yylval.string.text = "IS MODEL FOR THE CONCEPT";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccInheritsFrom{w} {
	            yylval.string.text = "INHERITS FROM";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccParameters{w} {
	            yylval.string.text = "PARAMETERS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccConstants{w} {
	            yylval.string.text = "CONSTANTS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccTypes{w} {
	            yylval.string.text = "TYPES";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccCreation{w}  {
	            yylval.string.text = "CREATION";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccOperations{w} {
	            yylval.string.text = "OPERATIONS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccAccessFunctions{w} {
	            yylval.string.text = "ACCESS FUNCTIONS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccPredicates{w} {
	            yylval.string.text = "PREDICATES";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccModifiers{w} {
	            yylval.string.text = "MODIFIERS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccHasModels{w} {
	            yylval.string.text = "HAS MODELS";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccSeeAlso{w} {
	            yylval.string.text = "SEE ALSO";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccImplementation{w} {
	            yylval.string.text = "IMPLEMENTATION";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccExample{w}   {
	            yylval.string.text = "EXAMPLE";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]begin[\{]ccAdvanced[\}]   {
	            yylval.string.text = "(begin of an advanced section)";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]end[\{]ccAdvanced[\}]   {
	            yylval.string.text = "(end of an advanced section)";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccInclude{w}   {
		    return INCLUDE;
                 }
[\\]ccHeading{w}   {
		    return HEADING;
                 }
[\\]ccPrecond/{noletter}      {
	            yylval.string.text = "Precondition:";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccPostcond/{noletter}     {
	            yylval.string.text = "Postcondition:";
		    yylval.string.len  = -1;
		    return STRING;
                 }
[\\]ccCommentHeading{w}/{noletter}   {
		    return COMMENTHEADING;
                 }
[\\]le[q]?/{noletter}       {
	            yylval.string.text = "<=";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]ge[q]?/{noletter}       {
	            yylval.string.text = ">=";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]neq/{noletter}          {
	            yylval.string.text = "!=";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]equiv/{noletter}        {
	            yylval.string.text = "==";
		    yylval.string.len  = 2;
		    return STRING;
                 }
[\\]ccSetTwoOfThreeColumns{w} { return GOBBLETWOPARAMS; }
[\\]ccSetThreeColumns{w}      { return GOBBLETHREEPARAMS; }
[\\]ccThree/{noletter}        { return GOBBLETHREEPARAMS; }
[\\]ccSetOneOfTwoColumns{w}   { return GOBBLEONEPARAM; }
[\\]ccSetTwoColumns{w}        { return GOBBLETWOPARAMS; }
[\\]ccTwo/{noletter}          { return GOBBLETWOPARAMS; }
[\\]cc((PropagateThreeToTwoColumns)(ThreeToTwo)){w}  {
		    yylval.string.text = " ";
		    yylval.string.len  = 0;
	  	    return SPACE;
}
[\\]cc((Glue((Begin)|(End)|(Declarations))?)|(ParDims)){w}  {
		    yylval.string.text = " ";
		    yylval.string.len  = 0;
	  	    return SPACE;
}
[\\]g?def{w}[\\]{idfier}{w}    {
		    return GOBBLEONEPARAM;
}
[\\](re)?newcommand{w}    {
		    return NEWCOMMAND;
}
[\\]ccTag((Defaults)|(FullDeclarations)){w}  {
		    yylval.string.text = " ";
		    yylval.string.len  = 0;
	  	    return SPACE;
}
[\\]ccChapter((Author)|(Release)|(SubTitle)){w}   {
		    return GOBBLEONEPARAM;
}

 /* HTML support                        */
 /* ----------------------------------- */
[\\]ccAnchor{w}  {
		    return GOBBLEONEPARAM;
}
[\\]ccTexHtml{w}  {
		    return GOBBLEAFTERONEPARAM;
}
[\\]begin{w}[\{]ccTexOnly[\}]{w}   { break;}
[\\]end{w}[\{]ccTexOnly[\}]        { break;}

[\\]begin{w}[\{]verbatim[\}]{w}      |
[\\]begin{w}[\{]ccExampleCode[\}]{w} |
[\\]begin{w}[\{]cprog[\}]{w}         |
[\\]begin{w}[\{]alltt[\}]{w}         |
[\\]begin{w}[\{]ccHtmlOnly[\}]{w}    {
		    BEGIN( VerbatimMode);
		    break;
		}
<VerbatimMode>[\n]	{
		    line_number++;
		    break;
}
<VerbatimMode>.  	{ /* ignore */ 
		    break;
}

<VerbatimMode>[\\]end{w}[\{]verbatim[\}]{w}      |
<VerbatimMode>[\\]end{w}[\{]ccExampleCode[\}]{w} |
<VerbatimMode>[\\]end{w}[\{]cprog[\}]{w}         |
<VerbatimMode>[\\]end{w}[\{]alltt[\}]{w}         |
<VerbatimMode>[\\]end{w}[\{]ccHtmlOnly[\}]       {
	            BEGIN(INITIAL);
		    break;
		}

[\\]begin{w}[\{]ccHtmlBlock[\}]{w}   {
		    BEGIN( HtmlBlockMode);
		    break;
		}
<HtmlBlockMode>[\n]	{
		    line_number++;
		    break;
}
<HtmlBlockMode>.  	{ /* ignore */ 
		    break;
}

<HtmlBlockMode>[\\]end{w}[\{]ccHtmlBlock[\}]      {
	            BEGIN(INITIAL);
		    break;
		}

 /* Flexibility for HTML class files. */
 /* -------------------------------------------------------------- */
<INITIAL,NestingMode>[\\]ccHtmlNoClassFile/{noletter}    {
		    skipspaces();
		    break;
}
<INITIAL,NestingMode>[\\]ccHtmlNo((Class)?)Links/{noletter}   {
		    skipspaces();
		    break;
}
<INITIAL,NestingMode>[\\]ccHtmlNo((Class)?)Index/{noletter}   {
		    skipspaces();
		    break;
}
[\\]begin{w}[\{]ccHtmlClassFile[\}]{w}   {
		    return GOBBLETWOPARAMS;
}
[\\]end{w}[\{]ccHtmlClassFile[\}]   {
		    skipspaces();
		    break;
}
[\\]ccHtmlIndex/{noletter}                                     {
		    skipspaces();
		    return GOBBLEONEPARAM;
}
[\\]ccHtmlIndex[\[][^\]][\]]/{noletter}                        {
		    skipspaces();
		    return GOBBLEONEPARAM;
}
[\\]ccHtmlIndexC/{noletter}                                    {
		    skipspaces();
		    return GOBBLEONEPARAM;
}
[\\]ccHtmlIndexC[\[][^\]][\]]/{noletter}                       {
		    skipspaces();
		    return GOBBLEONEPARAM;
}
[\\]ccHtmlCrossLink/{noletter}                                    {
		    skipspaces();
		    return GOBBLEONEPARAM;
}

 /* make the $ delimiters for math mode disappear: */
 /* -------------------------------------------------------------- */
[$]              {
		    break;
}
[\\]"&"          {
		     yylval.character   = '&';
		     return CHAR;
		 }
[\\]"_"          {
		     yylval.character   = '_';
		     return CHAR;
		 }
[\\]"^"          {
		     yylval.character   = '^';
		     return CHAR;
		 }
[\\]"#"          {
		     yylval.character   = '#';
		     return CHAR;
		 }
[~]|([\\]" ")|([\\][\n])         {
	            yylval.string.text = " ";
		    yylval.string.len  = 1;
	  	    return SPACE;
		 }
[\\][\\]         {
	  	    return NEWLINE;
		 }


 /* Grouping symbols */
 /* ---------------- */
[\\][\{]        {
	            yylval.character = '{';
	  	    return CHAR;
		}
[\\][\}]        {
	            yylval.character = '}';
	  	    return CHAR;
		}

<INITIAL,CCMode,NestingMode>[\{]    {
		    return '{';
		}
<INITIAL,CCMode,NestingMode>[\}]    {
	  	    return '}';
		}

{blockintro}    {  /* A couple of TeX styles like {\tt ... */
		    return BLOCKINTRO;
		}

<NestingMode>[\[]    {
		    return '[';
		}
<NestingMode>[\]]    {
	  	    return ']';
		}

<NestingMode>[\(]    {
		    return '(';
		}
<NestingMode>[\)]    {
	  	    return ')';
		}
<NestingMode>[\\]"&"    {
		     yylval.character   = '&';
		     return CHAR;
		 }
<NestingMode>[\\]"_"    {
		     yylval.character   = '_';
		     return CHAR;
		 }
<NestingMode>[\\]"^"    {
		     yylval.character   = '^';
		     return CHAR;
		 }
<NestingMode>[\\]"#"    {
		     yylval.character   = '#';
		     return CHAR;
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
	  	    return SPACE;
		}
<CCMode,NestingMode>{ws}	{
	            yylval.string.text = yytext;
		    yylval.string.len  = yyleng;
	  	    return SPACE;
		}
<INITIAL,CCMode,NestingMode>.	{
	            yylval.character = yytext[0];
	  	    return CHAR;
		}
%%

void skipspaces( void) {
    int c = yyinput();
    while( c && c <= ' ')
        c = yyinput();
    unput( c);
}

void init_scanner( FILE* in){
    line_number  = 1;
    set_CCMode   = 0;
    set_NestingMode   = 0;
    set_INITIAL  = 0;
    unchecked_tag = 0;
    yyrestart( in);
}



