/**************************************************************************
 
  html_error.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Error messages and protocol output.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#ifndef HTML_ERROR_H
#define HTML_ERROR_H 1

#include <iostream.h>

// Error messages and Numbers
// ==========================
enum ErrorNumber {
    NoError,
    ParseError,         // must stay at position 1
    VariableUsedError,
    ClassnameUsedError,
    TemplateParamExpectedError,
    MalformedTemplateParamError,
    MalformedFunctionDeclaration,
    SemicolonMissingError,
    ChapterStructureError,
    UnknownIndexCategoryError,
    EmptyClassNameError,
    EmptyCrossLinkError,
    ParameterOptionError,
    NParamRangeError,
    ParamIndexError,
    MacroUndefinedError,
    MacroDefUnknownError,
    MacroParamNumberError,
    EOFInIncludeFilenameError,
    EOFInMacroExpansionError,
    UnknownKeyError,
    ParseCCError,
    OutputStackEmptyError,
    OutputStackKeyError,
    MacroInCModeError,
    ParsingStateError,
    MacroStackUnderflowError,
    IncludeStackUnderflowError,
    FileReadOpenError,
    UserDefinedError
};


// Functions belonging to the Error messages
// -----------------------------------------

const char* errorMessage( ErrorNumber n);
void        printErrorMessage( ErrorNumber n);
int         yyerror( char *s);  // needed for yyparse()

#endif // HTML_ERROR_H 1 //
// EOF //

