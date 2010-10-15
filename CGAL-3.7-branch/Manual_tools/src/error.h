/**************************************************************************
 
  html_error.h
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Error messages and protocol output.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#ifndef ERROR_H
#define ERROR_H 1

#include <iostream>
#include <string>

using namespace std;

// Color attributes
// ----------------
extern string BlueColor;     // blue boldface
extern string BoldColor;     // boldface
extern string OkColor;       // green
extern string ErrorColor;    // red boldface
extern string WarnColor;     // magenta
extern string ResetColor;    // black, reset attribute

// Error and warning lead texts, terminate with ResetColor
// -------------------------------------------------------
extern string ErrorText;
extern string WarnText;

void enableColor(); // switch to Ansi Color codes for above strings, default
                    // is no color codes.

// Error messages and Numbers
// ==========================
enum ErrorNumber {
    NoError = 0,
    ParseError,         // must stay at position 1
    VariableUsedError,
    ClassnameUsedError,
    TemplateParamExpectedError,
    MalformedTemplateParamError,
    MalformedFunctionDeclaration,
    SemicolonMissingError,
    ChapterStructureError,
    ChapterAbsolutePathError,
    ClassAbsolutePathError,
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
    RomansOutOfBoundsError,
    AlphaOutOfBoundsError,
    RefPageNotClosedError,
    UserDefinedError
};


// Functions belonging to the Error messages
// -----------------------------------------

// Returns the first error number reported since program start or since
// the last call to firstError(). Returns 'NoError' if no error was
// reported yet. Resets firstError() to 'NoError'.
ErrorNumber firstError();

// Returns the formatted error message, where sequences of %<i> are 
// replaced with the argument text arg<i> for 1 <= i <= 6.
// Includes the leading ERROR: (incl. color), filename and linenumber.
string formattedErrorMessage( ErrorNumber n,
                              const string& arg1, const string& arg2,
                              const string& arg3, const string& arg4,
                              const string& arg5, const string& arg6);

// The raw error message incl. control codes for %<i> expansion.
const char* errorMessage( ErrorNumber n);

// Prints the error message and updates the firstError() status.
// Depending the error message, optional arguments arg<i> are 
// used and expaned where the raw message contains %<i> control 
// sequences.
void        printErrorMessage( ErrorNumber n,
                               string arg1 = string(),
                               string arg2 = string(),
                               string arg3 = string(),
                               string arg4 = string(),
                               string arg5 = string(),
                               string arg6 = string());

int         yyerror( char *s);  // needed for yyparse()

#endif // ERROR_H 1 //
// EOF //

