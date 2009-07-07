/**************************************************************************
 
  error.cpp
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Error messages and protocol output.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <error.h>
#include <config.h>
#include <input.h>
#include <string_conversion.h>

// external lex declaration
void printScannerState( ostream& out);


// Color attributes
// ----------------
string BlueColor;     // blue boldface
string BoldColor;     // boldface
string OkColor;       // green
string ErrorColor;    // red boldface
string WarnColor;     // magenta
string ResetColor;    // black, reset attribute

// Error and warning lead texts, terminate with ResetColor
// -------------------------------------------------------
string ErrorText( "ERROR: ");
string WarnText( "WARNING: ");

void enableColor() {
    BlueColor  = string( "[1m[34m");    // blue boldface
    BoldColor  = string( "[1m");          // boldface
    OkColor    = string( "[32m");         // green
    ErrorColor = string( "[1m[31m");    // red boldface
    WarnColor  = string( "[35m");         // magenta
    ResetColor = string( "[30m[00m");   // black, reset attribute

    ErrorText  = ErrorColor + string( "ERROR: ");
    WarnText   = WarnColor  + string( "WARNING: ");
}

// Error number reporting
// ----------------------

static ErrorNumber firstErrorVar = NoError;

// Returns the first error number reported since program start or since
// the last call to firstError(). Returns 'NoError' if no error was
// reported yet. Resets firstError() to 'NoError'.
ErrorNumber firstError() {
    ErrorNumber result = firstErrorVar;
    firstErrorVar = NoError;
    return result;
}

// Functions belonging to the Error messages
// -----------------------------------------
const char* errorMessage( ErrorNumber n) {
    if ( firstErrorVar == NoError)
        firstErrorVar = n;
    switch ( n) {
    case NoError:
	return "NO ERROR";
    case ParseError:
	return "parse error";
    case VariableUsedError:
	return "The creationvariable was used but not defined";
    case ClassnameUsedError:
	return "The classname was used out of scope of any class";
    case TemplateParamExpectedError:
        return "A template parameter is missing";
    case MalformedTemplateParamError:
        return "The template parameter is malformed (<> nesting ..)";
    case MalformedFunctionDeclaration:
        return "The function declaration is malformed";
    case SemicolonMissingError:
        return "The declaration does not end in a semicolon";
    case ChapterStructureError:
        return "Malformed chapter structure: one chapter per file";
    case ChapterAbsolutePathError:
        return "Input filename must be given with a relative path for files containing chapters";
    case ClassAbsolutePathError:
        return "Input filename must be given with a relative path for files containing class or reference page environments";
    case UnknownIndexCategoryError:
        return "Unknown index category in optional argument of \\ccHtmlIndex";
    case EmptyClassNameError:
	return "The classname was empty";
    case EmptyCrossLinkError:
	return "The key for a cross link was empty";
    case ParameterOptionError:
	return "Mixed up nesting of {} and [] for optional macro parameter";
    case NParamRangeError:
	return "Parameter number in macro definition out of range 1..9, a..z";
    case ParamIndexError:
	return "Index too large or illegal character behind # encountered";
    case MacroUndefinedError:
	return "Macro '%1' is undefined";
    case MacroDefUnknownError:
	return "TeX style macro definition is not understandable";
    case MacroParamNumberError:
	return "The number of (optional) parameters does not match macro definition";
    case EOFInIncludeFilenameError:
	return "EOF in input/include file name detected";
    case EOFInMacroExpansionError:
	return "EOF detected while expanding a macro or removing its trailing spaces";
    case ParseCCError:
	return "Need an opening parantheses to parse C++ expression";
    case UnknownKeyError:
	return "Unknown key in C++ declaration formatting (internal error)";
    case OutputStackEmptyError:
	return "Output stack is empty while doing a pop";
    case OutputStackKeyError:
	return "Unknown key for output stack push operation";
    case MacroInCModeError:
	return "Unknown macro occured in C++ text (only restricted macro expansion)";
    case ParsingStateError:
	return "Parsing teminates in wrong state (should be INITIAL)";
    case MacroStackUnderflowError:
	return "Nesting error: Underflow of the macro scope stack";
    case IncludeStackUnderflowError:
	return "Internal error: Underflow of the include file and macro stack";
    case FileReadOpenError:
	return "Open file for read failed";
    case RomansOutOfBoundsError:
	return "Conversion-to-Roman-digits parameter out of bounds";
    case AlphaOutOfBoundsError:
	return "Conversion-to-Alpha-digit parameter out of bounds";
    case RefPageNotClosedError:
	return "Previous reference page is not closed";
    case UserDefinedError:
	return "User defined error message";
    }
    return "UNKNOWN ERROR MESSAGE NUMBER";
}

// Returns the formatted error message, where sequences of %<i> are 
// replaced with the argument text arg<i> for 1 <= i <= 6.
// Includes the leading ERROR: (incl. color), filename and linenumber.
string formattedErrorMessage( ErrorNumber n,
                              const string& arg1, const string& arg2,
                              const string& arg3, const string& arg4,
                              const string& arg5, const string& arg6) {
    string err = errorMessage( n);
    size_t i = 0;
    while ( i+1 < err.size()) {
        if ( err[i] == '%' && err[i+1] == '1') {
            err.replace( i, 2, arg1);
            i += arg1.size() - 1;
        }
        ++i;
    }
    string result = ErrorText;
    if ( in_file)
        result = result + string("'") + in_file->name() + string( "' line ") 
                 + int_to_string(in_file->line()) + string(": ");
    else
	result = result + string( " at top level (no file): ");
    result = result + err + string(".") + ResetColor;
    return result;
}

void  printErrorMessage( ErrorNumber n, 
                         string arg1, string arg2, string arg3,
                         string arg4, string arg5, string arg6){
    if ( firstErrorVar == NoError)
        firstErrorVar = n;
    cerr << endl  << formattedErrorMessage(n, arg1, arg2, arg3, arg4, 
                                           arg5, arg6);
    cerr << ErrorColor;
    if ( in_file != in_string)
	cerr << "\n    while expanding macro '" << in_string->name() 
             << "' line " << in_string->line() << '.';
    if ( verbose_switch) {
        cerr << "\n    (parser state: ";
        printScannerState( cerr);
        cerr << ')';
    }
    if ( stack_trace_switch)
	cerr << '\n' << include_stack;
    cerr << ResetColor << endl;
}

int yyerror( char *s) {
    printErrorMessage( ParseError);
    return 0;
}


// EOF //


