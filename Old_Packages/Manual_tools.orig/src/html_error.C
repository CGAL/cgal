/**************************************************************************
 
  html_error.C
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Error messages and protocol output.
  System    : bison, flex, C++ (g++)
  Author    : (c) 1998 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#include <html_error.h>
#include <html_config.h>
#include <lex_include.h>

// external lex declaration
void printScannerState( ostream& out);


// Functions belonging to the Error messages
// -----------------------------------------
const char* errorMessage( ErrorNumber n) {
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
	return "Macro is undefined";
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
    case UserDefinedError:
	return "User defined error message";
    }
    return "UNKNOWN ERROR MESSAGE NUMBER";
}

void  printErrorMessage( ErrorNumber n){
    cerr << endl;
    if ( in_file)
	cerr << "*** Error " << int(n) << " in line " << in_file->line() 
	     << " in `"  << in_file->name() << "': " << errorMessage( n) 
	     << '.' << endl;
    else
	cerr << "*** Error " << int(n) << " at top level (no file): " 
	     << errorMessage( n) << '.' << endl;
    if ( in_file != in_string)
	cerr << "    while expanding macro `" << in_string->name() 
	     << "' in line " << in_string->line() << '.' << endl;
    if ( stack_trace_switch)
	cerr << include_stack;
    cerr << "Parser state: ";
    printScannerState( cerr);
    cerr << endl;
}

int yyerror( char *s) {
    printErrorMessage( ParseError);
    return 0;
}


// EOF //


