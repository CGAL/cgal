/**************************************************************************
 
  config.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Configuration constants and variables
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_CONFIG)
#define MODULE_CONFIG 1

#include <buffer.h>

// Global declarations that are implemented in the main module.
// There they can be taylored to the specific application, i.e.
// extraction or checker.

void handleComment( const Buffer_list& T);
void handleConstructorComment( const Buffer_list& T);
void handleMainComment( const Buffer_list& T);

void handleClass( const char* classname);
void handleClassEnd( void);
void handleRefPage( const char* token);
void handleRefPageEnd( void);

void handleDeclaration( const char* decl);
void handleNestedType( const char* decl);
void handleMethodDeclaration( const char* decl);
void handleConstructorDeclaration( const char* decl);
void handleFunctionDeclaration( const char* decl);
void handleFunctionTemplateDeclaration( const char* templ, const char* decl);

// Error messages and Numbers
// ==========================
enum ErrorNumber {
    NoError = 0,
    ParseError,         // must stay at position 1
    VariableUsedError,
    ClassnameUsedError,
    RefNameUsedError,
    TemplateParamExpectedError,
    MalformedTemplateParamError,
    MalformedFunctionDeclaration,
    SemicolonMissingError
};


// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in syntax.y

const char* errorMessage( ErrorNumber n);
void  printErrorMessage( ErrorNumber n);

#endif // MODULE_CONFIG //
