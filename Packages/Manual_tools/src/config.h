/**************************************************************************
 
  config.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Configuration constants and variables
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_CONFIG)
#define MODULE_CONFIG 1


// Global declarations that are implemented in the main module.
// There they can be taylored to the specific application, i.e.
// extraction or checker.

void handleComment( const Text& T);
void handleConstructorComment( const Text& T);
void handleMainComment( const Text& T);

void handleClass( const char* classname);
void handleClassEnd( void);
void handleClassTemplate( const char* classname);
void handleClassTemplateEnd( void);

void handleDeclaration( const char* decl);
void handleNestedType( const char* decl);
void handleMethodDeclaration( const char* decl);
void handleConstructorDeclaration( const char* decl);
void handleFunctionDeclaration( const char* decl);
void handleFunctionTemplateDeclaration( const char* templ, const char* decl);

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
    SemicolonMissingError
};


// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in syntax.y

const char* errorMessage( ErrorNumber n);
void  printErrorMessage( ErrorNumber n);

#endif // MODULE_CONFIG //
