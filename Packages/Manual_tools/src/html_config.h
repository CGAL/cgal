/**************************************************************************
 
  html_config.h
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : Configuration constants and variables
  System    : C++ (g++)
  Author    : (c) 1997 Lutz Kettner
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/

#if ! defined( MODULE_CONFIG)
#define MODULE_CONFIG 1


// Global declarations that are implemented in the main module.
// There they can be taylored to the specific application, i.e.
// extraction or checker.

/* Customization tags for the style */
/* ================================ */
extern bool tag_chapter_author;
extern bool tag_replace_prefix;
extern bool tag_replace_include;
extern bool tag_long_param_layout;

extern bool tag_rm_const_ref_pair;
extern bool tag_rm_eigen_class_name;
extern bool tag_operator_layout;
extern bool tag_rm_trailing_const;

void tag_defaults();
void tag_full_declarations();


/* read a file into a buffer */
/* The name has a trailing '}' */
Buffer* readFileInBuffer(const char* name);

/* An empty List as empty comment for global declarations */
/* ====================================================== */
extern Text empty_comment;

// void handleComment( const Text& T);
// void handleConstructorComment( const Text& T);
void handleMainComment( const Text& T);

void handleChapter(  const Text& T);
void handleSection(  const Text& T);
void handleLabel(    const char* l);

void handleText(       const Text&      T, bool check_nlnl = false);
void handleBuffer(     const Buffer&    B);
void handleTextToken(  const TextToken& TT);
void handleString(     const char*      s);
void handleChar(       char             c);

void handleBiblio(  const Text& T);
Buffer* handleCite( const char* l);
// for an empty item name use the key name as item name
Buffer* handleBibItem( const char* key_name, const char* item_name = 0);

void handleClass( const char* classname);
void handleClassEnd( void);
void handleClassTemplate( const char* classname);
void handleClassTemplateEnd( void);

void handleDeclaration( const char* decl);

void handleMethodDeclaration( const char* decl,
			      const Text& T = empty_comment);
void handleConstructorDeclaration( const char* decl,
				   const Text& T = empty_comment);
void handleFunctionDeclaration( const char* decl,
				const Text& T = empty_comment);
void handleFunctionTemplateDeclaration( const char* templ,
					const char* decl,
					const Text& T = empty_comment);
void handleVariableDeclaration( const char* decl,
				const Text& T = empty_comment);
void handleTypedefDeclaration( const char* decl,
			       const Text& T = empty_comment);
void handleNestedTypeDeclaration( const char* decl,
			          const Text& T = empty_comment);
void handleEnumDeclaration( const char* decl,
			    const Text& T = empty_comment);
void handleStructDeclaration( const char* decl,
			      const Text& T = empty_comment);

void handleNewCommand( char* idfier, char* body);


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
    IncludeNestingTooDeepError,
    IncludeOpenError,
    ChapterStructureError
};


// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in syntax.y

const char* errorMessage( ErrorNumber n);
void  printErrorMessage( ErrorNumber n);

#endif // MODULE_CONFIG //
