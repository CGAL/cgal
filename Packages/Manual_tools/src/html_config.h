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

extern const char* prog_name;
extern char* current_filename;


extern char* cgal_lib_dir;

/* An object storing the current font */
/* ================================== */
/* It is only used within CCMode at the moment */
enum Font { unknown_font = -1, 
	     rm_font, 
	     tt_font, 
	     bf_font, 
	     it_font, 
	     sl_font,
             sc_font,
             sf_font,
             var_font,
             math_font,
             end_font_array};

extern Font current_font;
const char* font_changing_tags( Font old_font, Font new_font);
const char* new_font_tags( Font new_font);
const char* lazy_new_font_tags( Font new_font);

// Index sorting.
// ==============
// sort_keys are used to sort the index according to different sections.
// A sort_key is followed by a 0 to indicate the section title.
// A sort_key is followed by a 1 to indicate a normal entry.

extern const char* sort_key_class;
extern const char* sort_key_nested_type;
extern const char* sort_key_struct;
extern const char* sort_key_enum;
extern const char* sort_key_enum_tags;
extern const char* sort_key_typedef;
extern const char* sort_key_variable;
extern const char* sort_key_function;
extern const char* sort_key_member_function;

const char* find_sort_key( const char* txt);

/* Flexibility for HTML class files. */
/* ================================= */
extern bool html_no_class_links;
extern bool html_no_class_file;
extern bool html_no_class_index;

extern bool html_inline_classes;

void handleHtmlClassFile( const char* filename, const Text& T);
void handleHtmlClassFileEnd();

const char* handleHtmlIndexC( const char* category, const char* item);
const char* handleHtmlIndex( const char* category, const char* item);
const char* handleHtmlCrossLink( const char* key, bool tmpl_class = false);

/* Functions to manage footnotes. */
/* ============================== */
void insertFootnote( char* s);
// increments counter and returns current value.
int nextFootnoteCounter();
// format footnote reference and hyperlink based on actual counter 
// into a Buffer.
Buffer* formattedFootnoteNumber();
// prints footnote reference and hyperlink based on actual counter.
void printFootnoteCounter( ostream& out);
// prints footnotes and resets counter.
void printFootnotes( ostream& out);

/* A map for user defined macros     */
/* ================================= */
void insertMacro( const char* key, const char* value);
const char* fetchMacro( const char* key);
char* extractRCS( char* s);
char* extractRCSDate( char* s);


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

extern bool tag_rm_template;
extern bool tag_template_inline;

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
void handleLabel(    const char* l, size_t len); // trust only len!

void handleText(       const Text&      T, bool check_nlnl = false);
void handleBuffer(     const Buffer&    B);
void handleTextToken(  const TextToken& TT);
void handleString(     const char*      s);
void handleChar(       char             c);

void handleBiblio(  const Text& T);
Buffer* handleCite( const char* cite_keys, const char* option = 0);
// Defines a mapping between cite keys and the visible item for a cite.
Buffer* handleBibCite( const char* key, const char* item);
// for an empty item name use the key name as item name
Buffer* handleBibItem( const char* key_name, const char* item_name = 0);

void handleClass( const char* classname);
void handleClassEnd( void);
void handleClassTemplate( const char* classname);
void handleClassTemplateEnd( void);

void handleClassNameEnd( void);
void handleClassFileEnd( void);

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
void handleTexDef( const char* idfier, char* body);


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
    ChapterStructureError,
    UnknownIndexCategoryError,
    EmptyClassNameError,
    EmptyCrossLinkError
};


// Functions belonging to the Error messages
// -----------------------------------------
// See their implementations in syntax.y

const char* errorMessage( ErrorNumber n);
void  printErrorMessage( ErrorNumber n);

#endif // MODULE_CONFIG //
