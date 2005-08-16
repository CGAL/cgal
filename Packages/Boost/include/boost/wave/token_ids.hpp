/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library
    The definition of a default set of token identifiers and related 
    functions.
    
    http://www.boost.org/

    Copyright (c) 2001-2005 Hartmut Kaiser. Distributed under the Boost
    Software License, Version 1.0. (See accompanying file
    LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(TOKEN_IDS_HPP_414E9A58_F079_4789_8AFF_513815CE475B_INCLUDED)
#define TOKEN_IDS_HPP_414E9A58_F079_4789_8AFF_513815CE475B_INCLUDED

#include <string>
#include <boost/assert.hpp>

#include <boost/wave/wave_config.hpp>

///////////////////////////////////////////////////////////////////////////////
//  Allow external redefinition of the token identifiers to use
#if !defined(BOOST_WAVE_TOKEN_IDS_DEFINED)
#define BOOST_WAVE_TOKEN_IDS_DEFINED

#if defined (__FreeBSD__) && defined (T_DIVIDE)
#undef T_DIVIDE
#endif

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {

///////////////////////////////////////////////////////////////////////////////
//  assemble tokenid's
#define TOKEN_FROM_ID(id, cat)   ((id) | (cat))
#define ID_FROM_TOKEN(tok)       ((tok) & ~TokenTypeMask)
#define BASEID_FROM_TOKEN(tok)   ((tok) & ~ExtTokenTypeMask)
#define CATEGORY_FROM_TOKEN(tok) ((tok) & TokenTypeMask)
#define EXTCATEGORY_FROM_TOKEN(tok) ((tok) & ExtTokenTypeMask)
#define IS_CATEGORY(tok, cat)    \
    ((CATEGORY_FROM_TOKEN(tok) == (cat)) ? true : false) \
    /**/
#define IS_EXTCATEGORY(tok, cat)    \
    ((EXTCATEGORY_FROM_TOKEN(tok) == (cat)) ? true : false) \
    /**/

///////////////////////////////////////////////////////////////////////////////
//  the token_category helps to classify the different token types 
enum token_category {
    IdentifierTokenType         = 0x10000000,
    ParameterTokenType          = 0x11000000,
    ExtParameterTokenType       = 0x11100000,
    KeywordTokenType            = 0x20000000,
    OperatorTokenType           = 0x30000000,
    LiteralTokenType            = 0x40000000,
    IntegerLiteralTokenType     = 0x41000000,
    FloatingLiteralTokenType    = 0x42000000,
    StringLiteralTokenType      = 0x43000000,
    CharacterLiteralTokenType   = 0x44000000,
    BoolLiteralTokenType        = 0x45000000,
    PPTokenType                 = 0x50000000,
    PPConditionalTokenType      = 0x50800000,

    UnknownTokenType            = 0xA0000000,
    EOLTokenType                = 0xB0000000,
    EOFTokenType                = 0xC0000000,
    WhiteSpaceTokenType         = 0xD0000000,
    InternalTokenType           = 0xE0000000,
    
    TokenTypeMask               = 0xFF000000,
    AltTokenType                = 0x00100000,
    TriGraphTokenType           = 0x00200000,
    AltExtTokenType             = 0x00500000,   // and, bit_and etc.
    ExtTokenTypeMask            = 0xFFF00000,
    ExtTokenOnlyMask            = 0x00F00000,
    TokenValueMask              = 0x000FFFFF,
    MainTokenMask               = 0xFF0FFFFF    // TokenTypeMask|TokenValueMask
};

///////////////////////////////////////////////////////////////////////////////
//  the token_id assigns unique numbers to the different C++ lexemes
enum token_id {
    T_FIRST_TOKEN  = 256,
    T_AND          = TOKEN_FROM_ID(T_FIRST_TOKEN, OperatorTokenType),
    T_AND_ALT      = TOKEN_FROM_ID(T_FIRST_TOKEN, OperatorTokenType|AltExtTokenType),
    T_ANDAND       = TOKEN_FROM_ID(257, OperatorTokenType),
    T_ANDAND_ALT   = TOKEN_FROM_ID(257, OperatorTokenType|AltExtTokenType),
    T_ASSIGN       = TOKEN_FROM_ID(258, OperatorTokenType),
    T_ANDASSIGN    = TOKEN_FROM_ID(259, OperatorTokenType),
    T_ANDASSIGN_ALT     = TOKEN_FROM_ID(259, OperatorTokenType|AltExtTokenType),
    T_OR           = TOKEN_FROM_ID(260, OperatorTokenType),
    T_OR_ALT       = TOKEN_FROM_ID(260, OperatorTokenType|AltExtTokenType),
    T_OR_TRIGRAPH  = TOKEN_FROM_ID(260, OperatorTokenType|TriGraphTokenType),
    T_ORASSIGN     = TOKEN_FROM_ID(261, OperatorTokenType),
    T_ORASSIGN_ALT          = TOKEN_FROM_ID(261, OperatorTokenType|AltExtTokenType),
    T_ORASSIGN_TRIGRAPH     = TOKEN_FROM_ID(261, OperatorTokenType|TriGraphTokenType),
    T_XOR          = TOKEN_FROM_ID(262, OperatorTokenType),
    T_XOR_ALT      = TOKEN_FROM_ID(262, OperatorTokenType|AltExtTokenType),
    T_XOR_TRIGRAPH = TOKEN_FROM_ID(262, OperatorTokenType|TriGraphTokenType),
    T_XORASSIGN    = TOKEN_FROM_ID(263, OperatorTokenType),
    T_XORASSIGN_ALT         = TOKEN_FROM_ID(263, OperatorTokenType|AltExtTokenType),
    T_XORASSIGN_TRIGRAPH    = TOKEN_FROM_ID(263, OperatorTokenType|TriGraphTokenType),
    T_COMMA        = TOKEN_FROM_ID(264, OperatorTokenType),
    T_COLON        = TOKEN_FROM_ID(265, OperatorTokenType),
    T_DIVIDE       = TOKEN_FROM_ID(266, OperatorTokenType),
    T_DIVIDEASSIGN = TOKEN_FROM_ID(267, OperatorTokenType),
    T_DOT          = TOKEN_FROM_ID(268, OperatorTokenType),
    T_DOTSTAR      = TOKEN_FROM_ID(269, OperatorTokenType),
    T_ELLIPSIS     = TOKEN_FROM_ID(270, OperatorTokenType),
    T_EQUAL        = TOKEN_FROM_ID(271, OperatorTokenType),
    T_GREATER      = TOKEN_FROM_ID(272, OperatorTokenType),
    T_GREATEREQUAL = TOKEN_FROM_ID(273, OperatorTokenType),
    T_LEFTBRACE    = TOKEN_FROM_ID(274, OperatorTokenType),
    T_LEFTBRACE_ALT         = TOKEN_FROM_ID(274, OperatorTokenType|AltTokenType),
    T_LEFTBRACE_TRIGRAPH    = TOKEN_FROM_ID(274, OperatorTokenType|TriGraphTokenType),
    T_LESS         = TOKEN_FROM_ID(275, OperatorTokenType),
    T_LESSEQUAL    = TOKEN_FROM_ID(276, OperatorTokenType),
    T_LEFTPAREN    = TOKEN_FROM_ID(277, OperatorTokenType),
    T_LEFTBRACKET  = TOKEN_FROM_ID(278, OperatorTokenType),
    T_LEFTBRACKET_ALT       = TOKEN_FROM_ID(278, OperatorTokenType|AltTokenType),
    T_LEFTBRACKET_TRIGRAPH  = TOKEN_FROM_ID(278, OperatorTokenType|TriGraphTokenType),
    T_MINUS        = TOKEN_FROM_ID(279, OperatorTokenType),
    T_MINUSASSIGN  = TOKEN_FROM_ID(280, OperatorTokenType),
    T_MINUSMINUS   = TOKEN_FROM_ID(281, OperatorTokenType),
    T_PERCENT      = TOKEN_FROM_ID(282, OperatorTokenType),
    T_PERCENTASSIGN = TOKEN_FROM_ID(283, OperatorTokenType),
    T_NOT          = TOKEN_FROM_ID(284, OperatorTokenType),
    T_NOT_ALT      = TOKEN_FROM_ID(284, OperatorTokenType|AltExtTokenType),
    T_NOTEQUAL     = TOKEN_FROM_ID(285, OperatorTokenType),
    T_NOTEQUAL_ALT      = TOKEN_FROM_ID(285, OperatorTokenType|AltExtTokenType),
    T_OROR         = TOKEN_FROM_ID(286, OperatorTokenType),
    T_OROR_ALT     = TOKEN_FROM_ID(286, OperatorTokenType|AltExtTokenType),
    T_OROR_TRIGRAPH     = TOKEN_FROM_ID(286, OperatorTokenType|TriGraphTokenType),
    T_PLUS         = TOKEN_FROM_ID(287, OperatorTokenType),
    T_PLUSASSIGN   = TOKEN_FROM_ID(288, OperatorTokenType),
    T_PLUSPLUS     = TOKEN_FROM_ID(289, OperatorTokenType),
    T_ARROW        = TOKEN_FROM_ID(290, OperatorTokenType),
    T_ARROWSTAR    = TOKEN_FROM_ID(291, OperatorTokenType),
    T_QUESTION_MARK = TOKEN_FROM_ID(292, OperatorTokenType),
    T_RIGHTBRACE   = TOKEN_FROM_ID(293, OperatorTokenType),
    T_RIGHTBRACE_ALT        = TOKEN_FROM_ID(293, OperatorTokenType|AltTokenType),
    T_RIGHTBRACE_TRIGRAPH   = TOKEN_FROM_ID(293, OperatorTokenType|TriGraphTokenType),
    T_RIGHTPAREN   = TOKEN_FROM_ID(294, OperatorTokenType),
    T_RIGHTBRACKET = TOKEN_FROM_ID(295, OperatorTokenType),
    T_RIGHTBRACKET_ALT      = TOKEN_FROM_ID(295, OperatorTokenType|AltTokenType),
    T_RIGHTBRACKET_TRIGRAPH = TOKEN_FROM_ID(295, OperatorTokenType|TriGraphTokenType),
    T_COLON_COLON  = TOKEN_FROM_ID(296, OperatorTokenType),
    T_SEMICOLON    = TOKEN_FROM_ID(297, OperatorTokenType),
    T_SHIFTLEFT    = TOKEN_FROM_ID(298, OperatorTokenType),
    T_SHIFTLEFTASSIGN = TOKEN_FROM_ID(299, OperatorTokenType),
    T_SHIFTRIGHT   = TOKEN_FROM_ID(300, OperatorTokenType),
    T_SHIFTRIGHTASSIGN = TOKEN_FROM_ID(301, OperatorTokenType),
    T_STAR         = TOKEN_FROM_ID(302, OperatorTokenType),
    T_COMPL        = TOKEN_FROM_ID(303, OperatorTokenType),
    T_COMPL_ALT         = TOKEN_FROM_ID(303, OperatorTokenType|AltExtTokenType),
    T_COMPL_TRIGRAPH    = TOKEN_FROM_ID(303, OperatorTokenType|TriGraphTokenType),
    T_STARASSIGN   = TOKEN_FROM_ID(304, OperatorTokenType),
    T_ASM          = TOKEN_FROM_ID(305, KeywordTokenType),
    T_AUTO         = TOKEN_FROM_ID(306, KeywordTokenType),
    T_BOOL         = TOKEN_FROM_ID(307, KeywordTokenType),
    T_FALSE        = TOKEN_FROM_ID(308, BoolLiteralTokenType),
    T_TRUE         = TOKEN_FROM_ID(309, BoolLiteralTokenType),
    T_BREAK        = TOKEN_FROM_ID(310, KeywordTokenType),
    T_CASE         = TOKEN_FROM_ID(311, KeywordTokenType),
    T_CATCH        = TOKEN_FROM_ID(312, KeywordTokenType),
    T_CHAR         = TOKEN_FROM_ID(313, KeywordTokenType),
    T_CLASS        = TOKEN_FROM_ID(314, KeywordTokenType),
    T_CONST        = TOKEN_FROM_ID(315, KeywordTokenType),
    T_CONSTCAST    = TOKEN_FROM_ID(316, KeywordTokenType),
    T_CONTINUE     = TOKEN_FROM_ID(317, KeywordTokenType),
    T_DEFAULT      = TOKEN_FROM_ID(318, KeywordTokenType),
    T_DEFINED      = TOKEN_FROM_ID(319, KeywordTokenType),
    T_DELETE       = TOKEN_FROM_ID(320, KeywordTokenType),
    T_DO           = TOKEN_FROM_ID(321, KeywordTokenType),
    T_DOUBLE       = TOKEN_FROM_ID(322, KeywordTokenType),
    T_DYNAMICCAST  = TOKEN_FROM_ID(323, KeywordTokenType),
    T_ELSE         = TOKEN_FROM_ID(324, KeywordTokenType),
    T_ENUM         = TOKEN_FROM_ID(325, KeywordTokenType),
    T_EXPLICIT     = TOKEN_FROM_ID(326, KeywordTokenType),
    T_EXPORT       = TOKEN_FROM_ID(327, KeywordTokenType),
    T_EXTERN       = TOKEN_FROM_ID(328, KeywordTokenType),
    T_FLOAT        = TOKEN_FROM_ID(329, KeywordTokenType),
    T_FOR          = TOKEN_FROM_ID(330, KeywordTokenType),
    T_FRIEND       = TOKEN_FROM_ID(331, KeywordTokenType),
    T_GOTO         = TOKEN_FROM_ID(332, KeywordTokenType),
    T_IF           = TOKEN_FROM_ID(333, KeywordTokenType),
    T_INLINE       = TOKEN_FROM_ID(334, KeywordTokenType),
    T_INT          = TOKEN_FROM_ID(335, KeywordTokenType),
    T_LONG         = TOKEN_FROM_ID(336, KeywordTokenType),
    T_MUTABLE      = TOKEN_FROM_ID(337, KeywordTokenType),
    T_NAMESPACE    = TOKEN_FROM_ID(338, KeywordTokenType),
    T_NEW          = TOKEN_FROM_ID(339, KeywordTokenType),
    T_OPERATOR     = TOKEN_FROM_ID(340, KeywordTokenType),
    T_PRIVATE      = TOKEN_FROM_ID(341, KeywordTokenType),
    T_PROTECTED    = TOKEN_FROM_ID(342, KeywordTokenType),
    T_PUBLIC       = TOKEN_FROM_ID(343, KeywordTokenType),
    T_REGISTER     = TOKEN_FROM_ID(344, KeywordTokenType),
    T_REINTERPRETCAST = TOKEN_FROM_ID(345, KeywordTokenType),
    T_RETURN       = TOKEN_FROM_ID(346, KeywordTokenType),
    T_SHORT        = TOKEN_FROM_ID(347, KeywordTokenType),
    T_SIGNED       = TOKEN_FROM_ID(348, KeywordTokenType),
    T_SIZEOF       = TOKEN_FROM_ID(349, KeywordTokenType),
    T_STATIC       = TOKEN_FROM_ID(350, KeywordTokenType),
    T_STATICCAST   = TOKEN_FROM_ID(351, KeywordTokenType),
    T_STRUCT       = TOKEN_FROM_ID(352, KeywordTokenType),
    T_SWITCH       = TOKEN_FROM_ID(353, KeywordTokenType),
    T_TEMPLATE     = TOKEN_FROM_ID(354, KeywordTokenType),
    T_THIS         = TOKEN_FROM_ID(355, KeywordTokenType),
    T_THROW        = TOKEN_FROM_ID(356, KeywordTokenType),
    T_TRY          = TOKEN_FROM_ID(357, KeywordTokenType),
    T_TYPEDEF      = TOKEN_FROM_ID(358, KeywordTokenType),
    T_TYPEID       = TOKEN_FROM_ID(359, KeywordTokenType),
    T_TYPENAME     = TOKEN_FROM_ID(360, KeywordTokenType),
    T_UNION        = TOKEN_FROM_ID(361, KeywordTokenType),
    T_UNSIGNED     = TOKEN_FROM_ID(362, KeywordTokenType),
    T_USING        = TOKEN_FROM_ID(363, KeywordTokenType),
    T_VIRTUAL      = TOKEN_FROM_ID(364, KeywordTokenType),
    T_VOID         = TOKEN_FROM_ID(365, KeywordTokenType),
    T_VOLATILE     = TOKEN_FROM_ID(366, KeywordTokenType),
    T_WCHART       = TOKEN_FROM_ID(367, KeywordTokenType),
    T_WHILE        = TOKEN_FROM_ID(368, KeywordTokenType),
    T_PP_DEFINE    = TOKEN_FROM_ID(369, PPTokenType),
    T_PP_IF        = TOKEN_FROM_ID(370, PPConditionalTokenType),
    T_PP_IFDEF     = TOKEN_FROM_ID(371, PPConditionalTokenType),
    T_PP_IFNDEF    = TOKEN_FROM_ID(372, PPConditionalTokenType),
    T_PP_ELSE      = TOKEN_FROM_ID(373, PPConditionalTokenType),
    T_PP_ELIF      = TOKEN_FROM_ID(374, PPConditionalTokenType),
    T_PP_ENDIF     = TOKEN_FROM_ID(375, PPConditionalTokenType),
    T_PP_ERROR     = TOKEN_FROM_ID(376, PPTokenType),
    T_PP_LINE      = TOKEN_FROM_ID(377, PPTokenType),
    T_PP_PRAGMA    = TOKEN_FROM_ID(378, PPTokenType),
    T_PP_UNDEF     = TOKEN_FROM_ID(379, PPTokenType),
    T_PP_WARNING   = TOKEN_FROM_ID(380, PPTokenType),
    T_IDENTIFIER   = TOKEN_FROM_ID(381, IdentifierTokenType),
    T_OCTALINT     = TOKEN_FROM_ID(382, IntegerLiteralTokenType),
    T_DECIMALINT   = TOKEN_FROM_ID(383, IntegerLiteralTokenType),
    T_HEXAINT      = TOKEN_FROM_ID(384, IntegerLiteralTokenType),
    T_INTLIT       = TOKEN_FROM_ID(385, IntegerLiteralTokenType),
    T_LONGINTLIT   = TOKEN_FROM_ID(386, IntegerLiteralTokenType),
    T_FLOATLIT     = TOKEN_FROM_ID(387, FloatingLiteralTokenType),
    T_FIXEDPOINTLIT = TOKEN_FROM_ID(387, FloatingLiteralTokenType|AltTokenType),  // IDL specific
    T_CCOMMENT     = TOKEN_FROM_ID(388, WhiteSpaceTokenType),
    T_CPPCOMMENT   = TOKEN_FROM_ID(389, WhiteSpaceTokenType),
    T_CHARLIT      = TOKEN_FROM_ID(390, CharacterLiteralTokenType),
    T_STRINGLIT    = TOKEN_FROM_ID(391, StringLiteralTokenType),
    T_CONTLINE     = TOKEN_FROM_ID(392, EOLTokenType),
    T_SPACE        = TOKEN_FROM_ID(393, WhiteSpaceTokenType),
    T_SPACE2       = TOKEN_FROM_ID(394, WhiteSpaceTokenType),
    T_NEWLINE      = TOKEN_FROM_ID(395, EOLTokenType),
    T_POUND_POUND           = TOKEN_FROM_ID(396, OperatorTokenType),
    T_POUND_POUND_ALT       = TOKEN_FROM_ID(396, OperatorTokenType|AltTokenType),
    T_POUND_POUND_TRIGRAPH  = TOKEN_FROM_ID(396, OperatorTokenType|TriGraphTokenType),
    T_POUND                 = TOKEN_FROM_ID(397, OperatorTokenType),
    T_POUND_ALT             = TOKEN_FROM_ID(397, OperatorTokenType|AltTokenType),
    T_POUND_TRIGRAPH        = TOKEN_FROM_ID(397, OperatorTokenType|TriGraphTokenType),
    T_ANY          = TOKEN_FROM_ID(398, UnknownTokenType),
    T_ANY_TRIGRAPH = TOKEN_FROM_ID(398, UnknownTokenType|TriGraphTokenType),
    T_PP_INCLUDE   = TOKEN_FROM_ID(399, PPTokenType),
    T_PP_QHEADER   = TOKEN_FROM_ID(400, PPTokenType),
    T_PP_HHEADER   = TOKEN_FROM_ID(401, PPTokenType),
    T_PP_INCLUDE_NEXT   = TOKEN_FROM_ID(399, PPTokenType|AltTokenType),
    T_PP_QHEADER_NEXT   = TOKEN_FROM_ID(400, PPTokenType|AltTokenType),
    T_PP_HHEADER_NEXT   = TOKEN_FROM_ID(401, PPTokenType|AltTokenType),
    T_EOF          = TOKEN_FROM_ID(402, EOFTokenType),      // end of file reached
    T_EOI          = TOKEN_FROM_ID(403, EOFTokenType),      // end of input reached
    
// MS extensions
    T_MSEXT_INT8   = TOKEN_FROM_ID(404, KeywordTokenType),
    T_MSEXT_INT16  = TOKEN_FROM_ID(405, KeywordTokenType),
    T_MSEXT_INT32  = TOKEN_FROM_ID(406, KeywordTokenType),
    T_MSEXT_INT64  = TOKEN_FROM_ID(407, KeywordTokenType),
    T_MSEXT_BASED  = TOKEN_FROM_ID(408, KeywordTokenType),
    T_MSEXT_DECLSPEC = TOKEN_FROM_ID(409, KeywordTokenType),
    T_MSEXT_CDECL  = TOKEN_FROM_ID(410, KeywordTokenType),
    T_MSEXT_FASTCALL = TOKEN_FROM_ID(411, KeywordTokenType),
    T_MSEXT_STDCALL = TOKEN_FROM_ID(412, KeywordTokenType),
    T_MSEXT_TRY    = TOKEN_FROM_ID(413, KeywordTokenType),
    T_MSEXT_EXCEPT = TOKEN_FROM_ID(414, KeywordTokenType),
    T_MSEXT_FINALLY = TOKEN_FROM_ID(415, KeywordTokenType),
    T_MSEXT_LEAVE  = TOKEN_FROM_ID(416, KeywordTokenType),
    T_MSEXT_INLINE = TOKEN_FROM_ID(417, KeywordTokenType),
    T_MSEXT_ASM    = TOKEN_FROM_ID(418, KeywordTokenType),

    T_MSEXT_PP_REGION    = TOKEN_FROM_ID(419, PPTokenType),
    T_MSEXT_PP_ENDREGION = TOKEN_FROM_ID(420, PPTokenType),
    
    T_LAST_TOKEN_ID,
    T_LAST_TOKEN = ID_FROM_TOKEN(T_LAST_TOKEN_ID),
    
// pseudo tokens to help streamlining macro replacement, these should not 
// returned from the lexer nor should these be returned from the pp-iterator
    T_NONREPLACABLE_IDENTIFIER = TOKEN_FROM_ID(T_LAST_TOKEN+1, IdentifierTokenType),
    T_PLACEHOLDER = TOKEN_FROM_ID(T_LAST_TOKEN+2, WhiteSpaceTokenType),
    T_PLACEMARKER = TOKEN_FROM_ID(T_LAST_TOKEN+3, InternalTokenType),
    T_PARAMETERBASE = TOKEN_FROM_ID(T_LAST_TOKEN+4, ParameterTokenType),
    T_EXTPARAMETERBASE = TOKEN_FROM_ID(T_LAST_TOKEN+5, ExtParameterTokenType)
};

///////////////////////////////////////////////////////////////////////////////
//  redefine the TOKEN_FROM_ID macro to be more type safe
#undef TOKEN_FROM_ID
#define TOKEN_FROM_ID(id, cat)   boost::wave::token_id((id) | (cat))
#define BASE_TOKEN(tok)                                                       \
    boost::wave::token_id((tok) & MainTokenMask)                              \
  /**/
#undef BASEID_FROM_TOKEN
#define BASEID_FROM_TOKEN(tok)                                                \
    boost::wave::token_id(((tok) & ~ExtTokenTypeMask))                        \
  /**/

///////////////////////////////////////////////////////////////////////////////
//  return a token name
inline BOOST_WAVE_STRINGTYPE
get_token_name(token_id tokid)
{
//  Table of token names
//
//      Please note that the sequence of token names must match the sequence of
//      token id's defined in then enum token_id above.
static char const *tok_names[] = {
    /* 256 */   "AND",
    /* 257 */   "ANDAND",
    /* 258 */   "ASSIGN",
    /* 259 */   "ANDASSIGN",
    /* 260 */   "OR",
    /* 261 */   "ORASSIGN",
    /* 262 */   "XOR",
    /* 263 */   "XORASSIGN",
    /* 264 */   "COMMA",
    /* 265 */   "COLON",
    /* 266 */   "DIVIDE",
    /* 267 */   "DIVIDEASSIGN",
    /* 268 */   "DOT",
    /* 269 */   "DOTSTAR",
    /* 270 */   "ELLIPSIS",
    /* 271 */   "EQUAL",
    /* 272 */   "GREATER",
    /* 273 */   "GREATEREQUAL",
    /* 274 */   "LEFTBRACE",
    /* 275 */   "LESS",
    /* 276 */   "LESSEQUAL",
    /* 277 */   "LEFTPAREN",
    /* 278 */   "LEFTBRACKET",
    /* 279 */   "MINUS",
    /* 280 */   "MINUSASSIGN",
    /* 281 */   "MINUSMINUS",
    /* 282 */   "PERCENT",
    /* 283 */   "PERCENTASSIGN",
    /* 284 */   "NOT",
    /* 285 */   "NOTEQUAL",
    /* 286 */   "OROR",
    /* 287 */   "PLUS",
    /* 288 */   "PLUSASSIGN",
    /* 289 */   "PLUSPLUS",
    /* 290 */   "ARROW",
    /* 291 */   "ARROWSTAR",
    /* 292 */   "QUESTION_MARK",
    /* 293 */   "RIGHTBRACE",
    /* 294 */   "RIGHTPAREN",
    /* 295 */   "RIGHTBRACKET",
    /* 296 */   "COLON_COLON",
    /* 297 */   "SEMICOLON",
    /* 298 */   "SHIFTLEFT",
    /* 299 */   "SHIFTLEFTASSIGN",
    /* 300 */   "SHIFTRIGHT",
    /* 301 */   "SHIFTRIGHTASSIGN",
    /* 302 */   "STAR",
    /* 303 */   "COMPL",
    /* 304 */   "STARASSIGN",
    /* 305 */   "ASM",
    /* 306 */   "AUTO",
    /* 307 */   "BOOL",
    /* 308 */   "FALSE",
    /* 309 */   "TRUE",
    /* 310 */   "BREAK",
    /* 311 */   "CASE",
    /* 312 */   "CATCH",
    /* 313 */   "CHAR",
    /* 314 */   "CLASS",
    /* 315 */   "CONST",
    /* 316 */   "CONSTCAST",
    /* 317 */   "CONTINUE",
    /* 318 */   "DEFAULT",
    /* 319 */   "DEFINED",
    /* 320 */   "DELETE",
    /* 321 */   "DO",
    /* 322 */   "DOUBLE",
    /* 323 */   "DYNAMICCAST",
    /* 324 */   "ELSE",
    /* 325 */   "ENUM",
    /* 326 */   "EXPLICIT",
    /* 327 */   "EXPORT",
    /* 328 */   "EXTERN",
    /* 329 */   "FLOAT",
    /* 330 */   "FOR",
    /* 331 */   "FRIEND",
    /* 332 */   "GOTO",
    /* 333 */   "IF",
    /* 334 */   "INLINE",
    /* 335 */   "INT",
    /* 336 */   "LONG",
    /* 337 */   "MUTABLE",
    /* 338 */   "NAMESPACE",
    /* 339 */   "NEW",
    /* 340 */   "OPERATOR",
    /* 341 */   "PRIVATE",
    /* 342 */   "PROTECTED",
    /* 343 */   "PUBLIC",
    /* 344 */   "REGISTER",
    /* 345 */   "REINTERPRETCAST",
    /* 346 */   "RETURN",
    /* 347 */   "SHORT",
    /* 348 */   "SIGNED",
    /* 349 */   "SIZEOF",
    /* 350 */   "STATIC",
    /* 351 */   "STATICCAST",
    /* 352 */   "STRUCT",
    /* 353 */   "SWITCH",
    /* 354 */   "TEMPLATE",
    /* 355 */   "THIS",
    /* 356 */   "THROW",
    /* 357 */   "TRY",
    /* 358 */   "TYPEDEF",
    /* 359 */   "TYPEID",
    /* 360 */   "TYPENAME",
    /* 361 */   "UNION",
    /* 362 */   "UNSIGNED",
    /* 363 */   "USING",
    /* 364 */   "VIRTUAL",
    /* 365 */   "VOID",
    /* 366 */   "VOLATILE",
    /* 367 */   "WCHART",
    /* 368 */   "WHILE",
    /* 369 */   "PP_DEFINE",
    /* 370 */   "PP_IF",
    /* 371 */   "PP_IFDEF",
    /* 372 */   "PP_IFNDEF",
    /* 373 */   "PP_ELSE",
    /* 374 */   "PP_ELIF",
    /* 375 */   "PP_ENDIF",
    /* 376 */   "PP_ERROR",
    /* 377 */   "PP_LINE",
    /* 378 */   "PP_PRAGMA",
    /* 379 */   "PP_UNDEF",
    /* 380 */   "PP_WARNING",
    /* 381 */   "IDENTIFIER",
    /* 382 */   "OCTALINT",
    /* 383 */   "DECIMALINT",
    /* 384 */   "HEXAINT",
    /* 385 */   "INTLIT",
    /* 386 */   "LONGINTLIT",
    /* 387 */   "FLOATLIT",
    /* 388 */   "CCOMMENT",
    /* 389 */   "CPPCOMMENT",
    /* 390 */   "CHARLIT",
    /* 391 */   "STRINGLIT",
    /* 392 */   "CONTLINE",
    /* 393 */   "SPACE",
    /* 394 */   "SPACE2",
    /* 395 */   "NEWLINE",
    /* 396 */   "POUND_POUND",
    /* 397 */   "POUND",
    /* 398 */   "ANY",
    /* 399 */   "PP_INCLUDE",
    /* 400 */   "PP_QHEADER",
    /* 401 */   "PP_HHEADER",
    /* 402 */   "EOF",
    /* 403 */   "EOI",

              // MS extensions
    /* 404 */   "MSEXT_INT8",
    /* 405 */   "MSEXT_INT16",
    /* 406 */   "MSEXT_INT32",
    /* 407 */   "MSEXT_INT64",
    /* 408 */   "MSEXT_BASED",
    /* 409 */   "MSEXT_DECLSPEC",
    /* 410 */   "MSEXT_CDECL",
    /* 411 */   "MSEXT_FASTCALL",
    /* 412 */   "MSEXT_STDCALL",
    /* 413 */   "MSEXT_TRY",
    /* 414 */   "MSEXT_EXCEPT",
    /* 415 */   "MSEXT_FINALLY",
    /* 416 */   "MSEXT_LEAVE",
    /* 417 */   "MSEXT_INLINE",
    /* 418 */   "MSEXT_ASM",
    /* 419 */   "MSEXT_REGION",
    /* 420 */   "MSEXT_ENDREGION",
    };   
     
    unsigned int id = BASEID_FROM_TOKEN(tokid)-T_FIRST_TOKEN;

    BOOST_ASSERT(id < T_LAST_TOKEN-T_FIRST_TOKEN);
    return tok_names[id];
}

///////////////////////////////////////////////////////////////////////////////
//  return a token name
inline char const *
get_token_value(token_id tokid)
{
//  Table of token values
//
//      Please note that the sequence of token names must match the sequence of
//      token id's defined in then enum token_id above.
static char const *tok_names[] = {
    /* 256 */   "&",
    /* 257 */   "&&",
    /* 258 */   "=",
    /* 259 */   "&=",
    /* 260 */   "|",
    /* 261 */   "|=",
    /* 262 */   "^",
    /* 263 */   "^=",
    /* 264 */   ",",
    /* 265 */   ":",
    /* 266 */   "/",
    /* 267 */   "/=",
    /* 268 */   ".",
    /* 269 */   ".*",
    /* 270 */   "...",
    /* 271 */   "==",
    /* 272 */   ">",
    /* 273 */   ">=",
    /* 274 */   "{",
    /* 275 */   "<",
    /* 276 */   "<=",
    /* 277 */   "(",
    /* 278 */   "[",
    /* 279 */   "-",
    /* 280 */   "-=",
    /* 281 */   "--",
    /* 282 */   "%",
    /* 283 */   "%=",
    /* 284 */   "!",
    /* 285 */   "!=",
    /* 286 */   "||",
    /* 287 */   "+",
    /* 288 */   "+=",
    /* 289 */   "++",
    /* 290 */   "->",
    /* 291 */   "->*",
    /* 292 */   "?",
    /* 293 */   "}",
    /* 294 */   ")",
    /* 295 */   "]",
    /* 296 */   "::",
    /* 297 */   ";",
    /* 298 */   "<<",
    /* 299 */   "<<=",
    /* 300 */   ">>",
    /* 301 */   ">>=",
    /* 302 */   "*",
    /* 303 */   "~",
    /* 304 */   "*=",
    /* 305 */   "asm",
    /* 306 */   "auto",
    /* 307 */   "bool",
    /* 308 */   "false",
    /* 309 */   "true",
    /* 310 */   "break",
    /* 311 */   "case",
    /* 312 */   "catch",
    /* 313 */   "char",
    /* 314 */   "class",
    /* 315 */   "const",
    /* 316 */   "const_cast",
    /* 317 */   "continue",
    /* 318 */   "default",
    /* 319 */   "defined",
    /* 320 */   "delete",
    /* 321 */   "do",
    /* 322 */   "double",
    /* 323 */   "dynamic_cast",
    /* 324 */   "else",
    /* 325 */   "enum",
    /* 326 */   "explicit",
    /* 327 */   "export",
    /* 328 */   "extern",
    /* 329 */   "float",
    /* 330 */   "for",
    /* 331 */   "friend",
    /* 332 */   "goto",
    /* 333 */   "if",
    /* 334 */   "inline",
    /* 335 */   "int",
    /* 336 */   "long",
    /* 337 */   "mutable",
    /* 338 */   "namespace",
    /* 339 */   "new",
    /* 340 */   "operator",
    /* 341 */   "private",
    /* 342 */   "protected",
    /* 343 */   "public",
    /* 344 */   "register",
    /* 345 */   "reinterpret_cast",
    /* 346 */   "return",
    /* 347 */   "short",
    /* 348 */   "signed",
    /* 349 */   "sizeof",
    /* 350 */   "static",
    /* 351 */   "static_cast",
    /* 352 */   "struct",
    /* 353 */   "switch",
    /* 354 */   "template",
    /* 355 */   "this",
    /* 356 */   "throw",
    /* 357 */   "try",
    /* 358 */   "typedef",
    /* 359 */   "typeid",
    /* 360 */   "typename",
    /* 361 */   "union",
    /* 362 */   "unsigned",
    /* 363 */   "using",
    /* 364 */   "virtual",
    /* 365 */   "void",
    /* 366 */   "volatile",
    /* 367 */   "wchar_t",
    /* 368 */   "while",
    /* 369 */   "#define",
    /* 370 */   "#if",
    /* 371 */   "#ifdef",
    /* 372 */   "#ifndef",
    /* 373 */   "#else",
    /* 374 */   "#elif",
    /* 375 */   "#endif",
    /* 376 */   "#error",
    /* 377 */   "#line",
    /* 378 */   "#pragma ",
    /* 379 */   "#undef ",
    /* 380 */   "#warning",
    /* 381 */   "",   // identifier
    /* 382 */   "",   // octalint
    /* 383 */   "",   // decimalint
    /* 384 */   "",   // hexlit
    /* 385 */   "",   // intlit
    /* 386 */   "",   // longintlit
    /* 387 */   "",   // floatlit
    /* 388 */   "",   // ccomment
    /* 389 */   "",   // cppcomment
    /* 390 */   "",   // charlit
    /* 391 */   "",   // stringlit
    /* 392 */   "",   // contline
    /* 393 */   "",   // space
    /* 394 */   "",   // space2
    /* 395 */   "\n",
    /* 396 */   "##",
    /* 397 */   "#",
    /* 398 */   "",   // any
    /* 399 */   "#include",
    /* 400 */   "#include",
    /* 401 */   "#include",
    /* 402 */   "",   // eof
    /* 403 */   "",   // eoi

              // MS extensions
    /* 404 */   "__int8",
    /* 405 */   "__int16",
    /* 406 */   "__int32",
    /* 407 */   "__int64",
    /* 408 */   "__based",
    /* 409 */   "__declspec",
    /* 410 */   "__cdecl",
    /* 411 */   "__fastcall",
    /* 412 */   "__stdcall",
    /* 413 */   "__try",
    /* 414 */   "__except",
    /* 415 */   "__finally",
    /* 416 */   "__leave",
    /* 417 */   "__inline",
    /* 418 */   "__asm",
    /* 419 */   "#region",
    /* 420 */   "#endregion",
    };   
     
    unsigned int id = BASEID_FROM_TOKEN(tokid)-T_FIRST_TOKEN;

    BOOST_ASSERT(id < T_LAST_TOKEN-T_FIRST_TOKEN);
    return tok_names[id];
}

///////////////////////////////////////////////////////////////////////////////
}   // namespace wave
}   // namespace boost

#endif // #if !defined(BOOST_WAVE_TOKEN_IDS_DEFINED)

#endif // !defined(TOKEN_IDS_HPP_414E9A58_F079_4789_8AFF_513815CE475B_INCLUDED)

