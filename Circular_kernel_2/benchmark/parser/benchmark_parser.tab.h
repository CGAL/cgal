/* A Bison parser, made by GNU Bison 1.875c.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

   SPDX-License-Identifier: GPL-2.0-or-later

*/

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     INTEGER = 258,
     FNUMBER = 259,
     STRING = 260,
     ERROR = 261,
     UNKNOWN_TOKEN = 262,
     MINUS_INFTY = 263,
     PLUS_INFTY = 264,
     COUNTERCLOCKWISE = 265,
     CLOCKWISE = 266,
     VOID = 267,
     FileFormat = 268,
     BenchmarkName = 269,
     Classification = 270,
     List = 271,
     Rational = 272,
     Polynomial_1 = 273,
     Point_2 = 274,
     AlgebraicReal = 275,
     ConicPoint_2 = 276,
     LineSegment_2 = 277,
     Conic_2 = 278,
     CircularArc_2 = 279,
     LineArc_2 = 280,
     CircularPoint_2 = 281,
     ConicArc_2 = 282,
     Circle_2 = 283,
     Cubic_2 = 284,
     Quadric_3 = 285
   };
#endif
#define INTEGER 258
#define FNUMBER 259
#define STRING 260
#define ERROR 261
#define UNKNOWN_TOKEN 262
#define MINUS_INFTY 263
#define PLUS_INFTY 264
#define COUNTERCLOCKWISE 265
#define CLOCKWISE 266
#define VOID 267
#define FileFormat 268
#define BenchmarkName 269
#define Classification 270
#define List 271
#define Rational 272
#define Polynomial_1 273
#define Point_2 274
#define AlgebraicReal 275
#define ConicPoint_2 276
#define LineSegment_2 277
#define Conic_2 278
#define CircularArc_2 279
#define LineArc_2 280
#define CircularPoint_2 281
#define ConicArc_2 282
#define Circle_2 283
#define Cubic_2 284
#define Quadric_3 285




#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE yylval;



