/* A Bison parser, made by GNU Bison 1.875c.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003 Free Software Foundation, Inc.

   SPDX-License-Identifier: GPL-2.0-or-later

*/

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



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




/* Copy the first part of user declarations.  */
#line 21 "benchmark_parser.y"

/* C/C++ declaration section */
/* ========================= */
//#include <stdlib.h>  /* for atoi */
#include <string>    /* for std::string */
#include <fstream>   /* for std::ifstream */
#include <benchmark_parser.h>
#include <benchmark_visitor.h>
#include <iostream>
static Benchmark_visitor* visitor; /* global visitor used during parsing */

/* declaration for flex parser call yylex */
int yylex( void);

/* error function called for parse errors */
void yyerror( char *s) { visitor->parse_error( std::string(s)); }

/* Use C++ std::string as semantic value to communicate with lexer */
#define YYSTYPE std::string

// Public parser interface, check with decl. in benchmark_parser.h
// ---------------------------------------------------------------

// Opens file 'name' and parses it. Uses visitor 'v' while parsing.
// Returns false if something went wrong. See the visitor for details.
bool benchmark_parse_file( std::string name, Benchmark_visitor* v);

// Starts parsing from stream 'in' with the associated filename 'name'
// (or analogous meaning for different streams) counting linenumbers
// starting from 'n'. Uses visitor 'v' while parsing. Returns false if
// something went wrong. See the visitor for details of the error reporting.
bool benchmark_parse_stream( std::istream&  in, std::string name,
                             Benchmark_visitor* v, int n);



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 214 of yacc.c.  */
#line 183 "benchmark_parser.tab.c"

#if ! defined (yyoverflow) || YYERROR_VERBOSE

# ifndef YYFREE
#  define YYFREE free
# endif
# ifndef YYMALLOC
#  define YYMALLOC malloc
# endif

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   define YYSTACK_ALLOC alloca
#  endif
# else
#  if defined (alloca) || defined (_ALLOCA_H)
#   define YYSTACK_ALLOC alloca
#  else
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
         || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))                                \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)                \
      do                                        \
        {                                        \
          YYSIZE_T yyi;                \
          for (yyi = 0; yyi < (Count); yyi++)        \
            (To)[yyi] = (From)[yyi];                \
        }                                        \
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)                                        \
    do                                                                        \
      {                                                                        \
        YYSIZE_T yynewbytes;                                                \
        YYCOPY (&yyptr->Stack, Stack, yysize);                                \
        Stack = &yyptr->Stack;                                                \
        yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / sizeof (*yyptr);                                \
      }                                                                        \
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   322

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  34
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  43
/* YYNRULES -- Number of rules. */
#define YYNRULES  92
/* YYNRULES -- Number of states. */
#define YYNSTATES  284

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   285

#define YYTRANSLATE(YYX)                                                 \
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      31,    33,     2,     2,    32,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned short yyprhs[] =
{
       0,     0,     3,     4,     6,     9,    14,    16,    25,    36,
      37,    40,    45,    47,    49,    52,    67,    68,    71,    72,
      74,    76,    80,    81,    85,    86,    92,    93,    99,   102,
     103,   111,   114,   137,   160,   161,   167,   168,   174,   175,
     181,   182,   184,   186,   190,   192,   194,   198,   208,   210,
     215,   230,   236,   238,   242,   246,   248,   262,   264,   265,
     266,   281,   282,   288,   290,   292,   296,   298,   305,   314,
     331,   333,   334,   342,   344,   348,   352,   354,   355,   356,
     357,   358,   376,   378,   385,   387,   389,   391,   393,   395,
     397,   399,   401
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      35,     0,    -1,    -1,     1,    -1,    35,    36,    -1,    37,
      38,    40,    42,    -1,    76,    -1,    13,    31,     5,    32,
       3,    32,     3,    33,    -1,    13,    31,     5,    32,     3,
      32,     3,    32,     5,    33,    -1,    -1,    38,    39,    -1,
      14,    31,     5,    33,    -1,    76,    -1,    41,    -1,    40,
      41,    -1,    15,    31,     5,    32,     5,    32,     5,    32,
       5,    32,     5,    32,     5,    33,    -1,    -1,    42,    45,
      -1,    -1,    44,    -1,    45,    -1,    44,    32,    45,    -1,
      -1,     1,    46,    76,    -1,    -1,    16,    47,    31,    43,
      33,    -1,    -1,    28,    48,    31,    56,    33,    -1,    23,
      57,    -1,    -1,    22,    49,    31,    63,    32,    63,    33,
      -1,    27,    58,    -1,    29,    31,     3,    32,     3,    32,
       3,    32,     3,    32,     3,    32,     3,    32,     3,    32,
       3,    32,     3,    32,     3,    33,    -1,    30,    31,     3,
      32,     3,    32,     3,    32,     3,    32,     3,    32,     3,
      32,     3,    32,     3,    32,     3,    32,     3,    33,    -1,
      -1,    26,    50,    31,    53,    33,    -1,    -1,    25,    51,
      31,    54,    33,    -1,    -1,    24,    52,    31,    55,    33,
      -1,    -1,    76,    -1,    63,    -1,    20,    32,    20,    -1,
      76,    -1,    22,    -1,    63,    32,    63,    -1,    26,    31,
      53,    33,    32,    26,    31,    53,    33,    -1,    76,    -1,
      28,    31,    56,    33,    -1,    28,    31,    56,    33,    32,
      26,    31,    53,    33,    32,    26,    31,    53,    33,    -1,
      63,    32,    63,    32,    72,    -1,    76,    -1,    63,    32,
      72,    -1,    63,    32,     3,    -1,    76,    -1,    31,     3,
      32,     3,    32,     3,    32,     3,    32,     3,    32,     3,
      33,    -1,    76,    -1,    -1,    -1,    31,    23,    59,    57,
      32,    21,    64,    32,    21,    64,    32,    75,    60,    33,
      -1,    -1,    31,    21,    61,    64,    33,    -1,    76,    -1,
       3,    -1,    62,    32,     3,    -1,    76,    -1,    19,    31,
       3,    32,     3,    33,    -1,    19,    31,     3,    32,     3,
      32,     3,    33,    -1,    19,    31,    17,    31,     3,    32,
       3,    33,    32,    17,    31,     3,    32,     3,    33,    33,
      -1,    76,    -1,    -1,    31,    23,    65,    57,    32,    66,
      33,    -1,    76,    -1,    67,    32,    73,    -1,    74,    32,
       3,    -1,    76,    -1,    -1,    -1,    -1,    -1,    20,    68,
      31,    18,    69,    31,    62,    33,    70,    32,    72,    32,
      72,    32,     3,    71,    33,    -1,    76,    -1,    17,    31,
       3,    32,     3,    33,    -1,     3,    -1,    74,    -1,     8,
      -1,     9,    -1,    10,    -1,    11,    -1,    12,    -1,     6,
      -1,     7,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short yyrline[] =
{
       0,   103,   103,   105,   106,   110,   114,   115,   119,   125,
     127,   131,   135,   136,   137,   141,   147,   149,   152,   154,
     158,   159,   163,   163,   165,   165,   167,   167,   169,   170,
     170,   173,   175,   180,   186,   186,   189,   189,   192,   192,
     196,   198,   199,   200,   203,   204,   205,   206,   210,   211,
     212,   214,   217,   218,   219,   223,   225,   233,   234,   238,
     234,   243,   243,   250,   251,   252,   257,   258,   260,   262,
     269,   270,   270,   278,   279,   280,   286,   287,   288,   290,
     294,   287,   299,   300,   305,   306,   310,   311,   315,   316,
     317,   326,   327
};
#endif

#if YYDEBUG || YYERROR_VERBOSE
/* YYTNME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "INTEGER", "FNUMBER", "STRING", "ERROR",
  "UNKNOWN_TOKEN", "MINUS_INFTY", "PLUS_INFTY", "COUNTERCLOCKWISE",
  "CLOCKWISE", "VOID", "FileFormat", "BenchmarkName", "Classification",
  "List", "Rational", "Polynomial_1", "Point_2", "AlgebraicReal",
  "ConicPoint_2", "LineSegment_2", "Conic_2", "CircularArc_2", "LineArc_2",
  "CircularPoint_2", "ConicArc_2", "Circle_2", "Cubic_2", "Quadric_3",
  "'('", "','", "')'", "$accept", "input", "file", "file_format",
  "file_header_options", "file_header_option", "file_classification",
  "classificat", "file_body", "stmt_sequence", "stmt_sequence_non_empty",
  "stmt", "@1", "@2", "@3", "@4", "@5", "@6", "@7", "circular_arc_point",
  "line_arc_2", "circular_arc_2", "circle_2", "conic_2", "conic_arc_2",
  "@8", "@9", "@10", "integer_sequence1", "point_2", "conic_point_2",
  "@11", "algorint", "algebraic_real", "@12", "@13", "@14", "@15",
  "rational", "inti", "infty", "orientation", "error_rules", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,    40,    44,    41
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    34,    35,    35,    35,    36,    37,    37,    37,    38,
      38,    39,    40,    40,    40,    41,    42,    42,    43,    43,
      44,    44,    46,    45,    47,    45,    48,    45,    45,    49,
      45,    45,    45,    45,    50,    45,    51,    45,    52,    45,
      53,    53,    53,    53,    54,    54,    54,    54,    55,    55,
      55,    55,    56,    56,    56,    57,    57,    58,    59,    60,
      58,    61,    58,    62,    62,    62,    63,    63,    63,    63,
      64,    65,    64,    66,    66,    66,    67,    68,    69,    70,
      71,    67,    72,    72,    73,    73,    74,    74,    75,    75,
      75,    76,    76
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     1,     2,     4,     1,     8,    10,     0,
       2,     4,     1,     1,     2,    14,     0,     2,     0,     1,
       1,     3,     0,     3,     0,     5,     0,     5,     2,     0,
       7,     2,    22,    22,     0,     5,     0,     5,     0,     5,
       0,     1,     1,     3,     1,     1,     3,     9,     1,     4,
      14,     5,     1,     3,     3,     1,    13,     1,     0,     0,
      14,     0,     5,     1,     1,     3,     1,     6,     8,    16,
       1,     0,     7,     1,     3,     3,     1,     0,     0,     0,
       0,    17,     1,     6,     1,     1,     1,     1,     1,     1,
       1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       0,     3,     0,     1,    91,    92,     0,     4,     9,     6,
       0,     0,     0,     0,     0,    10,    16,    13,    12,     0,
       0,     0,    14,     0,     0,     0,     0,    22,    24,    29,
       0,    38,    36,    34,     0,    26,     0,     0,    17,     0,
      11,     0,     0,     0,     0,     0,    28,    55,     0,     0,
       0,     0,    31,    57,     0,     0,     0,     0,     0,    23,
       0,     0,     0,     0,     0,    40,    61,    58,     0,     0,
       0,     0,     7,     0,     0,    19,    20,     0,     0,    66,
       0,     0,     0,     0,    48,    45,     0,     0,     0,    44,
       0,     0,    42,    41,     0,     0,     0,     0,    52,     0,
       0,     0,     0,    25,     0,     0,     0,     0,     0,    39,
       0,    40,    37,     0,     0,    35,     0,     0,    70,     0,
      27,     0,     0,     0,     8,     0,    21,     0,     0,     0,
       0,     0,     0,     0,    46,    43,    71,    62,     0,    54,
       0,    53,    82,     0,     0,     0,     0,     0,    30,     0,
      49,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    51,     0,     0,     0,     0,     0,     0,
       0,     0,    67,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    40,    40,    86,    87,
      77,     0,     0,     0,    73,     0,     0,     0,     0,     0,
      68,     0,     0,     0,     0,     0,    72,     0,     0,     0,
      83,     0,     0,    15,     0,     0,     0,    47,     0,    84,
      74,    85,    75,     0,     0,     0,     0,     0,     0,    78,
      88,    89,    90,    59,     0,     0,     0,    56,     0,     0,
       0,     0,     0,     0,    40,     0,    60,     0,     0,     0,
       0,    64,     0,    63,     0,     0,     0,    50,     0,    79,
       0,     0,     0,    65,     0,     0,     0,    69,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    32,    33,
       0,    80,     0,    81
};

/* YYDEFGOTO[NTERM-NUM]. */
static const short yydefgoto[] =
{
      -1,     2,     7,     8,    11,    15,    16,    17,    23,    74,
      75,    38,    42,    43,    54,    44,    50,    49,    48,    91,
      87,    82,    96,    46,    52,    95,   240,    94,   252,    92,
     117,   153,   191,   192,   205,   239,   264,   282,   141,   220,
     193,   233,    93
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -150
static const short yypact[] =
{
     117,  -150,   147,  -150,  -150,  -150,   -21,  -150,  -150,  -150,
      25,   149,     3,    19,    40,  -150,    30,  -150,  -150,   102,
      98,   105,  -150,   115,    75,    81,    93,  -150,  -150,  -150,
       6,  -150,  -150,  -150,     8,  -150,    96,   119,  -150,   126,
    -150,   143,    20,   128,   130,   159,  -150,  -150,   137,   138,
     139,    42,  -150,  -150,   140,   169,   170,    58,   142,  -150,
      56,    68,   144,    14,    12,    48,  -150,  -150,    68,   145,
     146,   174,  -150,   175,   148,   150,  -150,   152,   153,  -150,
     181,   156,   157,   154,   160,  -150,   158,   161,   163,   160,
     164,   165,  -150,  -150,    10,     6,   166,   168,   160,   185,
     188,   171,   173,  -150,    72,    26,    68,   176,    68,  -150,
      68,    48,  -150,    68,   177,  -150,   178,   179,  -150,   182,
    -150,    41,   183,   184,  -150,   197,  -150,   186,   172,   180,
     203,   187,   189,   190,  -150,  -150,  -150,  -150,   196,  -150,
     191,  -150,  -150,   204,   206,   192,   207,   208,  -150,   193,
     194,   129,   195,     6,    10,   216,   198,   199,   223,   100,
     200,   226,   209,  -150,   210,   201,   202,   205,   235,   236,
     212,   237,  -150,   238,   213,   211,   215,    16,   227,   244,
     217,   218,   246,   219,   220,   251,    48,    48,  -150,  -150,
    -150,   222,   224,   225,   228,    10,   229,   255,   256,   230,
    -150,   232,   233,   234,   240,   239,  -150,    61,   258,   242,
    -150,   243,   245,  -150,   252,   265,   247,  -150,   253,  -150,
    -150,  -150,  -150,   155,   273,   275,   249,   248,   257,  -150,
    -150,  -150,  -150,  -150,   250,   254,   281,  -150,   259,   260,
     261,   282,   284,   263,    48,   106,  -150,   264,   266,   285,
     267,  -150,   125,  -150,   286,   289,   268,  -150,   290,  -150,
     270,   271,   272,  -150,   274,   294,   296,  -150,   129,   276,
     277,   278,   301,   304,   129,   279,   280,   283,  -150,  -150,
     308,  -150,   287,  -150
};

/* YYPGOTO[NTERM-NUM].  */
static const short yypgoto[] =
{
    -150,  -150,  -150,  -150,  -150,  -150,  -150,   298,  -150,  -150,
    -150,   -52,  -150,  -150,  -150,  -150,  -150,  -150,  -150,  -110,
    -150,  -150,   214,   -93,  -150,  -150,  -150,  -150,  -150,   -57,
    -149,  -150,  -150,  -150,  -150,  -150,  -150,  -150,  -148,  -150,
     109,  -150,    -2
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -77
static const short yytable[] =
{
       9,   133,   119,   163,    78,   166,    83,    88,    76,    18,
      10,    97,     4,     5,     4,     5,     4,     5,     4,     5,
       4,     5,     4,     5,   188,   189,     4,     5,    47,   127,
      12,    77,    53,    77,    85,    19,   190,    45,    86,    51,
      59,   116,    81,   128,   139,    14,   209,     4,     5,   129,
      20,    97,   126,   132,     4,     5,   134,    27,   140,    79,
     165,    84,    89,    66,   219,    67,    98,    77,    90,   188,
     189,    21,    28,    27,     4,     5,   203,   204,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    77,    28,   -18,
      71,    72,   118,    47,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    25,    79,    24,    98,    39,    79,   251,
      26,    79,     4,     5,    40,    -5,    27,    -2,     1,   142,
     271,    -5,    -5,    -2,    -2,    41,   277,    55,    -5,    57,
      -2,    28,   171,   172,   250,     4,     5,    29,    30,    31,
      32,    33,    34,    35,    36,    37,   140,     3,    58,   142,
      56,    47,   118,     4,     5,     4,     5,   258,   259,    60,
       6,    61,    62,    13,    14,   230,   231,   232,    63,    64,
      65,    68,    69,    70,    73,   194,    80,    99,   100,   101,
     102,   103,   104,   105,   107,   106,   110,   108,   122,   111,
     109,   123,   -66,   118,   112,   113,   114,   135,   115,   120,
     121,   136,   145,   147,   124,   125,   149,   156,   130,   157,
     159,   160,   137,   148,   138,   143,   144,   154,   146,   167,
     150,   151,   155,   152,   158,   161,   162,   164,   170,   174,
     168,   169,   173,   177,   178,   175,   176,   179,   180,   181,
     183,   184,   186,   253,   182,   185,   187,   196,   195,   197,
     198,   199,   200,   201,   202,   206,   207,   208,   211,   212,
     -76,   222,   210,   213,   214,   215,   142,   216,   227,   226,
     218,   229,   142,   217,   223,   224,   234,   225,   235,   228,
     236,   237,   241,   238,   243,   247,   242,   248,   256,   260,
     244,   245,   261,   263,   246,   249,   254,   269,   255,   270,
     257,   262,   265,   266,   275,   267,   268,   276,   272,   273,
     274,   281,   278,   279,    22,   280,   221,     0,     0,     0,
     283,     0,   131
};

static const short yycheck[] =
{
       2,   111,    95,   151,    61,   154,    63,    64,    60,    11,
      31,    68,     6,     7,     6,     7,     6,     7,     6,     7,
       6,     7,     6,     7,     8,     9,     6,     7,    30,     3,
       5,    19,    34,    19,    22,    32,    20,    31,    26,    31,
      42,    31,    28,    17,     3,    15,   195,     6,     7,   106,
      31,   108,   104,   110,     6,     7,   113,     1,    17,    61,
     153,    63,    64,    21,     3,    23,    68,    19,    20,     8,
       9,    31,    16,     1,     6,     7,   186,   187,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    19,    16,    33,
      32,    33,    94,    95,    22,    23,    24,    25,    26,    27,
      28,    29,    30,     5,   106,     3,   108,    32,   110,     3,
       5,   113,     6,     7,    33,     0,     1,     0,     1,   121,
     268,     6,     7,     6,     7,    32,   274,    31,    13,     3,
      13,    16,    32,    33,   244,     6,     7,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    17,     0,     5,   151,
      31,   153,   154,     6,     7,     6,     7,    32,    33,    31,
      13,    31,     3,    14,    15,    10,    11,    12,    31,    31,
      31,    31,     3,     3,    32,   177,    32,    32,    32,     5,
       5,    33,    32,    31,     3,    32,    32,    31,     3,    31,
      33,     3,    32,   195,    33,    32,    32,    20,    33,    33,
      32,    23,     5,    31,    33,    32,     3,     3,    32,     3,
       3,     3,    33,    33,    32,    32,    32,    21,    32,     3,
      33,    32,    31,    33,    32,    32,    32,    32,     5,     3,
      32,    32,    32,    32,    32,    26,    26,    32,     3,     3,
       3,     3,    31,   245,    32,    32,    31,     3,    21,    32,
      32,     5,    33,    33,     3,    33,    32,    32,     3,     3,
      32,     3,    33,    33,    32,    32,   268,    33,     3,    17,
      31,    18,   274,    33,    32,    32,     3,    32,     3,    32,
      31,    33,    32,    26,     3,     3,    32,     3,     3,     3,
      31,    31,     3,     3,    33,    32,    32,     3,    32,     3,
      33,    33,    32,    32,     3,    33,    32,     3,    32,    32,
      32,     3,    33,    33,    16,    32,   207,    -1,    -1,    -1,
      33,    -1,   108
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,     1,    35,     0,     6,     7,    13,    36,    37,    76,
      31,    38,     5,    14,    15,    39,    40,    41,    76,    32,
      31,    31,    41,    42,     3,     5,     5,     1,    16,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    45,    32,
      33,    32,    46,    47,    49,    31,    57,    76,    52,    51,
      50,    31,    58,    76,    48,    31,    31,     3,     5,    76,
      31,    31,     3,    31,    31,    31,    21,    23,    31,     3,
       3,    32,    33,    32,    43,    44,    45,    19,    63,    76,
      32,    28,    55,    63,    76,    22,    26,    54,    63,    76,
      20,    53,    63,    76,    61,    59,    56,    63,    76,    32,
      32,     5,     5,    33,    32,    31,    32,     3,    31,    33,
      32,    31,    33,    32,    32,    33,    31,    64,    76,    57,
      33,    32,     3,     3,    33,    32,    45,     3,    17,    63,
      32,    56,    63,    53,    63,    20,    23,    33,    32,     3,
      17,    72,    76,    32,    32,     5,    32,    31,    33,     3,
      33,    32,    33,    65,    21,    31,     3,     3,    32,     3,
       3,    32,    32,    72,    32,    57,    64,     3,    32,    32,
       5,    32,    33,    32,     3,    26,    26,    32,    32,    32,
       3,     3,    32,     3,     3,    32,    31,    31,     8,     9,
      20,    66,    67,    74,    76,    21,     3,    32,    32,     5,
      33,    33,     3,    53,    53,    68,    33,    32,    32,    64,
      33,     3,     3,    33,    32,    32,    33,    33,    31,     3,
      73,    74,     3,    32,    32,    32,    17,     3,    32,    18,
      10,    11,    12,    75,     3,     3,    31,    33,    26,    69,
      60,    32,    32,     3,    31,    31,    33,     3,     3,    32,
      53,     3,    62,    76,    32,    32,     3,    33,    32,    33,
       3,     3,    33,     3,    70,    32,    32,    33,    32,     3,
       3,    72,    32,    32,    32,     3,     3,    72,    33,    33,
      32,     3,    71,    33
};

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok                (yyerrstatus = 0)
#define yyclearin        (yychar = YYEMPTY)
#define YYEMPTY                (-2)
#define YYEOF                0

#define YYACCEPT        goto yyacceptlab
#define YYABORT                goto yyabortlab
#define YYERROR                goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL                goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                        \
do                                                                \
  if (yychar == YYEMPTY && yylen == 1)                                \
    {                                                                \
      yychar = (Token);                                                \
      yylval = (Value);                                                \
      yytoken = YYTRANSLATE (yychar);                                \
      YYPOPSTACK;                                                \
      goto yybackup;                                                \
    }                                                                \
  else                                                                \
    {                                                                 \
      yyerror ("syntax error: cannot back up");\
      YYERROR;                                                        \
    }                                                                \
while (0)

#define YYTERROR        1
#define YYERRCODE        256

/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)                \
   ((Current).first_line   = (Rhs)[1].first_line,        \
    (Current).first_column = (Rhs)[1].first_column,        \
    (Current).last_line    = (Rhs)[N].last_line,        \
    (Current).last_column  = (Rhs)[N].last_column)
#endif

/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                                \
  if (yydebug)                                        \
    YYFPRINTF Args;                                \
} while (0)

# define YYDSYMPRINT(Args)                        \
do {                                                \
  if (yydebug)                                        \
    yysymprint Args;                                \
} while (0)

# define YYDSYMPRINTF(Title, Token, Value, Location)                \
do {                                                                \
  if (yydebug)                                                        \
    {                                                                \
      YYFPRINTF (stderr, "%s ", Title);                                \
      yysymprint (stderr,                                         \
                  Token, Value);        \
      YYFPRINTF (stderr, "\n");                                        \
    }                                                                \
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short *bottom, short *top)
#else
static void
yy_stack_print (bottom, top)
    short *bottom;
    short *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                                \
do {                                                                \
  if (yydebug)                                                        \
    yy_stack_print ((Bottom), (Top));                                \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %u), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname [yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname [yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)                \
do {                                        \
  if (yydebug)                                \
    yy_reduce_print (Rule);                \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YYDSYMPRINT(Args)
# define YYDSYMPRINTF(Title, Token, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef        YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if defined (YYMAXDEPTH) && YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

#endif /* !YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    {
      YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
# ifdef YYPRINT
      YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
    }
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yydestruct (int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yytype, yyvaluep)
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  switch (yytype)
    {

      default:
        break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{

  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short        yyssa[YYINITDEPTH];
  short *yyss = yyssa;
  short *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;                /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
        /* Give user a chance to reallocate the stack. Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        YYSTYPE *yyvs1 = yyvs;
        short *yyss1 = yyss;


        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow ("parser stack overflow",
                    &yyss1, yysize * sizeof (*yyssp),
                    &yyvs1, yysize * sizeof (*yyvsp),

                    &yystacksize);

        yyss = yyss1;
        yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyoverflowlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyoverflowlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        short *yyss1 = yyss;
        union yyalloc *yyptr =
          (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
        if (! yyptr)
          goto yyoverflowlab;
        YYSTACK_RELOCATE (yyss);
        YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YYDSYMPRINTF ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */
  YYDPRINTF ((stderr, "Shifting token %s, ", yytname[yytoken]));

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
#line 105 "benchmark_parser.y"
    { /* parse error restart here */ ;}
    break;

  case 6:
#line 114 "benchmark_parser.y"
    {;}
    break;

  case 7:
#line 116 "benchmark_parser.y"
    { visitor->accept_file_format( yyvsp[-5],
                                      atoi( yyvsp[-3].c_str()), atoi( yyvsp[-1].c_str()),
                                      std::string("")); ;}
    break;

  case 8:
#line 120 "benchmark_parser.y"
    { visitor->accept_file_format( yyvsp[-7],
                                     atoi( yyvsp[-5].c_str()), atoi( yyvsp[-3].c_str()),
                                                              yyvsp[-1]); ;}
    break;

  case 11:
#line 131 "benchmark_parser.y"
    { visitor->accept_benchmark_name( yyvsp[-1]); ;}
    break;

  case 12:
#line 135 "benchmark_parser.y"
    {;}
    break;

  case 15:
#line 143 "benchmark_parser.y"
    { visitor->accept_classification( yyvsp[-11], yyvsp[-9], yyvsp[-7],
                                                              yyvsp[-5], yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 22:
#line 163 "benchmark_parser.y"
    { /* parse error restart here */ ;}
    break;

  case 23:
#line 164 "benchmark_parser.y"
    {;}
    break;

  case 24:
#line 165 "benchmark_parser.y"
    { visitor->begin_list(); ;}
    break;

  case 25:
#line 166 "benchmark_parser.y"
    { visitor->end_list(); ;}
    break;

  case 26:
#line 167 "benchmark_parser.y"
    { visitor->begin_circle_2(); ;}
    break;

  case 27:
#line 168 "benchmark_parser.y"
    { visitor->end_circle_2(); ;}
    break;

  case 29:
#line 170 "benchmark_parser.y"
    { visitor->begin_line_segment_2(); ;}
    break;

  case 30:
#line 172 "benchmark_parser.y"
    { visitor->end_line_segment_2(); ;}
    break;

  case 32:
#line 177 "benchmark_parser.y"
    { visitor->accept_cubic_2( yyvsp[-19], yyvsp[-17], yyvsp[-15], yyvsp[-13],
                                                            yyvsp[-11], yyvsp[-9], yyvsp[-7],
                                                            yyvsp[-5], yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 33:
#line 182 "benchmark_parser.y"
    { visitor->accept_quadric_3( yyvsp[-19], yyvsp[-17], yyvsp[-15], yyvsp[-13],
                                                              yyvsp[-11], yyvsp[-9], yyvsp[-7],
                                                              yyvsp[-5], yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 34:
#line 186 "benchmark_parser.y"
    {visitor->begin_CircularPoint_2();;}
    break;

  case 35:
#line 188 "benchmark_parser.y"
    {visitor->end_CircularPoint_2();;}
    break;

  case 36:
#line 189 "benchmark_parser.y"
    {visitor->begin_LineArc_2();;}
    break;

  case 37:
#line 191 "benchmark_parser.y"
    {visitor->end_LineArc_2();;}
    break;

  case 38:
#line 192 "benchmark_parser.y"
    { visitor->begin_CircularArc_2();;}
    break;

  case 39:
#line 194 "benchmark_parser.y"
    { visitor->end_CircularArc_2();;}
    break;

  case 41:
#line 198 "benchmark_parser.y"
    {;}
    break;

  case 44:
#line 203 "benchmark_parser.y"
    {;}
    break;

  case 48:
#line 210 "benchmark_parser.y"
    {;}
    break;

  case 52:
#line 217 "benchmark_parser.y"
    {;}
    break;

  case 54:
#line 219 "benchmark_parser.y"
    { visitor->accept_integer(yyvsp[0]); ;}
    break;

  case 55:
#line 223 "benchmark_parser.y"
    {;}
    break;

  case 56:
#line 227 "benchmark_parser.y"
    { visitor->accept_conic_2( yyvsp[-11], yyvsp[-9], yyvsp[-7],
                                                            yyvsp[-5], yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 57:
#line 233 "benchmark_parser.y"
    {;}
    break;

  case 58:
#line 234 "benchmark_parser.y"
    { visitor->begin_conic_arc_2(); ;}
    break;

  case 59:
#line 238 "benchmark_parser.y"
    { visitor->accept_orientation( yyvsp[0]); ;}
    break;

  case 60:
#line 240 "benchmark_parser.y"
    { visitor->end_conic_arc_2(); ;}
    break;

  case 61:
#line 243 "benchmark_parser.y"
    { visitor->begin_conic_arc_2(); ;}
    break;

  case 62:
#line 245 "benchmark_parser.y"
    { visitor->end_conic_arc_2(); ;}
    break;

  case 63:
#line 250 "benchmark_parser.y"
    {;}
    break;

  case 64:
#line 251 "benchmark_parser.y"
    { visitor->accept_integer( yyvsp[0]); ;}
    break;

  case 65:
#line 252 "benchmark_parser.y"
    { visitor->accept_integer( yyvsp[0]); ;}
    break;

  case 66:
#line 257 "benchmark_parser.y"
    {;}
    break;

  case 67:
#line 259 "benchmark_parser.y"
    { visitor->accept_point_2( yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 68:
#line 261 "benchmark_parser.y"
    { visitor->accept_point_2( yyvsp[-5], yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 69:
#line 264 "benchmark_parser.y"
    { visitor->accept_point_2( yyvsp[-11], yyvsp[-9],
                                                             yyvsp[-4], yyvsp[-2]); ;}
    break;

  case 70:
#line 269 "benchmark_parser.y"
    {;}
    break;

  case 71:
#line 270 "benchmark_parser.y"
    { visitor->begin_conic_point_2();;}
    break;

  case 72:
#line 273 "benchmark_parser.y"
    { visitor->end_conic_point_2(); ;}
    break;

  case 73:
#line 278 "benchmark_parser.y"
    {;}
    break;

  case 75:
#line 280 "benchmark_parser.y"
    { visitor->accept_infty( yyvsp[-2]);
                                    visitor->accept_integer( yyvsp[0]); ;}
    break;

  case 76:
#line 286 "benchmark_parser.y"
    {;}
    break;

  case 77:
#line 287 "benchmark_parser.y"
    { visitor->begin_algebraic_real(); ;}
    break;

  case 78:
#line 288 "benchmark_parser.y"
    { visitor->begin_polynomial_1(); ;}
    break;

  case 79:
#line 290 "benchmark_parser.y"
    { visitor->end_polynomial_1(); ;}
    break;

  case 80:
#line 294 "benchmark_parser.y"
    { visitor->accept_integer(yyvsp[0]); ;}
    break;

  case 81:
#line 295 "benchmark_parser.y"
    { visitor->end_algebraic_real(); ;}
    break;

  case 82:
#line 299 "benchmark_parser.y"
    {;}
    break;

  case 83:
#line 301 "benchmark_parser.y"
    { visitor->accept_rational( yyvsp[-3], yyvsp[-1]); ;}
    break;

  case 84:
#line 305 "benchmark_parser.y"
    { visitor->accept_integer( yyvsp[0]); ;}
    break;

  case 85:
#line 306 "benchmark_parser.y"
    { visitor->accept_infty( yyvsp[0]); ;}
    break;

  case 91:
#line 326 "benchmark_parser.y"
    { visitor->parse_error( yyvsp[0]); ;}
    break;

  case 92:
#line 327 "benchmark_parser.y"
    { visitor->unknown_token( yyvsp[0]); ;}
    break;


    }

/* Line 1000 of yacc.c.  */
#line 1613 "benchmark_parser.tab.c"

  yyvsp -= yylen;
  yyssp -= yylen;


  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
        {
          YYSIZE_T yysize = 0;
          int yytype = YYTRANSLATE (yychar);
          const char* yyprefix;
          char *yymsg;
          int yyx;

          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;

          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yycount = 0;

          yyprefix = ", expecting ";
          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
              {
                yysize += yystrlen (yyprefix) + yystrlen (yytname [yyx]);
                yycount += 1;
                if (yycount == 5)
                  {
                    yysize = 0;
                    break;
                  }
              }
          yysize += (sizeof ("syntax error, unexpected ")
                     + yystrlen (yytname[yytype]));
          yymsg = (char *) YYSTACK_ALLOC (yysize);
          if (yymsg != 0)
            {
              char *yyp = yystpcpy (yymsg, "syntax error, unexpected ");
              yyp = yystpcpy (yyp, yytname[yytype]);

              if (yycount < 5)
                {
                  yyprefix = ", expecting ";
                  for (yyx = yyxbegin; yyx < yyxend; ++yyx)
                    if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
                      {
                        yyp = yystpcpy (yyp, yyprefix);
                        yyp = yystpcpy (yyp, yytname[yyx]);
                        yyprefix = " or ";
                      }
                }
              yyerror (yymsg);
              YYSTACK_FREE (yymsg);
            }
          else
            yyerror ("syntax error; also virtual memory exhausted");
        }
      else
#endif /* YYERROR_VERBOSE */
        yyerror ("syntax error");
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* If at end of input, pop the error token,
             then the rest of the stack, then return failure.  */
          if (yychar == YYEOF)
             for (;;)
               {
                 YYPOPSTACK;
                 if (yyssp == yyss)
                   YYABORT;
                 YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
                 yydestruct (yystos[*yyssp], yyvsp);
               }
        }
      else
        {
          YYDSYMPRINTF ("Error: discarding", yytoken, &yylval, &yylloc);
          yydestruct (yytoken, &yylval);
          yychar = YYEMPTY;

        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

#ifdef __GNUC__
  /* Pacify GCC when the user code never invokes YYERROR and the label
     yyerrorlab therefore never appears in user code.  */
  if (0)
     goto yyerrorlab;
#endif

  yyvsp -= yylen;
  yyssp -= yylen;
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;        /* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
        {
          yyn += YYTERROR;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;

      YYDSYMPRINTF ("Error: popping", yystos[*yyssp], yyvsp, yylsp);
      yydestruct (yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  YYDPRINTF ((stderr, "Shifting error token, "));

  *++yyvsp = yylval;


  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*----------------------------------------------.
| yyoverflowlab -- parser overflow comes here.  |
`----------------------------------------------*/
yyoverflowlab:
  yyerror ("parser stack overflow");
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 332 "benchmark_parser.y"


// Opens file 'name' and parses it. Uses visitor 'v' while parsing.
// Returns false if something went wrong. See the visitor for details.
bool benchmark_parse_file( std::string name, Benchmark_visitor* v) {
    std::ifstream in( name.c_str());
    if ( ! in) {
        v->parse_error( std::string( "cannot open file '") + name +
                        std::string( "'."));
        return false;
    }
    visitor = v;
    benchmark_init_lexer( in, name, 1);
    yyparse();
    return ! v->error();
}

// Starts parsing from stream 'in' with the associated filename 'name'
// (or analogous meaning for different streams) counting linenumbers
// starting from 'n'. Uses visitor 'v' while parsing. Returns false if
// something went wrong. See the visitor for details of the error reporting.
bool benchmark_parse_stream( std::istream&  in, std::string name,
                             Benchmark_visitor* v, int n) {
    visitor = v;
    benchmark_init_lexer( in, name, n);
    yyparse();
    return ! v->error();
}

/* EOF */

