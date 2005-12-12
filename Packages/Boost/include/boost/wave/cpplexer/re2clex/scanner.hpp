/*=============================================================================
    Boost.Wave: A Standard compliant C++ preprocessor library

    http://www.boost.org/
    
    Copyright (c) 2001 Daniel C. Nuffer.
    Copyright (c) 2001-2005 Hartmut Kaiser. 
    Distributed under the Boost Software License, Version 1.0. (See accompanying 
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/

#if !defined(SCANNER_HPP_F4FB01EB_E75C_4537_A146_D34B9895EF37_INCLUDED)
#define SCANNER_HPP_F4FB01EB_E75C_4537_A146_D34B9895EF37_INCLUDED

#include <boost/wave/cpplexer/re2clex/aq.hpp>

///////////////////////////////////////////////////////////////////////////////
namespace boost {
namespace wave {
namespace cpplexer {
namespace re2clex {

struct Scanner;
typedef unsigned char uchar;
typedef int (* ReportErrorProc)(struct Scanner const *, char const *, ...);

typedef struct Scanner {
    int    fd;  /* file descriptor */
    uchar* first;   /* start of input buffer (if fd == -1) */
    uchar* act;     /* act position of input buffer (if fd == -1) */
    uchar* last;    /* end (one past last char) of input buffer (if fd == -1) */
    uchar* bot; /* beginning of the current buffer */
    uchar* top; /* top of the current buffer */
    uchar* eof; /* when we read in the last buffer, will point 1 past the 
                   end of the file, otherwise 0 */
    uchar* tok; /* points to the beginning of the current token */
    uchar* ptr; /* used for YYMARKER - saves backtracking info */
    uchar* cur; /* saves the cursor (maybe is redundant with tok?) */
    uchar* lim; /* used for YYLIMIT - points to the end of the buffer */
                /* (lim == top) except for the last buffer, it points to
                   the end of the input (lim == eof - 1) */
    unsigned int line;    /* current line being lexed */
    unsigned int column;        /* current token start column position */
    unsigned int curr_column;   /* current column position */
    ReportErrorProc error_proc;     /* if != 0 this function is called to 
                report an error */
    char const *file_name;  /* name of the lexed file */
    aq_queue eol_offsets;
    int enable_ms_extensions;   /* enable MS extensions */
    int act_in_c99_mode;        /* lexer works in C99 mode */
    int act_in_cpp0x_mode;      /* lexer works in C++0x mode */
} Scanner;

///////////////////////////////////////////////////////////////////////////////
}   // namespace re2clex
}   // namespace cpplexer
}   // namespace wave
}   // namespace boost

#endif // !defined(SCANNER_HPP_F4FB01EB_E75C_4537_A146_D34B9895EF37_INCLUDED)
