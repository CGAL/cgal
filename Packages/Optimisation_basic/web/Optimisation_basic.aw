@! ============================================================================
@! The CGAL Library
@! Implementation: basic stuff for optimisation algorithms
@! ----------------------------------------------------------------------------
@! file  : web/Optimisation/optimisation_basic.aw
@! author: Sven Schönherr <sven@inf.fu-berlin.de>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Optimisation_basic WIP $
@! $Revision$
@! $Date$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amssymb}
@usepackage{cc_manual}
@article

\setlength{\parskip}{1ex}
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 5 mm
@t title titlefont centre "Basic Stuff for Optimisation Algorithms*"
@t vskip 1 mm
@t title smalltitlefont centre "Sven Schönherr"
\smallskip
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Chapter: Geometric Optimisation $+ \\
    \verb+$CGAL_Package: Optimisation_basic WIP+\cgalWIP\verb+$+ \\
  \end{tabular}
\end{center}
@t vskip 1 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Projects
  No.~21957 (CGAL) and No.~28155 (GALIA).}

@! ============================================================================
@! Introduction and Contents
@! ============================================================================

\section*{Introduction}

We provide macros for checking assertions, pre- and postconditions and
warnings, for encapsulating code for debugging, and a failure function for the
checking functions of the optimisation algorithms.

\tableofcontents

@! ============================================================================
@! Macros
@! ============================================================================

\clearpage
\section{Macros}

@! --------
@! Checking
@! --------

\subsection{Checking}

The following macros are used to perform several checks. We distinguish four
groups of checks, namely checking for assertions, for preconditions, for
postconditions and for warnings. Each group consists of macros for normal
checks, exactness checks, expensive checks, and expensive exactness checks.

@macro<check macros>(2) many = @begin
    @<normal checks>(@1,@2)

    @<exactness checks>(@1,@2)

    @<expensive checks>(@1,@2)

    @<expensive exactness checks>(@1,@2)

@end

@! Normal Checks
\subsubsection*{Normal Checks}

@macro<normal checks>(2) many = @begin
    #if (    defined( CGAL_OPTIMISATION_NO_@1) \
          || defined( CGAL_NO_@1) || defined( NDEBUG))
    #  define  CGAL_optimisation_@2(EX)         ((void)0)
    #  define  CGAL_optimisation_@2_msg(EX,MSG) ((void)0)
    #  define  CGAL_optimisation_@2_code(CODE)
    #else
    #  define  CGAL_optimisation_@2(EX) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,0))
    #  define  CGAL_optimisation_@2_msg(EX,MSG) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,MSG))
    #  define  CGAL_optimisation_@2_code(CODE) CODE
    #endif // optimisation @2s
@end

@! Exactness Checks
\subsubsection*{Exactness Checks}

@macro<exactness checks>(2) many = @begin
    #if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
                 || defined( CGAL_CHECK_EXACTNESS)              ) \
          || defined( CGAL_OPTIMISATION_NO_@1) \
          || defined( CGAL_NO_@1) || defined( NDEBUG))
    #  define  CGAL_optimisation_exactness_@2(EX)         ((void)0)
    #  define  CGAL_optimisation_exactness_@2_msg(EX,MSG) ((void)0)
    #  define  CGAL_optimisation_exactness_@2_code(CODE)
    #else
    #  define  CGAL_optimisation_exactness_@2(EX) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,0))
    #  define  CGAL_optimisation_exactness_@2_msg(EX,MSG) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,MSG))
    #  define  CGAL_optimisation_exactness_@2_code(CODE) CODE
    #endif // optimisation exactness @2s
@end

@! Expensive Checks
\subsubsection*{Expensive Checks}

@macro<expensive checks>(2) many = @begin
    #if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
                 || defined( CGAL_CHECK_EXPENSIVE)              ) \
          || defined( CGAL_OPTIMISATION_NO_@1) \
          || defined( CGAL_NO_@1) || defined( NDEBUG))
    #  define  CGAL_optimisation_expensive_@2(EX)         ((void)0)
    #  define  CGAL_optimisation_expensive_@2_msg(EX,MSG) ((void)0)
    #  define  CGAL_optimisation_expensive_@2_code(CODE)
    #else
    #  define  CGAL_optimisation_expensive_@2(EX) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,0))
    #  define  CGAL_optimisation_expensive_@2_msg(EX,MSG) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,MSG))
    #  define  CGAL_optimisation_expensive_@2_code(CODE) CODE
    #endif // optimisation expensive @2s
@end

@! Expensive Exactness Checks
\subsubsection*{Expensive Exactness Checks}

@macro<expensive exactness checks>(2) many = @begin
    #if (    ! (    defined( CGAL_OPTIMISATION_CHECK_EXACTNESS) \
                 || defined( CGAL_OPTIMISATION_CHECK_EXPENSIVE) \
                 || defined( CGAL_CHECK_EXACTNESS)              \
                 || defined( CGAL_CHECK_EXPENSIVE)              ) \
          || defined( CGAL_OPTIMISATION_NO_@1) \
          || defined( CGAL_NO_@1) || defined( NDEBUG))
    #  define  CGAL_optimisation_expensive_exactness_@2(EX) \
                                                                      ((void)0)
    #  define  CGAL_optimisation_expensive_exactness_@2_msg(EX,MSG) \
                                                                      ((void)0)
    #  define  CGAL_optimisation_expensive_exactness_@2_code(CODE)
    #else
    #  define  CGAL_optimisation_expensive_exactness_@2(EX) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,0))
    #  define  CGAL_optimisation_expensive_exactness_@2_msg(EX,MSG) \
         ((EX)?((void)0): ::CGAL::@2_fail( # EX ,__FILE__,__LINE__,MSG))
    #  define  CGAL_optimisation_expensive_exactness_@2_code(CODE) CODE
    #endif // optimisation expensive exactness @2s
@end

@! ---------
@! Debugging
@! ---------

\subsection{Debugging}

The \ccc{CGAL_optimisation_debug} macro can be used to encapsulate code, that
should not be compiled into the production version, e.g.~diagnose output. The
following piece of code

@macro<optimisation debug macro usage> zero = @begin
    CGAL_optimisation_debug {
        // do something
    }
@end

is expanded to 

@macro<optimisation debug macro on> zero = @begin
    if ( 1) {
        // do something
    }
@end

in development versions, and to

@macro<optimisation debug macro off> zero = @begin
    if ( 0) {
        // do something
    }
@end

in production versions of the program. The latter will be discarded if
compiled with the \ccc{-O} option.

@macro<optimisation debug macro> = @begin
    #if (    defined( CGAL_OPTIMISATION_NO_DEBUG) \
          || defined( CGAL_NO_DEGUG) || defined( NDEBUG))
    #  define  CGAL_optimisation_debug  if ( 0)
    #else
    #  define  CGAL_optimisation_debug  if ( 1)
    #endif // optimisation debug
@end


@! ============================================================================
@! Functions
@! ============================================================================

\section{Functions}

@! ----------------------------------------------------------------------------
@! Function _optimisation_is_valid_fail
@! ----------------------------------------------------------------------------

\subsection{Function \ccFont \_optimisation\_is\_valid\_fail}

This function is called from the checking functions of the
optimisation algorithms if a check fails.  First, we declare the
function.

@macro<Function _optimisation_is_valid_fail declaration> = @begin
    bool
    _optimisation_is_valid_fail( CGAL::Verbose_ostream& verr,
                                 const char*            message);
@end

The function prints a failure \ccc{message} to the \cgal\ verbose
stream \ccc{verr}.

@macro<Function _optimisation_is_valid_fail> = @begin
    bool
    _optimisation_is_valid_fail( CGAL::Verbose_ostream& verr,
                                 const char*            message)
    {
        verr << "FAILED." << std::endl;
        verr << "  --> " << message << std::endl;
        verr << "  object is NOT valid!" << std::endl;
        return( false);
    }
@end


@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@i ../namespace.awi

@! ----------------------------------------------------------------------------
@! assertions.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation/assertions.h}

@file <include/CGAL/Optimisation/assertions.h> = @begin
    @<file header with coauthor>(
        "include/CGAL/Optimisation/assertions.h",
        "assertion macros for optimisation algorithms",
        "Geert-Jan Giezeman")

    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #define CGAL_OPTIMISATION_ASSERTIONS_H

    // macro definitions
    // =================

    // assertions
    // ----------
    @<check macros>("ASSERTIONS","assertion")

    // preconditions
    // -------------
    @<check macros>("PRECONDITIONS","precondition")

    // postconditions
    // --------------
    @<check macros>("POSTCONDITIONS","postcondition")

    // warnings
    // --------
    @<check macros>("WARNINGS","warning")

    #endif // CGAL_OPTIMISATION_ASSERTIONS_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! debug.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation/debug.h}

@file <include/CGAL/Optimisation/debug.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation/debug.h",
        "debug macro for optimisation algorithms")

    #ifndef CGAL_OPTIMISATION_DEBUG_H
    #define CGAL_OPTIMISATION_DEBUG_H

    // macro definitions
    // =================

    // debug
    // -----
    @<optimisation debug macro>

    #endif // CGAL_OPTIMISATION_DEBUG_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! basic.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Optimisation/basic.h}

@file <include/CGAL/Optimisation/basic.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation/basic.h",
        "basic things for optimisation algorithms")

    #ifndef CGAL_OPTIMISATION_BASIC_H
    #define CGAL_OPTIMISATION_BASIC_H

    // includes
    #ifndef CGAL_BASIC_H
    #  include <CGAL/basic.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif
    #ifndef CGAL_OPTIMISATION_DEBUG_H
    #  include <CGAL/Optimisation/debug.h>
    #endif
    #ifndef CGAL_IO_VERBOSE_OSTREAM_H
    #  include <CGAL/IO/Verbose_ostream.h>
    #endif

    @<namespace begin>("CGAL")

    // Function declarations
    // =====================

    // is_valid failure function
    // -------------------------
    @<Function _optimisation_is_valid_fail declaration>

    @<namespace end>("CGAL")

    #endif // CGAL_OPTIMISATION_BASIC_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! basic.C
@! ----------------------------------------------------------------------------

\subsection{src/Optimisation/basic.C}

@file <src/Optimisation/basic.C> = @begin
    @<file header>(
        "src/Optimisation/basic.C",
        "basic things for optimisation algorithms")

    #include <CGAL/Optimisation/basic.h>

    @<namespace begin>("CGAL")

    // Function implementations
    // ========================

    // is_valid failure function
    // -------------------------
    @<Function _optimisation_is_valid_fail>

    @<namespace end>("CGAL")

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! optimisation_Window_stream.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/IO/optimisation\_Window\_stream.h}

@file <include/CGAL/IO/optimisation_Window_stream.h> = @begin
    @<file header>(
        "include/CGAL/IO/optimisation_Window_stream.h",
        "graphical output to `leda_window' for optimisation algorith.")

    #include <CGAL/IO/Min_circle_2_Window_stream.h>
    #include <CGAL/IO/Min_ellipse_2_Window_stream.h>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! File Header
@! ----------------------------------------------------------------------------

\subsection*{File Header}

@i ../file_header.awi
 
And here come two specific file headers for the product files of this
web file.

@macro <file header>(2) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Geometric Optimisation",
        "Optimisation_basic","Optimisation/optimisation_basic",
        "$Revision$","$Date$",
        "Sven Schönherr",
        "Sven Schönherr <sven@@inf.fu-berlin.de>",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@macro <file header with coauthor>(3) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Geometric Optimisation",
        "Optimisation_basic","Optimisation/optimisation_basic",
        "$Revision$","$Date$",
        "@3, Sven Schönherr",
        "Sven Schönherr <sven@@inf.fu-berlin.de>",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ===== EOF ==================================================================
