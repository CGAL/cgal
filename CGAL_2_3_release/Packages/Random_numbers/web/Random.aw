@! ============================================================================
@! The CGAL Library
@! Implementation: Random Numbers Generator
@! ----------------------------------------------------------------------------
@! file  : web/Random.aw
@! author: Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Random Numbers Generator $
@! $CGAL_Package: Random_numbers $
@! $Revision$
@! $Date$
@! ============================================================================

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{cc_manual,cc_manual_index}
@article

\setlength{\skip\footins}{3ex}

\pagestyle{headings}

@! LaTeX macros
\newcommand{\remark}[2]{[\textbf{#1:} \emph{#2}]}

\newcommand{\linebreakByHand}{\ccTexHtml{\linebreak[4]}{}}
\newcommand{  \newlineByHand}{\ccTexHtml{\\}{}}

\renewcommand{\sectionmark}[1]{\markboth{\uppercase{#1}}{}}

\newcommand{\subsectionRef}[2]{
  \addtocounter{subsection}{1}
  \addcontentsline{toc}{subsection}{\protect\numberline{\thesubsection}#1: #2}
  \markright{\thesubsection~~#1: #2}}

@! settings for `cc_manual.sty'
\ccDefGlobalScope{CGAL::}
\renewcommand{\ccRefPageEnd}{\clearpage}

\newcommand{\cgalColumnLayout}{%
  \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
  \ccPropagateThreeToTwoColumns}

\newcommand{\ccRequirements}{\ccHeading{Requirements}}
\newcommand{\ccRequire}{\ccCommentHeading{Requirements}}

@! ============================================================================
@! Title
@! ============================================================================

\thispagestyle{empty}

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 20 mm
@t title titlefont centre "Random Numbers Generator*"
@t vskip 15 mm
@t title smalltitlefont centre "Sven Schönherr"
\begin{center}
  \textbf{FU Berlin / ETH Z{\"u}rich}
\end{center}
@t vskip 10 mm
{\small
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Package: Random_numbers WIP+\cgalWIP\verb+$+ \\
    \verb+$CGAL_Chapter: Geometric Object Generators $+ \\
  \end{tabular}
\end{center}
}
@t vskip 30 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Projects
  No.~21957 (CGAL) and No.~28155 (GALIA).}

\renewcommand{\thefootnote}{\arabic{footnote}}

@! --------
@! Abstract
@! --------

\begin{abstract}
  We provide an implementation of a random numbers generator. It allows to
  generate uniformly distributed random @prg{bool}s, @prg{int}s, and
  @prg{double}s. The interface fulfills the requirements of an STL random
  number generating function object.
\end{abstract}

@! --------
@! Contents
@! --------

\clearpage

\newlength{\defaultparskip}
\setlength{\defaultparskip}{\parskip}
\setlength{\parskip}{1ex}

\tableofcontents

\setlength{\parskip}{\defaultparskip}


@! ============================================================================
@! Introduction
@! ============================================================================

\clearpage
\markright{\uppercase{Introduction}}
\section{Introduction}

We provide an implementation of a random numbers generator. It allows
to generate uniformly distributed random @prg{bool}s, @prg{int}s, and
@prg{double}s. The interface fulfills the requirements of an STL
random number generating function object, e.g.\ for the STL algorithm
@prg{random_shuffle}.

This document is organized as follows. Section~1 contains the
specification as it appears in the CGAL Reference Manual. Section~2
gives the implementation. In Section~3 we provide a test program which
performs some correctness checks. Finally the product files are
created in Section~4.


@! ============================================================================
@! Reference Pages
@! ============================================================================

\clearpage
\section{Reference Pages} \label{sec:reference_pages}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

@! ----------------------------------------------------------------------------
@! Class: Random
@! ----------------------------------------------------------------------------

\renewcommand{\ccSection}{\ccSubsection}
\input{../doc_tex/support/Generator/Random.tex}


@! ============================================================================
@! Implementation
@! ============================================================================

\clearpage
\section{Implementation}

This section describes the implementation of the random numbers
generator. We use the function \ccc{rand} from the standard C library
to generate the random numbers.

First, we declare the class \ccc{Random}.

@macro<Random declaration> = @begin
    class Random;
@end

The class interface looks as follows.

@macro <Random interface> = @begin
    class Random {
      public:
        @<Random public interface>
      private:
        @<Random data members>  
    };
@end   

@! ----------------------------------------------------------------------------
\subsection{Public Interface}

The functionality is described and documented in
Section~\ref{sec:random_numbers_generator}, so we do not comment on it
here.

@macro <Random public interface> = @begin
    // creation
    Random( );
    Random( unsigned int  seed);

    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);

    int     operator () ( int upper);
@end

@! ----------------------------------------------------------------------------
\subsection{Data Members}

We store the value of \ccc{RAND_MAX} plus one as a \ccc{double}
constant in \ccc{rand_max_plus_1}. This value is used in the
\ccc{get_...}  functions.

@macro <Random data members> = @begin
    // data members
    const double  rand_max_plus_1;
@end

@! ----------------------------------------------------------------------------
\subsection{Constructors}

In the default constructor the seed is set using the system time.

@macro <Random constructors> = @begin
    Random::
    Random( )
        : rand_max_plus_1( RAND_MAX+1.0)
    {
        // get system's time
        time_t s;
        time( &s);
        unsigned int  seed = s;

        // initialize random numbers generator
        srand( seed);
    }

    Random::
    Random( unsigned int  seed)
        : rand_max_plus_1( RAND_MAX+1.0)
    {
        // initialize random numbers generator
        srand( seed);
    }
@end

@! ----------------------------------------------------------------------------
\subsection{Operations}

The C library function \ccc{rand} returns a random \ccc{int},
uniformly chosen from the interval $[\ccc{0},\ccc{RAND_MAX}]$.
The result is converted to a number in the given range.

@macro <Random operations> = @begin
    inline
    bool
    Random::
    get_bool( )
    {
        return( static_cast< bool>( rand() & 1));
    }

    inline
    int
    Random::
    get_int( int lower, int upper)
    {
        return( lower + static_cast< int>(
          ( static_cast< double>( upper) - lower) * rand() / rand_max_plus_1));
    }

    inline
    double
    Random::
    get_double( double lower, double upper)
    {
        return( lower + ( ( upper-lower) * rand() / rand_max_plus_1));
    }

    inline
    int
    Random::
    operator () ( int upper)
    {
        return( get_int( 0, upper));
    }
@end

@! ----------------------------------------------------------------------------
\subsection{Global Variable}

The global variable \ccc{default_random} is the default random numbers
generator.

@macro <Random global variable declaration> = @begin
    extern  Random  default_random;
@end

@macro <Random global variable definition> = @begin
    Random  default_random;
@end


@! ============================================================================
@! Test Program
@! ============================================================================

\clearpage
\section{Test Program}

We call each function of class \ccc{Random} at least once to
ensure code coverage. In addition, we check if the generated random
numbers lie in the given ranges.

@macro <Random tests> = @begin
    // test get_bool
    {
        bool b = CGAL::default_random.get_bool();
        assert( ! b || b);
    }

    // test get_int
    {
        int  l = CGAL::default_random.get_int( -100, 0);
        int  u = CGAL::default_random.get_int( 0, 1000);
        int  i = CGAL::default_random.get_int( l, u);
        assert( ( l <= i) && ( i < u));
    }

    // test get_double
    {
        double  l = CGAL::default_random.get_double( -123.45, -0.99);
        double  u = CGAL::default_random.get_double( 22.0/7.0, 33.3);
        double  d = CGAL::default_random.get_double( l, u);
        assert( ( l <= d) && ( d < u));
    }

    // test operator()
    {
        int  i = CGAL::default_random( 5555);
        assert( ( 0 <= i) && ( i < 5555));
    }
@end


@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Random.h
@! ----------------------------------------------------------------------------

\subsection{include/CGAL/Random.h}

@file <include/CGAL/Random.h> = @begin
    @<file header>("include/CGAL/Random.h","Random Numbers Generator")

    #ifndef CGAL_RANDOM_H
    #define CGAL_RANDOM_H

    // includes
    // --------
    #ifndef CGAL_BASIC_H
    #  include <CGAL/basic.h>
    #endif

    @<namespace begin>("CGAL")

    // Class declaration
    // =================
    @<Random declaration>
    
    // Class interface
    // ===============
    @<Random interface>

    // Global variables
    // ================
    @<Random global variable declaration>

    @<namespace end>("CGAL")
    
    @<dividing line>

    // Class implementation (inline functions)
    // =======================================
    // includes
    #ifndef CGAL_PROTECT_CSTDLIB
    #  include <cstdlib>
    #  define CGAL_PROTECT_CSTDLIB
    #endif

    @<namespace begin>("CGAL")

    // operations
    @<Random operations>

    @<namespace end>("CGAL")
    
    #endif // CGAL_RANDOM_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Random.C
@! ----------------------------------------------------------------------------

\subsection{src/Random.C}

@file <src/Random.C> = @begin
    @<file header>("src/Random.C","Random Numbers Generator")

    #include <CGAL/Random.h>

    // additional includes
    #ifndef CGAL_PROTECT_CTIME
    #  include <ctime>
    #  define CGAL_PROTECT_CTIME
    #endif

    @<namespace begin>("CGAL")

    // Class implementation (continued)
    // ================================

    // constructors
    @<Random constructors>

    // Global variables
    // ================
    @<Random global variable definition>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Random.C
@! ----------------------------------------------------------------------------

\subsection{test/Random\_numbers/test\_Random.C}

@file <test/Random_numbers/test_Random.C> = @begin
    @<file header>(
        "test/Random_numbers/test_Random.C",
        "test program for Random Numbers Generator")

    // includes
    #include <CGAL/Random.h>
    #include <cassert>

    int
    main( int, char**)
    {
        @<Random tests>

        return( 0);
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! File Header
@! ----------------------------------------------------------------------------

\subsection*{File Header}

@i share/file_header.awi

@macro <file header>(2) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Random_numbers",
        "Random Numbers Generator",
        "Random",
        "$Revision$","$Date$",
        "Sven Schönherr <sven@@inf.ethz.ch>",
        "INRIA Sophia-Antipolis",
        "@2")
@end    

@! ===== EOF ==================================================================
