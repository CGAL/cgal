@! ============================================================================
@! The CGAL Library
@! Implementation: Random Numbers Generator
@! ----------------------------------------------------------------------------
@! file  : web/Random/Random.aw
@! author: Sven Schönherr <sven@inf.fu-berlin.de>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Random Sources $
@! $CGAL_Package: Random WIP $
@! $Revision$
@! $Date$
@! ============================================================================

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{cc_manual}
@article

\setlength{\parskip}{1ex}
\addtolength{\textheight}{5ex}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsrevision}{$Revision$}
\RCSdefDate{\rcsdate}{$Date$}

@t vskip 5 mm
@t title titlefont centre "CGAL -- Random Numbers Generator*"
@t vskip 1 mm
@t title smalltitlefont centre "Implementation Documentation"
@t vskip 5 mm
@t title smalltitlefont centre "Sven Schönherr"
\smallskip
\centerline{\rcsrevision\ , \rcsdate}
@t vskip 1 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Project
  No.~21957 (CGAL).}

@! ============================================================================
@! Introduction and Contents
@! ============================================================================

\section*{Introduction}

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

\tableofcontents

@! ============================================================================
@! Specification
@! ============================================================================

\clearpage
\section{Specification}

\renewcommand{\ccSection}{\ccSubsection}
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{}
\ccSetThreeColumns{CGAL_Random}{random.restore_state( State state)}{}
\ccPropagateThreeToTwoColumns
\input{../../doc_tex/support/Random/Random.tex}

@! ============================================================================
@! Implementation
@! ============================================================================

\clearpage
\section{Implementation}

This section describes the implementation of the random numbers
generator. We use the C library function \ccc{erand48} to generate the
random numbers. It behaves like the well-known function \ccc{drand48}
but operates on a user supplied storage for the 48-Bit state. This
makes different instances of the random number generator independent.

First, we declare the class \ccc{CGAL_Random}.

@macro<Random declaration> = @begin
    class CGAL_Random;
@end

The class interface looks as follows.

@macro <Random interface> = @begin
    class CGAL_Random {
      public:
        @<Random public interface>

      private:
        @<Random private data members>
    };
@end   

\subsection{Public Interface}

The functionality is described and documented in the specification
section, so we do not comment on it here.

@macro <Random public interface> = @begin
    // types
    typedef  unsigned short  State[3];                  // 48 Bits

    // creation
    CGAL_Random( );
    CGAL_Random( long seed);
    CGAL_Random( State state);

    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);

    int     operator () ( int upper);

    // state functions
    void       save_state(       State& state) const;
    void    restore_state( const State& state);

    // equality test
    bool  operator == ( const CGAL_Random& rnd) const;
@end


%\pagebreak

\subsection{Global Variable}

The global variable \ccc{CGAL_random} is the default random numbers
generator.

@macro <Random global variable declaration> = @begin
    extern  CGAL_Random  CGAL_random;
@end

@macro <Random global variable definition> = @begin
    CGAL_Random  CGAL_random;
@end


\subsection{Private Data Members}

The state is stored in an array of three \ccc{unsigned short}s.

@macro <Random private data members> = @begin
    // data members
    unsigned short  _state[3];                          // 48 Bits
@end


\subsection{Constructors}

In the default constructor the seed is set using the system time.

@macro <Random constructors> = @begin
    #ifndef CGAL_PROTECT_CTIME
    #  include <ctime>
    #  define CGAL_PROTECT_CTIME
    #endif

    CGAL_Random::
    CGAL_Random( )
    {
        // get system's microseconds
        std::timeval tv;
        std::gettimeofday( &tv, NULL);
        unsigned long  ms = tv.tv_sec*1000000+tv.tv_usec;

        // initialize random numbers generator
        _state[ 0] = _state[ 2] = CGAL_static_cast( unsigned short, ms >> 16);
        _state[ 1] =            CGAL_static_cast( unsigned short, ms & 65535);
    }

    CGAL_Random::
    CGAL_Random( long seed)
    {
        // initialize random numbers generator
        _state[ 0] = _state[ 2] = CGAL_static_cast( unsigned short,seed >> 16);
        _state[ 1] =            CGAL_static_cast( unsigned short,seed & 65535);
    }

    CGAL_Random::
    CGAL_Random( State state)
    {
        // initialize random numbers generator
        _state[ 0] = state[ 0];
        _state[ 1] = state[ 1];
        _state[ 2] = state[ 2];
    }

@end


\subsection{Operations}

The C library function \ccc{erand48} returns a random \ccc{double},
uniformly chosen from the interval $[\ccc{0.0},\ccc{1.0})$.
The result is converted to a number in the given range.

@macro <Random operations> = @begin
    #ifndef CGAL_PROTECT_CSTDLIB
    #  include <cstdlib>
    #  define CGAL_PROTECT_CSTDLIB
    #endif

    using std::erand48;

    inline
    bool
    CGAL_Random::
    get_bool( )
    {
        return( CGAL_static_cast( bool, ( erand48( _state) >= 0.5)));
    }

    inline
    int
    CGAL_Random::
    get_int( int lower, int upper)
    {
        return( lower + CGAL_static_cast( int,
                  CGAL_static_cast( double, upper-lower) * erand48( _state)));
    }

    inline
    double
    CGAL_Random::
    get_double( double lower, double upper)
    {
        return( lower + ( upper-lower) * erand48( _state));
    }

    inline
    int
    CGAL_Random::
    operator () ( int upper)
    {
        return( get_int( 0, upper));
    }
@end


\subsection{State Functions}

The state functions just copy the internal state to or from the given
state variable, respectively.

@macro <Random state functions> = @begin
    void
    CGAL_Random::
    save_state( State& state) const
    {
        state[ 0] = _state[ 0];
        state[ 1] = _state[ 1];
        state[ 2] = _state[ 2];
    }

    void
    CGAL_Random::
    restore_state( const State& state)
    {
        _state[ 0] = state[ 0];
        _state[ 1] = state[ 1];
        _state[ 2] = state[ 2];
    }
@end


\subsection{Equality Test}

The equality test compares the internal states of the two operands.

@macro <Random equality test> = @begin
    inline
    bool    
    CGAL_Random::
    operator == ( const CGAL_Random& rnd) const
    {
        return( CGAL_static_cast( bool,
                    ( _state[ 0] == rnd._state[ 0]) &&
                    ( _state[ 1] == rnd._state[ 1]) &&
                    ( _state[ 2] == rnd._state[ 2]) ) );
    }
@end

@! ============================================================================
@! Test
@! ============================================================================

\clearpage
\section{Test}

We call each function of class \ccc{CGAL_Random} at least once to
ensure code coverage. In addition, we check if the generated random
numbers lie in the given ranges, and if two random numbers generators
initialized with the same seed generate the same sequence of random
numbers.

@macro <Random tests> = @begin
    CGAL_Random::State  state;
    CGAL_random.save_state( state);

    // test get_bool
    {
        bool b = CGAL_random.get_bool();
        assert( ! b || b);
    }

    // test get_int
    {
        int  l = CGAL_random.get_int( -100, 0);
        int  u = CGAL_random.get_int( 0, 1000);
        int  i = CGAL_random.get_int( l, u);
        assert( ( l <= i) && ( i < u));
    }

    // test get_double
    {
        double  l = CGAL_random.get_double( -123.45, -0.99);
        double  u = CGAL_random.get_double( 22.0/7.0, 33.3);
        double  d = CGAL_random.get_double( l, u);
        assert( ( l <= d) && ( d < u));
    }

    // test operator()
    {
        int  i = CGAL_random( 5555);
        assert( ( 0 <= i) && ( i < 5555));
    }

    // test state functions
    {
        CGAL_random.restore_state( state);      // `CGAL_Random' and `rnd'
        CGAL_Random rnd( state);                // have the same state now
        assert( CGAL_random.get_bool()         == rnd.get_bool()        );
        assert( CGAL_random.get_int( -100,100) == rnd.get_int( -100,100));
        assert( CGAL_random.get_double()       == rnd.get_double()      );
        assert( CGAL_random                    == rnd                   );

        long init = CGAL_random( 9999);
        CGAL_Random rnd1( init), rnd2( init);
        assert( rnd1.get_bool()         == rnd2.get_bool()        );
        assert( rnd1.get_int( -100,100) == rnd2.get_int( -100,100));
        assert( rnd1.get_double()       == rnd2.get_double()      );
        assert( rnd1                    == rnd2                   );
    }
@end

@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@! ----------------------------------------------------------------------------
@! Random.h
@! ----------------------------------------------------------------------------

@file <include/CGAL/Random.h> = @begin
    @<file header>("include/CGAL/Random.h","Random Numbers Generator")

    #ifndef CGAL_RANDOM_H
    #define CGAL_RANDOM_H

    // Class declaration
    // =================
    @<Random declaration>

    // Class interface
    // ===============
    // includes
    #ifndef CGAL_BASIC_H
    #  include <CGAL/basic.h>
    #endif

    @<Random interface>

    // Global variables
    // ================
    @<Random global variable declaration>

    @<dividing line>

    // Class implementation (inline functions)
    // =======================================
    // operations
    @<Random operations>

    @<Random equality test>

    #endif // CGAL_RANDOM_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Random.C
@! ----------------------------------------------------------------------------

@file <src/Random.C> = @begin
    @<file header>("src/Random.C","Random Numbers Generator")

    #include <CGAL/Random.h>

    // Class implementation (continued)
    // ================================
    // constructors
    @<Random constructors>

    // state functions
    @<Random state functions>

    // Global variables
    // ================
    @<Random global variable definition>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Random.C
@! ----------------------------------------------------------------------------

@file <test/Random/test_Random.C> = @begin
    @<file header>(
        "test/Random/test_Random.C",
        "test program for Random Numbers Generator")

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

@i ../file_header.awi

@macro <file header>(2) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Random Sources and Geometric Object Genera.",
        "Random","Random/Random",
        "$Revision$","$Date$",
        "Sven Schönherr <sven@@inf.fu-berlin.de>","N.N.",
        "INRIA Sophia-Antipolis (<Herve.Bronnimann@@sophia.inria.fr>)",
        "@2")
@end    

@! ===== EOF ==================================================================
