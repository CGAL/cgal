@! ============================================================================
@! The CGAL Library
@! Implementation: Random Numbers Generator
@! ----------------------------------------------------------------------------
@! file  : web/Random_numbers/Random.aw
@! author: Sven Schönherr <sven@inf.fu-berlin.de>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Random Numbers Generator $
@! $CGAL_Package: Random_numbers $
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
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Projects
  No.~21957 (CGAL) and No.~28155 (GALIA).}

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
\ccSetThreeColumns{Random}{random.restore_state( State state)}{}
\ccPropagateThreeToTwoColumns
\input{../../doc_tex/general/Random.tex}

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
    Random( );
    Random( long seed);
    Random( State state);

    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);

    int     operator () ( int upper);

    // state functions
    void       save_state(       State& state) const;
    void    restore_state( const State& state);

    // equality test
    bool  operator == ( const Random& rnd) const;
@end


%\pagebreak

\subsection{Global Variable}

The global variable \ccc{default_random} is the default random numbers
generator.

@macro <Random global variable declaration> = @begin
    extern  Random  default_random;
@end

@macro <Random global variable definition> = @begin
    Random  default_random;
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
    Random::
    Random( )
    {
        // get system's microseconds
        timeval tv;
        gettimeofday( &tv, NULL);
        unsigned long  ms = tv.tv_sec*1000000+tv.tv_usec;

        // initialize random numbers generator
        _state[ 0] = _state[ 2] = CGAL_static_cast( unsigned short, ms >> 16);
        _state[ 1] =            CGAL_static_cast( unsigned short, ms & 65535);
    }

    Random::
    Random( long seed)
    {
        // initialize random numbers generator
        _state[ 0] = _state[ 2] = CGAL_static_cast( unsigned short,seed >> 16);
        _state[ 1] =            CGAL_static_cast( unsigned short,seed & 65535);
    }

    Random::
    Random( State state)
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
    inline
    bool
    Random::
    get_bool( )
    {
        return( CGAL_static_cast( bool, ( erand48( _state) >= 0.5)));
    }

    inline
    int
    Random::
    get_int( int lower, int upper)
    {
        return( lower + CGAL_static_cast( int,
            CGAL_static_cast( double, upper-lower) * erand48( _state)));
    }

    inline
    double
    Random::
    get_double( double lower, double upper)
    {
        return( lower + ( upper-lower) * erand48( _state));
    }

    inline
    int
    Random::
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
    Random::
    save_state( State& state) const
    {
        state[ 0] = _state[ 0];
        state[ 1] = _state[ 1];
        state[ 2] = _state[ 2];
    }

    void
    Random::
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
    Random::
    operator == ( const Random& rnd) const
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

We call each function of class \ccc{Random} at least once to
ensure code coverage. In addition, we check if the generated random
numbers lie in the given ranges, and if two random numbers generators
initialized with the same seed generate the same sequence of random
numbers.

@macro <Random tests> = @begin
    Random::State  state;
    default_random.save_state( state);

    // test get_bool
    {
        bool b = default_random.get_bool();
        assert( ! b || b);
    }

    // test get_int
    {
        int  l = default_random.get_int( -100, 0);
        int  u = default_random.get_int( 0, 1000);
        int  i = default_random.get_int( l, u);
        assert( ( l <= i) && ( i < u));
    }

    // test get_double
    {
        double  l = default_random.get_double( -123.45, -0.99);
        double  u = default_random.get_double( 22.0/7.0, 33.3);
        double  d = default_random.get_double( l, u);
        assert( ( l <= d) && ( d < u));
    }

    // test operator()
    {
        int  i = default_random( 5555);
        assert( ( 0 <= i) && ( i < 5555));
    }

    // test state functions
    {
        default_random.restore_state( state);     // `default_random' and `rnd'
        Random rnd( state);                       // have the same state now
        assert( default_random.get_bool()         == rnd.get_bool()        );
        assert( default_random.get_int( -100,100) == rnd.get_int( -100,100));
        assert( default_random.get_double()       == rnd.get_double()      );
        assert( default_random                    == rnd                   );

        long init = default_random( 9999);
        Random rnd1( init), rnd2( init);
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

@i ../namespace.awi

@! ----------------------------------------------------------------------------
@! Random.h
@! ----------------------------------------------------------------------------

\subsection{Random.h}

@file <include/CGAL/Random.h> = @begin
    @<file header>("include/CGAL/Random.h","Random Numbers Generator")

    #ifndef CGAL_RANDOM_H
    #define CGAL_RANDOM_H

    // includes
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

    @<Random equality test>

    @<namespace end>("CGAL")
    
    #endif // CGAL_RANDOM_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Random.C
@! ----------------------------------------------------------------------------

\subsection{Random.C}

@file <src/Random.C> = @begin
    @<file header>("src/Random.C","Random Numbers Generator")

    #include <CGAL/Random.h>

    // additional includes
    #ifndef CGAL_PROTECT_CTIME
    #  include <ctime>
    #  define CGAL_PROTECT_CTIME
    #endif
    #ifndef CGAL_PROTECT_SYS_TIME_H
    #  include <sys/time.h>
    #  define CGAL_PROTECT_SYS_TIME_H
    #endif

    @<namespace begin>("CGAL")

    // Class implementation (continued)
    // ================================

    // constructors
    @<Random constructors>

    // state functions
    @<Random state functions>

    // Global variables
    // ================
    @<Random global variable definition>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Random.C
@! ----------------------------------------------------------------------------

\subsection{test\_Random.C}

@file <test/Random_numbers/test_Random.C> = @begin
    @<file header>(
        "test/Random_numbers/test_Random.C",
        "test program for Random Numbers Generator")

    // includes
    #include <CGAL/Random.h>
    #include <cassert>

    using namespace CGAL;
    
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
        "Random Numbers Generator",
        "Random_numbers","Random_numbers/Random",
        "$Revision$","$Date$",
        "Sven Schönherr <sven@@inf.fu-berlin.de>","N.N.",
        "INRIA Sophia-Antipolis (<Herve.Bronnimann@@sophia.inria.fr>)",
        "@2")
@end    

@! ===== EOF ==================================================================
