@! ============================================================================
@! The CGAL Project
@! Implementation: Random Numbers Generator
@! ----------------------------------------------------------------------------
@! file  : Kernel/web/Random.aw
@! author: Sven Schönherr (sven@inf.fu-berlin.de)
@! ----------------------------------------------------------------------------
@! $Revision$
@! $Date$
@! ============================================================================

@p maximum_input_line_length = 120
@p maximum_output_line_length = 120

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{cc_manual}
@article

\setlength{\parskip}{1ex}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsrevision}{$Revision$}
\RCSdefDate{\rcsdate}{$Date$}

@!thickline

@t vskip 5 mm
@t title titlefont centre "CGAL -- Random Numbers Generator*"
@t vskip 1 mm
@t title smalltitlefont centre "Implementation Documentation"
@t vskip 5 mm
@t title smalltitlefont centre "Sven Schönherr"
\smallskip
\centerline{\rcsrevision\ , \rcsdate}
@t vskip 1 mm
@!thickline

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
random number generating function object, e.g. for the STL algorithm
@prg{random_shuffle}.

This document is organized as follows. Section~1 contains the
specification as it appears in the CGAL Reference Manual. Section~2
gives the implementation. In Section~3 we provide a test program which
performs some correctness checks. Finally the product files are
created in Section~4.

\tableofcontents
\clearpage

@! ============================================================================
@! Specification
@! ============================================================================

\section{Specification}

\renewcommand{\ccSection}{\ccSubsection}
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{}
\input{../spec/Random.tex}
\clearpage

@! ============================================================================
@! Implementation
@! ============================================================================

\section{Implementation}

This section describes the implementation of the random numbers
generator. We use the C library function @prg{erand48} to generate the
random numbers. It behaves like the well-known function @prg{drand48} but
operates on a user supplied storage for the 48-Bit seed. This makes different
instances of the random number generator independent.

First, we declare the class @prg{CGAL_Random}.

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
    typedef typename  unsigned short  Seed[3];		// 48 Bits

    // creation
    CGAL_Random( );
    CGAL_Random( Seed seed);

    // operations
    bool    get_bool  ( );
    int     get_int   ( int lower, int upper);
    double  get_double( double lower = 0.0, double upper = 1.0);

    int     operator () ( int upper);

    void    get_seed( Seed& seed) const;
    void    set_seed( Seed const& seed);
@end

  
\subsection{Private Data Members}

The seed is stored in an array of three @prg{unsigned short}s.

@macro <Random private data members> = @begin
    // data members
    Seed  _seed;			
@end


\subsection{Constructors}

In the default constructor te seed is initialized using the system
time.

@macro <Random constructors> = @begin
    #include <stdlib.h>
    #include <sys/time.h>

    CGAL_Random::
    CGAL_Random( )
    {
	// get system's microseconds
	timeval tv;
	gettimeofday( &tv, NULL);
	unsigned long  ms = tv.tv_sec*1000000+tv.tv_usec;

	// initialize random numbers generator
	_seed[ 0] = _seed[ 2] = CGAL_static_cast( unsigned short, ms >> 16);
	_seed[ 1] =             CGAL_static_cast( unsigned short, ms & 65535);
    }

    CGAL_Random::
    CGAL_Random( Seed seed)
    {
	// initialize random numbers generator
	_seed[ 0] = seed[ 0];
	_seed[ 1] = seed[ 1];
	_seed[ 2] = seed[ 2];
    }
@end


\subsection{Operations}

The C library function @prg{erand48} returns a random @prg{double},
uniformly chosen from the interval $[@prg{0.0},@prg{1.0})$.
The result is converted to a number in the given range.

@macro <Random operations> = @begin
    inline
    bool
    CGAL_Random::
    get_bool( )
    {
	return( CGAL_static_cast( bool, ( erand48( _seed) >= 0.5)));
    }

    inline
    int
    CGAL_Random::
    get_int( int lower, int upper)
    {
	return( lower + CGAL_static_cast( int,
	          CGAL_static_cast( double, upper-lower) * erand48( _seed)));
    }

    inline
    double
    CGAL_Random::
    get_double( double lower, double upper)
    {
	return( lower + ( upper-lower) * erand48( _seed));
    }

    inline
    int
    CGAL_Random::
    operator () ( int upper)
    {
	return( get_int( 0, upper));
    }
@end

The seed operations just copy the internal seed to or from the given
seed variable, respectively.

@macro <Random seed operations> = @begin
    void
    CGAL_Random::
    get_seed( Seed& seed) const
    {
	seed[ 0] = _seed[ 0];
	seed[ 1] = _seed[ 1];
	seed[ 2] = _seed[ 2];
    }

    void
    CGAL_Random::
    set_seed( Seed const& seed)
    {
	_seed[ 0] = seed[ 0];
	_seed[ 1] = seed[ 1];
	_seed[ 2] = seed[ 2];
    }
@end

@! ============================================================================
@! Test
@! ============================================================================

\section{Test}

We call each function of class @prg{CGAL_Random} at least once to
ensure code coverage. In addition we check if the generated random
numbers lie in the given ranges, and if two random numbers generators
initialized with the same seed generate the same sequence of random
numbers.

@macro <Random tests> = @begin
    CGAL_Random        rnd;
    CGAL_Random::Seed  seed;
    rnd.get_seed( seed);

    // test get_bool
    {
	bool b = rnd.get_bool();
	CGAL_assertion( ! b || b);
    }

    // test get_int
    {
        int  l = rnd.get_int( -100, 0);
	int  u = rnd.get_int( 0, 1000);
	int  i = rnd.get_int( l, u);
	CGAL_assertion( ( l <= i) && ( i < u));
    }

    // test get_double
    {
        double  l = rnd.get_double( -123.45, -0.99);
	double  u = rnd.get_double( 22.0/7.0, 33.3);
	double  d = rnd.get_double( l, u);
	CGAL_assertion( ( l <= d) && ( d < u));
    }

    // test operator()
    {
	int  i = rnd( 5555);
	CGAL_assertion( ( 0 <= i) && ( i < 5555));
    }

    // test seed funtions
    {
	rnd.set_seed( seed);			// now `rnd' and `rnd2'
	CGAL_Random rnd2( seed);		// have the same seed

	CGAL_assertion( rnd.get_bool()          == rnd2.get_bool());
	CGAL_assertion( rnd.get_int( -100, 100) == rnd2.get_int( -100, 100));
	CGAL_assertion( rnd.get_double()        == rnd2.get_double());
    }
@end

@! ==========================================================================
@! Files
@! ==========================================================================

\section{Files}

@file <Random.h> = @begin
    @<Random header>("include/CGAL/Random.h")

    #ifndef CGAL_RANDOM_H
    #define CGAL_RANDOM_H

    #define typename

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

    @<dividing line>

    // Class implementation (inline functions)
    // =======================================
    // operations
    @<Random operations>

    #endif // CGAL_RANDOM_H

    @<end of file line>
@end

@file <Random.C> = @begin
    @<Random header>("src/Random.C")

    #include <CGAL/Random.h>

    // Class implementation (continued)
    // ================================
    // constructors
    @<Random constructors>

    // ssed operations
    @<Random seed operations>

    @<end of file line>
@end

@file <test_Random.C> = @begin
    @<Random header>("test/test_Random.C")

    #include <CGAL/Random.h>

    #define CGAL_assertion CGAL_kernel_assertion

    int
    main( int, char**)
    {
	@<Random tests>
    }

    @<end of file line>
@end

@i file_header.awlib

@macro <Random header>(1) many = @begin
    @<file header>("Random Numbers Generator",@1,"Random",
		   "Sven Schönherr (sven@@inf.fu-berlin.de)",
		   "$Revision$","$Date$")
@end

@! ===== EOF ==================================================================
