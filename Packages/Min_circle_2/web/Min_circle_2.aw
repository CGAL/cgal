@! ============================================================================
@! The CGAL Library
@! Implementation: 2D Smallest Enclosing Circle
@! ----------------------------------------------------------------------------
@! file  : web/Min_circle_2.aw
@! author: Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Min_circle_2 WIP $
@! $Revision$
@! $Date$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amssymb}
@usepackage{path}
@usepackage{cc_manual,cc_manual_index}
@article

\input{cprog.sty}
\setlength{\parskip}{1ex}

@! LaTeX macros
\newenvironment{pseudocode}[1]%
  {\vspace*{-0.5\baselineskip} \upshape \begin{tabbing}
     99 \= bla \= bla \= bla \= bla \= bla \= bla \= bla \= bla \kill
     #1 \+ \\}%
  {\end{tabbing}}
\newcommand{\keyword}[1]{\texttt{#1}}
\newcommand{\IF}{\keyword{IF} }
\newcommand{\THEN}{\keyword{THEN} \+ \\}
\newcommand{\ELSE}{\< \keyword{ELSE} \\}
\newcommand{\END}{\< \keyword{END} \- \\ }
\newcommand{\OR}{\keyword{OR} }
\newcommand{\FOR}{\keyword{FOR} }
\newcommand{\TO}{\keyword{TO} }
\newcommand{\DO}{\keyword{DO} \+ \\}
\newcommand{\RETURN}{\keyword{RETURN} }

\newcommand{\mc}{\texttt{mc}}

\newcommand{\linebreakByHand}{\ccTexHtml{\linebreak[4]}{}}
\newcommand{  \newlineByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 5 mm
@t title titlefont centre "2D Smallest Enclosing Circle*"
@t vskip 1 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
\smallskip
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Chapter: Geometric Optimisation $+ \\
    \verb+$CGAL_Package: Min_circle_2 WIP+\cgalWIP\verb+$+ \\
  \end{tabular}
\end{center}
@t vskip 1 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Project
  No.~21957 (CGAL).}

@! ============================================================================
@! Introduction and Contents
@! ============================================================================

\section*{Introduction}

We provide an implementation of an optimisation algorithm for computing
the smallest (w.r.t.\ area) enclosing circle of a finite point set $P$
in the plane. The class template \ccc{Min_circle_2} is implemented
as a semi-dynamic data structure, thus allowing to insert points while
maintaining the smallest enclosing circle. It is parameterized with a
traits class, that defines the abstract interface between the
optimisation algorithm and the primitives it uses. For ease of use, we
provide traits class adapters that interface the optimisation algorithm
with user supplied point classes.

This document is organized as follows. The algorithm is described in
Section~1. Section~2 contains the specifications as they appear in the
CGAL Reference Manual. Section~3 gives the implementations. In
Section~4 we provide a test program which performs some correctness
checks. Finally the product files are created in Section~5.

\tableofcontents

@! ============================================================================
@! The Algorithm
@! ============================================================================

\clearpage
\section{The Algorithm} \label{sec:algo}

The implementation is based on an algorithm by Welzl~\cite{w-sedbe-91a},
which we shortly describe now. The smallest (w.r.t.\ area) enclosing
circle of a finite point set $P$ in the plane, denoted by $mc(P)$, is
built up incrementally, adding one point after another. Assume $mc(P)$
has been constructed, and we would like to obtain $mc(P \cup \{p\})$, $p$
some new point. There are two cases: if $p$ already lies inside $mc(P)$,
then $mc(P \cup \{p\}) = mc(P)$. Otherwise $p$ must lie on the boundary
of $mc(P \cup \{p\})$ (this is proved in~\cite{w-sedbe-91a} and not hard
to see), so we need to compute $mc(P,\{p\})$, the smallest circle
enclosing $P$ with $p$ on the boundary. This is recursively done in the
same manner. In general, for point sets $P$,$B$, define $mc(P,B)$ as the
smallest circle enclosing $P$ that has the points of $B$ on the boundary
(if defined). Although the algorithm finally delivers a circle
$mc(P,\emptyset)$, it internally deals with circles that have a possibly
nonempty set $B$. Here is the pseudo-code of Welzl's method.  To compute
$mc(P)$, it is called with the pair $(P,\emptyset)$, assuming that
$P=\{p_1,\ldots,p_n\}$ is stored in a linked list.

\begin{pseudocode}{$\mc(P,B)$:}
  $mc := \mc(\emptyset,B)$ \\
  \IF $|B| = 3$ \keyword{THEN} \RETURN $mc$ \\
  \FOR $i := 1$ \TO $n$ \DO
      \IF $p_i \not\in mc$ \THEN
          $mc := \mc(\{p_1,\ldots,p_{i-1}\}, B \cup \{p_i\})$ \\
          move $p_i$ to the front of $P$ \\
      \END
  \END
  \RETURN $mc$ \\
\end{pseudocode}

Note the following: (a) $|B|$ is always bounded by 3, thus the
computation of $mc(\emptyset,B)$ is easy. In our implementation, it is
done by the private member function \ccc{compute_circle}. (b) One can
check that the method maintains the invariant `$mc(P,B)$ exists'. This
justifies termination if $|B| = 3$, because then $mc(P,B)$ must be the
unique circle with the points of $B$ on the boundary, and $mc(P,B)$
exists if and only if this circle contains the points of $P$. Thus, no
subsequent in-circle tests are necessary anymore (for details
see~\cite{w-sedbe-91a}). (c) points which are found to lie outside the
current circle $mc$ are considered `important' and are moved to the front
of the linked list that stores $P$. This is crucial for the method's
efficiency.

It can also be advisable to bring $P$ into random order before
computation starts. There are `bad' insertion orders which cause the
method to be very slow -- random shuffling gives these orders a very
small probability.

@! ============================================================================
@! Specifications
@! ============================================================================

@! \clearpage
@! \section{Specifications}
@! 
@! \emph{Note:} Below some references are undefined, they refer to sections
@! in the \cgal\ Reference Manual.
@! 
@! \renewcommand{\ccFont}{\tt}
@! \renewcommand{\ccEndFont}{}
@! \newcommand{\cgalColumnLayout}{
@!   \ccSetThreeColumns{Oriented_side}{}{\hspace*{10cm}}
@!   \ccPropagateThreeToTwoColumns}
@! \newcommand{\cgalSetMinCircleLayout}{%
@!   \ccSetThreeColumns{Support_point_iterator}{}{returns
@!     \ccc{ON_BOUNDED_SIDE}, \ccc{ON_BOUNDARY}, or \ccc{ON_UNBOUNDED_SIDE}}
@!     \ccPropagateThreeToTwoColumns}
@! \input{../../doc_tex/basic/Optimisation/Min_circle_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Optimisation_circle_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_circle_2_traits_2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_circle_2_adapterC2.tex}
@! \input{../../doc_tex/basic/Optimisation/Min_circle_2_adapterH2.tex}

@! ============================================================================
@! Implementations
@! ============================================================================

\clearpage
\section{Implementations}

@! ----------------------------------------------------------------------------
@! Class template Min_circle_2<Traits>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_circle\_2<Traits>}

First, we declare the class template \ccc{Min_circle_2}.

@macro<Min_circle_2 declaration> = @begin
    template < class Traits_ >
    class Min_circle_2;
@end

The actual work of the algorithm is done in the private member
functions \ccc{mc} and \ccc{compute_circle}. The former directly
realizes the pseudo-code of $\mc(P,B)$, the latter solves the basic
case $\mc(\emptyset,B)$, see Section~\ref{sec:algo}.

The class interface looks as follows.

@macro <Min_circle_2 interface> = @begin
    template < class Traits_ >
    class Min_circle_2 {
      public:
        @<Min_circle_2 public interface>

      private:
        // private data members
        @<Min_circle_2 private data members>

        // copying and assignment not allowed!
        Min_circle_2( const Min_circle_2<Traits_>&);
        Min_circle_2<Traits_>&
            operator = ( const Min_circle_2<Traits_>&);

    @<dividing line>

    // Class implementation
    // ====================

      public:
        // Access functions and predicates
        // -------------------------------
        @<Min_circle_2 access functions `number_of_...'>

        @<Min_circle_2 predicates `is_...'>

        @<Min_circle_2 access functions>

        @<Min_circle_2 predicates>

      private:
        // Private member functions
        // ------------------------
        @<Min_circle_2 private member function `compute_circle'>

        @<Min_circle_2 private member function `mc'>

      public:
        // Constructors
        // ------------
        @<Min_circle_2 constructors>

        // Destructor
        // ----------
        @<Min_circle_2 destructor>

        // Modifiers
        // ---------
        @<Min_circle_2 modifiers>

        // Validity check
        // --------------
        @<Min_circle_2 validity check>

        // Miscellaneous
        // -------------
        @<Min_circle_2 miscellaneous>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Public Interface}

The functionality is described and documented in the specification
section, so we do not comment on it here.

@macro <Min_circle_2 public interface> = @begin
    // types
    typedef           Traits_                           Traits;
    typedef typename  Traits_::Point                    Point;
    typedef typename  Traits_::Circle                   Circle;
    typedef typename  std::list<Point>::const_iterator  Point_iterator;
    typedef           const Point *                     Support_point_iterator;

    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.

    // creation
    template < class InputIterator >
    Min_circle_2( InputIterator first,
                  InputIterator last,
                  bool          randomize = false,
                  Random&       random    = default_random,
                  const Traits& traits    = Traits());
    
    Min_circle_2( const Traits& traits = Traits());
    Min_circle_2( const Point&  p,
                  const Traits& traits = Traits());
    Min_circle_2( const Point&  p, const Point&  q,
                  const Traits& traits = Traits());
    Min_circle_2( const Point&  p1, const Point&  p2, const Point&  p3,
                  const Traits& traits = Traits());
    ~Min_circle_2( );

    // access functions
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;

    Point_iterator  points_begin( ) const;
    Point_iterator  points_end  ( ) const;

    Support_point_iterator  support_points_begin( ) const;
    Support_point_iterator  support_points_end  ( ) const;

    const Point&  support_point( int i) const;

    const Circle&  circle( ) const;

    // predicates
    Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;

    // modifiers
    void  insert( const Point& p);
    void  insert( const Point* first,
                  const Point* last );
    void  insert( std::list<Point>::const_iterator first,
                  std::list<Point>::const_iterator last );
    void  insert( std::istream_iterator<Point,std::ptrdiff_t> first,
                  std::istream_iterator<Point,std::ptrdiff_t> last );
    void  clear( );

    // validity check
    bool  is_valid( bool verbose = false, int level = 0) const;

    // miscellaneous
    const Traits&  traits( ) const;
    **************************************************************************/
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Data Members}

First, the traits class object is stored.

@macro <Min_circle_2 private data members> += @begin
    Traits       tco;                           // traits class object
@end

The points of $P$ are internally stored as a linked list that allows to
bring points to the front of the list in constant time. We use the
sequence container \ccc{list} from STL~\cite{sl-stl-95}.

@macro <Min_circle_2 private data members> += @begin
    std::list<Point>  points;                   // doubly linked list of points
@end

The support set $S$ of at most three support points is stored in an
array \ccc{support_points}, the actual number of support points is
given by \ccc{n_support_points}. During the computations, the set of
support points coincides with the set $B$ appearing in the pseudo-code
for $\mc(P,B)$, see Section~\ref{sec:algo}.

\emph{Workaround:} The array of support points is allocated dynamically,
because the SGI compiler (mipspro CC 7.1) does not accept a static
array here.

@macro <Min_circle_2 private data members> += @begin
    int          n_support_points;              // number of support points
    Point*       support_points;                // array of support points
@end

Finally, the actual circle is stored in a variable \ccc{circle}
provided by the traits class object, by the end of computation equal to
$mc(P)$. During computation, \ccc{tco.circle} equals the circle $mc$
appearing in the pseudo-code for $\mc(P,B)$, see Section~\ref{sec:algo}.

@! ----------------------------------------------------------------------------
\subsubsection{Constructors and Destructor}

We provide several different constructors, which can be put into two
groups. The constructors in the first group, i.e. the more important
ones, build the smallest enclosing circle $mc(P)$ from a point set
$P$, given by a begin iterator and a past-the-end iterator, realized
as a member template.

All constructors of the first group copy the points into the internal
list \ccc{points}. If randomization is demanded, the points are copied to
a vector and shuffled at random, before being copied to \ccc{points}.
Finally the private member function $mc$ is called to compute
$mc(P)=mc(P,\emptyset)$.

@macro <Min_circle_2 constructors> += @begin
    // STL-like constructor (member template)
    template < class InputIterator >
    Min_circle_2( InputIterator first,
                  InputIterator last,
                  bool          randomize
#if !defined(__BORLANDC__) && (!defined(_MSC_VER) || _MSC_VER > 1200)
                                          = false
#endif
                                                 ,
                  Random&       random    = default_random,
                  const Traits& traits    = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 3];

        // range of points not empty?
        if ( first != last) {

            // store points
            if ( randomize) {

                // shuffle points at random
                std::vector<Point> v( first, last);
                std::random_shuffle( v.begin(), v.end(), random);
                std::copy( v.begin(), v.end(),
                           std::back_inserter( points)); }
            else
                std::copy( first, last, std::back_inserter( points)); }

        // compute mc
        mc( points.end(), 0);
    }
@end

The remaining constructors are actually specializations of the
previous ones, building the smallest enclosing circle for up to three
points. The idea is the following: recall that for any point set $P$
there exists $S \subseteq P$, $|S| \leq 3$ with $mc(S) = mc(P)$ (in
fact, such a set $S$ is determined by the algorithm). Once $S$ has
been computed (or given otherwise), $mc(P)$ can easily be
reconstructed from $S$ in constant time. To make this reconstruction
more convenient, a constructor is available for each size of $|S|$,
ranging from 0 to 3. For $|S|=0$, we get the default constructor,
building $mc(\emptyset)$.

@macro <Min_circle_2 constructors> += @begin

    // default constructor
    inline
    Min_circle_2( const Traits& traits = Traits())
        : tco( traits), n_support_points( 0)
    {
        // allocate support points' array
        support_points = new Point[ 3];

        // initialize circle
        tco.circle.set();

        CGAL_optimisation_postcondition( is_empty());
    }

    // constructor for one point
    inline
    Min_circle_2( const Point& p, const Traits& traits = Traits())
        : tco( traits), points( 1, p), n_support_points( 1)
    {
        // allocate support points' array
        support_points = new Point[ 3];

        // initialize circle
        support_points[ 0] = p;
        tco.circle.set( p);

        CGAL_optimisation_postcondition( is_degenerate());
    }

    // constructor for two points
    inline
    Min_circle_2( const Point& p1, const Point& p2,
                  const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 3];

        // store points
        points.push_back( p1);
        points.push_back( p2);

        // compute mc
        mc( points.end(), 0);
    }

    // constructor for three points
    inline
    Min_circle_2( const Point& p1, const Point& p2, const Point& p3,
                  const Traits& traits = Traits())
        : tco( traits)
    {
        // allocate support points' array
        support_points = new Point[ 3];

        // store points
        points.push_back( p1);
        points.push_back( p2);
        points.push_back( p3);

        // compute mc
        mc( points.end(), 0);
    }
@end

The destructor only frees the memory of the support points' array.

@macro <Min_circle_2 destructor> = @begin
    inline
    ~Min_circle_2( )
    {
        // free support points' array
        delete[] support_points;
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Access Functions}

These functions are used to retrieve information about the current
status of the \ccc{Min_circle_2<Traits>} object. They are all very
simple (and therefore \ccc{inline}) and mostly rely on corresponding
access functions of the data members of \ccc{Min_circle_2<Traits>}.

First, we define the \ccc{number_of_...} methods.

@macro <Min_circle_2 access functions `number_of_...'> = @begin
    // #points and #support points
    inline
    int
    number_of_points( ) const
    {
        return( points.size());
    }

    inline
    int
    number_of_support_points( ) const
    {
        return( n_support_points);
    }
@end

Then, we have the access functions for points and support points.

@macro <Min_circle_2 access functions> += @begin
    // access to points and support points
    inline
    Point_iterator
    points_begin( ) const
    {
        return( points.begin());
    }

    inline
    Point_iterator
    points_end( ) const
    {
        return( points.end());
    }

    inline
    Support_point_iterator
    support_points_begin( ) const
    {
        return( support_points);
    }

    inline
    Support_point_iterator
    support_points_end( ) const
    {
        return( support_points+n_support_points);
    }

    // random access for support points    
    inline
    const Point&
    support_point( int i) const
    {
        CGAL_optimisation_precondition( (i >= 0) &&
                                        (i <  number_of_support_points()));
        return( support_points[ i]);
    }
@end

Finally, the access function \ccc{circle}.

@macro <Min_circle_2 access functions> += @begin
    // circle
    inline
    const Circle&
    circle( ) const
    {
        return( tco.circle);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The predicates \ccc{is_empty} and \ccc{is_degenerate} are used in
preconditions and postconditions of some member functions. Therefore we
define them \ccc{inline} and put them in a separate macro.

@macro <Min_circle_2 predicates `is_...'> = @begin
    // is_... predicates
    inline
    bool
    is_empty( ) const
    {
        return( number_of_support_points() == 0);
    }

    inline
    bool
    is_degenerate( ) const
    {
        return( number_of_support_points() <  2);
    }
@end

The remaining predicates perform in-circle tests, based on the
corresponding predicates of class \ccc{Circle}.

@macro <Min_circle_2 predicates> = @begin
    // in-circle test predicates
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( tco.circle.bounded_side( p));
    }

    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( tco.circle.has_on_bounded_side( p));
    }

    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( tco.circle.has_on_boundary( p));
    }

    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( tco.circle.has_on_unbounded_side( p));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Modifiers}

There is another way to build up $mc(P)$, other than by supplying
the point set $P$ at once. Namely, $mc(P)$ can be built up
incrementally, adding one point after another. If you look at the
pseudo-code in the introduction, this comes quite naturally. The
modifying method \ccc{insert}, applied with point $p$ to a
\ccc{Min_circle_2<Traits>} object representing $mc(P)$,
computes $mc(P \cup \{p\})$, where work has to be done only if $p$
lies outside $mc(P)$. In this case, $mc(P \cup \{p\}) = mc(P,\{p\})$
holds, so the private member function \ccc{mc} is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudo-code in
Section~\ref{sec:algo}.

@macro <Min_circle_2 modifiers> += @begin
    void
    insert( const Point& p)
    {
        // p not in current circle?
        if ( has_on_unbounded_side( p)) {

            // p new support point
            support_points[ 0] = p;

            // recompute mc
            mc( points.end(), 1);

            // store p as the first point in list
            points.push_front( p); }
        else

            // append p to the end of the list
            points.push_back( p);
    }
@end

Inserting a range of points is done by a single member template. In case
a compiler does not support this yet, we provide specialized \ccc{insert}
functions for C~arrays (using pointers as iterators), for STL sequence
containers \ccc{vector<Point>} and \ccc{list<Point>} and for the STL
input stream iterator \ccc{istream_iterator<Point>}. Actually, the
\ccc{insert} function for a C~array and a \ccc{vector<point>} are the
same, since the random access iterator of \ccc{vector<Point>} is
implemented as \ccc{Point*}.

The following \ccc{insert} functions perform a call \ccc{insert(p)} for
each point \ccc{p} in the range $[\mbox{\ccc{first}},\mbox{\ccc{last}})$.

@macro <Min_circle_2 modifiers> += @begin
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES

    template < class InputIterator >
    void
    insert( InputIterator first, InputIterator last)
    {
        for ( ; first != last; ++first)
            insert( *first);
    }

#else

    inline
    void
    insert( const Point* first, const Point* last)
    {
        for ( ; first != last; ++first)
            insert( *first);
    }

    inline
    void
    insert( std::list<Point>::const_iterator first,
            std::list<Point>::const_iterator last )
    {
        for ( ; first != last; ++first)
            insert( *first);
    }

    inline
    void
    insert( std::istream_iterator<Point,std::ptrdiff_t>  first,
            std::istream_iterator<Point,std::ptrdiff_t>  last )
    {
        for ( ; first != last; ++first)
            insert( *first);
    }

#endif // CGAL_CFG_NO_MEMBER_TEMPLATES
@end

The member function \ccc{clear} deletes all points from a
\ccc{Min_circle_2<Traits>} object and resets it to the
empty circle.

@macro <Min_circle_2 modifiers> += @begin
    void
    clear( )
    {
        points.erase( points.begin(), points.end());
        n_support_points = 0;
        tco.circle.set();
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Validity Check}

A \ccc{Min_circle_2<Traits>} object can be checked for validity.  This
means, it is checked whether (a) the circle contains all points of its
defining set $P$, (b) the circle is the smallest circle spanned by its
support set, and (c) the support set is minimal, i.e.\ no support point is
redundant. If \ccc{verbose} is \ccc{true}, some messages concerning the
performed checks are written to standard error stream. The second parameter
\ccc{level} is not used, we provide it only for consistency with interfaces
of other classes.

@macro <Min_circle_2 validity check> = @begin
    bool
    is_valid( bool verbose = false, int level = 0) const
    {
        CGAL_USING_NAMESPACE_STD
        
        CGAL::Verbose_ostream verr( verbose);
        verr << endl;
        verr << "CGAL::Min_circle_2<Traits>::" << endl;
        verr << "is_valid( true, " << level << "):" << endl;
        verr << "  |P| = " << number_of_points()
             << ", |S| = " << number_of_support_points() << endl;

        // containment check (a)
        @<Min_circle_2 containment check>

        // support set checks (b)+(c)
        @<Min_circle_2 support set checks>

        verr << "  object is valid!" << endl;
        return( true);
    }
@end

The containment check (a) is easy to perform, just a loop over all
points in \ccc{points}.

@macro <Min_circle_2 containment check> = @begin
    verr << "  a) containment check..." << flush;
    Point_iterator point_iter;
    for ( point_iter  = points_begin();
          point_iter != points_end();
          ++point_iter)
        if ( has_on_unbounded_side( *point_iter)) 
            return( CGAL::_optimisation_is_valid_fail( verr,
                        "circle does not contain all points"));
    verr << "passed." << endl;
@end

To check the support set (b) and (c), we distinguish four cases with
respect to the number of support points, which may range from 0 to 3.

@macro <Min_circle_2 support set checks> = @begin
    verr << "  b)+c) support set checks..." << flush;
    switch( number_of_support_points()) {

      case 0:
        @<Min_circle_2 check no support point>
        break;

      case 1:
        @<Min_circle_2 check one support point>
        break;

      case 2: {
        @<Min_circle_2 check two support points> }
        break;

      case 3: {
        @<Min_circle_2 check three support points> }
        break;

      default:
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "illegal number of support points, \
                     not between 0 and 3."));
    };
    verr << "passed." << endl;
@end

The case of no support point happens if and only if the defining
point set $P$ is empty.

@macro <Min_circle_2 check no support point> = @begin
    if ( ! is_empty())
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "P is nonempty, \
                     but there are no support points."));
@end

If the smallest enclosing circle has one support point $p$, it must
be equal to that point, i.e.\ its center must be $p$ and its radius
$0$.

@macro <Min_circle_2 check one support point> = @begin
    if ( ( circle().center() != support_point( 0)    ) ||
         ( ! CGAL_NTS is_zero( circle().squared_radius())) )
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "circle differs from the circle \
                     spanned by its single support point."));
@end

In case of two support points $p,q$, these points must form a diameter
of the circle. The support set $\{p,q\}$ is minimal if and only if
$p,q$ are distinct.

The diameter property is checked as follows. If $p$ and $q$ both lie
on the circle's boundary and if $p$, $q$ (knowing they are distinct)
and the circle's center are collinear, then $p$ and $q$ form a
diameter of the circle.

@macro <Min_circle_2 check two support points> = @begin
    const Point& p = support_point( 0),
                 q = support_point( 1);

    // p equals q?
    if ( p == q)
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "the two support points are equal."));

    // segment(p,q) is not diameter?
    if ( ( ! has_on_boundary( p)                                ) ||
         ( ! has_on_boundary( q)                                ) ||
         ( tco.orientation( p, q,
                            circle().center()) != CGAL::COLLINEAR) )
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "circle does not have its \
                     two support points as diameter."));
@end

If the number of support points is three (and they are distinct and
not collinear), the circle through them is unique, and must therefore
equal the current circle stored in \ccc{circle}. It is the smallest
one containing the three points if and only if the center of the
circle lies inside or on the boundary of the triangle defined by the
three points. The support set is minimal only if the center lies
properly inside the triangle.

Both triangle properties are checked by comparing the orientations of
three point triples, each containing two of the support points and the
center of the current circle, resp. If one of these orientations equals
\ccc{CGAL::COLLINEAR}, the center lies on the boundary of the triangle.
Otherwise, if two triples have opposite orientations, the center is not
contained in the triangle.

@macro <Min_circle_2 check three support points> = @begin
    const Point& p = support_point( 0),
                 q = support_point( 1),
                 r = support_point( 2);

    // p, q, r not pairwise distinct?
    if ( ( p == q) || ( q == r) || ( r == p))
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "at least two of the three \
                     support points are equal."));

    // p, q, r collinear?
    if ( tco.orientation( p, q, r) == CGAL::COLLINEAR)
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "the three support points are collinear."));

    // current circle not equal the unique circle through p,q,r ?
    Circle c = circle();
    c.set( p, q, r);
    if ( circle() != c)
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "circle is not the unique circle \
                     through its three support points."));

    // circle's center on boundary of triangle(p,q,r)?
    const Point& center = circle().center();
    CGAL::Orientation o_pqz = tco.orientation( p, q, center);
    CGAL::Orientation o_qrz = tco.orientation( q, r, center);
    CGAL::Orientation o_rpz = tco.orientation( r, p, center);
    if ( ( o_pqz == CGAL::COLLINEAR) ||
         ( o_qrz == CGAL::COLLINEAR) ||
         ( o_rpz == CGAL::COLLINEAR) )
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "one of the three support points is redundant."));

    // circle's center not inside triangle(p,q,r)?
    if ( ( o_pqz != o_qrz) || ( o_qrz != o_rpz) || ( o_rpz != o_pqz))
        return( CGAL::_optimisation_is_valid_fail( verr,
                    "circle's center is not in the \
                     convex hull of its three support points."));
@end

@! ----------------------------------------------------------------------------
\subsubsection{Miscellaneous}

The member function \ccc{traits} returns a const reference to the
traits class object.

@macro <Min_circle_2 miscellaneous> = @begin
    inline
    const Traits&
    traits( ) const
    {
        return( tco);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Min_circle_2 I/O operators declaration> = @begin
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os, const Min_circle_2<Traits_>& mc);

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is,       Min_circle_2<Traits_>& mc);
@end

@macro <Min_circle_2 I/O operators> = @begin
    template < class Traits_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const Min_circle_2<Traits_>& min_circle)
    {
        CGAL_USING_NAMESPACE_STD

        typedef  Min_circle_2<Traits_>::Point  Point;
        typedef  ostream_iterator<Point>       Os_it;

        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << endl;
            os << "CGAL::Min_circle_2( |P| = " << min_circle.number_of_points()
               << ", |S| = " << min_circle.number_of_support_points() << endl;
            os << "  P = {" << endl;
            os << "    ";
            copy( min_circle.points_begin(), min_circle.points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  S = {" << endl;
            os << "    ";
            copy( min_circle.support_points_begin(),
                  min_circle.support_points_end(),
                  Os_it( os, ",\n    "));
            os << "}" << endl;
            os << "  circle = " << min_circle.circle() << endl;
            os << ")" << endl;
            break;

          case CGAL::IO::ASCII:
            copy( min_circle.points_begin(), min_circle.points_end(),
                  Os_it( os, "\n"));
            break;

          case CGAL::IO::BINARY:
            copy( min_circle.points_begin(), min_circle.points_end(),
                  Os_it( os));
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class Traits_ >
    std::istream&
    operator >> ( std::istream& is, CGAL::Min_circle_2<Traits_>& min_circle)
    {
        CGAL_USING_NAMESPACE_STD
       
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << endl;
            cerr << "Stream must be in ascii or binary mode" << endl;
            break;

          case CGAL::IO::ASCII:
          case CGAL::IO::BINARY:
            typedef  CGAL::Min_circle_2<Traits_>::Point  Point;
            typedef  istream_iterator<Point>            Is_it;
            min_circle.clear();
            min_circle.insert( Is_it( is), Is_it());
            break;

          default:
            CGAL_optimisation_assertion_msg( false, "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Graphical Output}

@macro<Min_circle_2 graphical output operator> = @begin
    #ifdef CGAL_MIN_CIRCLE_2_H
    #ifndef CGAL_IO_WINDOW_STREAM_MIN_CIRCLE_2
    #define CGAL_IO_WINDOW_STREAM_MIN_CIRCLE_2

    template< class Traits_ >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream &ws,
                  const CGAL::Min_circle_2<Traits_>& min_circle)
    {
        typedef  CGAL::Min_circle_2<Traits_>::Point_iterator  Point_iterator;

        Point_iterator  first( min_circle.points_begin());
        Point_iterator  last ( min_circle.points_end());
        for ( ; first != last; ++first)
            ws << *first;
        return( ws << min_circle.circle());
    }

    #endif // CGAL_IO_WINDOW_STREAM_MIN_CIRCLE_2
    #endif // CGAL_MIN_CIRCLE_2_H
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Member Function {\ccFont compute\_circle}}

This is the method for computing the basic case $\mc(\emptyset,B)$,
the set $B$ given by the first \ccc{n_support_points} in the array
\ccc{support_points}. It is realized by a simple case analysis,
noting that $|B| \leq 3$.

@macro <Min_circle_2 private member function `compute_circle'> = @begin
    // compute_circle
    inline
    void
    compute_circle( )
    {
        switch ( n_support_points) {
          case 3:
            tco.circle.set( support_points[ 0],
                            support_points[ 1],
                            support_points[ 2]);
            break;
          case 2:
            tco.circle.set( support_points[ 0], support_points[ 1]);
            break;
          case 1:
            tco.circle.set( support_points[ 0]);
            break;
          case 0:
            tco.circle.set( );
            break;
          default:
            CGAL_optimisation_assertion( ( n_support_points >= 0) &&
                                         ( n_support_points <= 3) ); }
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Member Function {\ccFont mc}}

This function computes the general circle $mc(P,B)$, where $P$ contains
the points in the range $[$\ccc{points.begin()}$,$\ccc{last}$)$ and $B$
is given by the first \ccc{n_sp} support points in the array
\ccc{support_points}. The function is directly modeled after the
pseudo-code above.

@macro <Min_circle_2 private member function `mc'> = @begin
    void
    mc( const Point_iterator& last, int n_sp)
    {
        // compute circle through support points
        n_support_points = n_sp;
        compute_circle();
        if ( n_sp == 3) return;

        // test first n points
        typename std::list<Point>::iterator  point_iter = points.begin();
        for ( ; last != point_iter; ) {
            const Point& p = *point_iter;

            // p not in current circle?
            if ( has_on_unbounded_side( p)) {

                // recursive call with p as additional support point
                support_points[ n_sp] = p;
                mc( point_iter, n_sp+1);

                // move current point to front
                points.splice( points.begin(), points, point_iter++); }

            else
                ++point_iter; }
    }
@end

@! ----------------------------------------------------------------------------
@! Class template Optimisation_circle_2<K>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Optimisation\_circle\_2<K>}

First, we declare the class template \ccc{Optimisation_circle_2},

@macro<Optimisation_circle_2 declaration> = @begin
    template < class K_ >
    class Optimisation_circle_2;
@end

The class interface looks as follows.

@macro <Optimisation_circle_2 interface> = @begin
    template < class K_ >
    class Optimisation_circle_2 {
      public:
        @<Optimisation_circle_2 public interface>

      private:
        // private data members
        @<Optimisation_circle_2 private data members>

    @<dividing line>

    // Class implementation
    // ====================

      public:
        // Constructor
        // -----------
        @<Optimisation_circle_2 constructor>
        
        // Set functions
        // -------------
        @<Optimisation_circle_2 set functions>

        // Access functions
        // ----------------
        @<Optimisation_circle_2 access functions>

        // Equality tests
        // --------------
        @<Optimisation_circle_2 equality tests>

        // Predicates
        // ----------
        @<Optimisation_circle_2 predicates>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Public Interface}

The functionality is described and documented in the specification
section, so we do not comment on it here.

@macro <Optimisation_circle_2 public interface> = @begin
    // types
    typedef           K_                K;
    typedef           CGAL::Point_2<K>  Point;
    typedef typename  K_::FT            Distance;

    /**************************************************************************
    WORKAROUND: Some compilers are unable to match member functions defined
    outside the class template. Therefore, all member functions are implemented
    in the class interface.

    // creation
    Optimisation_circle_2( );

    void  set( );
    void  set( const Point& p);
    void  set( const Point& p, const Point& q);
    void  set( const Point& p, const Point& q, const Point& r);
    void  set( const Point& center, const Distance& squared_radius);

    // access functions    
    const Point&     center        ( ) const;
    const Distance&  squared_radius( ) const

    // equality tests
    bool  operator == ( const Optimisation_circle_2<K>& c) const;
    bool  operator != ( const Optimisation_circle_2<K>& c) const;

    // predicates
    CGAL::Bounded_side  bounded_side( const Point& p) const;
    bool  has_on_bounded_side      ( const Point& p) const;
    bool  has_on_boundary          ( const Point& p) const;
    bool  has_on_unbounded_side    ( const Point& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
    **************************************************************************/
@end

@! ----------------------------------------------------------------------------
\subsubsection{Private Data Members}

The circle is represented by its center and squared radius.

@macro <Optimisation_circle_2 private data members> = @begin
    Point     _center;
    Distance  _squared_radius;
@end

@! ----------------------------------------------------------------------------
\subsubsection{Constructor}

Only a default constructor is needed.

@macro <Optimisation_circle_2 constructor> = @begin
    inline
    Optimisation_circle_2( )
    { }
@end
    
@! ----------------------------------------------------------------------------
\subsubsection{Set Functions}

We provide set functions taking zero, one, two, or three boundary
points and another set function taking a center point and a squared
radius.

@macro <Optimisation_circle_2 set functions> = @begin
    inline
    void
    set( )
    {
        _center         =  Point( CGAL::ORIGIN);
        _squared_radius = -Distance( 1);
    }
    
    inline
    void
    set( const Point& p)
    {
        _center         = p;
        _squared_radius = Distance( 0);
    }
    
    inline
    void
    set( const Point& p, const Point& q)
    {
        _center         = CGAL::midpoint( p, q);
        _squared_radius = CGAL::squared_distance( p, _center);
    }
    
    inline
    void
    set( const Point& p, const Point& q, const Point& r)
    {
        _center         = CGAL::circumcenter( p, q, r);
        _squared_radius = CGAL::squared_distance( p, _center);
    }
    
    inline
    void
    set( const Point& center, const Distance& squared_radius)
    {
        _center         = center;
        _squared_radius = squared_radius;
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Access Functions}

These functions are used to get the current center point or
squared radius, resp.

@macro <Optimisation_circle_2 access functions> = @begin
    inline
    const Point&
    center( ) const
    {
        return( _center);
    }

    inline
    const Distance&
    squared_radius( ) const
    {
        return( _squared_radius);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Equality Tests}

@macro <Optimisation_circle_2 equality tests> = @begin
    bool
    operator == ( const Optimisation_circle_2<K>& c) const
    {
        return( ( _center          == c._center        ) &&
                ( _squared_radius  == c._squared_radius) );
    }
    
    bool
    operator != ( const Optimisation_circle_2<K>& c) const
    {
        return( ! operator==( c));
    }
@end
    
@! ----------------------------------------------------------------------------
\subsubsection{Predicates}

The following predicates perform in-circle tests and check for
emptiness and degeneracy, resp.

@macro <Optimisation_circle_2 predicates> = @begin
    inline
    CGAL::Bounded_side
    bounded_side( const Point& p) const
    {
        return( CGAL::Bounded_side( CGAL_NTS sign(
            _squared_radius - CGAL::squared_distance( p, _center))));
    }

    inline
    bool
    has_on_bounded_side( const Point& p) const
    {
        return( CGAL::squared_distance( p, _center) < _squared_radius);
    }

    inline
    bool
    has_on_boundary( const Point& p) const
    {
        return( CGAL::squared_distance( p, _center) == _squared_radius);
    }

    inline
    bool
    has_on_unbounded_side( const Point& p) const
    {
        return( _squared_radius < CGAL::squared_distance( p, _center));
    }

    inline
    bool
    is_empty( ) const
    {
        return( CGAL::is_negative( _squared_radius));
    }

    inline
    bool
    is_degenerate( ) const
    {
        return( ! CGAL::is_positive( _squared_radius));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{I/O}

@macro <Optimisation_circle_2 I/O operators declaration> = @begin
    template < class K_ >
    std::ostream&
    operator << ( std::ostream&, const CGAL::Optimisation_circle_2<K_>&);

    template < class K_ >
    std::istream&
    operator >> ( std::istream&, CGAL::Optimisation_circle_2<K_>&);
@end

@macro <Optimisation_circle_2 I/O operators> = @begin
    template < class K_ >
    std::ostream&
    operator << ( std::ostream& os, const CGAL::Optimisation_circle_2<K_>& c)
    {
        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Optimisation_circle_2( "
               << c.center() << ", "
               << c.squared_radius() << ')';
            break;

          case CGAL::IO::ASCII:
            os << c.center() << ' ' << c.squared_radius();
            break;

          case CGAL::IO::BINARY:
            os << c.center();
            CGAL::write( os, c.squared_radius());
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class K_ >
    std::istream&
    operator >> ( std::istream& is, CGAL::Optimisation_circle_2<K_>& c)
    {
        typedef  CGAL::Optimisation_circle_2<K_>::Point     Point;
        typedef  CGAL::Optimisation_circle_2<K_>::Distance  Distance;

        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII: {
            Point     center;
            Distance  squared_radius;
            is >> center >> squared_radius;
            c.set( center, squared_radius); }
            break;

          case CGAL::IO::BINARY: {
            Point     center;
            Distance  squared_radius;
            is >> center;
            CGAL::read( is, squared_radius);
            c.set( center, squared_radius); }
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( is) invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Graphical Output}

@macro<Optimisation_circle_2 graphical output operator> = @begin
    #ifdef CGAL_OPTIMISATION_CIRCLE_2_H
    #ifndef CGAL_IO_WINDOW_STREAM_OPTIMISATION_CIRCLE_2
    #define CGAL_IO_WINDOW_STREAM_OPTIMISATION_CIRCLE_2

    template< class K_ >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream &ws,
                  const CGAL::Optimisation_circle_2<K_>& oc)
    {
        double  cx( CGAL::to_double( oc.center().x()));
        double  cy( CGAL::to_double( oc.center().y()));
        double  sr( CGAL::to_double( oc.squared_radius()));

        if ( ! CGAL_NTS is_negative( sr))
            ws.draw_circle( cx, cy, CGAL::sqrt( sr));
        return( ws);
    }

    #endif // CGAL_IO_WINDOW_STREAM_OPTIMISATION_CIRCLE_2
    #endif // CGAL_OPTIMISATION_CIRCLE_2_H
@end


@! ----------------------------------------------------------------------------
@! Class template Min_circle_2_traits_2<K>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_circle\_2\_traits\_2<K>}

First, we declare the class templates \ccc{Min_circle_2} and
\ccc{Min_circle_2_traits_2}.

@macro<Min_circle_2_traits_2 declarations> = @begin
    template < class Traits_ >
    class Min_circle_2;

    template < class K_ >
    class Min_circle_2_traits_2;
@end

Since the actual work of the traits class is done in the nested type
\ccc{Circle}, we implement the whole class template in its interface.

The variable \ccc{circle} containing the current circle is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current circle,
it is declared \ccc{friend}.

@macro <Min_circle_2_traits_2 interface and implementation> = @begin
    template < class K_ >
    class Min_circle_2_traits_2 {
      public:
        // types
        typedef  K_                              K;
        typedef  CGAL::Point_2<K>                Point;
        typedef  CGAL::Optimisation_circle_2<K>  Circle;

    private:
        // data members
        Circle  circle;                                 // current circle

        // friends
        friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_traits_2<K> >;

      public:
        // creation (use default implementations)
        // CGAL::Min_circle_2_traits_2( );
        // CGAL::Min_circle_2_traits_2( CGAL::Min_circle_2_traits_2<K> const&);

        // operations
        inline
        CGAL::Orientation
        orientation( const Point& p, const Point& q, const Point& r) const
        {
            return( CGAL::orientation( p, q, r));
        }
    };
@end   

@! ----------------------------------------------------------------------------
@! Class template Min_circle_2_adapterC2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_circle\_2\_adapterC2<PT,DA>}

First, we declare the class templates \ccc{Min_circle_2},
\ccc{Min_circle_2_adapterC2} and
\ccc{_Min_circle_2_adapterC2__Circle}.

@macro<Min_circle_2_adapterC2 declarations> = @begin
    template < class Traits_ >
    class Min_circle_2;

    template < class PT_, class DA_ >
    class Min_circle_2_adapterC2;

    template < class PT_, class DA_ >
    class _Min_circle_2_adapterC2__Circle;
@end

The actual work of the adapter is done in the nested class
\ccc{Circle}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{circle} containing the current circle is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current circle,
it is declared \ccc{friend}.

@macro <Min_circle_2_adapterC2 interface and implementation> = @begin
    template < class PT_, class DA_ >
    class Min_circle_2_adapterC2 {
      public:
        // types
        typedef  PT_  PT;
        typedef  DA_  DA;

        // nested types
        typedef  PT                                            Point;
        typedef  CGAL::_Min_circle_2_adapterC2__Circle<PT,DA>  Circle;

      private:
        DA      dao;                                    // data accessor object
        Circle  circle;                                 // current circle
        friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_adapterC2<PT,DA> >;

      public:
        // creation
        @<Min_circle_2_adapterC2 constructors>

        // operations
        @<Min_circle_2_adapterC2 operations>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Constructors}

@macro <Min_circle_2_adapterC2 constructors> = @begin
    Min_circle_2_adapterC2( const DA& da = DA())
        : dao( da), circle( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_circle_2_adapterC2 operations> = @begin
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::FT  FT;

        FT  px;
        FT  py;
        FT  qx;
        FT  qy;
        FT  rx;
        FT  ry;
        
        dao.get( p, px, py);
        dao.get( q, qx, qy);
        dao.get( r, rx, ry);

        return( static_cast< CGAL::Orientation>(
                  CGAL_NTS sign( ( px-rx) * ( qy-ry) - ( py-ry) * ( qx-rx))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Circle}

@macro <Min_circle_2_adapterC2 nested type `Circle'> = @begin
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream&,
                  const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream&,
                  CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

    template < class PT_, class DA_ >
    class _Min_circle_2_adapterC2__Circle {
      public:
        // typedefs
        typedef  PT_  PT;
        typedef  DA_  DA;

        typedef  typename DA_::FT  FT;

      private:
        // data members
        DA  dao;                                // data accessor object

        FT  center_x;                           // center's x-coordinate
        FT  center_y;                           // center's y-coordinate
        FT  sqr_rad;                            // squared radius

        // private member functions
        FT
        sqr_dist( const FT& px, const FT& py, const FT& qx, const FT& qy) const
        {
            FT  dx( px - qx);
            FT  dy( py - qy);
            return( dx*dx + dy*dy);
        }

        friend  std::ostream&  operator << CGAL_NULL_TMPL_ARGS ( std::ostream&,
            const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

        friend  std::istream&  operator >> CGAL_NULL_TMPL_ARGS ( std::istream&,
            CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>&);

      public:
        // types
        typedef  PT  Point;
        typedef  FT  Distance;

        // creation
        _Min_circle_2_adapterC2__Circle( const DA& da) : dao( da) { }

        void  set( )
        {
            center_x =  FT( 0);
            center_y =  FT( 0);
            sqr_rad  = -FT( 1);
        }

        void  set( const Point& p)
        {
            dao.get( p, center_x, center_y);
            sqr_rad = FT( 0);
        }

        void  set( const Point& p, const Point& q)
        {
            FT  px;
            FT  py;
            FT  qx;
            FT  qy;

            dao.get( p, px, py);
            dao.get( q, qx, qy);

            center_x = ( px+qx) / FT( 2);
            center_y = ( py+qy) / FT( 2);
            sqr_rad  = sqr_dist( px, py, center_x, center_y);
        }

        void  set( const Point& p, const Point& q, const Point& r)
        {
            FT  px;
            FT  py;
            FT  qx;
            FT  qy;
            FT  rx;
            FT  ry;

            dao.get( p, px, py);
            dao.get( q, qx, qy);
            dao.get( r, rx, ry);

            FT  qx_px( qx - px);
            FT  qy_py( qy - py);
            FT  rx_px( rx - px);
            FT  ry_py( ry - py);
 
            FT  p2   ( px*px + py*py);
            FT  q2_p2( qx*qx + qy*qy - p2); 
            FT  r2_p2( rx*rx + ry*ry - p2); 
            FT  denom( ( qx_px*ry_py - rx_px*qy_py) * FT( 2));

            center_x = ( q2_p2*ry_py - r2_p2*qy_py) / denom;
            center_y = ( r2_p2*qx_px - q2_p2*rx_px) / denom;
            sqr_rad  = sqr_dist( px, py, center_x, center_y);
        }

        // predicates
        CGAL::Bounded_side
        bounded_side( const Point& p) const
        {
            FT  px;
            FT  py;
            dao.get( p, px, py);
            return( CGAL::Bounded_side( 
             CGAL_NTS sign( sqr_rad - sqr_dist( px, py, center_x, center_y))));
        }

        bool
        has_on_bounded_side( const Point& p) const
        {
            FT  px;
            FT  py;
            dao.get( p, px, py);
            return( sqr_dist( px, py, center_x, center_y) < sqr_rad);
        }

        bool
        has_on_boundary( const Point& p) const
        {
            FT  px;
            FT  py;
            dao.get( p, px, py);
            return( sqr_dist( px, py, center_x, center_y) == sqr_rad);
        }

        bool
        has_on_unbounded_side( const Point& p) const
        {
            FT  px;
            FT  py;
            dao.get( p, px, py);
            return( sqr_rad < sqr_dist( px, py, center_x, center_y));
        }

        bool
        is_empty( ) const
        {
            return( CGAL::is_negative( sqr_rad));
        }

        bool
        is_degenerate( ) const
        {
            return( ! CGAL::is_positive( sqr_rad));
        }

        // additional operations for checking
        bool
        operator == (
            const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c) const
        {
            return( ( center_x == c.center_x) &&
                    ( center_y == c.center_y) &&
                    ( sqr_rad  == c.sqr_rad ) );
        }

        bool
        operator != (
            const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c) const
        {
            return( ! ( *this == c));
        }

        Point
        center( ) const
        {
            Point  p;
            dao.set( p, center_x, center_y);
            return( p);
        }

        const Distance&
        squared_radius( ) const
        {
            return( sqr_rad);
        }
    };

    // I/O
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c)
    {
        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Min_circle_2_adapterC2::Circle( "
               << c.center_x << ", "
               << c.center_y << ", "
               << c.sqr_rad  << ')';
            break;

          case CGAL::IO::ASCII:
            os << c.center_x << ' ' << c.center_y << ' ' << c.sqr_rad;
            break;

          case CGAL::IO::BINARY:
            CGAL::write( os, c.center_x);
            CGAL::write( os, c.center_y);
            CGAL::write( os, c.sqr_rad);
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream& is,
                  CGAL::_Min_circle_2_adapterC2__Circle<PT_,DA_>& c)
    {
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII:
            is >> c.center_x >> c.center_y >> c.sqr_rad;
            break;

          case CGAL::IO::BINARY:
            CGAL::read( is, c.center_x);
            CGAL::read( is, c.center_y);
            CGAL::read( is, c.sqr_rad);
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ----------------------------------------------------------------------------
@! Class template Min_circle_2_adapterH2<PT,DA>
@! ----------------------------------------------------------------------------

\subsection{Class template \ccFont Min\_circle\_2\_adapterH2<PT,DA>}

First, we declare the class templates \ccc{Min_circle_2},
\ccc{Min_circle_2_adapterH2} and
\ccc{_Min_circle_2_adapterH2__Circle}.

@macro<Min_circle_2_adapterH2 declarations> = @begin
    template < class Traits_ >
    class Min_circle_2;

    template < class PT_, class DA_ >
    class Min_circle_2_adapterH2;

    template < class PT_, class DA_ >
    class _Min_circle_2_adapterH2__Circle;
@end

The actual work of the adapter is done in the nested class
\ccc{Circle}. Therefore, we implement the whole adapter in its
interface.

The variable \ccc{circle} containing the current circle is declared
\ccc{private} to disallow the user from directly accessing or modifying
it. Since the algorithm needs to access and modify the current circle,
it is declared \ccc{friend}.

@macro <Min_circle_2_adapterH2 interface and implementation> = @begin
    template < class PT_, class DA_ >
    class Min_circle_2_adapterH2 {
      public:
        // types
        typedef  PT_  PT;
        typedef  DA_  DA;

        // nested types
        typedef  PT                                            Point;
        typedef  CGAL::_Min_circle_2_adapterH2__Circle<PT,DA>  Circle;

      private:
        DA      dao;                                    // data accessor object
        Circle  circle;                                 // current circle
        friend  class CGAL::Min_circle_2< CGAL::Min_circle_2_adapterH2<PT,DA> >;

      public:
        // creation
        @<Min_circle_2_adapterH2 constructors>

        // operations
        @<Min_circle_2_adapterH2 operations>
    };
@end   

@! ----------------------------------------------------------------------------
\subsubsection{Constructors}

@macro <Min_circle_2_adapterH2 constructors> = @begin
    Min_circle_2_adapterH2( const DA& da = DA())
        : dao( da), circle( da)
    { }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Operations}

@macro <Min_circle_2_adapterH2 operations> = @begin
    CGAL::Orientation
    orientation( const Point& p, const Point& q, const Point& r) const
    {
        typedef  typename DA_::RT  RT;

        RT  phx;
        RT  phy;
        RT  phw;
        RT  qhx;
        RT  qhy;
        RT  qhw;
        RT  rhx;
        RT  rhy;
        RT  rhw;
        
        dao.get( p, phx, phy, phw);
        dao.get( q, qhx, qhy, qhw);
        dao.get( r, rhx, rhy, rhw);

        return( static_cast< CGAL::Orientation>(
                 CGAL_NTS sign( ( phx*rhw - rhx*phw) * ( qhy*rhw - rhy*qhw)
                              - ( phy*rhw - rhy*phw) * ( qhx*rhw - rhx*qhw))));
    }
@end

@! ----------------------------------------------------------------------------
\subsubsection{Nested Type \ccFont Circle}

@macro <Min_circle_2_adapterH2 nested type `Circle'> = @begin
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream&,
                  const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream&,
                  CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

    template < class PT_, class DA_ >
    class _Min_circle_2_adapterH2__Circle {
      public:
        // typedefs
        typedef  PT_  PT;
        typedef  DA_  DA;

        typedef  typename DA_::RT             RT;
        typedef           CGAL::Quotient<RT>  FT;

      private:
        // data members
        DA  dao;                                // data accessor object

        RT  center_hx;                          // center's hx-coordinate
        RT  center_hy;                          // center's hy-coordinate
        RT  center_hw;                          // center's hw-coordinate
        FT  sqr_rad;                            // squared radius

        // private member functions
        FT
        sqr_dist( const RT& phx, const RT& phy, const RT& phw,
                  const RT& qhx, const RT& qhy, const RT& qhw) const
        {
            RT  dhx( phx*qhw - qhx*phw);
            RT  dhy( phy*qhw - qhy*phw);
            RT  dhw( phw*qhw);
            return( FT( dhx*dhx + dhy*dhy, dhw*dhw));
        }

        friend  std::ostream&  operator << CGAL_NULL_TMPL_ARGS ( std::ostream&,
            const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

        friend  std::istream&  operator >> CGAL_NULL_TMPL_ARGS ( std::istream&,
            CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>&);

      public:
        // types
        typedef  PT  Point;
        typedef  FT  Distance;

        // creation
        _Min_circle_2_adapterH2__Circle( const DA& da) : dao( da) { }

        void  set( )
        {
            center_hx =  RT( 0);
            center_hy =  RT( 0);
            center_hw =  RT( 1);
            sqr_rad   = -FT( 1);
        }

        void  set( const Point& p)
        {
            dao.get( p, center_hx, center_hy, center_hw);
            sqr_rad = FT( 0);
        }

        void  set( const Point& p, const Point& q)
        {
            RT  phx;
            RT  phy;
            RT  phw;
            RT  qhx;
            RT  qhy;
            RT  qhw;

            dao.get( p, phx, phy, phw);
            dao.get( q, qhx, qhy, qhw);
            center_hx = ( phx*qhw + qhx*phw);
            center_hy = ( phy*qhw + qhy*phw);
            center_hw = ( phw*qhw * RT( 2));
            sqr_rad   = sqr_dist( phx, phy, phw,
                                  center_hx, center_hy, center_hw);
        }

        void  set( const Point& p, const Point& q, const Point& r)
        {
            RT  phx;
            RT  phy;
            RT  phw;
            RT  qhx;
            RT  qhy;
            RT  qhw;
            RT  rhx;
            RT  rhy;
            RT  rhw;

            dao.get( p, phx, phy, phw);
            dao.get( q, qhx, qhy, qhw);
            dao.get( r, rhx, rhy, rhw);

            RT  qhx_phx( qhx*phw - phx*qhw);
            RT  qhy_phy( qhy*phw - phy*qhw);    // denominator: qhw*phw

            RT  rhx_phx( rhx*phw - phx*rhw);
            RT  rhy_phy( rhy*phw - phy*rhw);    // denominator: rhw*phw
 
            RT  phw2( phw*phw);
            RT  qhw2( qhw*qhw);
            RT  rhw2( rhw*rhw);

            RT  p2( phx*phx + phy*phy);         // denominator: phw2

            RT  q2_p2( ( qhx*qhx + qhy*qhy) * phw2 - p2 * qhw2);
                                                // denominator: qhw2*phw2

            RT  r2_p2( ( rhx*rhx + rhy*rhy) * phw2 - p2 * rhw2);
                                                // denominator: rhw2*phw2

            center_hx = q2_p2*rhy_phy * rhw - r2_p2*qhy_phy * qhw;
            center_hy = r2_p2*qhx_phx * qhw - q2_p2*rhx_phx * rhw;
            center_hw = ( qhx_phx*rhy_phy - rhx_phx*qhy_phy)
                          * phw*qhw*rhw * RT( 2);
            sqr_rad   = sqr_dist( phx, phy, phw,
                                  center_hx, center_hy, center_hw);
        }

        // predicates
        CGAL::Bounded_side
        bounded_side( const Point& p) const
        {
            RT  phx;
            RT  phy;
            RT  phw;
            dao.get( p, phx, phy, phw);
            return( CGAL::Bounded_side( CGAL_NTS sign(
                sqr_rad - sqr_dist( phx, phy, phw,
                                    center_hx, center_hy, center_hw))));
        }

        bool
        has_on_bounded_side( const Point& p) const
        {
            RT  phx;
            RT  phy;
            RT  phw;
            dao.get( p, phx, phy, phw);
            return( sqr_dist( phx, phy, phw,
                              center_hx, center_hy, center_hw) < sqr_rad);
        }

        bool
        has_on_boundary( const Point& p) const
        {
            RT  phx;
            RT  phy;
            RT  phw;
            dao.get( p, phx, phy, phw);
            return( sqr_dist( phx, phy, phw,
                              center_hx, center_hy, center_hw) == sqr_rad);
        }

        bool
        has_on_unbounded_side( const Point& p) const
        {
            RT  phx;
            RT  phy;
            RT  phw;
            dao.get( p, phx, phy, phw);
            return( sqr_rad < sqr_dist( phx, phy, phw,
                                        center_hx, center_hy, center_hw));
        }

        bool
        is_empty( ) const
        {
            return( CGAL::is_negative( sqr_rad));
        }

        bool
        is_degenerate( ) const
        {
            return( ! CGAL::is_positive( sqr_rad));
        }

        // additional operations for checking
        bool
        operator == (
            const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c) const
        {
            return( ( center_hx*c.center_hw == c.center_hx*center_hw) &&
                    ( center_hy*c.center_hw == c.center_hy*center_hw) &&
                    ( sqr_rad  == c.sqr_rad ) );
        }

        bool
        operator != (
            const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c) const
        {
            return( ! ( *this == c));
        }

        Point
        center( ) const
        {
            Point  p;
            dao.set( p, center_hx, center_hy, center_hw);
            return( p);
        }

        const Distance&
        squared_radius( ) const
        {
            return( sqr_rad);
        }
    };

    // I/O
    template < class PT_, class DA_ >
    std::ostream&
    operator << ( std::ostream& os,
                  const CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c)
    {
        switch ( CGAL::get_mode( os)) {

          case CGAL::IO::PRETTY:
            os << "CGAL::Min_circle_2_adapterH2::Circle( "
               << c.center_hx << ", "
               << c.center_hy << ", "
               << c.center_hw << ", "
               << c.sqr_rad   << ')';
            break;

          case CGAL::IO::ASCII:
            os << c.center_hx << ' '
               << c.center_hy << ' '
               << c.center_hw << ' '
               << c.sqr_rad;
            break;

          case CGAL::IO::BINARY:
            CGAL::write( os, c.center_hx);
            CGAL::write( os, c.center_hy);
            CGAL::write( os, c.center_hw);
            CGAL::write( os, c.sqr_rad);
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                            "CGAL::get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class PT_, class DA_ >
    std::istream&
    operator >> ( std::istream& is,
                  CGAL::_Min_circle_2_adapterH2__Circle<PT_,DA_>& c)
    {
        switch ( CGAL::get_mode( is)) {

          case CGAL::IO::PRETTY:
            cerr << std::endl;
            cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case CGAL::IO::ASCII:
            is >> c.center_hx >> c.center_hy >> c.center_hw >> c.sqr_rad;
            break;

          case CGAL::IO::BINARY:
            CGAL::read( is, c.center_hx);
            CGAL::read( is, c.center_hy);
            CGAL::read( is, c.center_hw);
            CGAL::read( is, c.sqr_rad);
            break;

          default:
            CGAL_optimisation_assertion_msg( false,
                                             "CGAL::IO::mode invalid!");
            break; }

        return( is);
    }
@end

@! ============================================================================
@! Tests
@! ============================================================================

\clearpage
\section{Test}

We test \ccc{Min_circle_2} with the traits class implementation
for optimisation algorithms, using exact arithmetic, i.e.\ Cartesian
representation with number type \ccc{Quotient<Gmpz>} and
homogeneous representation with number type \ccc{Gmpz}.

@macro <Min_circle_2 test (includes and typedefs)> = @begin
    #include <CGAL/Cartesian.h>
    #include <CGAL/Homogeneous.h>
    #include <CGAL/Min_circle_2.h>
    #include <CGAL/Min_circle_2_traits_2.h>
    #include <CGAL/Min_circle_2_adapterC2.h>
    #include <CGAL/Min_circle_2_adapterH2.h>
    #include <CGAL/IO/Verbose_ostream.h>
    #include <cassert>
    #include <cstring>
    #include <fstream>

    #ifdef CGAL_USE_LEDA_FOR_OPTIMISATION_TEST
    #  include <CGAL/leda_integer.h>
       typedef  leda_integer                      Rt;
       typedef  CGAL::Quotient< leda_integer >    Ft;
    #else
    #  include <CGAL/Gmpz.h>
       typedef  CGAL::Gmpz                        Rt;
       typedef  CGAL::Quotient< CGAL::Gmpz >      Ft;
    #endif

    typedef  CGAL::Cartesian< Ft >                KerC;
    typedef  CGAL::Homogeneous< Rt >              KerH;
    typedef  CGAL::Min_circle_2_traits_2< KerC >  TraitsC;
    typedef  CGAL::Min_circle_2_traits_2< KerH >  TraitsH;
@end

The command line option \ccc{-verbose} enables verbose output.

@macro <Min_circle_2 test (verbose option)> = @begin
    bool  verbose = false;
    if ( ( argc > 1) && ( CGAL_CLIB_STD::strcmp( argv[ 1], "-verbose") == 0)) {
        verbose = true;
        --argc;
        ++argv; }
@end

@! ----------------------------------------------------------------------------
@! Code Coverage
@! ----------------------------------------------------------------------------

\subsection{Code Coverage}

We call each function of class \ccc{Min_circle_2<Traits>} at least
once to ensure code coverage.

@macro <Min_circle_2 test (code coverage)> = @begin
    cover_Min_circle_2( verbose, TraitsC(), Rt());
    cover_Min_circle_2( verbose, TraitsH(), Rt());
@end

@macro <Min_circle_2 test (code coverage test function)> = @begin
    template < class Traits, class RT >
    void
    cover_Min_circle_2( bool verbose, const Traits&, const RT&)
    {
        CGAL_USING_NAMESPACE_STD

        typedef  CGAL::Min_circle_2< Traits >  Min_circle;
        typedef  typename Min_circle::Point    Point;
        typedef  typename Min_circle::Circle   Circle;

        CGAL::Verbose_ostream verr( verbose);

        // generate `n' points at random
        const int     n = 20;
        CGAL::Random  random_x, random_y;
        Point         random_points[ n];
        int           i;
        verr << n << " random points from [0,128)^2:" << endl;
        for ( i = 0; i < n; ++i) {
            random_points[ i] = Point( RT( random_x( 128)),
                                       RT( random_y( 128)));
            verr << i << ": " << random_points[ i] << endl;
        }

        // cover code
        verr << endl << "default constructor...";
        {
            Min_circle  mc;
            bool  is_valid = mc.is_valid( verbose);
            bool  is_empty = mc.is_empty();
            assert( is_valid); 
            assert( is_empty);
        }

        verr << endl << "one point constructor...";
        {
            Min_circle  mc( random_points[ 0]);
            bool  is_valid      = mc.is_valid( verbose);
            bool  is_degenerate = mc.is_degenerate();
            assert( is_valid);
            assert( is_degenerate);
        }

        verr << endl << "two points constructor...";
        {
            Min_circle  mc( random_points[ 1],
                            random_points[ 2]);
            bool  is_valid = mc.is_valid( verbose);
            assert( is_valid);
            assert( mc.number_of_points() == 2);
        }

        verr << endl << "three points constructor...";
        {    
            Min_circle  mc( random_points[ 3],
                            random_points[ 4],
                            random_points[ 5]);
            bool  is_valid = mc.is_valid( verbose);
            assert( is_valid);
            assert( mc.number_of_points() == 3);
        }

        verr << endl << "Point* constructor...";
        Min_circle  mc( random_points, random_points+9, false);
        {
            Min_circle  mc2( random_points, random_points+9, true);
            bool  is_valid  = mc .is_valid( verbose);
            bool  is_valid2 = mc2.is_valid( verbose);
            assert( is_valid);
            assert( is_valid2);
            assert( mc .number_of_points() == 9);
            assert( mc2.number_of_points() == 9);
            assert( mc.circle() == mc2.circle());
        }

        verr << endl << "list<Point>::const_iterator constructor...";
        {
            Min_circle  mc1( mc.points_begin(), mc.points_end(), false);
            Min_circle  mc2( mc.points_begin(), mc.points_end(), true);
            bool  is_valid1 = mc1.is_valid( verbose);
            bool  is_valid2 = mc2.is_valid( verbose);
            assert( is_valid1);
            assert( is_valid2);
            assert( mc1.number_of_points() == 9);
            assert( mc2.number_of_points() == 9);
            assert( mc.circle() == mc1.circle());
            assert( mc.circle() == mc2.circle());
        }

        verr << endl << "#points already called above.";

        verr << endl << "points access already called above.";

        verr << endl << "support points access...";
        {
            typedef  typename Min_circle::Support_point_iterator
                                              Support_point_iterator;
            Point                   support_point;
            Support_point_iterator  iter( mc.support_points_begin());
            for ( i = 0; i < mc.number_of_support_points(); ++i, ++iter) {
                support_point = mc.support_point( i);
                assert( support_point == *iter); }
            Support_point_iterator  end_iter( mc.support_points_end());
            assert( iter == end_iter);
        }

        verr << endl << "circle access already called above...";

        verr << endl << "in-circle predicates...";
        {
            Point               p;
            CGAL::Bounded_side  bounded_side;
            bool                has_on_bounded_side;
            bool                has_on_boundary;
            bool                has_on_unbounded_side;
            for ( i = 0; i < 9; ++i) {
                p = random_points[ i];
                bounded_side          = mc.bounded_side( p);
                has_on_bounded_side   = mc.has_on_bounded_side( p);
                has_on_boundary       = mc.has_on_boundary( p);
                has_on_unbounded_side = mc.has_on_unbounded_side( p);
            assert( bounded_side != CGAL::ON_UNBOUNDED_SIDE);
            assert( has_on_bounded_side || has_on_boundary);
            assert( ! has_on_unbounded_side); }
        }

        verr << endl << "is_... predicates already called above.";

        verr << endl << "single point insert...";
        mc.insert( random_points[ 9]);
        {
            bool  is_valid = mc.is_valid( verbose);
            assert( is_valid);
            assert( mc.number_of_points() == 10);
        }

        verr << endl << "Point* insert...";
        mc.insert( random_points+10, random_points+n);
        {
            bool  is_valid = mc.is_valid( verbose);
            assert( is_valid);
            assert( mc.number_of_points() == n);
        }

        verr << endl << "list<Point>::const_iterator insert...";
        {
            Min_circle  mc2;
            mc2.insert( mc.points_begin(), mc.points_end());
            bool  is_valid = mc2.is_valid( verbose);
            assert( is_valid);
            assert( mc2.number_of_points() == n);
            
            verr << endl << "clear...";
            mc2.clear();        
                  is_valid = mc2.is_valid( verbose);
            bool  is_empty = mc2.is_empty();
            assert( is_valid); 
            assert( is_empty);
        }

        verr << endl << "validity check already called several times.";

        verr << endl << "traits class access...";
        {
            Traits  traits( mc.traits());
        }

        verr << endl << "I/O...";
        {
            verr << endl << "  writing `test_Min_circle_2.ascii'...";
            ofstream os( "test_Min_circle_2.ascii");
            CGAL::set_ascii_mode( os);
            os << mc;
        }
        {
            verr << endl << "  writing `test_Min_circle_2.pretty'...";
            ofstream os( "test_Min_circle_2.pretty");
            CGAL::set_pretty_mode( os);
            os << mc;
        }
        {
            verr << endl << "  writing `test_Min_circle_2.binary'...";
            ofstream os( "test_Min_circle_2.binary");
            CGAL::set_binary_mode( os);
            os << mc;
        }
        {
            verr << endl << "  reading `test_Min_circle_2.ascii'...";
            Min_circle mc_in;
            ifstream is( "test_Min_circle_2.ascii");
            CGAL::set_ascii_mode( is);
            is >> mc_in;
            bool    is_valid = mc_in.is_valid( verbose);
            assert( is_valid);
            assert( mc_in.number_of_points() == n);
            assert( mc_in.circle() == mc.circle());
        }
        verr << endl;
    }
@end

@! ----------------------------------------------------------------------------
@! Adapters
@! ----------------------------------------------------------------------------

\subsection{Traits Class Adapters}

We define two point classes (one with Cartesian, one with homogeneous
representation) and corresponding data accessors.

@macro <Min_circle_2 test (point classes)> = @begin
    // 2D Cartesian point class
    class MyPointC2;

    std::ostream&  operator << ( std::ostream&, const MyPointC2&);
    std::istream&  operator >> ( std::istream&,       MyPointC2&);
    
    class MyPointC2 {
      public:
        typedef  ::Ft  FT;
      private:
        FT x_;
        FT y_;
      public:
        MyPointC2( ) { }
        MyPointC2( const FT& x, const FT& y) : x_( x), y_( y) { }

        const FT&  x( ) const { return( x_); }
        const FT&  y( ) const { return( y_); }

        bool
        operator == ( const MyPointC2& p) const
        {
            return( ( x_ == p.x_) && ( y_ == p.y_));
        }

        bool
        operator != ( const MyPointC2& p) const
        {
            return( ( x_ != p.x_) || ( y_ != p.y_));
        }

        friend
        std::ostream&  operator << ( std::ostream& os, const MyPointC2& p);

        friend
        std::istream&  operator >> ( std::istream& is,       MyPointC2& p);
    };

    std::ostream&
    operator << ( std::ostream& os, const MyPointC2& p)
    {
        return( os << p.x_ << ' ' << p.y_);
    }

    std::istream&
    operator >> ( std::istream& is, MyPointC2& p)
    {
        return( is >> p.x_ >> p.y_);
    }
    
    // 2D Cartesian point class data accessor
    class MyPointC2DA {
      public:
        typedef  ::Ft  FT;

	MyPointC2DA( ) { }

        const FT&  get_x( const MyPointC2& p) const { return( p.x()); }
        const FT&  get_y( const MyPointC2& p) const { return( p.y()); }

        void
        get( const MyPointC2& p, FT& x, FT& y) const
        {
            x = get_x( p);
            y = get_y( p);
        }

        void
        set( MyPointC2& p, const FT& x, const FT& y) const
        {
            p = MyPointC2( x, y);
        }
    };


    // 2D homogeneous point class
    class MyPointH2;

    std::ostream&  operator << ( std::ostream&, const MyPointH2&);
    std::istream&  operator >> ( std::istream&,       MyPointH2&);
    
    class MyPointH2 {
      public:
        typedef  ::Rt  RT;
      private:
        RT hx_;
        RT hy_;
        RT hw_;
      public:
        MyPointH2( ) { }
        MyPointH2( const RT& hx, const RT& hy, const RT& hw = RT( 1))
            : hx_( hx), hy_( hy), hw_( hw) { }

        const RT&  hx( ) const { return( hx_); }
        const RT&  hy( ) const { return( hy_); }
        const RT&  hw( ) const { return( hw_); }

        bool
        operator == ( const MyPointH2& p) const
        {
            return( ( hx_*p.hw_ == p.hx_*hw_) && ( hy_*p.hw_ == p.hy_*hw_));
        }

        bool
        operator != ( const MyPointH2& p) const
        {
            return( ( hx_*p.hw_ != p.hx_*hw_) || ( hy_*p.hw_ != p.hy_*hw_));
        }

        friend
        std::ostream&  operator << ( std::ostream& os, const MyPointH2& p);

        friend
        std::istream&  operator >> ( std::istream& is,       MyPointH2& p);
    };

    std::ostream&
    operator << ( std::ostream& os, const MyPointH2& p)
    {
        return( os << p.hx_ << ' ' << p.hy_ << ' ' << p.hw_);
    }

    std::istream&
    operator >> ( std::istream& is, MyPointH2& p)
    {
        return( is >> p.hx_ >> p.hy_ >> p.hw_);
    }

    // 2D homogeneous point class data accessor
    class MyPointH2DA {
      public:
        typedef  ::Rt  RT;

	MyPointH2DA( ) { }

        const RT&  get_hx( const MyPointH2& p) const { return( p.hx()); }
        const RT&  get_hy( const MyPointH2& p) const { return( p.hy()); }
        const RT&  get_hw( const MyPointH2& p) const { return( p.hw()); }

        void
        get( const MyPointH2& p, RT& hx, RT& hy, RT& hw) const
        {
            hx = get_hx( p);
            hy = get_hy( p);
            hw = get_hw( p);
        }

        void
        set( MyPointH2& p, const RT& hx, const RT& hy, const RT& hw) const
        {
            p = MyPointH2( hx, hy, hw);
        }
    };

    CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( MyPointC2)
    CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( MyPointH2)
@end

To test the traits class adapters we use the code coverage test function.

@macro <Min_circle_2 test (adapters test)> = @begin
    typedef  CGAL::Min_circle_2_adapterC2< MyPointC2, MyPointC2DA >  AdapterC2;
    typedef  CGAL::Min_circle_2_adapterH2< MyPointH2, MyPointH2DA >  AdapterH2;
    cover_Min_circle_2( verbose, AdapterC2(), Rt());
    cover_Min_circle_2( verbose, AdapterH2(), Rt());
@end

@! ----------------------------------------------------------------------------
@! External Test Sets
@! ----------------------------------------------------------------------------

\subsection{External Test Sets}

In addition, some data files can be given as command line arguments.
A data file contains pairs of \ccc{int}s, namely the x- and
y-coordinates of a set of points. The first number in the file is the
number of points. A short description of the test set is given at the
end of each file.

@macro <Min_circle_2 test (external test sets)> = @begin
    while ( argc > 1) {

        typedef  CGAL::Min_circle_2< TraitsH >  Min_circle;
        typedef  Min_circle::Point              Point;
        typedef  Min_circle::Circle             Circle;

        CGAL::Verbose_ostream verr( verbose);

        // read points from file
        verr << std::endl << "input file: `" << argv[ 1] << "'" << std::flush;

        std::list<Point>  points;
        int               n, x, y;
        std::ifstream     in( argv[ 1]);
        in >> n;
        assert( in);
        for ( int i = 0; i < n; ++i) {
            in >> x >> y;
            assert( in);
            points.push_back( Point( x, y)); }

        // compute and check min_circle
        Min_circle  mc2( points.begin(), points.end(), false);
        bool  is_valid = mc2.is_valid( verbose);
        assert( is_valid);

        // next file
        --argc;
        ++argv; }
@end

@! ==========================================================================
@! Files
@! ==========================================================================

\clearpage
\section{Files}

@i share/namespace.awi

@! ----------------------------------------------------------------------------
@! Min_circle_2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2.h}

@file <include/CGAL/Min_circle_2.h> = @begin
    @<file header>(
        "include/CGAL/Min_circle_2.h",
        "2D Smallest Enclosing Circle")

    #ifndef CGAL_MIN_CIRCLE_2_H
    #define CGAL_MIN_CIRCLE_2_H

    // includes
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif
    #ifndef CGAL_RANDOM_H
    #  include <CGAL/Random.h>
    #endif
    #ifndef CGAL_PROTECT_LIST
    #  include <list>
    #  define CGAL_PROTECT_LIST
    #endif
    #ifndef CGAL_PROTECT_VECTOR
    #  include <vector>
    #  define CGAL_PROTECT_VECTOR
    #endif
    #ifndef CGAL_PROTECT_ALGORITHM
    #  include <algorithm>
    #  define CGAL_PROTECT_ALGORITHM
    #endif
    #ifndef CGAL_PROTECT_IOSTREAM
    #  include <iostream>
    #  define CGAL_PROTECT_IOSTREAM
    #endif

    @<namespace begin>("CGAL")
    
    // Class declaration
    // =================
    @<Min_circle_2 declaration>

    // Class interface
    // ===============
    @<Min_circle_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Min_circle_2 I/O operators declaration>

    @<namespace end>("CGAL")
    
    #ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #  include <CGAL/Min_circle_2.C>
    #endif

    #endif // CGAL_MIN_CIRCLE_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_circle_2.C
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2.C}

@file <include/CGAL/Min_circle_2.C> = @begin
    @<file header>(
        "include/CGAL/Min_circle_2.C",
        "2D Smallest Enclosing Circle")

    @<namespace begin>("CGAL")
    
    // Class implementation (continued)
    // ================================
    // I/O
    // ---
    @<Min_circle_2 I/O operators>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_circle_2.h
@! ----------------------------------------------------------------------------

\subsection{Optimisation\_circle\_2.h}

@file <include/CGAL/Optimisation_circle_2.h> = @begin
    @<file header>(
        "include/CGAL/Optimisation_circle_2.h",
        "2D Optimisation Circle")

    #ifndef CGAL_OPTIMISATION_CIRCLE_2_H
    #define CGAL_OPTIMISATION_CIRCLE_2_H

    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_BASIC_CONSTRUCTIONS_2_H
    #  include <CGAL/basic_constructions_2.h>
    #endif
    #ifndef CGAL_SQUARED_DISTANCE_2_H
    #  include <CGAL/squared_distance_2.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declaration
    // =================
    @<Optimisation_circle_2 declaration>

    // Class interface
    // ===============
    @<Optimisation_circle_2 interface>

    // Function declarations
    // =====================
    // I/O
    // ---
    @<Optimisation_circle_2 I/O operators declaration>

    @<namespace end>("CGAL")
    
    #ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #  include <CGAL/Optimisation_circle_2.C>
    #endif

    #endif // CGAL_OPTIMISATION_CIRCLE_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Optimisation_circle_2.C
@! ----------------------------------------------------------------------------

\subsection{Optimisation\_circle\_2.C}

@file <include/CGAL/Optimisation_circle_2.C> = @begin
    @<file header>(
        "include/CGAL/Optimisation_circle_2.C",
        "2D Optimisation Circle")

    // includes
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class implementation (continued)
    // ================================

    // I/O
    // ---
    @<Optimisation_circle_2 I/O operators>

    @<namespace end>("CGAL")
    
    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_circle_2_traits_2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2\_traits\_2.h}

@file <include/CGAL/Min_circle_2_traits_2.h> = @begin
    @<file header>(
        "include/CGAL/Min_circle_2_traits_2.h",
        "default traits class for 2D Smallest Enclosing Circle")

    #ifndef CGAL_MIN_CIRCLE_2_TRAITS_2_H
    #define CGAL_MIN_CIRCLE_2_TRAITS_2_H

    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_OPTIMISATION_CIRCLE_2_H
    #  include <CGAL/Optimisation_circle_2.h>
    #endif
    #ifndef CGAL_PREDICATES_ON_POINTS_2_H
    #  include <CGAL/predicates_on_points_2.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_circle_2_traits_2 declarations>

    // Class interface and implementation
    // ==================================
    @<Min_circle_2_traits_2 interface and implementation>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_CIRCLE_2_TRAITS_2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_circle_2_adapterC2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2\_adapterC2.h}

@file <include/CGAL/Min_circle_2_adapterC2.h> = @begin
    @<file header>(
        "include/CGAL/Min_circle_2_adapterC2.h",
        "traits class adapter for 2D Smallest Enclosing Circle")

    #ifndef CGAL_MIN_CIRCLE_2_ADAPTERC2_H
    #define CGAL_MIN_CIRCLE_2_ADAPTERC2_H

    // includes
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_circle_2_adapterC2 declarations>

    // Class interface and implementation
    // ==================================
    @<Min_circle_2_adapterC2 interface and implementation>

    // Nested type `Circle'
    @<Min_circle_2_adapterC2 nested type `Circle'>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_CIRCLE_2_ADAPTERC2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_circle_2_adapterH2.h
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2\_adapterH2.h}

@file <include/CGAL/Min_circle_2_adapterH2.h> = @begin
    @<file header>(
        "include/CGAL/Min_circle_2_adapterH2.h",
        "traits class adapter for 2D Smallest Enclosing Circle")

    #ifndef CGAL_MIN_CIRCLE_2_ADAPTERH2_H
    #define CGAL_MIN_CIRCLE_2_ADAPTERH2_H

    // includes
    #ifndef CGAL_OPTIMISATION_BASIC_H
    #  include <CGAL/Optimisation/basic.h>
    #endif

    @<namespace begin>("CGAL")
    
    // Class declarations
    // ==================
    @<Min_circle_2_adapterH2 declarations>

    // Class interface and implementation
    // ==================================
    @<Min_circle_2_adapterH2 interface and implementation>

    // Nested type `Circle'
    @<Min_circle_2_adapterH2 nested type `Circle'>

    @<namespace end>("CGAL")
    
    #endif // CGAL_MIN_CIRCLE_2_ADAPTERH2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Min_circle_2_Window_stream.h
@! ----------------------------------------------------------------------------

\subsection{Min\_circle\_2\_Window\_stream.h}

@file <include/CGAL/IO/Min_circle_2_Window_stream.h> = @begin
    @<file header>(
        "include/CGAL/IO/Min_circle_2_Window_stream.h",
        "graphical output to `leda_window' for Min_circle_2 algorith.")

    // Each of the following operators is individually 
    // protected against multiple inclusion.

    // Window_stream I/O operators
    // ===========================

    // Optimisation_circle_2
    // ---------------------
    @<Optimisation_circle_2 graphical output operator>

    // Min_circle_2
    // ------------
    @<Min_circle_2 graphical output operator>

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! test_Min_circle_2.C
@! ----------------------------------------------------------------------------

\subsection{test\_Min\_circle\_2.C}

@file <test/Min_circle_2/test_Min_circle_2.C> = @begin
    @<file header>(
        "test/Min_circle_2/test_Min_circle_2.C",
        "test program for 2D Smallest Enclosing Circle")

    @<Min_circle_2 test (includes and typedefs)>

    // code coverage test function
    // ---------------------------
    @<Min_circle_2 test (code coverage test function)>

    // point classes for adapters test
    // -------------------------------
    @<Min_circle_2 test (point classes)>

    // main
    // ----
    int
    main( int argc, char* argv[])
    {
        // command line options
        // --------------------
        // option `-verbose'
        @<Min_circle_2 test (verbose option)>

        // code coverage
        // -------------
        @<Min_circle_2 test (code coverage)>

        // adapters test
        // -------------
        @<Min_circle_2 test (adapters test)>

        // external test sets
        // -------------------
        @<Min_circle_2 test (external test sets)>

        return( 0);
    }

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! File Header
@! ----------------------------------------------------------------------------

\subsection*{File Header}

@i share/file_header.awi
 
And here comes the specific file header for the product files of this
web file.

@macro <file header>(2) many = @begin
    @<copyright notice>
    @<file name>(@1)
    @<file description>(
        "Min_circle_2",
        "Geometric Optimisation",
        "Min_circle_2",
        "$Revision$","$Date$",
        "Sven Schönherr <sven@@inf.ethz.ch>, Bernd Gärtner",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ============================================================================
@! Bibliography
@! ============================================================================

\clearpage
\bibliographystyle{plain}
\bibliography{geom,../doc_tex/basic/Optimisation/cgal}

@! ===== EOF ==================================================================
