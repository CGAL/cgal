@! ============================================================================
@! The CGAL Project
@! Implementation: 2D Smallest Enclosing Circle
@! ----------------------------------------------------------------------------
@! file  : Library/web/Min_circle_2.aw
@! author: Bernd Gärtner, Sven Schönherr (sven@inf.fu-berlin.de)
@! $Id$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage{a4wide2}

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

@article
@p maximum_input_line_length = 180
@p maximum_output_line_length = 180

@thickline

@t vskip 5 mm
@t title titlefont centre "CGAL: 2D Smallest Enclosing Circle"
@t vskip 1 mm
@t title smalltitlefont centre "Implementation Documentation"
@t vskip 8 mm
@t title smalltitlefont centre "Bernd Gärtner, Sven Schönherr"
@t title normalfont centre "$Revision$, $Date$"
@t vskip 1 mm

@thickline
@thinline
@t table_of_contents

@thinline

@! ============================================================================
@section{ Introduction}
@! ============================================================================

We define a class template @prg{CGAL_Min_circle_2<R>}. An object of
this class represents the smallest (w.r.t. area) enclosing circle of a
finite point set $P$ in the plane, denoted by $\mc(P)$. The template
parameter of the representation class $R$ (the domain of the point
coordinates) is any number type that fulfills the CGAL number type
requirements, but correct results are in this version only guaranteed
if the number type is an exact one (like LEDA's @prg{integer}). In
particular, using the number type @prg{double} might lead to wrong
results. A correct @prg{double} implementation is planned for the next
release.

The implementation is based on an algorithm by Welzl \cite{Wel}, which
we shortly describe now. $\mc(P)$ is built up incrementally, adding
one point after another. Assume $\mc(P)$ has been constructed, and we
would like to obtain $\mc(P \cup \{p\})$, $p$ some new point. There
are two cases: if $p$ already lies inside $\mc(P)$, then $\mc(P \cup
\{p\}) = \mc(P)$. Otherwise $p$ must lie on the boundary of $\mc(P
\cup \{p\})$ (this is proved in \cite{Wel} and not hard to see), so we
need to compute $\mc(P,\{p\})$, the smallest circle enclosing $P$ with
$p$ on the boundary. This is recursively done in the same manner. In
general, for point sets $P$,$B$, define $\mc(P,B)$ as the smallest
circle enclosing $P$ that has the points of $B$ on the boundary (if
defined). Although the algorithm finally delivers a circle
$\mc(P,\emptyset)$, it internally deals with circles that have a
possibly nonempty set $B$. Here is the pseudocode of Welzl's method.
To compute $\mc(P)$, it is called with the pair $(P,\emptyset)$,
assuming that $P=\{p_1,\ldots,p_n\}$ is stored in a linked list.

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
computation of $\mc(\emptyset,B)$ is easy. In our implementation, it
is done by the private member function @prg{compute_circle}. (b) One
can check that the method maintains the invariant `$\mc(P,B)$ exists'.
This justifies termination if $|B| = 3$, because then $\mc(P,B)$ must
be the unique circle with the points of $B$ on the boundary, and
$\mc(P,B)$ exists if and only if this circle contains the points of
$P$. Thus, no subsequent in-circle tests are necessary anymore (for
details see \cite{Wel}). (c) points which are found to lie outside the
current circle $mc$ are considered `important' and are moved to the
front of the linked list that stores $P$. This is crucial for the
method's efficiency.

It can also be advisable to bring $P$ into random order before
computation starts. There are `bad' insertion orders which cause the
method to be very slow -- random shuffling gives these orders a very
small probability.

The private member function @prg{mc} directly realizes the pseudocode
above.

@! ============================================================================
@section{ Class \texttt{CGAL\_Min\_circle\_2}: Interface}
@! ============================================================================

The class interface looks as follows.

@macro <Min_circle_2 interface> = @begin
    template < class R >
    class CGAL_Min_circle_2 {
      public:
	@<public interface>

      private:
	@<private data members>
	@<private member functions declaration>
    };
@end   

The public interface is described and documented in the CGAL Reference
Manual, so we do not comment on it here.

@macro <public interface> = @begin
    // creation
    CGAL_Min_circle_2( );
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p);
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2);
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2,
                       const CGAL_Point_2<R>& p3);
    CGAL_Min_circle_2( const CGAL_Point_2<R>* first,
                       const CGAL_Point_2<R>* last,
                       bool randomize = false);
    ~CGAL_Min_circle_2( );

    // access
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;

    const CGAL_Point_2<R>&  point        ( int i) const;
    const CGAL_Point_2<R>&  support_point( int i) const;
    const CGAL_Point_2<R>&  operator []  ( int i) const;

    const CGAL_Circle_2<R>&  circle( ) const;

    CGAL_Bbox_2  bbox( ) const;

    // updates
    void  insert ( const CGAL_Point_2<R>& p);
    void  reserve( int n);

    // check
    bool  check( bool verbose = false) const;

    // predicates
    CGAL_Bounded_side  bounded_side( const CGAL_Point_2<R>& p) const;
    bool  has_on_bounded_side      ( const CGAL_Point_2<R>& p) const;
    bool  has_on_boundary          ( const CGAL_Point_2<R>& p) const;
    bool  has_on_unbounded_side    ( const CGAL_Point_2<R>& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
@end

The private member functions -- the important ones being
@prg{compute_circle} and @prg{mc} -- are in detail described in the
implementation section below.

@macro <private member functions declaration> = @begin
    // copying and assignment not allowed!
    CGAL_Min_circle_2( const CGAL_Min_circle_2<R>&);
    CGAL_Min_circle_2<R>& operator = ( const CGAL_Min_circle_2<R>&);

    // private member functions
    void mc( int n, int n_sp);
    void compute_circle( );
    void check_error_msg( char* error_msg) const;
@end

Now let us focus on the private data members.

The points of $P$ are internally stored as a linked list that allows
to bring points to the front of the list in constant time. The linked
list is realized as an STL-@prg{vector} whose entries are points with
successor and predecessor information, realized in the class
@prg{_Point}.

@macro <private data members> += @begin
    // doubly linked list element, used for move-to-front heuristic
    class _Point {
      public:
        CGAL_Point_2<R>  point;
        int              pred, succ;
        _Point( ) { }
        _Point( const CGAL_Point_2<R>& _point, int _pred = -1, int _succ = -1)
            : point( _point), pred( _pred), succ( _succ)
        { }
    };
@end

Here is the vector @prg{points} itself, along with an index
@prg{i_first_point} in @prg{points} referring to the current front
element of the linked list.

@macro <private data members> += @begin
    vector< _Point >  points;
    int               i_first_point;
@end

A @prg{CGAL_Min_circle_2} object maintains an array
@prg{support_points} of at most three \emph{support points} (the
actual number is given by @prg{n_support_points}, which at the end of
the computation contains the indices of a minimal subset $S \subseteq
P$ with $\mc(P)= \mc(S)$. (Note: that subset is not necessarily
\emph{minimum}). During the computations, the set of support points
coincides with the set $B$ appearing in the pseudocode for $\mc(P,B)$,
see the introduction above.

@macro <private data members> += @begin
    int               n_support_points;
    CGAL_Point_2<R>   support_points[ 3];
@end

Finally, an actual circle is stored in a @prg{CGAL_Min_circle_2}
object, by the end of computation equal to $\mc(P)$. During
computation, this circle equals the circle $mc$ appearing in the
pseudocode for $\mc(P,B)$, see the introduction above.

@macro <private data members> += @begin
    CGAL_Circle_2<R>  min_circle;
@end  


@! ============================================================================
@section{ Class \texttt{CGAL\_Min\_circle\_2}: Implementation}
@! ============================================================================

The implementation consists of several parts, each of which is
described in the sequel. The actual work is hidden in the functions
@prg{mc} and @prg{compute_circle}.

@macro <Min_circle_2 implementation> = @begin
    @<access operations and predicates>

    @<private member functions>

    @<constructors>
    @<update operations>
    @<check operation>
@end


@! ----------------------------------------------------------------------------
@subsection{ Constructors}
@! ----------------------------------------------------------------------------

In this version, @prg{CGAL_Min_circle_2} has five different
construction methods, where the most important one builds the
@prg{CGAL_Min_circle_2} $\mc(P)$ from a point set $P$, stored in a
piece of memory delimited by two pointers (As soon as member templates
are available, the pointers will become iterators). The constructor
copies the points into the internal array @prg{points}, performs a
random reordering if the corresponding flag requires that, then builds
up the linked list over the point set, and finally calls the private
method $\mc$ to compute $\mc(P)=\mc(P,\emptyset)$.

@macro <constructors> += @begin
    // constructors
    // ------------
    template < class R >
    CGAL_Min_circle_2<R>::
    CGAL_Min_circle_2( const CGAL_Point_2<R>* first,
		       const CGAL_Point_2<R>* last,
		       bool randomize)
    {
	// store points
	int n = 0;
	distance( first, last, n);
	if ( n > 0) {
	    points.reserve( n);
	    copy( first, last, back_inserter( points));

	    // shuffle points at random
	    if ( randomize)
		random_shuffle( points.begin(), points.end());

	    // link points
	    for ( int i = 0; i < number_of_points(); ++i) {
		points[ i].pred = i-1;
		points[ i].succ = i+1; }
	    points[ 0].pred = points[ number_of_points()-1].succ = -1;
	    i_first_point = 0; }
	else
	    i_first_point = -1;

	// compute mc
	mc( points.size(), 0);
    }
@end

The remaining constructors are actually specializations of the
previous one, building a @prg{CGAL_Min_circle_2} for up to three
points. The idea is the following: recall that for any point set $P$
there exists $S \subseteq P$, $|S| \leq 3$ with $\mc(S) = \mc(P)$ (in
fact, such a set $S$ is determined during a @prg{CGAL_Min_circle_2}
computation). Once $S$ has been computed (or given otherwise),
$\mc(P)$ can easily be reconstructed from $S$ in constant time. To
make this reconstruction more convenient, a method is available for
each size of $|S|$, ranging from 0 to 3. For $|S| = 0$, we get the
default constructor, building the empty @prg{CGAL_Min_circle_2}
$\mc(\emptyset)$.

@macro <constructors> += @begin

    // default constructor
    template < class R >
    CGAL_Min_circle_2<R>::
    CGAL_Min_circle_2( )
      : i_first_point( -1),
	n_support_points( 0),
	min_circle( CGAL_Point_2<R>( ), R::FT( 0))
    { }

    // constructor for one point
    template < class R >
    CGAL_Min_circle_2<R>::
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p)
      : points( 1, _Point( p)),
	i_first_point( 0),
	n_support_points( 1),
	min_circle( p, R::FT( 0))
    {
	support_points[ 0] = p;
    }

    // constructor for two points
    template < class R >
    CGAL_Min_circle_2<R>::
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2)
    {
	// store points
	points.reserve( 2);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0, -1));
	i_first_point = 0;

	// compute mc
	mc( 2, 0);
    }

    // constructor for three points
    template < class R >
    CGAL_Min_circle_2<R>::
    CGAL_Min_circle_2( const CGAL_Point_2<R>& p1,
		       const CGAL_Point_2<R>& p2,
		       const CGAL_Point_2<R>& p3)
    {
	// store points
	points.reserve( 3);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0,  2));
	points.push_back( _Point( p3,  1, -1));
	i_first_point = 0;

	// compute mc
	mc( 3, 0);
    }
@end

Finally, we have a (default) destructor.

@macro <constructors> += @begin

    // destructor
    template < class R >
    CGAL_Min_circle_2<R>::
    ~CGAL_Min_circle_2( )
    { }
@end


@! ----------------------------------------------------------------------------
@subsection{ Update operations}
@! ----------------------------------------------------------------------------

There is another way to build up $\mc(P)$, other than by supplying the
point set $P$ at once. Namely, $\mc(P)$ can be built up incrementally,
adding one point after another. If you look at the pseudocode in the
introduction, this comes quite naturally. The method @prg{insert},
applied with point $p$ to a @prg{CGAL_Min_circle_2} object
representing $\mc(P)$, computes $\mc(P \cup \{p\})$, where work has to
be done only if $p$ lies outside $\mc(P)$.  In this case, $\mc(P \cup
\{p\}) = \mc(P,\{p\})$ holds, so the method $\mc$ is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudocode and
the `main' constructor above.

@macro <update operations> += @begin
    // update operations
    // -----------------
    template < class R >
    void
    CGAL_Min_circle_2<R>::
    insert( const CGAL_Point_2<R>& p)
    {
	// store point
	int old_n = number_of_points();
	points.push_back( _Point( p));  // NOTE: p is not linked with list yet!

	// p not in current circle?
	if ( has_on_unbounded_side( p)) {

	    // p new support point
	    support_points[ 0] = p;

	    // recompute mc
	    mc( old_n, 1); }

	// make p the first point in list
	if ( old_n > 0)                                  // old list not empty?
	    points[ i_first_point].pred = old_n;
	points[ old_n].succ = i_first_point;
	i_first_point = old_n;
    }
@end

The operation @prg{reserve} does nothing but tell the
@prg{CGAL_Min_circle_2} that some number $n$ of points might
eventually be inserted, allowing the object to allocate storage for
them at once. Inserting the points without doing this might lead to
overhead caused by moving the array @prg{points} around in memory
while it grows.

@macro <update operations> += @begin

    template < class R > inline
    void
    CGAL_Min_circle_2<R>::
    reserve( int n)
    {
	points.reserve( n);
    }
@end


@! ----------------------------------------------------------------------------
@subsection{ Access Operations and Predicates}
@! ----------------------------------------------------------------------------

These operations are used to retrieve information about the current
status of the @prg{CGAL_Min_circle_2} object. They are all very simple
and mostly rely on corresponding access functions and predicates for
the data members of @prg{CGAL_Min_circle_2}.

@macro <access operations and predicates> = @begin
    // 'number_of_' operations
    // -----------------------
    template < class R > inline
    int
    CGAL_Min_circle_2<R>::
    number_of_points( ) const
    {
	return( points.size());
    }

    template < class R > inline
    int
    CGAL_Min_circle_2<R>::
    number_of_support_points( ) const
    {
	return( n_support_points);
    }

    // 'is_' predicates
    // ----------------
    template < class R > inline
    bool
    CGAL_Min_circle_2<R>::
    is_empty( ) const
    {
	return( number_of_support_points() == 0);
    }

    template < class R > inline
    bool
    CGAL_Min_circle_2<R>::
    is_degenerate( ) const
    {
	return( number_of_support_points() < 2);
    }

    // access operations
    // -----------------
    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_circle_2<R>::
    point( int i) const
    {
	CGAL_Min_circle_2_precondition( (i >= 0) && (i < number_of_points()));
	return( points[ i].point);
    }

    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_circle_2<R>::
    support_point( int i) const
    {
	CGAL_Min_circle_2_precondition( (i >= 0) &&
					(i <  number_of_support_points()));
	return( support_points[ i]);
    }

    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_circle_2<R>::
    operator [] ( int i) const
    {
	return( point( i));
    }

    template < class R > inline
    const CGAL_Circle_2<R>&
    CGAL_Min_circle_2<R>::
    circle( ) const
    {
        CGAL_Min_circle_2_precondition( ! is_empty());
	// ensure positive orientation
	if ( min_circle.orientation() == CGAL_NEGATIVE)
	    const_cast( CGAL_Circle_2<R>&, min_circle) = min_circle.opposite();
	CGAL_Min_circle_2_assertion(min_circle.orientation() == CGAL_POSITIVE);
	return( min_circle);
    }

    template < class R > inline
    CGAL_Bbox_2
    CGAL_Min_circle_2<R>::
    bbox( ) const
    {
        CGAL_Min_circle_2_precondition( ! is_empty());
	return( min_circle.bbox());
    }

    // incircle predicates
    // -------------------
    template < class R > inline
    CGAL_Bounded_side
    CGAL_Min_circle_2<R>::
    bounded_side( const CGAL_Point_2<R>& p) const
    {
	return( is_empty() ? CGAL_ON_UNBOUNDED_SIDE 
			   : min_circle.bounded_side( p));
    }

    template < class R > inline
    bool
    CGAL_Min_circle_2<R>::
    has_on_bounded_side( const CGAL_Point_2<R>& p) const
    {
	return( ( ! is_empty()) && ( min_circle.has_on_bounded_side( p)));
    }

    template < class R > inline
    bool
    CGAL_Min_circle_2<R>::
    has_on_boundary( const CGAL_Point_2<R>& p) const
    {
	return( ( ! is_empty()) && ( min_circle.has_on_boundary( p)));
    }

    template < class R > inline
    bool
    CGAL_Min_circle_2<R>::
    has_on_unbounded_side( const CGAL_Point_2<R>& p) const
    {
	return( ( is_empty()) || ( min_circle.has_on_unbounded_side( p)));
    }
@end


@! ----------------------------------------------------------------------------
@subsection{ Private Member Function \texttt{compute\_circle}}
@! ----------------------------------------------------------------------------

This is the method for computing $mc(\emptyset,B)$ the set $B$ given
by the first @prg{n_support_points} indices in the array
@prg{i_support_points}. It is realized by a simple case analysis,
noting that $|B| \leq 3$.

@macro <private member functions> += @begin
    // private member functions
    // ------------------------
    template < class R >
    void
    CGAL_Min_circle_2<R>::
    compute_circle( )
    {
	switch ( n_support_points) {
	  case 3:
	    min_circle = CGAL_Circle_2<R>( support_points[0],
					   support_points[1],
					   support_points[2]);
            break;
	  case 2: {
	    const CGAL_Point_2<R>& p0( support_points[ 0]);
	    const CGAL_Point_2<R>& p1( support_points[ 1]);
	    min_circle = CGAL_Circle_2<R>( p0 + (p1-p0)/R::RT(2), p0); }
	    break;
	  case 1:
	    min_circle = CGAL_Circle_2<R>( support_points[ 0], R::FT( 0));
	    break;
	  case 0:
	    min_circle = CGAL_Circle_2<R>( CGAL_Point_2<R>( ), R::FT( 0));
	    break;
	  default:
	    CGAL_Min_circle_2_assertion( ( n_support_points >= 0) &&
					 ( n_support_points <= 3) ); }
    }
@end


@! ----------------------------------------------------------------------------
@subsection{ Private Member Function \texttt{mc}}
@! ----------------------------------------------------------------------------

This function computes the general circle $\mc(P,B)$, where $P$
contains the first @prg{n} points stored in the array @prg{points} and
$B$ is given by the first @prg{n_support_points} indices in the array
@prg{i_support_points}. The function is directly modelled after the
pseudocode above.

@macro <private member functions> += @begin

    template < class R >
    void
    CGAL_Min_circle_2<R>::
    mc( int n, int n_sp)
    {
	// compute circle through support points
	n_support_points = n_sp;
	compute_circle( );
	if ( n_sp == 3) return;

	// test first n points
	int index = i_first_point, succ;
	for ( int i = 0; i < n; ++i) {
	    _Point& p = points[ index];
	    succ = p.succ;

	    // p not in current circle?
	    if ( has_on_unbounded_side( p.point)) {

		// recursive call with p as additional support point
		support_points[ n_sp] = p.point;
		mc( i, n_sp+1);

		// move current point to front
		if ( index != i_first_point) {               // p not first?
		    points[ p.pred].succ = succ;
		    if ( succ != -1)                         // p not last?
			points[ succ].pred = p.pred;
		    points[ i_first_point].pred = index;
		    p.pred = -1;
		    p.succ = i_first_point;
		    i_first_point = index; } }

	    // next point
	    index = succ; }
    }
@end


@! ----------------------------------------------------------------------------
@subsection{Checking}
@! ----------------------------------------------------------------------------

A @prg{CGAL_Min_circle_2} object can be checked for consistency. This
means, it is checked whether (a) the circle contains all points of its
defining set $P$, (b) the circle is the smallest circle spanned by its
support set, (c) the support set is minimal, i.e. no support point is
redundant. Ideally, the checking should be able to cope with a
malicious implementor whose intention is not to provide correct code
but to make the checker believe a wrong result.  However, since the
malicious implementor in this case also provides the checker, this
philosophy is not appropriate here. Instead, the checker is meant as a
\emph{simple} device to detect problems the implementor overlooked
although he/she tried the best to write correct code. Routines which
are seen to be correct by a simple code inspection are not checked
(for example, one can trust the in-circle test, because it is directly
mapped to the corresponding test of the kernel class @prg{Circle_2},
which me must assume to be correct here). This also applies to the
check routine itself, which we keep simple in order to allow its
correctness to be checked by just inspecting the code.

The checking may be run in verbose mode where error messages are
written to standard error stream, using the following auxiliary
function.

@macro <check operation> += @begin
    // check operation
    // ---------------
    template < class R >
    void CGAL_Min_circle_2<R>::
    check_error_msg (char* error_msg) const
    {
	cerr << "CGAL_Min_circle_2 check error: " << error_msg << endl;
	cerr.flush();
    }
@end  

Here is the actual check function, returning @prg{true} if and only if
neither the containment check (a) nor the checks (b) and (c) -- which
we perform for each number of support points separately -- return
@prg{false} before.

@macro <check operation> += @begin 

    template < class R >
    bool
    CGAL_Min_circle_2<R>::
    check( bool verbose) const
    {
	// containment check
	@<containment check>

	// support set check
	@<support set check>

	return (true);
    }
@end   

The containment check (a) is easy to perform, without reference to the
support set.

@macro <containment check> = @begin
    for ( int i = 0; i < number_of_points(); ++i) 
        if ( has_on_unbounded_side( point( i))) {
            if ( verbose)
                check_error_msg( "Min_circle_2 does not enclose all points");
            return (false); }
@end

Checks (b) and (c) are performed by distingushing between the number
of support points which may range between 0 and 3.

@macro <support set check> = @begin
    switch( n_support_points) {

      case 0: {
	@<check 0 support points> }
	break;

      case 1: {
	@<check 1 support point> }
	break;

      case 2: {
	@<check 2 support points> }
	break;

      case 3: {
	@<check 3 support points> }
	break;

      default:
	if (verbose)
	    check_error_msg( "Min_circle_2 has illegal number of \
			      support points, not between 0 and 3");
	return (false); }
@end

The case of 0 support points happens if and only if the defining point
set $P$ is empty.

@macro <check 0 support points> = @begin
    if ( ! is_empty()) {
	if (verbose)
	    check_error_msg( "Min_circle is nonempty \
			      but has no support points");
	return (false); }
@end

If the circle has one support point $p$, it must be equal to that
point, i.e. its center must be $p$ and its radius $0$.

@macro <check 1 support point> = @begin
    if ( ( circle().center()         != support_point(0)) ||
	 ( circle().squared_radius() != R::FT(0))         ) {
	if (verbose)
	    check_error_msg( "Min_circle differs from the circle spanned \
                              by its single support point");  
	return (false); }
@end

In case of two support points $p,q$, these points must form a diameter
of the circle. The support set $\{p,q\}$ is minimal if and only if
$p,q$ are distinct.

The diameter property is checked by comparing the midpoint of the
segment $\overline{pq}$ with the circle's center and the squared
distance of $p$ to the midpoint with the circle's squared radius.

@macro <check 2 support points> = @begin
    CGAL_Point_2<R> p( support_point(0)),
		    q( support_point(1));
    if ( p == q) {
	if (verbose)
	    check_error_msg( "Min_circle has two support points \
			      which are equal");
	return (false); }
    CGAL_Point_2<R> c( p + (q-p)/R::RT( 2));       // should be circle's center
    R::FT           sqr_r( (p-c)*(p-c));   // should be circle's squared radius
    if ( ( c != circle().center()) || ( sqr_r != circle().squared_radius())) {
	if (verbose)
            check_error_msg( "Min_circle does not have its two support \
			      points as a diameter");
	return (false); }
@end

If the number of support points is three (and they are distinct and
not collinear), the circle through them is unique, and must therefore
equal the @prg{min_circle}. It is the smallest one defined by the
three points if and only if the center of the circle lies inside or on
the boundary of the triangle defined by the three points. The support
set is minimal only if the center lies properly inside the triangle.

@macro <check 3 support points> = @begin
    CGAL_Point_2<R> p( support_point(0)),
		    q( support_point(1)),
		    r( support_point(2));
    if ( ( p == q) || ( q == r) || ( r == p)) {
	if (verbose)
	    check_error_msg( "Min_circle has three support points, \
			      two of which are equal");
	return (false); }
    CGAL_Circle_2<R> circ( p, q, r);
    if ( circ.orientation() == CGAL_COLLINEAR) {
	if (verbose)
	    check_error_msg( "Min_circle has three support points, \
			      which are collinear");
	return (false); }
    // circle() should equal circ up to its orientation
    if ( ( circle().center()         != circ.center())         ||
	 ( circle().squared_radius() != circ.squared_radius()) ) {
	if (verbose)
	    check_error_msg( "Min_circle is not the unique circle \
			      through its three support points");
	return (false); }
    CGAL_Triangle_2<R> t( p, q, r);
    if ( t.has_on_unbounded_side( circle().center())) {
	if (verbose)
	    check_error_msg( "Min_circle has three support points \
			      whose convex hull does not contain \
			      the circle's center");
	return (false); }
    if ( t.has_on_boundary( circle().center())) {
	if (verbose)
	    check_error_msg( "Min_circle has three support points, \
			      one of which is redundant");
	return (false); }
@end


@! ==========================================================================
@section{ File Organisation}
@! ==========================================================================

@file <Min_circle_2.h> = @begin
    @<file header>("2D Smallest Enclosing Circle",
		   "include/CGAL/Min_circle_2.h","Min_circle_2",
		   "Bernd Gärtner, Sven Schönherr (sven@@inf.fu-berlin.de)")

    #ifndef CGAL_MIN_CIRCLE_2_H
    #define CGAL_MIN_CIRCLE_2_H

    // Class declaration
    // =================
    template < class R >
    class CGAL_Min_circle_2;

    class CGAL_Bbox_2;

    // Class interface
    // ===============
    // includes
    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif
    #ifndef CGAL_CIRCLE_2_H
    #  include <CGAL/Circle_2.h>
    #endif
    #include <vector.h>

    @<Min_circle_2 interface>

    @<dividing line>

    // Class definition
    // ================
    // includes
    // --------
    #ifndef CGAL_UTILS_H
    #include <CGAL/utils.h>
    #endif
    #ifndef CGAL_TRIANGLE_2_H
    #include <CGAL/Triangle_2.h>
    #endif
    #include <algo.h>

    // check macros
    // ------------
    #ifdef CGAL_CHECK_ASSERTIONS
    #define CGAL_Min_circle_2_assertion(EX) \
     ((EX) ? ((void)0) : cgal_assertion_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_circle_2_assertion_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_assertion_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_circle_2_assertion(EX) ((void)0)
    #define CGAL_Min_circle_2_assertion_msg(EX,MSG) ((void)0)
    #endif // CGAL_CHECK_ASSERTIONS

    #ifdef CGAL_CHECK_PRECONDITIONS
    #define CGAL_Min_circle_2_precondition(EX) \
     ((EX) ? ((void)0) : cgal_precondition_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_circle_2_precondition_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_precondition_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_circle_2_precondition(EX) ((void)0)
    #define CGAL_Min_circle_2_precondition_msg(EX,MSG) ((void)0)
    #endif // CGAL_CHECK_PRECONDITIONS

    #ifdef CGAL_CHECK_POSTCONDITIONS
    #define CGAL_Min_circle_2_postcondition(EX) \
     ((EX) ? ((void)0) : cgal_postcondition_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_circle_2_postcondition_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_postcondition_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_circle_2_postcondition(EX) ((void)0)
    #define CGAL_Min_circle_2_postcondition_msg(EX,MSG) ((void)0)
    #endif // CGAL_CHECK_POSTCONDITIONS

    // workaround for new C++-style casts
    // ----------------------------------
    #if (__SUNPRO_CC)
    #define      static_cast(type,expr) (type)( expr)
    #define       const_cast(type,expr) (type)( expr)
    #define reinterpret_cast(type,expr) (type)( expr)
    #define     dynamic_cast(type,expr) (type)( expr)
    #else
    #define      static_cast(type,expr)      static_cast< type>( expr)
    #define       const_cast(type,expr)       const_cast< type>( expr)
    #define reinterpret_cast(type,expr) reinterpret_cast< type>( expr)
    #define     dynamic_cast(type,expr)     dynamic_cast< type>( expr)
    #endif // (__SUNPRO_CC)

    @<Min_circle_2 implementation>

    // specializations
    // ---------------
    #if ( defined( CGAL_INTEGER_H) && defined( CGAL_HOMOGENEOUS_H))
    #include <CGAL/Min_circle_2__Homogeneous_integer.h>
    #endif

    #endif // CGAL_MIN_CIRLCE_2_H

    @<end of file line>
@end

@i file_header.awlib
 

\begin{thebibliography}{Wel} 
\bibitem{Wel}
E.~Welzl.
\newblock Smallest enclosing disks (balls and ellipsoids).
\newblock In H.~Maurer, editor, {\em New Results and New Trends in Computer
  Science}, volume 555 of {\em Lecture Notes in Computer Science}, pages
  359--370. Springer-Verlag, 1991.
\end{thebibliography}
