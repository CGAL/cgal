@! ============================================================================
@! The CGAL Project
@! Implementation: 2D Smallest Enclosing Ellipse
@! ----------------------------------------------------------------------------
@! file  : Library/web/Min_ellipse_2.aw
@! author: Bernd Gärtner, Sven Schönherr (sven@inf.fu-berlin.de)
@! $Id$
@! ============================================================================
 
@documentclass[twoside]{article}
@usepackage{a4wide2}
@usepackage{amssymb}

@! LaTeX macros
\newcommand{\R}{\mathbb{R}}
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

\newcommand{\me}{\texttt{me}}

@article
@p maximum_input_line_length = 180
@p maximum_output_line_length = 180

@thickline

@t vskip 5 mm
@t title titlefont centre "CGAL: 2D Smallest Enclosing Ellipse"
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

We define a class template @prg{CGAL_Min_ellipse_2<R>}. An object of
this class represents the smallest (w.r.t. area) enclosing ellipse of a
finite point set $P$ in the plane, denoted by $\me(P)$. The template
parameter of the representation class $R$ (the domain of the point
coordinates) is any number type that fulfills the CGAL number type
requirements, but correct results are in this version only guaranteed
if the number type is an exact one (like LEDA's @prg{integer}). In
particular, using the number type @prg{double} might lead to wrong
results. A correct @prg{double} implementation is planned for the next
release.

The implementation is based on an algorithm by Welzl \cite{Wel}, which
we shortly describe now. $\me(P)$ is built up incrementally, adding
one point after another. Assume $\me(P)$ has been constructed, and we
would like to obtain $\me(P \cup \{p\})$, $p$ some new point. There
are two cases: if $p$ already lies inside $\me(P)$, then $\me(P \cup
\{p\}) = \me(P)$. Otherwise $p$ must lie on the boundary of $\me(P
\cup \{p\})$ (this is proved in \cite{Wel} and not hard to see), so we
need to compute $\me(P,\{p\})$, the smallest ellipse enclosing $P$ with
$p$ on the boundary. This is recursively done in the same manner. In
general, for point sets $P$,$B$, define $\me(P,B)$ as the smallest
ellipse enclosing $P$ that has the points of $B$ on the boundary (if
defined). Although the algorithm finally delivers an ellipse
$\me(P,\emptyset)$, it internally deals with ellipses that have a
possibly nonempty set $B$. Here is the pseudocode of Welzl's method.
To compute $\me(P)$, it is called with the pair $(P,\emptyset)$,
assuming that $P=\{p_1,\ldots,p_n\}$ is stored in a linked list.

\begin{pseudocode}{$\me(P,B)$:}
  $me := \me(\emptyset,B)$ \\
  \IF $|B| = 5$ \keyword{THEN} \RETURN $me$ \\
  \FOR $i := 1$ \TO $n$ \DO
      \IF $p_i \not\in me$ \THEN
          $me := \me(\{p_1,\ldots,p_{i-1}\}, B \cup \{p_i\})$ \\
          move $p_i$ to the front of $P$ \\
      \END
  \END
  \RETURN $me$ \\
\end{pseudocode}

Note the following: (a) $|B|$ is always bounded by 5, thus the
computation of $\me(\emptyset,B)$ is a constant-time operation. 
In our implementation, it
is done by the private member function @prg{compute_ellipse}. (b) One
can check that the method maintains the invariant `$\me(P,B)$ exists'.
This justifies termination if $|B| = 5$, because then $\me(P,B)$ must
be the unique ellipse with the points of $B$ on the boundary, and
$\me(P,B)$ exists if and only if this ellipse contains the points of
$P$. Thus, no subsequent in-ellipse tests are necessary anymore (for
details see \cite{Wel}). (c) points which are found to lie outside the
current ellipse $me$ are considered `important' and are moved to the
front of the linked list that stores $P$. This is crucial for the
method's efficiency.

It can also be advisable to bring $P$ into random order before
computation starts. There are `bad' insertion orders which cause the
method to be very slow -- random shuffling gives these orders a very
small probability.

The private member function @prg{me} directly realizes the pseudocode
above.

@! ============================================================================
@section{ Class \texttt{CGAL\_Min\_ellipse\_2}: Interface}
@! ============================================================================

The class interface looks as follows.

@macro <Min_ellipse_2 interface> = @begin
    template < class R >
    class CGAL_Min_ellipse_2 {

      public:
	@<public interface>

      private:
        @<private classes>
	@<private data members>
	@<private member functions declaration>
    };
@end   

@! ============================================================================
@subsection{ Public Interface}
@! ============================================================================

The public interface is described and documented in the CGAL Reference
Manual, so we do not comment on it here. The interface is very similar
to the one of \texttt{CGAL\_Min\_circle\_2}. The major difference is
that no method exists to convert a \texttt{CGAL\_Min\_ellipse\_2} to
a \texttt{CGAL\_Ellipse\_2}, because this type does not exist yet.
Thus, the usage of \texttt{CGAL\_Min\_ellipse\_2} is restricted; in this
version, it can merely be drawn. Also, a bounding box and a check function 
do not exist
in this version. To check, we would need to test the sign of a root of
a cubic equation; future releases of the LEDA-type \texttt{real} will
provide this feature, and we add a checker as soon as a stable version
is available. Thus, checking is done visually, by looking at the drawing.

@macro <public interface> = @begin
    // creation
    CGAL_Min_ellipse_2( );
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p);
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2);
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2,
                       const CGAL_Point_2<R>& p3);
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2,
                       const CGAL_Point_2<R>& p3,
		       const CGAL_Point_2<R>& p4);
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
                       const CGAL_Point_2<R>& p2,
                       const CGAL_Point_2<R>& p3,
		       const CGAL_Point_2<R>& p4,
	               const CGAL_Point_2<R>& p5);
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>* first,
                       const CGAL_Point_2<R>* last,
                       bool randomize = false);
    ~CGAL_Min_ellipse_2( );

    // access
    int  number_of_points        ( ) const;
    int  number_of_support_points( ) const;

    const CGAL_Point_2<R>&  point        ( int i) const;
    const CGAL_Point_2<R>&  support_point( int i) const;
    const CGAL_Point_2<R>&  operator []  ( int i) const;

    // bounding box not supported in this version
    // CGAL_Bbox_2  bbox( ) const;

    // updates
    void  insert ( const CGAL_Point_2<R>& p);
    void  reserve( int n);

    // predicates; note: in-ellipse tests are non-const functions here
    CGAL_Bounded_side  bounded_side( const CGAL_Point_2<R>& p) const; 
    bool  has_on_bounded_side      ( const CGAL_Point_2<R>& p) const;
    bool  has_on_boundary          ( const CGAL_Point_2<R>& p) const;
    bool  has_on_unbounded_side    ( const CGAL_Point_2<R>& p) const;

    bool  is_empty     ( ) const;
    bool  is_degenerate( ) const;
@end

@! ----------------------------------------------------------------------------
@subsection{ Private Interface}
@! ----------------------------------------------------------------------------

The private member functions @prg{compute_ellipse} and @prg{me}  
are in detail described in the implementation section below. 

@macro <private member functions declaration> = @begin
    // copying and assignment not allowed!
    CGAL_Min_ellipse_2( const CGAL_Min_ellipse_2<R>&);
    CGAL_Min_ellipse_2<R>& operator = ( const CGAL_Min_ellipse_2<R>&);

    // functions for actual computation
    void me( int n, int n_sp);
    void compute_ellipse( );
@end

Now let us focus on the private classes and data members, where the 
difference to the
class @prg{CGAL_Min_circle_2} is most apparant. The linked list
to store the points, however, is still the same.

The points of $P$ are internally stored as a linked list that allows
to bring points to the front of the list in constant time. The linked
list is realized as an STL-@prg{vector} whose entries are points with
successor and predecessor information, realized in the class
@prg{_Point}.

@macro <private classes> += @begin
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

A @prg{CGAL_Min_ellipse_2} object maintains an array
@prg{support_points} of at most five \emph{support points} (the
actual number is given by @prg{n_support_points}, which at the end of
the computation contains the indices of a minimal subset $S \subseteq
P$ with $\me(P)= \me(S)$. (Note: that subset is not necessarily
\emph{minimum}). During the computations, the set of support points
coincides with the set $B$ appearing in the pseudocode for $\me(P,B)$,
see the introduction above.

@macro <private data members> += @begin
    int               n_support_points;
    CGAL_Point_2<R>   support_points[ 5];
@end

Most important, a @prg{CGAL_Min_ellipse_2} object at any time
keeps the actual ellipse $\me(\emptyset,B)$. Unlike in the class
@prg{CGAL_Min_circle_2}, this is not a @prg{CGAL_Ellipse_2}, on
the one hand, because this type does not exist yet, on the other
hand because the representation of $\me(\emptyset,B)$ is different
for any number of support points. For $|B|=3,5$, we store an explicit
ellipse (although in different formats, see below), for $|B|=4$, the
representation is implicit, because an explicit one might contain
irrational coordinates. For $|B|\leq 2$, we have no representation
at all. In this case, the ellipse is degenerate, and the in-ellipse test
(which is the only routine we need the ellipse for) can be decided directly
from the set $B$. 

Although at any point in time, only one representation is
valid and will be accessed, all three representations (for $|B|=3,4,5$)
are simultaneously stored in a @prg{CGAL_Min_ellipse_2} object. Since
none of them is really space-consuming, this is not a problem. 
(Unions, which would be just the right concept in this situation, are not
allowed to be used with classes as members, even if the respective
constructors do absolutely nothing). Each representation has its own
class. 

@macro <private classes> += @begin
    @<Ellipse_3 class>
    @<Ellipse_4 class>
    @<Ellipse_5 class>
@end

Here are the actual representations. 

@macro <private data members> += @begin
    Ellipse_3 ellipse_3;
    Ellipse_4 ellipse_4;
    Ellipse_5 ellipse_5;
@end   
    

@! ----------------------------------------------------------------------------
@subsubsection{ Ellipse representation classes}
@! ----------------------------------------------------------------------------
\label{ellipse_rep}
We have three classes,
@prg{Ellipse_3}, @prg{Ellipse_4} and @prg{Ellipse_5}, being in charge
of sets $B$ with 3,4, or 5 points. All classes have a method @prg{set} to
compute some representation of $\me(\emptyset,B)$, suitable to do in-ellipse
tests `$p\in \me(\emptyset,B)$'. To perform these tests, each class has a
@prg{bounded_side} method.

@macro <Ellipse_3 class> = @begin
    class Ellipse_3 {
	public:
	    @<Ellipse_3 data members>

	    // default constructor
	    // -------------------
	    Ellipse_3 () {}

	    // set method
	    // ----------
	    void set (const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2,
	              const CGAL_Point_2<R>& p3)
	    {
	    	@<Ellipse_3::set body>
	    }

	    // bounded_side method
	    // -------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p) const
	    {
	    	@<Ellipse_3::bounded_side body>
	    }
    };
@end 

@macro <Ellipse_4 class> = @begin
    class Ellipse_4 {
	public:
	    @<Ellipse_4 data members>

	    // default constructor
	    // -------------------
	    Ellipse_4 () {}

	    // set method
	    // ----------
	    void set (const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2,
	              const CGAL_Point_2<R>& p3, const CGAL_Point_2<R>& p4)
	    {
	    	@<Ellipse_4::set body>
	    }

	    // bounded_side method, non-const
	    // ------------------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p)
	    {
	    	@<Ellipse_4::bounded_side body>
	    }
    };
@end 

The class @prg{Ellipse_5} does not have a method to compute the ellipse 
from its five support points. The reason is that by the time an 
@prg{Ellipse_5} object is set up, the ellipse through the five points 
has already been computed and stored in the @prg{Ellipse_4} representation, 
see implementation section below. Thus it suffices to `steal' this ellipse.
To this end, a reference to the @prg{Ellipse_4} object is passed to the
@prg{set} method. Because the ellipse in question is computed by the
@prg{bounded_side} method of the @prg{Ellipse_4} class, the latter method
is not a @prg{const} method, unlike all other @prg{bounded_side} methods.

@macro <Ellipse_5 class> = @begin
    class Ellipse_5 {
	public:
	    @<Ellipse_5 data members>

 	    // default constructor
	    // -------------------
	    Ellipse_5 () {}

	    // set method
	    // ----------
 	    void set (Ellipse_4& ellipse_4)
	    {
	    	@<Ellipse_5::set body>	
	    }    
	
	    // bounded_side method
	    // -------------------
	    CGAL_Bounded_side bounded_side (const CGAL_Point_2<R>& p) const
	    {
	 	@<Ellipse_5::bounded_side body>
	    }	   
    };
@end 	    

@! ----------------------------------------------------------------------------
@section{Class \texttt{CGAL\_Min\_ellipse\_2}: Implementation}
@! ----------------------------------------------------------------------------

The implementation consists of several parts, each of which is
described in the sequel. The actual work is hidden in the member
functions @prg{set} and @prg{bounded_side} of the ellipse
representation classes. The remaining functions -- in particular
the private member functions @prg{me} and @prg{compute_ellipse} -- are
implemented completely similar to the corresponding functions of the
class @prg{CGAL_Min_circle_2}. An exception are the predicates for
the in-ellipse tests. In @prg{CGAL_Min_circle_2}, the in-circle predicates
were directly mapped to the corresponding predicates over the 
@prg{CGAL_circle_2} object stored in the class. In our case,
the @prg{bounded_side} method is mapped to the corresponding one
of the respective ellipse representation class (and evaluated directly
if the current ellipse $me$ is defined by less than three points); the predicates (@prg{has_on_bounded_side}, \ldots) are then implemented by 
directly referring to the @prg{bounded_side} method.  

@macro <Min_ellipse_2 implementation> = @begin
    @<access operations and predicates>
    @<private member functions>
    @<constructors>
    @<update operations>
@end

@! ----------------------------------------------------------------------------
@subsection{ Constructors}
@! ----------------------------------------------------------------------------

In this version, @prg{CGAL_Min_ellipse_2} has seven different
construction methods, where the most important one builds the
@prg{CGAL_Min_ellipse_2} $\me(P)$ from a point set $P$, stored in a
piece of memory delimited by two pointers (As soon as member templates
are available, the pointers will become iterators). The constructor
copies the points into the internal array @prg{points}, performs a
random reordering if the corresponding flag requires that, then builds
up the linked list over the point set, and finally calls the private
method $\me$ to compute $\me(P)=\me(P,\emptyset)$.

@macro <constructors> += @begin
    // constructors
    // ------------
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>* first,
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

	// compute me
	me( points.size(), 0);
    }
@end

The remaining constructors are actually specializations of the
previous one, building a @prg{CGAL_Min_ellipse_2} for up to five
points. The idea is the following: recall that for any point set $P$
there exists $S \subseteq P$, $|S| \leq 5$ with $\me(S) = \me(P)$ (in
fact, such a set $S$ is determined during a @prg{CGAL_Min_ellipse_2}
computation). Once $S$ has been computed (or given otherwise),
$\me(P)$ can easily be reconstructed from $S$ in constant time. To
make this reconstruction more convenient, specialized construction 
methods are available for point sets $P$ of size at most 5. For $|P| = 0$, 
we get the default constructor, building the empty @prg{CGAL_Min_ellipse_2}
$\me(\emptyset)$.

@macro <constructors> += @begin

    // default constructor
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( )
      : i_first_point( -1),
	n_support_points( 0)
    { }

    // constructor for one point
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p)
      : points( 1, _Point( p)),
	i_first_point( 0),
	n_support_points( 1)
    {
	support_points[ 0] = p;
    }

    // constructor for two points
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1, const CGAL_Point_2<R>& p2)
    {
	// store points
	points.reserve( 2);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0, -1));
	i_first_point = 0;

	// compute me
	me( 2, 0);
    }

    // constructor for three points
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
		        const CGAL_Point_2<R>& p2,
		        const CGAL_Point_2<R>& p3)
    {
	// store points
	points.reserve( 3);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0,  2));
	points.push_back( _Point( p3,  1, -1));
	i_first_point = 0;

	// compute me
	me( 3, 0);
    }

    // constructor for four points
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
		        const CGAL_Point_2<R>& p2,
		        const CGAL_Point_2<R>& p3,
			const CGAL_Point_2<R>& p4)
    {
	// store points
	points.reserve( 4);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0,  2));
	points.push_back( _Point( p3,  1,  3));
	points.push_back( _Point( p3,  2, -1));
	i_first_point = 0;

	// compute me
	me( 4, 0);
    }

    // constructor for five points
    template < class R >
    CGAL_Min_ellipse_2<R>::
    CGAL_Min_ellipse_2( const CGAL_Point_2<R>& p1,
		        const CGAL_Point_2<R>& p2,
		        const CGAL_Point_2<R>& p3,
			const CGAL_Point_2<R>& p4,
			const CGAL_Point_2<R>& p5)
    {
	// store points
	points.reserve( 5);
	points.push_back( _Point( p1, -1,  1));
	points.push_back( _Point( p2,  0,  2));
	points.push_back( _Point( p3,  1,  3));
	points.push_back( _Point( p3,  2,  4));
	points.push_back( _Point( p3,  3, -1));
	i_first_point = 0;

	// compute me
	me( 5, 0);
    }
@end

Finally, we have a (default) destructor.

@macro <constructors> += @begin
    // destructor
    template < class R >
    CGAL_Min_ellipse_2<R>::
    ~CGAL_Min_ellipse_2( )
    { }
@end

@! ----------------------------------------------------------------------------
@subsection{ Update operations}
@! ----------------------------------------------------------------------------

There is another way to build up $\me(P)$, other than by supplying the
point set $P$ at once. Namely, $\me(P)$ can be built up incrementally,
adding one point after another. If you look at the pseudocode in the
introduction, this comes quite naturally. The method @prg{insert},
applied with point $p$ to a @prg{CGAL_Min_ellipse_2} object
representing $\me(P)$, computes $\me(P \cup \{p\})$, where work has to
be done only if $p$ lies outside $\me(P)$.  In this case, $\me(P \cup
\{p\}) = \me(P,\{p\})$ holds, so the method $\me$ is called with
support set $\{p\}$. After the insertion has been performed, $p$ is
moved to the front of the point list, just like in the pseudocode and
the `main' constructor above.

@macro <update operations> += @begin
    // update operations
    // -----------------
    template < class R >
    void
    CGAL_Min_ellipse_2<R>::
    insert( const CGAL_Point_2<R>& p)
    {
	// store point
	int old_n = number_of_points();
	points.push_back( _Point( p));  // NOTE: p is not linked with list yet!

	// p not in current ellipse?
	if ( has_on_unbounded_side( p)) {

	    // p new support point
	    support_points[ 0] = p;

	    // recompute mc
	    me( old_n, 1); }

	// make p the first point in list
	if ( old_n > 0)                                  // old list not empty?
	    points[ i_first_point].pred = old_n;
	points[ old_n].succ = i_first_point;
	i_first_point = old_n;
    }
@end

The operation @prg{reserve} does nothing but tell the
@prg{CGAL_Min_ellipse_2} that some number $n$ of points might
eventually be inserted, allowing the object to allocate storage for
them at once. Inserting the points without doing this might lead to
overhead caused by moving the array @prg{points} around in memory
while it grows.

@macro <update operations> += @begin

    template < class R > inline
    void
    CGAL_Min_ellipse_2<R>::
    reserve( int n)
    {
	points.reserve( n);
    }
@end

@! ----------------------------------------------------------------------------
@subsection{ Access Operations and Predicates}
@! ----------------------------------------------------------------------------

These operations are used to retrieve information about the current
status of the @prg{CGAL_Min_ellipse_2} object. We first have 
operations to get information about points or support points and 
their numbers.

@macro <access operations and predicates> += @begin
    // 'number_of_' operations
    // -----------------------
    template < class R > inline
    int
    CGAL_Min_ellipse_2<R>::
    number_of_points( ) const
    {
	return( points.size());
    }

    template < class R > inline
    int
    CGAL_Min_ellipse_2<R>::
    number_of_support_points( ) const
    {
	return( n_support_points);
    }

    // 'is_' predicates
    // ----------------
    template < class R > inline
    bool
    CGAL_Min_ellipse_2<R>::
    is_empty( ) const
    {
	return( number_of_support_points() == 0);
    }

    template < class R > inline
    bool
    CGAL_Min_ellipse_2<R>::
    is_degenerate( ) const
    {
	return( number_of_support_points() < 3);
    }

    // access operations
    // -----------------
    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_ellipse_2<R>::
    point( int i) const
    {
	CGAL_Min_ellipse_2_precondition( (i >= 0) && (i < number_of_points()));
	return( points[ i].point);
    }

    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_ellipse_2<R>::
    support_point( int i) const
    {
	CGAL_Min_ellipse_2_precondition( (i >= 0) &&
					(i <  number_of_support_points()));
	return( support_points[ i]);
    }

    template < class R > inline
    const CGAL_Point_2<R>&
    CGAL_Min_ellipse_2<R>::
    operator [] ( int i) const
    {
	return( point( i));
    }
@end

Next we have predicates testing the relative position of a point
w.r.t. the ellipse. They rely on the @prg{bounded_side} method
which for less than 3 support points is evaluated directly, while
for at least three support points, the corresponding methods
of the ellipse representation classes are used. A degenerate ellipse 
defined by less than two support points has no bounded side. In
this case, the ellipse itself is the boundary, everything else
belongs to the unbounded side. Note that the @prg{bounded_side} method
is a @prg{const} method although the underlying method of the class
@prg{Ellipse_4} is not. The difference is that, `from the outside', the
@prg{CGAL_Min_ellipse_2} object appears unchanged after calling
@prg{bounded_side} (logical constness), while internally, 
we explicitly make use of the changes to an @prg{Ellipse_4} object
caused by a call to @prg{bounded_side}, see the discussion in 
subsection \ref{ellipse_rep}.

@macro <access operations and predicates> += @begin
    // in-ellipse predicates
    // -------------------
    template < class R > inline
    CGAL_Bounded_side
    CGAL_Min_ellipse_2<R>::
    bounded_side( const CGAL_Point_2<R>& p) const
    {
	switch (n_support_points) {
	  case 5: 
	    return ellipse_5.bounded_side (p);
	    break;
	  case 4:
	    // cast away ellipse_4's constness 
	    return ( const_cast( Ellipse_4&, ellipse_4)).bounded_side (p);
	    break;
	  case 3:
	    return ellipse_3.bounded_side (p);
	    break;
	  case 2: // Min_ellipse is line segment spanned by support points
	    if (CGAL_between (support_points [0], p, support_points [1]))
		return CGAL_ON_BOUNDARY;
	    else
		return CGAL_ON_UNBOUNDED_SIDE;
	    break;
	  case 1: // Min_ellipse is equal to its support point
	    if (p==support_points [0])
		return CGAL_ON_BOUNDARY;
	    else
		return CGAL_ON_UNBOUNDED_SIDE;
	    break;
	  case 0: // Min_ellipse is empty
	    return CGAL_ON_UNBOUNDED_SIDE;
	  default:
 	    // we should never get here ...
	    CGAL_Min_ellipse_2_assertion 
		( (number_of_support_points() >= 0) &&
		  (number_of_support_points() <= 5) );
	    // ... and we definitely never get here
	    return CGAL_ON_BOUNDARY; 	
	    break; 
	}
    }
 
    template < class R > inline
    bool
    CGAL_Min_ellipse_2<R>::
    has_on_bounded_side( const CGAL_Point_2<R>& p) const
    {
	return( bounded_side( p) == CGAL_ON_BOUNDED_SIDE);
    }

    template < class R > inline
    bool
    CGAL_Min_ellipse_2<R>::
    has_on_boundary( const CGAL_Point_2<R>& p) const
    {
	return( bounded_side( p) == CGAL_ON_BOUNDARY);
    }

    template < class R > inline
    bool
    CGAL_Min_ellipse_2<R>::
    has_on_unbounded_side( const CGAL_Point_2<R>& p) const
    {
	return( bounded_side( p) == CGAL_ON_UNBOUNDED_SIDE);
    }
@end

@! ----------------------------------------------------------------------------
@subsection{ Private Member Function \texttt{compute\_ellipse}}
@! ----------------------------------------------------------------------------

This is the method for computing $me(\emptyset,B)$ the set $B$ given
by the first @prg{n_support_points} indices in the array
@prg{i_support_points}. It is realized by a case analysis,
noting that $|B| \leq 5$. If $|B|\leq 2$, nothing is done,
in the other cases the @prg{set} methods for @prg{Ellipse_3}, 
@prg{Ellipse_4}, or @prg{Ellipse_5} are called. 

@macro <private member functions> += @begin
    template < class R >
    void
    CGAL_Min_ellipse_2<R>::
    compute_ellipse( )
    {
	switch ( n_support_points) {
	  case 5:
	    // `steal' ellipse from Ellipse_4 represemtation
	    ellipse_5.set (ellipse_4);	    
            break;
	  case 4:
	    ellipse_4.set (support_points[0], support_points[1],
			   support_points[2], support_points[3]);
	    break;
	  case 3:
	    ellipse_3.set (support_points[0], support_points[1],
			   support_points[2]); 
	    break;
	  default:
	    // do not compute any representation of a degenerate ellipse
	    break; 
	}
    }
@end


@! ----------------------------------------------------------------------------
@subsection{ Private Member Function \texttt{me}}
@! ----------------------------------------------------------------------------

This function computes the general ellipse $\me(P,B)$, where $P$
contains the first @prg{n} points stored in the array @prg{points} and
$B$ is given by the first @prg{n_support_points} indices in the array
@prg{i_support_points}. The function is directly modelled after the
pseudocode in the introduction.

@macro <private member functions> += @begin

    template < class R >
    void
    CGAL_Min_ellipse_2<R>::
    me( int n, int n_sp)
    {
	// compute ellipse through support points
	n_support_points = n_sp;
	compute_ellipse( );
	if ( n_sp == 5) return;

	// test first n points
	int index = i_first_point, succ;
	for ( int i = 0; i < n; ++i) {
	    _Point& p = points[ index];
	    succ = p.succ;

	    // p not in current ellipse?
	    if ( has_on_unbounded_side( p.point)) {

		// recursive call with p as additional support point
		support_points[ n_sp] = p.point;
		me( i, n_sp+1);

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
@subsection{ Class \texttt{Conic<R>}}
@! ----------------------------------------------------------------------------

In their implementations, the ellipse representation classes rely on a
concept more general than ellipses, namely on {\em conics}. Ellipses are
special conics, in addition there are {\em hyperbolas} and {\em parabolas}. 
Conics play a particularly important role in the @prg{Ellipse_4} class.

A conic (in linear form) is the set of points $q=(x,y)^T$ satisfying
\begin{equation}
{\cal C}(q) := rx^2 + sy^2 + 2txy + 2ux + 2vy + w = 0. 
\end{equation}
$r,s,t,u,v,w$ being real parameters. 

Alternatively, ${\cal C}$ can be written as 
$${\cal C}(q) = (x,y)\left(\begin{array}{cc}r & t \\ t & s \end{array}
\right)(x,y)^T + 2(u,v)(x,y)^T + w,$$
and the determinant of the matrix 
$$M=\left(\begin{array}{cc}r & t \\ t & s \end{array}
\right)$$ tells the type of ${\cal C}$. We also denote $\det(M)$ 
by $\det({\cal C})$. 

If $\det({\cal C})>0$, ${\cal C}$ is an ellipse, for $\det({\cal C})<0$ 
we get a hyperbola, and if $\det({\cal C})=0$, ${\cal C}$ is a parabola. 

The {\em linear combination} $\lambda{\cal C}_1 +\mu{\cal C}_2$ of two
conics ${\cal C}_1$ and ${\cal C}_2$ is defined as the conic ${\cal C}$
with $${\cal C}(q) := \lambda{\cal C}_1(q) +\mu{\cal C}_2(q)=0.$$

We introduce a conic class that supports 
\begin{itemize}
\item definition of a conic directly from values $r,s,t,u,v,w$,
\item computation of the linear combination of two conics, 
\item conic evaluation (i.e. computation of the value 
${\cal C}(p)$ for given $p\in\R^2$), 
\item computation of $\det({\cal C})$,
\item normalization, i.e. scaling of the conic by $-1$ -- if 
necessary  -- to guarantee $r\geq 0$. 
If $\det(M) > 0$, normalization guarantees that $M$ 
is positive definite afterwards (which is important
for doing in-ellipse tests, see subsequent subsections).
\end{itemize}

@macro <conic class> = @begin
    template <class R>
    class Conic {
	public:
	    R::FT  r;
	    R::FT  s;
	    R::FT  t;
	    R::FT  u;
	    R::FT  v;
	    R::FT  w;

	    // default constructor
	    // -------------------
	    Conic () {} 
	    
	    // member functions
	    // ----------------
	    // definition of conic by values r, s, t, u, v, w
	    void set (R::FT _r, R::FT _s, R::FT _t, R::FT _u,
		      R::FT _v, R::FT _w) 
	    {
		r=_r; s=_s; t=_t; u=_u; v=_v; w=_w;
	    }

	    // definition of conic by linear combination of conics
	    void set (R::FT lambda, Conic& conic1, 
		      R::FT mu, Conic& conic2)
	    {
		r = lambda*conic1.r + mu*conic2.r;
		s = lambda*conic1.s + mu*conic2.s;
		t = lambda*conic1.t + mu*conic2.t;
		u = lambda*conic1.u + mu*conic2.u;
		v = lambda*conic1.v + mu*conic2.v;
		w = lambda*conic1.w + mu*conic2.w;
	    }

	    // evaluation of conic at point p
	    R::FT eval (const CGAL_Point_2<R>& p) const
	    {
		R::FT x = p.x(), y = p.y();
		return r*x*x + s*y*y + R::FT(2)*(t*x*y + u*x + v*y) + w;
	    }

	    // determinant of matrix M = [[r,t],[t,s]]
	    R::FT det () const
	    {
		return r*s-t*t;
	    }

	    // scaling of conic such that r becomes nonnegative
	    void normalize ()
	    {
	   	if (r<R::FT(0)) {
		   r=-r; s=-s; t=-t; u=-u; v=-v; w=-w;
	  	}
	    } 
    };
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_3 Class}
@! ----------------------------------------------------------------------------

This class stores $\me(\emptyset,B)$ directly as an ellipse 
in {\em center form}, i.e. in the form $$(q-c)^TM(q-c)-z = 0,$$
$c\in \R^2$ the center. 
For support points $p_1,p_2,p_3$, the ellipse is given by
$$c=\frac{1}{3}\sum_{i=1}^3p_i,\quad M^{-1} = 
\frac{2}{3}\sum_{i=1}^3(p_i-c)(p_i-c)^T, \quad z = 1,$$
see \cite{GS}.

Let $$M^{-1} = \frac{2}{3}
\left(\begin{array}{cc} m_{11} & m_{12} \\ m_{12} & m_{22}
\end{array}\right).$$ Then 
$$M = \frac{3}{2}\left(\begin{array}{cc} m_{11} & m_{12} \\ m_{12} & m_{22}
\end{array}\right)^{-1} = \frac{3}{2(m_{11}m_{22}-m_{12}^2)}
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right).$$

This implies 
\begin{eqnarray*}
(q-c)^TM(q-c)-z 
&=& \frac{3}{2(m_{11}m_{22}-m_{12}^2)}(q-c)^T
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right)(q-c)-1,
\end{eqnarray*}
which is equal in sign to 
\begin{equation}
\label{ellipse_3_test}
(q-c)^T
\left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right)(q-c)-\frac{2}{3}(m_{11}m_{22}-m_{12}^2).
\end{equation}
(For this, we need
the fact that $M$ is positive definite, hence $\det(M),\det(M')$ and
$\det(M'^{-1})$ are positive.)
If $z'$ denotes $2(m_{11}m_{22}-m_{12}^2)/3$, an
@prg{Ellipse_3} object stores the values 
$c, z', m_{11},2m_{12},m_{22}$.
This enables the subsequent
@prg{bounded_side} routine to perform the in-ellipse test by 
just evaluating the sign of (\ref{ellipse_3_test}).

@macro <Ellipse_3 data members> = @begin
    R::FT  cx;
    R::FT  cy;
    R::FT  z;
    R::FT  m11;
    R::FT  m12;
    R::FT  m22;
@end
  
@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{set} method}
@! ----------------------------------------------------------------------------

Using the formulas above, the private data members are
straightforward to set.

@macro <Ellipse_3::set body> = @begin
    R::FT  x1 = p1.x();
    R::FT  y1 = p1.y();
    R::FT  x2 = p2.x(); 
    R::FT  y2 = p2.y();
    R::FT  x3 = p3.x(); 
    R::FT  y3 = p3.y();

    // center c
    cx = (x1+x2+x3)/R::FT(3);
    cy = (y1+y2+y3)/R::FT(3);

    R::FT  x1_cx = x1-cx;
    R::FT  x2_cx = x2-cx;
    R::FT  x3_cx = x3-cx;
    R::FT  y1_cy = y1-cy;
    R::FT  y2_cy = y2-cy;
    R::FT  y3_cy = y3-cy;

    // matrix M
    m11 = (x1_cx)*(x1_cx) + (x2_cx)*(x2_cx) + (x3_cx)*(x3_cx);
    m12 = (x1_cx)*(y1_cy) + (x2_cx)*(y2_cy) + (x3_cx)*(y3_cy);
    m22 = (y1_cy)*(y1_cy) + (y2_cy)*(y2_cy) + (y3_cy)*(y3_cy);

    // z'
    z = R::FT(2)*(m11*m22-m12*m12)/R::FT(3);

    // assert positive definiteness
    CGAL_Min_ellipse_2_assertion( (z > R::FT(0)) && (m11 > R::FT(0)) );
@end

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{bounded\_side} method}
@! ----------------------------------------------------------------------------

This is quite obvious now by (\ref{ellipse_3_test}). 
Since $M$ is positive definite, we have
$p=(x,y)^T$ in/on/outside the ellipse iff
$$(x-c_x,y-c_y)^T \left(\begin{array}{cc} m_{22} & -m_{12} \\ -m_{12} & m_{11}
\end{array}\right) (x-c_x,y-c_y) - z' \left\{\begin{array}{c}<\\=\\>
\end{array}\right. 0.$$

@macro <Ellipse_3::bounded_side body> = @begin
    R::FT  x = p.x();
    R::FT  y = p.y();
    R::FT  x_cx = x-cx;
    R::FT  y_cy = y-cy;
    R::FT  discr = m22*(x_cx)*(x_cx) + m11*(y_cy)*(y_cy)
		   - R::FT(2)*m12*(x_cx)*(y_cy) - z;	    
    return static_cast( CGAL_Bounded_side, CGAL_sign( discr));
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_4 Class}
@! ----------------------------------------------------------------------------

This is by far the most complicated one among the three representation 
classes. The @prg{set} method computes some implicit representation of 
the ellipse, derived from the four support points. The @prg{bounded_side}
predicate works over this implicit representation and is not
straightforward anymore.  

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{set} method}
@! ----------------------------------------------------------------------------

Subsequent computations assume that $p_1,\ldots p_4$ are in convex position 
and enumerated in clockwise or counterclockwise order (we would like 
to work without the latter assumption but don't know how). While convex
position is guaranteed by the ambient algorithm, a (counter)clockwise
odering can be achieved in a preprocessing step, as follows.

We first check whether $p_1$ and $p_3$ lie on different sides of the
line through $p_2,p_4$. If so, the segment $\overline{p_2p_4}$ is a
diagonal of the convex quadrilateral defined by the four points, so
the given order is ok. Otherwise, we check whether $\overline{p_3p_2}$
defines a diagonal. From the previous test we already know the relative
position of $p_4$ w.r.t. this segment, so it remains to check the position
of $p_1$. If $p_1,p_4$ lie on different sides, $\overline{p_3p_2}$ is a
diagonal. To achieve correct enumeration of the points, we just need to
swap points $p_1$ and $p_2$. If $p_1,p_4$ lie on the same side of
$\overline{p_3p_2}$, then $\overline{p_1p_2}$ must be a diagonal. For
the correct order, we swap $p_2$ and $p_3$.

@macro <Ellipse_4::set body> = @begin
    // copy p_i's to allow swapping
    CGAL_Point_2<R>  q1 = p1;
    CGAL_Point_2<R>  q2 = p2;
    CGAL_Point_2<R>  q3 = p3;
    CGAL_Point_2<R>  q4 = p4;

    // preprocessing: bring q_i's into (counter)clockwise order
    CGAL_Orientation side1_24 = CGAL_orientation (q2, q4, q1),
	     	     side3_24 = CGAL_orientation (q2, q4, q3);
    if (side1_24 == side3_24) { 
	CGAL_Orientation side1_32 = CGAL_orientation (q3, q2, q1);
	if (side1_32 != side3_24)  // side3_24 == side4_32
	    swap (q1, q2); 	// swap from STL algobase
	else 
	    swap (q2, q3);
    }

    // assert correct ordering
    CGAL_Min_ellipse_2_assertion(
        ( CGAL_orientation( q1, q2, q3) == CGAL_orientation( q2, q3, q4)) &&
	( CGAL_orientation( q2, q3, q4) == CGAL_orientation( q3, q4, q1)) &&
	( CGAL_orientation( q3, q4, q1) == CGAL_orientation( q4, q1, q2)) &&
 	( CGAL_orientation( q1, q2, q3) != CGAL_COLLINEAR               ) );

    // do the actual computations
    @<Ellipse_4 setup>
@end 

In all terminology that follows, we refer to \cite{GS}. 
The implicit representation of
$\me(\emptyset,\{p_1,p_2,p_3,p_4\})$ first of all consists of two
special conics ${\cal C}_1,{\cal C}_2$ through the for points. 
Assuming that $q_i=(x_i,y_i)$, these conics are given by parameters 
$r_1,s_1,t_1,u_1,v_1,w_1$ and $r_2,s_2,t_2,u_2,v_2,w_2$, where
\begin{eqnarray*}
  r_1 & = & (y_1-y_2)(y_3-y_4), \\
  s_1 & = & (x_1-x_2)(x_3-x_4), \\
  2t_1 & = & -(x_1-x_2)(y_3-y_4)-(y_1-y_2)(x_3-x_4), \\
  2u_1 & = & (x_1y_2-x_2y_1)(y_3-y_4) + (x_3y_4-x_4y_3)(y_1-y_2), \\
  2v_1 & = & (x_1y_2-x_2y_1)(x_4-x_3) + (x_3y_4-x_4y_3)(x_2-x_1), \\
  w_1 & = & (x_1y_2-x_2y_1)(x_3y_4-x_4y_3)
\end{eqnarray*}
and
\begin{eqnarray*}
  r_2 & = & (y_2-y_3)(y_4-y_1), \\
  s_2 & = & (x_2-x_3)(x_4-x_1), \\
  2t_2 & = & -(x_2-x_3)(y_4-y_1)-(y_2-y_3)(x_4-x_1), \\
  2u_2 & = & (x_2y_3 - x_3y_2)(y_4-y_1) + (x_4y_1-x_1y_4)(y_2-y_3), \\
  2v_2 & = & (x_2y_3 - x_3y_2)(x_1-x_4) + (x_4y_1-x_1y_4)(x_3-x_2), \\
  w_2 & = & (x_2y_3 - x_3y_2)(x_4y_1-x_1y_4).
\end{eqnarray*}

${\cal C}_1$ and ${\cal C}_2$ are stored as data members in an
@prg{Ellipse_4} object.

@macro <Ellipse_4 data members> += @begin
    Conic<R>  conic1;
    Conic<R>  conic2;
@end

@macro <Ellipse_4 setup> += @begin
    // point coordinates
    R::FT  x1 = q1.x();
    R::FT  y1 = q1.y();
    R::FT  x2 = q2.x();
    R::FT  y2 = q2.y();
    R::FT  x3 = q3.x();
    R::FT  y3 = q3.y();
    R::FT  x4 = q4.x();
    R::FT  y4 = q4.y();

    // auxiliary values
    R::FT  y1_y2 = y1-y2;
    R::FT  y3_y4 = y3-y4;
    R::FT  x1_x2 = x1-x2;
    R::FT  x3_x4 = x3-x4;
    R::FT  y2_y3 = y2-y3;
    R::FT  y4_y1 = y4-y1;
    R::FT  x2_x3 = x2-x3;
    R::FT  x4_x1 = x4-x1;
    R::FT  x1y2_x2y1 = x1*y2 - x2*y1;
    R::FT  x3y4_x4y3 = x3*y4 - x4*y3;
    R::FT  x2y3_x3y2 = x2*y3 - x3*y2;
    R::FT  x4y1_x1y4 = x4*y1 - x1*y4;

    // define conics
    conic1.set (  y1_y2 * y3_y4, 
		  x1_x2 * x3_x4, 
	        -(x1_x2 * y3_y4 + y1_y2 * x3_x4)/R::FT(2),
		 (x1y2_x2y1 * y3_y4 + x3y4_x4y3 * y1_y2)/R::FT(2),
		-(x1y2_x2y1 * x3_x4 + x3y4_x4y3 * x1_x2)/R::FT(2),
    		  x1y2_x2y1 * x3y4_x4y3 ); 
    conic2.set (  y2_y3 * y4_y1,
		  x2_x3 * x4_x1,
		-(x2_x3 * y4_y1 + y2_y3 * x4_x1)/R::FT(2),
		 (x2y3_x3y2 * y4_y1 + x4y1_x1y4 * y2_y3)/R::FT(2),
		-(x2y3_x3y2 * x4_x1 + x4y1_x1y4 * x2_x3)/R::FT(2),
		  x2y3_x3y2 * x4y1_x1y4 );
@end

In addition, we store values that are needed to perform the in-ellipse
test w.r.t. $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, given a fifth point $p$.
Some of these values are independent of $p$, they are computed here. 

In \cite{GS} it is shown that the unique conic ${\cal C}_0$ through $\{p_1,p_2,p_3,p_4,p\}$ is given by
$$\lambda_0 {\cal C}_1 + \mu_0 {\cal C}_2, $$ where
$\lambda_0,\mu_0$ are given by
$\lambda_0 = {\cal C}_2(p), \mu_0 = -{\cal C}_1(p).$

Furthermore, to do the in-ellipse test with $p$ over $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, the sign of a certain 
derivative needs to be computed. For this, we define
\begin{equation}
\left(\begin{array}{c}\lambda(\tau)\\ \mu(\tau)\end{array}\right) = 
\left(\begin{array}{c}\lambda_0\\ \mu_0\end{array}\right) + {\tau}
\left(\begin{array}{r} -r2 \\ r1 \end{array}\right)
\end{equation}
and then consider the conic
$$
{\cal C}(\tau) := \lambda(\tau) {\cal C}_1 + \mu(\tau) {\cal C}_2,
$$
given by parameters $r(\tau),s(\tau), t(\tau), u(\tau),v(\tau),w(\tau)$.
To evaluate the derivative in question, we need the values
$$r(0)=r_0, s(0)=s_0, t(0)=t_0, u(0)=u_0, v(0)=v_0, w(0)=w_0$$
and 
$$r'(0), s'(0), t'(0), u'(0), v'(0), w'(0),$$
primed functions being derivatives w.r.t. $\tau$
(see subsection \ref{ellipse_4::bounded_side}).
While the former values depend on the fifth point $p$, the primed values 
do not, so they can be computed and stored already here. We get
\begin{eqnarray*}
r'(0) &=& r_1 r_2 - r_2 r_1 = 0, \\
s'(0) &=& r_1 s_2 - r_2 s_1, \\
t'(0) &=& r_1 t_2 - r_2 t_1, \\
u'(0) &=& r_1 u_2 - r_2 u_1, \\
v'(0) &=& r_1 v_2 - r_2 v_1, \\
w'(0) &=& r_1 w_2 - r_2 w_1.
\end{eqnarray*}

@macro <Ellipse_4 data members> += @begin
    R::FT  ds;	  
    R::FT  dt;
    R::FT  du;
    R::FT  dv;	
    R::FT  dw;
@end

@macro <Ellipse_4 setup> += @begin
    ds = conic1.r*conic2.s - conic2.r*conic1.s;
    dt = conic1.r*conic2.t - conic2.r*conic1.t;
    du = conic1.r*conic2.u - conic2.r*conic1.u;
    dv = conic1.r*conic2.v - conic2.r*conic1.v;
    dw = conic1.r*conic2.w - conic2.r*conic1.w;
@end

Finally, we need some ellipse through the four points. Such an ellipse
can be obtained as $$E = \lambda {\cal C}_1 + \mu {\cal C}_2,$$
with \begin{eqnarray*}
\lambda &=& 2\gamma-\beta, \\
\mu &=& 2\alpha-\beta,
\end{eqnarray*}
$\alpha = \det({\cal C}_1), \gamma = \det({\cal C}_2), 
\beta = r_1s_2+r_2s_1-2t_1t_2$.
This ellipse is stored as another conic in @prg{Ellipse_4}.

@macro <Ellipse_4 data members> += @begin
    Conic<R> ellipse;
@end

@macro <Ellipse_4 setup> += @begin
    R::FT  beta   = conic1.r*conic2.s + conic2.r*conic1.s 
		    - R::FT(2)*conic1.t*conic2.t;
    R::FT  lambda = R::FT(2)*conic2.det()-beta;
    R::FT  mu     = R::FT(2)*conic1.det()-beta;
    ellipse.set (lambda, conic1, mu, conic2);
    ellipse.normalize();
    CGAL_Min_ellipse_2_assertion( ellipse.det() > R::FT(0));
@end

@! ----------------------------------------------------------------------------
@subsubsection{ The \texttt{bounded\_side} method}
@! ----------------------------------------------------------------------------
\label{ellipse_4::bounded_side}
Given a fifth point $p$, the first action is to compute the unique conic
${\cal C}_0$ through $p_1,p_2,p_3,p_4,p$, as already described in the 
previous subsubsection. Since ${\cal C}_0$ is possibly being
recycled for usage in the @prg{Ellipse_5} class later, we need to 
store it. Thus, ${\cal C}_0$ is the only data member of @prg{Ellipse_4}
that is not independent of the query point $p$. (For this reason, the
@prg{bounded_side} method is not a @prg{const} method; moreover, the change 
to the @prg{Ellipse_4} object caused by calling the method does not fall 
under logical constness, because the change becomes quite visible when 
${\cal C}_0$ is later passed to the @prg{set} method of @prg{Ellipse_5}.)

@macro <Ellipse_4 data members> += @begin
    Conic<R> conic0;
@end

Depending on the type of ${\cal C}_0$ (given by its determinant), 
different actions are taken.

@macro <Ellipse_4::bounded_side body> = @begin
    R::FT  lambda_0 =  conic2.eval(p);
    R::FT  mu_0     = -conic1.eval(p);
    conic0.set (lambda_0, conic1, mu_0, conic2);
    R::FT d = conic0.det();
    if (d > R::FT(0)) {
	@<handle ellipse case>
    }
    else {
	@<handle hyperbola/parabola case>
    }
@end

In the hyperbola/parabola case, we test $p$ against the ellipse
$E$ precomputed by the @prg{Ellipse_4} @prg{set} method, and just return
the result (which is the same for all ellipses, in particular for
$\me(\emptyset,B)$).

@macro <handle hyperbola/parabola case> = @begin
    R::FT discr = ellipse.eval(p); 
    return static_cast (CGAL_Bounded_side, CGAL_sign (discr));
@end

To handle the ellipse case, we need to evaluate two derivatives and 
combine the results. Let us review these two tasks.

(i) compute $$\rho := \frac{\partial}{\partial\tau}E^{\tau}(p),$$
where 
\begin{eqnarray*}
E^{\tau}(p) &=&  
(\lambda_0-\tau r_2){\cal C}_1(p) + (\mu_0+\tau r_1){\cal C}_2(p) \\
&=& \tau (r_1{\cal C}_2(p) - r_2{\cal C}_1(p)).
\end{eqnarray*}
This means, $$\rho = \lambda_0r_1 + \mu_0r_2 = r_0.$$ 

@macro <handle ellipse case> += @begin
    R::FT rho = conic0.r;
@end

(ii) in the terminology of the previous subsubsection, 
the second task is to decide the sign of 
$$\frac{\partial}{\partial\tau}\det
\left(\frac{M(\tau)}{z(\tau)}\right)\mid_{\tau=0},$$
where 
$$M(\tau) = \left(\begin{array}{cc} r(\tau) & t(\tau) \\ t(\tau) & s(\tau)
\end{array}\right)$$ and 

$$z(\tau)=\left(u(\tau), v(\tau)\right)M(\tau)^{-1}
\left(u(\tau), v(\tau)\right)^T - w(\tau).$$

For the sake of readability, we omit the parameter $\tau$ 
in the sequel. Noting that 
$$M^{-1} = \frac{1}{\det(M)}\left(\begin{array}{cc}s & - t \\ 
- t & r\end{array}\right),$$
we get
$$z = \frac{1}{\det(M)} (u^2s-2uvt+v^2r) - w.$$

Let us introduce the following abbreviations.
\begin{eqnarray*}
d &:=& \det(M), \\
Z &:=& u^2s-2uvt+v^2r.
\end{eqnarray*}
With primes ($d',Z'$ etc.) we denote derivatives w.r.t. $\tau$.
Now we can write $$\frac{\partial}{\partial \tau}\det(M/z) = (d/z^2)' =
\frac{d'z^2 - 2dzz'}{z^4} = \frac{d'z -2dz'}{z^3}.$$

Since $d(0),z(0)>0$ (recall that ${\cal C}_0$ is an ellipse), this is
equal in sign to $$f := d(d'z-2dz'),$$
at least for $\tau=0$ which is the value we are interested in. 
Now we have
\begin{eqnarray*}
d'z &=& d'(\frac{1}{d}Z-w) = \frac{d'}{d}Z-d'w, \\
dz' &=& d(\frac{Z'd - Zd'}{d^2} - w') = \frac{Z'd - Zd'}{d} - dw'.
\end{eqnarray*}
Hence 
\begin{eqnarray*} 
f    &=& d'Z-dd'w - 2(Z'd -Zd' - d^2w') \\
     &=& 3d'Z + d(2dw'-d'w - 2Z').
\end{eqnarray*}
Now it's time to look at the explicit expressions for $d,d',Z,Z'$. 
Rewriting $Z$ as
\begin{eqnarray*} 
Z &=& u (us - vt) + v (vr-ut) \\
  &=& uZ_1+vZ_2,
\end{eqnarray*} we get 
\begin{eqnarray*}
d   &=& rs-t^2, \\
d'  &=& r's + rs' - 2tt', \\
Z_1' &=& u's + us' - v't -vt', \\
Z_2' &=& v'r + vr'  -u't -ut', \\
Z'   &=& u'Z_1 + uZ_1' + v'Z_2 + v Z_2'.
\end{eqnarray*}
Evaluated for $\tau=0$, all these values can be computed from 
$r(0),s(0),t(0),u(0),v(0),w(0)$ (the defining values of the conic
${\cal C}_0$) and their corresponding primed values which have already
been computed by the @prg{set} method. $d(0)=\det({\cal C}_0)$ is already
defined, the other values follow (recall that $r'(0)=0$). For all this,
${\cal C}_0$ is assumed to be normalized.

@macro <handle ellipse case> += @begin
    conic0.normalize();
    R::FT  dd  = conic0.r*ds - R::FT(2)*conic0.t*dt;
    R::FT  Z1  = conic0.u*conic0.s - conic0.v*conic0.t;
    R::FT  Z2  = conic0.v*conic0.r - conic0.u*conic0.t;
    R::FT  dZ1 = du*conic0.s + conic0.u*ds - dv*conic0.t - conic0.v*dt;
    R::FT  dZ2 = dv*conic0.r - du*conic0.t - conic0.u*dt;
    R::FT  Z   = conic0.u*Z1 + conic0.v*Z2;
    R::FT  dZ  = du*Z1 + conic0.u*dZ1 +dv*Z2 + conic0.v*dZ2;
    R::FT  f   = R::FT(3)*dd*Z + d*(R::FT(2)*d*dw-dd*conic0.w-R::FT(2)*dZ);
@end

The sign of $f(0)$ is the one we are interested in. We deduce 
from \cite{GS} that
$$p \left\{\begin{array}{c} {\rm inside} \\ {\rm on} \\ {\rm outside}
\end{array}\right\} E^*
\Leftrightarrow \rho~f(0) \left\{\begin{array}{c}< 0 \\ = 0 \\ > 0
\end{array}\right.,$$
$E^{*}$ the smallest ellipse through $p_1,p_2,p_3,p_4$. 

@macro <handle ellipse case> += @begin
    return static_cast (CGAL_Bounded_side, CGAL_sign (rho)*CGAL_sign (f));
@end

@! ----------------------------------------------------------------------------
@subsection{ Ellipse\_5 Class}
@! ----------------------------------------------------------------------------

This class is the simplest one among the three representation classes, since
it does not need to compute its ellipse. Recall that the current support set
$B$ attains cardinality five only if previously, a point
$p$ was found to lie outside $\me(\emptyset,\{p_1,p_2,p_3,p_4\})$, 
$\{p_1,p_2,p_3,p_4\}$ the support set at that time. However, during this test
with $p$, the unique conic (which is then an ellipse)  through 
$\{p_1,p_2,p_3,p_4,p\}$ has already been computed by the @prg{Ellipse_4}
representation. In fact, it is the conic ${\cal C}_0$ addressed in the
previous subsection, and it suffices to store a pointer to it in the
@prg{Ellipse_5} object, which is then initialized by the @prg{set} method.

@macro <Ellipse_5 data members> = @begin
    Conic<R>* ellipse;
@end

@macro <Ellipse_5::set body> = @begin
    ellipse = &(ellipse_4.conic0);
@end

The @prg{bounded_side} method is then straightforward, noting that the
ellipse has already been normalized by @prg{ellipse_4}.

@macro <Ellipse_5::bounded_side body> = @begin
    R::FT discr = (*ellipse).eval(p);
    return static_cast (CGAL_Bounded_side, CGAL_sign (discr));
@end 
    
@! ==========================================================================
@section{ File Organisation}
@! ==========================================================================

@file <Min_ellipse_2.h> = @begin
    @<file header>("2D Smallest Enclosing Ellipse",
		   "include/CGAL/Min_ellipse_2.h","Min_ellipse_2",
		   "Bernd Gärtner, Sven Schönherr (sven@@inf.fu-berlin.de)")
    #ifndef CGAL_MIN_ELLIPSE_2_H
    #define CGAL_MIN_ELLIPSE_2_H

    // check macros
    // ------------
    #ifdef CGAL_CHECK_ASSERTIONS
    #define CGAL_Min_ellipse_2_assertion(EX) \
     ((EX) ? ((void)0) : cgal_assertion_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_ellipse_2_assertion_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_assertion_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_ellipse_2_assertion(EX) ((void)0)
    #define CGAL_Min_ellipse_2_assertion_msg(EX,MSG) ((void)0)
    #endif // CGAL_CHECK_ASSERTIONS

    #ifdef CGAL_CHECK_PRECONDITIONS
    #define CGAL_Min_ellipse_2_precondition(EX) \
     ((EX) ? ((void)0) : cgal_precondition_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_ellipse_2_precondition_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_precondition_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_ellipse_2_precondition(EX) ((void)0)
    #define CGAL_Min_ellipse_2_precondition_msg(EX,MSG) ((void)0)
    #endif // CGAL_CHECK_PRECONDITIONS

    #ifdef CGAL_CHECK_POSTCONDITIONS
    #define CGAL_Min_ellipse_2_postcondition(EX) \
     ((EX) ? ((void)0) : cgal_postcondition_fail( #EX , __FILE__, __LINE__, NULL))
    #define CGAL_Min_ellipse_2_postcondition_msg(EX,MSG) \
     ((EX) ? ((void)0) : cgal_postcondition_fail( #EX , __FILE__, __LINE__, MSG))
    #else
    #define CGAL_Min_ellipse_2_postcondition(EX) ((void)0)
    #define CGAL_Min_ellipse_2_postcondition_msg(EX,MSG) ((void)0)
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


    // Class declaration
    // =================
    template < class R >
    class CGAL_Min_ellipse_2;

    // class CGAL_Bbox_2;

    // Class interface
    // ===============
    // includes
    #ifndef CGAL_POINT_2_H
    #include <CGAL/Point_2.h>
    #endif
    #include <vector.h>

    @<Min_ellipse_2 interface>

    @<dividing line>

    // Class definition
    // ================
    // includes
    // --------
    #ifndef CGAL_PREDICATES_ON_POINTS_2_H
    #include <CGAL/predicates_on_points_2.h>
    #endif
    #ifndef CGAL_UTILS_H
    #include <CGAL/utils.h>
    #endif
    #include <algo.h>
    #include <algobase.h>

    @<Min_ellipse_2 implementation>

    // specialization does not exist yet
    // ---------------------------------
    // #if ( defined( CGAL_INTEGER_H) && defined( CGAL_HOMOGENEOUS_H))
    // #include <CGAL/Min_ellipse_2__Homogeneous_integer.h>
    // #endif

    // Conic class definition
    @<conic class>

    #endif // CGAL_MIN_ELLIPSE_2_H

    @<end of file line>
@end
   
@i file_header.awlib

\begin{thebibliography}{Wel} 
\bibitem{GS}
B.~G{\"a}rtner and S.~Sch{\"o}nherr.
\newblock Smallest enclosing ellipses -- fast end exact.
\newblock Manuscript, 1997

\bibitem{Wel}
E.~Welzl.
\newblock Smallest enclosing disks (balls and ellipsoids).
\newblock In H.~Maurer, editor, {\em New Results and New Trends in Computer
  Science}, volume 555 of {\em Lecture Notes in Computer Science}, pages
  359--370. Springer-Verlag, 1991.
\end{thebibliography}

