@! ============================================================================
@! The CGAL Library
@! Implementation: 2D Conic
@! ----------------------------------------------------------------------------
@! file  : web/Optimisation/Conic_2.aw
@! author: Bernd Gärtner, Sven Schönherr <gaertner@inf.ethz.ch>
@! ----------------------------------------------------------------------------
@! $CGAL_Chapter: Geometric Optimisation $
@! $CGAL_Package: Min_ellipse_2 WIP $
@! $Revision$
@! $Date$
@! ============================================================================

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{epsf}
@usepackage{cc_manual}
@article

\setlength{\parskip}{1ex}
\setcounter{secnumdepth}{4}
\setcounter{tocdepth}{4}

@! ============================================================================
@! Title
@! ============================================================================

\RCSdef{\rcsRevision}{$Revision$}
\RCSdefDate{\rcsDate}{$Date$}
\newcommand{\cgalWIP}{{\footnotesize{} (\rcsRevision{} , \rcsDate) }}

@t vskip 5 mm
@t title titlefont centre "2D Conic*"
@t vskip 1 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
\smallskip
\begin{center}
  \begin{tabular}{l}
    \verb+$CGAL_Chapter: Geometric Optimisation $+ \\
    \verb+$CGAL_Package: Min_ellipse_2 WIP+\cgalWIP\verb+$+ \\
  \end{tabular}
\end{center}
@t vskip 1 mm

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\footnotetext[1]{This work was supported by the ESPRIT IV LTR Project
  No.~21957 (CGAL).}

@t table_of_contents

\renewcommand{\R}{I\!\!R}
\newcommand{\C}{{\cal C}}
\renewcommand{\E}{{\cal E}}
\newcommand{\Vol}{\mathop{\rm Vol}}
\renewcommand{\r}{{\cal R}}
\renewcommand{\d}{\partial}
\newcommand{\sgn}{\mathop{\rm sgn}}
\newcommand{\I}{{\bf i}}
\newcommand{\acos}{\mathop{\rm acos}}
\newcommand{\asin}{\mathop{\rm asin}}
\newtheorem{Definition}{Definition}[section]
\newtheorem{Lemma}[Definition]{Lemma}

@! ============================================================================
@section{Introduction}
@! ============================================================================

We define a class template @prg{Conic_2<R>} to store, access and
manipulate @em{conics} in the plane, a.k.a. @em{second order curves}.

@! ---------------------------------------------------------------------------
@subsection{Conics}
@! ---------------------------------------------------------------------------

For a given real vector $\r=(r,s,t,u,v,w)$, a conic 
$\C=\C(\r)$ is the set of @em{homogeneous} points
$p = (x,y,h), h\neq 0$ satisfying

\begin{equation}
\label{conic_def_hom}
\r(p) := rx^2 + sy^2 + txy + uxh + vyh + wh^2 = 0,
\end{equation}
equivalently

\begin{equation}
\label{conic_hom}
(x, y, h) 
\left( \begin{array}{ccc}
2r & t & u \\
t & 2s & v \\
u & v & 2w
\end{array} \right)
\left( \begin{array}{c}
x \\
y \\
h
\end{array} \right) = 0.
\end{equation}

$\r$ is called a @em{representation} of $\C$.
Note that the homogeneous point $(x,y,h)$ corresponds 
to the Cartesian point $(x/h, y/h)$ in the plane. Also, any Cartesian 
point $p=(x,y)$ can be identified with the homogeneous point $(x,y,1)$,
in which case (\ref{conic_def_hom}) assumes the form
\begin{equation}
\label{conic_def_cart}
\r(p) = rx^2 + sy^2 + txy + ux + vy + w = 0,
\end{equation}
equivalently
\begin{equation}
\label{conic_cart}
(x,y) 
\left( \begin{array}{ccc}
2r & t \\
t & 2s 
\end{array} \right)
\left( \begin{array}{c}
x \\
y
\end{array} \right)
+ (2u, 2v) \left( \begin{array}{c}
x \\
y
\end{array} \right)
+ 2w = 0.
\end{equation}

Thus, under the condition $h\neq 0$, homogeneous and Cartesian 
representation are equivalent, and we frequently switch between them.
\footnote{The reason for scaling equations (\ref{conic_hom}) 
and (\ref{conic_cart}) by a factor of 2 is purely technical -- we 
want to argue with division-free terms.}

$\C$ is called @em{trivial} if $\C=\R^2$ which is equivalent to 
$\r=0$. $\C$ is @em{empty} if $\C=\emptyset$, i.e. if
(\ref{conic_def_hom}) has no real solutions $x,y,h$ with $h\neq 0$
(which happens e.g. in case of $\r = (1, 1, 0, 0, 0, 1)$, or
$\r=(0,0,0,0,0,1)$).

$\C$ is invariant under scaling its representation $\r$
by any nonzero factor. This means, a conic has five degrees of freedom,
and in fact it holds that any five points uniquely determine a nontrivial
conic passing through the points. Some care is in place: this does
@em{not} mean that five points uniquely determine the conic's 
representation, up to scaling. For example, the $x$-axis is a conic
uniquely determined by any five points on it, but as a representation
we may choose $\{y=0\}$ or $\{y^2=0\}$. Recall that we have already seen
two representations of the empty conic which are not multiples of each
other.

@! ---------------------------------------------------------------------------
@subsection{Conic types}
@! ---------------------------------------------------------------------------
\label{types_sec}
The number
\begin{equation}
\label{det}
\det(\r) := \det \left( \begin{array}{cc} 2r & t \\ t & 2s \end{array}\right)
\end{equation}
determines the type of $\C(\r)$. If $\det(\r)>0$,
$\C$ is an @em{ellipse}, if $\det(\r)<0$, we get a @em{hyperbola}, and
for $\det(\r)=0$, a @em{parabola} is obtained. While the trivial conic
is a degenerate parabola equal to the whole plane, any nontrivial one 
consists of at most two simple curves. As a special case, there is a 
conic $\C=\{p\}$ for any point $p=(x_0,y_0)$. It can be specified
as 
\[
\C = \{(x,y)\mid (x-x_0)^2 + (y-y_0)^2 = 0\},
\]
and a possible representation is 
$\r=(1, 1, 0, -2x_0, -2y_0, x_0^2+y_0^2)$. This implies
$\det(\r)=4$, so $\C$ is a degenerate ellipse, see 
subsection \ref{orientation_sec}.

Note that $\det(\r)> 0$ implies $r,s > 0$ or $r,s<0$ which is equivalent
to 
\begin{equation}
M :=  \left( \begin{array}{cc} 2r & t \\ t & 2s \end{array}\right)
\end{equation}
being positive definite ($x^TMx>0$ for $x\neq 0$) or negative definite
($x^TMx<0$ for $x\neq 0$). In case of $\det(\r)<0$, $M$ is indefinite,
meaning that $x^TMx$ assumes positive and negative values.     

@! ---------------------------------------------------------------------------
@subsection{Symmetry properties}
@! ---------------------------------------------------------------------------
\label{symmetry_sec}

If the conic $\C(\r)$ is not a parabola, the matrix 
\[
M = \left(\begin{array}{cc} 2r & t \\ t & 2s \end{array}\right)
\]
is regular, and $\C$ has a unique center of symmetry $c$, given as   
\[ 
c = -M^{-1}\left(\begin{array}{c} u \\ v \end{array}\right). 
\]

With this definition, (\ref{conic_cart}) can 
alternatively been written as 
\begin{equation}
\label{center_form}
(p-c)^T M (p-c) + 2w-c^TMc = 0,
\end{equation}
$p=(x,y)^T$, from which the symmetry is obvious. In case of a 
parabola, we get an axis of symmetry. 

@! ---------------------------------------------------------------------------
@subsection{Orientation and degeneracy}
@! ---------------------------------------------------------------------------
\label{orientation_sec}

An @em{oriented conic} is a pair $\C_{\r} = (\C(\r),\r)$, i.e. a conic with a
particular representation. $\C_{\r}$ subdivides $\R^2\setminus \C(\r)$
into a @em{positive} side, formed by the set of points such that $\r(p)>0$,
and a @em{negative} side ($\r(p)<0$). Replacing $\r$ with $-\r$ leads to
an oriented conic with positive and negative sides interchanged. This
concept of assigning positive and negative sides is purely algebraic.

In addition, there is a geometric way of assigning sides to an oriented
conic which does not depend on $\r$ (like an oriented circle has a bounded 
side, independent from its orientation). To this end, we define the 
@em{convex side} of a conic $\C$ as the the union of the convex 
connected components of $\R^2\setminus \C$. The @em{non-convex side} is 
then just the union of the non-convex components. Figure \ref{orientations}
depicts the convex sides of an ellipse, hyperbola and parabola, labeled
with the letter `$c$'. 

\begin{figure}
\begin{center}
\leavevmode
\epsfxsize=12cm
%\epsfbox{orientations.eps}
\end{center}
\caption{Ellipse, hyperbola and parabola with convex sides}
\label{orientations}
\end{figure}

A conic is defined to be @em{degenerate} if either its convex side or
its non-convex side is empty. 
Let us discuss the possibilities for this. First, the trivial conic
is degenerate, with both sides being empty. The empty conic is degenerate,
with the non-convex side being empty. An ellipse is
degenerate if and only if it consists of just one point (and so the
convex side is empty). A hyperbola is degenerate if and only if it 
contains its center of symmetry $c$. In this case, the hyperbola is a
pair of lines crossing at $c$, and so the non-convex side is empty.
A degenerate parabola is either a pair of parallel lines 
or just one line. In both cases, the non-convex side is empty.

In the non-degenerate case, the classification of points by positive 
and negative side coincides with the one by convex and non-convex side.
In the degenerate case, this exactly holds if positive or
negative side disappear, like for a degenerate ellipse (but not for a
degenerate hyperbola).

An oriented conic $\C_{\r}$ is said to have positive (negative)
orientation, if and only if the convex side coincides with the 
positive (negative) side. If neither is the case, the orientation 
is zero. Thus, a degenerate ellipse has nonzero orientation,
but a degenerate hyperbola has not. 

While it is clear that positive and negative side of a conic are only
defined with respect to some representation, it is interesting to
note that even the partition of $\R^2\setminus\C_{\r}$ into 
positive and negative side does in general depend on $\r$. Coming
back to the conic $\C = \{y=0\} = \{y^2=0\}$, the first representation
of it leads to nonempty positive and negative side, while in the second 
one, the negative side is empty. 

@! ---------------------------------------------------------------------------
@subsection{Ellipses and the volume formula}
@! ---------------------------------------------------------------------------

The volume of an ellipse $\E$, $\Vol(\E)$, is defined as the area of 
its convex side. If $\E$ is non-degenerate and presented in center form 
(\ref{center_form}), consider the matrix 
$$A := M / (2w - c^TMc).$$ It is easy to see that $A$ is 
invariant under scaling $\r$ by any nonzero factor. The following
holds.
\begin{Lemma}
\label{ellipse_volume}
$$\Vol(\E) = \frac{\pi}{\sqrt{\det(A)}}.$$
\end{Lemma}
For this note that $\det(M)>0$ implies $\det(A)>0$. 

@! ---------------------------------------------------------------------------
@section{Class Template {\tt CGAL\_Conic\_2<R>}}
@! ---------------------------------------------------------------------------

An object of the class @prg{Conic_2<R>} stores an oriented
conic $\C_{\r} = (\C(\r),\r)$ by its representation $\r$. The template 
parameter @prg{R} is the representation type (homogeneous or Cartesian). In
the sequel, a conic always means an oriented conic.

@macro <Conic_2 declaration> = @begin
    template < class R>
    class Conic_2;

@end

The class @prg{Conic_2<R>} has several member functions, some of which 
can be considered general purpose (and are declared public), others are 
more specialized and needed only by the class template 
@prg{Optimisation_ellipse_2<R>} -- they appear as private methods.

@macro <Optimisation_ellipse_2 declaration> = @begin
    template < class _R >
    class Optimisation_ellipse_2;
@end

@macro <Optimisation_ellipse_2 I/O operator declaration> = @begin
    template < class _R >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream&,
                  const CGAL::Optimisation_ellipse_2<_R>&);
@end
 
As a general rule, a method only qualifies for the public domain if 
it can be specified without referring to particular values of $\r$, in
other words, if its behavior only depends on the equivalence class of
all positive multiples of $\r$. For example, this holds for the sign of 
$\r(p)$, $p$ some point, but not for the value of $\r(p)$. As usual, there
is an exception to this rule: we have public methods to retrieve the 
components of the representation $\r$.   
 
All calls to member functions are 
directly mapped to calls of corresponding methods declared by 
the conic class @prg{R::Conic_2} of the representation type @prg{R}. 

@macro <Conic_2 interface and implementation> = @begin
    template < class _R>
    class Conic_2 : public _R::Conic_2 {

        friend  class Optimisation_ellipse_2<_R>;

      public:

        // types
        typedef  _R                    R;
        typedef  typename _R::RT       RT;
        typedef  typename _R::FT       FT;
        typedef  typename _R::Conic_2  _Conic_2;

        // construction
        @<Conic_2 constructors> 

        // general access
        @<Conic_2 general access methods>

        // type related access
        @<Conic_2 type access methods>

        // orientation related access 
        @<Conic_2 orientation access methods>

        // comparisons
        @<Conic_2 comparison methods>

        // set methods
        @<Conic_2 set methods>

      private:
        @<Conic_2 private methods>
    };

@end


@! ---------------------------------------------------------------------------
@subsection{Construction}
@! ---------------------------------------------------------------------------

We have a default constructor and a constructor to initialize a conic
at coordinate level, providing the components of a representation $\r$.  
The orientation of the conic is determined by $\r$, as described in 
subsection \ref{orientation_sec}.  

@macro <Conic_2 constructors> = @begin
    Conic_2 ()     
    {}

    Conic_2 (RT r, RT s, RT t, RT u, RT v, RT w)
        : R::Conic_2 (r, s, t, u, v, w)
    {}

@end

@! ---------------------------------------------------------------------------
@subsection{General access}
@! ---------------------------------------------------------------------------

We provide access to the coordinates of the representation $\r$.  
Even if $\C_{\r}$ was constructed from some specific characteristic
vector, the return coordinates do not necessarily belong to this 
vector (but at least to a nonnegative multiple of it). For example,
if the coordinates are from some integer domain, the implementation might 
decide to divide all of them by their gcd to get smaller numbers. (In this
implementation, this does not happen.)

@macro <Conic_2 general access methods> += @begin
    RT r () const 
    {
        return _Conic_2::r();
    }

    RT s () const 
    {
        return _Conic_2::s();
    }

    RT t () const 
    {
        return _Conic_2::t();
    }

    RT u () const 
    {
        return _Conic_2::u();
    }

    RT v () const 
    {
        return _Conic_2::v();
    }

    RT w () const 
    {
        return _Conic_2::w();
    }

@end

We can obtain the center of symmetry of $\C_{\r}$. Precondition is that
$\C$ is not a parabola. 

@macro <Conic_2 general access methods> += @begin
    CGAL::Point_2<R> center () const
    {
        return _Conic_2::center();
    }

@end

The symmetry axis of a parabola and the two axes of an ellipse resp. a
hyperbola require irrational coordinates in general, therefore they cannot
be added here without further assumptions about the representation type
@prg{R}. 

@! ---------------------------------------------------------------------------
@subsection{Type related access}
@! ---------------------------------------------------------------------------

Here we have access methods and predicates related to the type of 
$\C_{\r}$. 
We can either directly retrieve the type or ask whether $\C_{\r}$ is of 
some specific type. To realize the former, we provide a suitable 
enumeration type.

@macro <Conic_type declaration> = @begin
    enum Conic_type 
    {
        HYPERBOLA = -1,
        PARABOLA,
        ELLIPSE
    };

@end

@macro <Conic_2 type access methods> += @begin
    Conic_type conic_type () const
    {
        return _Conic_2::conic_type();
    }

    bool is_hyperbola () const
    {
        return _Conic_2::is_hyperbola();
    }

    bool is_parabola () const
    {
        return _Conic_2::is_parabola();
    }

    bool is_ellipse () const
    {
        return _Conic_2::is_ellipse();
    }

@end

We have three more predicates that -- in combination with the above -- 
yield a finer access to the type, namely a test for emptiness, triviality
and degeneracy, where both the empty and the trivial conic are also
degenerate. 

@macro <Conic_2 type access methods> += @begin
    bool is_empty () const
    {
        return _Conic_2::is_empty();
    }

    bool is_trivial () const
    {
        return _Conic_2::is_trivial();
    }

    bool is_degenerate () const
    {
        return _Conic_2::is_degenerate();
    }

@end

@! ---------------------------------------------------------------------------
@subsection{Orientation related access}
@! ---------------------------------------------------------------------------

We can retrieve the orientation of a conic as defined in subsection
\ref{orientation_sec}. A notable difference to the class 
@prg{Circle_2<R>} is that the orientation can be zero. 

@macro <Conic_2 orientation access methods> += @begin
    CGAL::Orientation orientation () const
    {
        return _Conic_2::orientation ();
    }

@end

The following methods classify points according to the positive
and negative side of $\C_{\r}$.  We can either retrieve the oriented 
side directly or ask whether the point lies on some specific side. 

@macro <Conic_2 orientation access methods> += @begin   
    CGAL::Oriented_side oriented_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::oriented_side (p);
    }

    bool has_on_positive_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_positive_side (p);
    }

    bool has_on_negative_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_negative_side (p);
    }

    bool has_on_boundary (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_boundary (p);
    }
        
    bool has_on (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on (p);
    }

@end

Then we have the classification according to convex and non-convex side. 
Because this is a natural generalization of the classification according
to bounded and unbounded side, we `extend' the type @prg{Bounded_side}.

@macro <Convex_side declaration> = @begin
    typedef CGAL::Bounded_side Convex_side;
    const Convex_side ON_CONVEX_SIDE    = CGAL::ON_BOUNDED_SIDE;
    const Convex_side ON_NONCONVEX_SIDE = CGAL::ON_UNBOUNDED_SIDE; 

@end 

@macro <Conic_2 orientation access methods> += @begin   
    Convex_side convex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::convex_side (p);
    }

    bool has_on_convex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_convex_side (p);
    }

    bool has_on_nonconvex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_nonconvex_side (p);
    }

@end

@! ---------------------------------------------------------------------------
@subsection{Comparison methods}
@! ---------------------------------------------------------------------------

We provide tests for equality and inequality of two conics.

@macro <Conic_2 comparison methods> = @begin
    bool operator == ( const Conic_2<_R>& c) const
    {
        return _Conic_2::operator == ( (Conic_2)c);
    }

    bool operator != ( const Conic_2<_R>& c) const
    {
        return( ! operator == ( c));
    }
@end

@! ---------------------------------------------------------------------------
@subsection{Public set methods}
@! ---------------------------------------------------------------------------

Apart from the most basic way of constructing an oriented 
conic from the coordinates of
a representation, there are methods at a higher level, building a conic
from points or other conics. Such methods appear here; they are not given 
as constructors, but as `set' functions that modify an already existing
conic. But of course, there is also a set function that works on the
coordinate level.

@macro <Conic_2 set methods> += @begin
    void set (RT r, RT s, RT t,
              RT u, RT v, RT w)
    {
        _Conic_2::set (r, s, t, u, v, w);
    }

@end


@! ---------------------------------------------------------------------------
@subsubsection{Opposite conic}
@! ---------------------------------------------------------------------------

We can obtain the conic $\C_{-\r}$ of opposite orientation, having positive 
and negative sides interchanged (convex and non-convex side stay the same,
of course). Note that if $\C_{\r}$ has zero orientation, so has $\C_{-\r}$. 

@macro <Conic_2 set methods> += @begin
    void set_opposite ()
    {
        _Conic_2::set_opposite();
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Pair of lines through four points}
@! ---------------------------------------------------------------------------

The method @prg{set_linepair} builds the conic equal to union of the
lines $\overline{p_1p_2}$ and $\overline{p_3p_4}$. This is either a
degenerate hyperbola in case the two lines are not parallel, or a
degenerate parabola.  The precondition is that $p_1\neq p_2$ and
$p_3\neq p_4$. The positive side is the region to the left resp. to
the right of both oriented lines $\overline{p_1p_2}, \overline{p_3p_4}$.

@macro <Conic_2 set methods> += @begin    
    void set_linepair (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                       const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4)
    {
        _Conic_2::set_linepair (p1, p2, p3, p4);
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Smallest ellipse through three points}
@! ---------------------------------------------------------------------------

The following method @prg{set_ellipse} generates the ellipse of smallest 
volume through $p_1,p_2$
and $p_3$. This is at the same time the unique ellipse through these points
whose center is equal to the center of gravity of the points. 
Precondition is that $p_1,p_2,p_3$ are not collinear. The orientation of
the ellipse is the orientation of the point triple. 

@macro <Conic_2 set methods> += @begin
    void set_ellipse (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                      const CGAL::Point_2<R>& p3)
    {
        _Conic_2::set_ellipse (p1, p2, p3);
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Some ellipse through four points in convex position}
@! ---------------------------------------------------------------------------

We have another @prg{set_ellipse} method, constructing an ellipse from
four given points, assumed to be in convex position (if not, the result
will still be some conic through the points, but not an ellipse). 
The result will be just some ellipse through the points; in particular, it
will usually @em{not} be the smallest one passing through the points. The
orientation of the ellipse can be specified and defaults to positive. 

@macro <Conic_2 set methods> += @begin
    void set_ellipse (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                      const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4,
                      CGAL::Orientation o = POSITIVE)
    {
        _Conic_2::set_ellipse (p1, p2, p3, p4, o);
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Unique conic through five points}
@! ---------------------------------------------------------------------------

The method @prg{set} generates the unique nontrivial conic containing
the points $p_1,p_2,p_3,p_4,p_5$. Precondition is that all points are 
distinct. The orientation can be specified but is automatically set to
zero if the resulting conic has zero orientation.  

@macro <Conic_2 set methods> += @begin         
    void set (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
              const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4,
              const CGAL::Point_2<R>& p5, 
              CGAL::Orientation o = POSITIVE)
    {
        _Conic_2::set (p1, p2, p3, p4, p5, o);
    }

@end


@! ---------------------------------------------------------------------------
@subsection{Private methods}
@! --------------------------------------------------------------------------- 

@! ---------------------------------------------------------------------------
@subsubsection{Linear combination of conics}
@! ---------------------------------------------------------------------------

The linear combination of two conics $\C_{\r_1}$ and $\C_{\r_2}$ with
coefficients $a_1,a_2$ is the conic with representation 
$\r := a_1\r_1+a_2\r_2$. This is not a geometric operation, 
because the result depends on $\r_1,\r_2$ and not only on 
$\C_{\r_1}$ and $\C_{\r_2}$.  

@macro <Conic_2 private methods> += @begin
    void set_linear_combination (
        const RT& a1, const Conic_2<R>& c1,
        const RT& a2, const Conic_2<R>& c2)
    {
        _Conic_2::set_linear_combination (a1, c1, a2, c2);
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Two pairs of lines through four points}
@! ---------------------------------------------------------------------------

The following method constructs two line-pairs through four
given points $p_1,p_2,p_3,p_4$, assumed to be in convex position. If
the points are enumerated in clockwise or counterclockwise order, these
will be the line-pairs $\overline{p_1p_2}\cup\overline{p_3p_4}$ and
$\overline{p_2p_3}\cup\overline{p_4p_1}$. Otherwise, points are renamed
first to achieve (counter)clockwise orientation. The two resulting
conics are passed to the method by reference, and their orientations 
depend on the points as described in the method @prg{set_linepair}.
If $p_1,p_2,p_3,p_4$ are not in convex position, the method still
computes two line-pairs through the points, but it is not specified
which ones. This is a static method.

@macro <Conic_2 private methods> += @begin
    static void set_two_linepairs (const CGAL::Point_2<R>& p1, 
                                   const CGAL::Point_2<R>& p2,
                                   const CGAL::Point_2<R>& p3,
                                   const CGAL::Point_2<R>& p4,
                                   Conic_2<R>& pair1, 
                                   Conic_2<R>& pair2)
    {
        _Conic_2::set_two_linepairs (p1, p2, p3, p4, pair1, pair2);
    }

@end 

@! ---------------------------------------------------------------------------
@subsubsection{Some ellipse from two pairs of lines}
@! ---------------------------------------------------------------------------

We have a method to construct an ellipse through $p_1,p_2,p_3,p_4$ in convex
position, using two line-pairs through the points as  
obtained from a call to the method @prg{set_two_linepairs} with parameters
$p_1,p_2,p_3,p_4$. The orientation of the ellipse is not specified. In case 
the argument conics have not been constructed by the method 
@prg{set_two_linepairs} in the described way, the resulting 
conic is unspecified.   

@macro <Conic_2 private methods> += @begin
    void set_ellipse (const Conic_2<R>& pair1, 
                      const Conic_2<R>& pair2)
    {
        _Conic_2::set_ellipse (pair1, pair2);
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Conic from two conics and a point}
@!----------------------------------------------------------------------------

Assuming that we already have two conics intersecting in four points, 
the following alternative @prg{set} method is more efficient in constructing
the unique nontrivial conic through them and another point $p$. It accepts 
two conics $\C_{\r_1}$ and $\C_{\r_2}$ and a point $p$, and computes some 
conic $\C_{\r}$ that goes through $\C_{\r_1}\cap \C_{\r_2} \cup\{p\}$. 
The method may construct the trivial conic. 
This happens if one of $\C_{\r_1}$ and $\C_{\r_2}$ is 
already trivial, $\C(\r_1)=\C(\r_2)$ holds, or if 
$p\in \C_{\r_1}\cap \C_{\r_2}$. In case of nonzero orientation, the 
actual orientation of $\C_{\r}$ is unspecified. 

@macro <Conic_2 private methods> += @begin
    void set (const Conic_2<R>& c1, const Conic_2<R>& c2,
              const CGAL::Point_2<R>& p)
    {
        _Conic_2::set( c1, c2, p);  analyse();
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Volume derivative of an ellipse}
@! ---------------------------------------------------------------------------

Given an ellipse $\E$ with representation $\r$ and some vector
$\d \r := (\d r, \d s, \d t, \d u, \d v, \d w)$, define
\[
\E(\tau) = \C_{\r + \tau \d \r}.
\]

For small values of $\tau$, $\E(\tau)$ is still an ellipse. The method
@prg{vol_derivative} computes the sign of 
\[
\frac{\partial}{\partial \tau} \Vol(\E(\tau))\left|_{\tau=0}\right.,
\]
i.e. decides how the volume develops when going from $\E(0)$ `in direction'
$\d\r$. If $\E$ is not an ellipse, the result is meaningless. 

@macro <Conic_2 private methods> += @begin
    CGAL::Sign vol_derivative (RT dr, RT ds,
                               RT dt, RT du,
                               RT dv, RT dw) const
    {
        return _Conic_2::vol_derivative (dr, ds, dt, du, dv, dw);
    }

@end

A related method computes the value $\tau^*$ such that 
$\Vol(\E(\tau))$ assumes its minimum over the set 
$T=\{\tau \mid \E(\tau) \mbox{~is an ellipse}\}$. Precondition is 
that this minimum exists. If so,  
$\tau^*$ is a local extremum of the volume function. $\tau^*$ might be
irrational, but because the value is only used for drawing the ellipse 
$\E(\tau^*)$, a @prg{double}-approximation suffices.
As before, if $\E$ is not an ellipse, the result is meaningless.

@macro <Conic_2 private methods> += @begin
    double vol_minimum (RT dr, RT ds,
                        RT dt, RT du,
                        RT dv, RT dw) const
    {
        return _Conic_2::vol_minimum (dr, ds, dt, du, dv, dw);
    }

@end

    
@! ---------------------------------------------------------------------------
@subsection{IO routines}
@! ---------------------------------------------------------------------------

@macro <Conic_2 I/O routines> = @begin
    template< class _R>
    std::ostream& operator << ( std::ostream& os, const Conic_2<_R>& c)
    {
        return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
                   << c.u() << ' ' << c.v() << ' ' << c.w());
    }
@end

@! ---------------------------------------------------------------------------
@subsubsection{Output to {\tt CGAL\_Window\_stream}}
@! ---------------------------------------------------------------------------

We provide an operator to write a conic to a
@prg{Window_stream}.  The function is not extraordinarily
efficient but simple (and works without `understanding' the conic;
other methods -- like the one by Maxwell \& Baker -- need to
determine the conic type first, compute start values etc.). The
method works in image space, proceeding in two phases.

Phase 1 draws the conic in $x$-direction. This means that the width
of the output window is scanned pixel-wise, for any $x$-value
computing the at most two corresponding values $y_1,y_2$ such that
$(x,y_1),(x,y_2)\in \C$. (This is done by solving a quadratic
equation). The resulting pixels are stored for output, which is
triggered after all $x$-values have been processed in this way.

Phase 2 draws the conic in $y$-direction, proceeding
similarly. Phases 1 and 2 together ensure that there are no gaps in
the drawn curve(s).

@macro<Conic_2 graphical output operator> = @begin
    #ifdef CGAL_CONIC_2_H
    #ifndef CGAL_IO_WINDOW_STREAM_CONIC_2
    #define CGAL_IO_WINDOW_STREAM_CONIC_2

    template< class R >
    CGAL::Window_stream&
    operator << ( CGAL::Window_stream& ws, const CGAL::Conic_2<R>& c)
    {
        // length of a pixel in window-coordinates
        double pixel = 1/ws.scale(); 

        // pixel dimensions of window
        int width  = (int)((ws.xmax() - ws.xmin()) * ws.scale()) + 1,
            height = (int)((ws.ymax() - ws.ymin()) * ws.scale()) + 1,
            dim    = std::max( width, height);

        // pixel coordinates, stored for faster output 
        double *X = new double [2*dim];
        double *Y = new double [2*dim];

        // actual number of pixels to be drawn
        int pixels;

        // conic coordinates
        double r = CGAL::to_double (c.r()),
               s = CGAL::to_double (c.s()),
               t = CGAL::to_double (c.t()),
               u = CGAL::to_double (c.u()),
               v = CGAL::to_double (c.v()),
               w = CGAL::to_double (c.w());

        // Phase I (drawing in x-direction)
        pixels = 0;
        // solve conic equation for y
        if (s != 0.0) 
            for (double x = ws.xmin(); x <= ws.xmax(); x+=pixel) {
                double discr = (t*t-4.0*r*s)*(x*x) + (2.0*t*v-4.0*s*u)*x + 
                                 v*v - 4.0*s*w;
                if (discr >= 0.0) {
                    double y1 = (-t*x - v - CGAL::sqrt(discr))/(2.0*s);
                    double y2 = (-t*x - v + CGAL::sqrt(discr))/(2.0*s);
                    X[pixels] = x; Y[pixels++] = y1;
                    X[pixels] = x; Y[pixels++] = y2; } }
        else
            for (double x = ws.xmin(); x <= ws.xmax(); x+=pixel) {
                double denom = t*x + v;
                if (denom != 0.0) {
                    double y = -(r*x*x + u*x + w)/denom;
                    X[pixels] = x; Y[pixels++] = y; } }
        ws.draw_pixels (pixels, X, Y);

        // Phase II (drawing in y-direction)
        pixels = 0;
        // solve conic equation for x
        if (r != 0.0) 
            for (double y = ws.ymin(); y <= ws.ymax(); y+=pixel) {
                double discr = (t*t-4.0*r*s)*(y*y) + (2.0*t*u-4.0*r*v)*y + 
                                 u*u - 4.0*r*w;
                if (discr >= 0.0) {
                    double x1 = (-t*y - u - CGAL::sqrt(discr))/(2.0*r);
                    double x2 = (-t*y - u + CGAL::sqrt(discr))/(2.0*r);
                    X[pixels] = x1; Y[pixels++] = y;
                    X[pixels] = x2; Y[pixels++] = y; } }
        else
            for (double y = ws.ymin(); y <= ws.ymax(); y+=pixel) {
                double denom = t*y + u;
                if (denom != 0.0) {
                    double x = -(s*y*y + v*y + w)/denom;
                    X[pixels] = x; Y[pixels++] = y; } }
        ws.draw_pixels (pixels, X, Y);

        // free memory
        delete[] Y;
        delete[] X;

        return( ws);
    }

    #endif // CGAL_IO_WINDOW_STREAM_CONIC_2
    #endif // CGAL_CONIC_2_H
@end


@! ---------------------------------------------------------------------------
@section{Class Templates {\tt CGAL\_ConicHPA2<PT,DA>} 
and {\tt CGAL\_ConicCPA2<PT,DA>}}
@! ---------------------------------------------------------------------------

The two classes described in this section realize the functionality of
the general class @prg{Conic_2<R>}, distinguished between Cartesian and
homogeneous representation. Unlike it is practice in the kernel, the
template parameters are not field and/or ring type but point and data
accessor type. The reason is that the classes described here are also
used to build traits class adapters for the 
class @prg{Min_ellipse_2<Traits>},
where the adapters themselves are templatized with a 
point and data accessor type. 

To realize this scheme, the representation type @prg{R} of the 
class @prg{Conic_2<R>} is enhanced with a data accessor class -- it will 
be the canonical one, realizing access to the coordinates of a 
@prg{Point_2<R>} (see files @prg{homogeneous_rep.h} and 
@prg{cartesian_rep.h}). 

@macro <ConicHPA2 declaration> = @begin
    template < class PT, class DA>
    class ConicHPA2;

    template < class PT, class DA>
    class _Min_ellipse_2_adapterH2__Ellipse;

@end


@macro <ConicCPA2 declaration> = @begin
    template < class PT, class DA>
    class ConicCPA2;

    template < class PT, class DA>
    class _Min_ellipse_2_adapterC2__Ellipse;

@end

@! ---------------------------------------------------------------------------
@subsection{Requirements for template parameters {\tt PT} and {\tt DA}}
@! ---------------------------------------------------------------------------

For @prg{ConicHPA2}, 
the data accessor type @prg{DA} must define a
ring type @prg{RT} and a function @prg{get} to access $x$-, $y$- and 
$h$-coordinate of a point $p=(x,y,h)$ of type @prg{PT}. The type @prg{RT} is
the coordinate type of @prg{PT}, and we expect @prg{PT} to have a constructor
with three arguments $x,y,h$ of type @prg{RT}. In addition, we need an
assignment constructor (for the @prg{center} method). A minimal interface 
of a type @prg{DA} for the homogeneous representation would look as follows.

@macro <ConicHPA2 DA requirements> zero = @begin
    class DA 
    {
        public: 
            // ring type
            typedef /* some CGAL-conform number type */ RT;
            // coordinate access
            void get (const PT& p, RT& x, RT& y, RT& h);
    };

@end

In the Cartesian case, @prg{DA} must declare a field type @prg{FT} and
access via @prg{get} to the $x$- and $y$-coordinate of a point $p=(x,y)$
of type @prg{FT}. @prg{PT} is expected to declare a constructor with two 
arguments of type @prg{FT}. In addition, we need an assignment constructor. 
A minimal interface of a type @prg{DA} in the Cartesian case would
look as follows.

@macro <ConicCPA2 DA requirements> zero = @begin
    class DA 
    {
        public: 
            // field type
            typedef /* some CGAL-conform number type */ FT;
            // coordinate access
            void get (const PT& p, FT& x, FT& y);
    };

@end

In both cases, class @prg{DA} needs to provide a default constructor. 

@! ---------------------------------------------------------------------------
@subsection{Interfaces}
@! ---------------------------------------------------------------------------

The interfaces look similar to the interface of the class
@prg{Conic_2<R>} in the sense that for any public (private)
member function of @prg{Conic_2<R>}, we have a corresponding
public (private) member function here. These are the @em{high-level}
member functions. 

In addition, there are
private data members to store the conic representation and a 
couple of additional @em{low-level} protected member functions. The 
reason for making them protected is that we do not want them to show
up explicitly in the interface of the class @prg{Conic_2<R>}, but
on the other hand, friends of @prg{Conic_2<R>} should be able to
use them.

@! ---------------------------------------------------------------------------
@subsubsection{{\tt CGAL\_ConicHPA2<PT,DA>}}
@! ---------------------------------------------------------------------------

@macro <ConicHPA2 interface and implementation> = @begin
    template < class _PT, class _DA>
    class ConicHPA2
    {
      public:   
        // types
        typedef           _PT      PT;
        typedef           _DA      DA;
        typedef  typename _DA::RT  RT;

      private:
        friend class Conic_2< CGAL::Homogeneous<RT> >;
        friend class _Min_ellipse_2_adapterH2__Ellipse<PT,DA>;

        @<ConicHPA2 private data members>
        @<ConicHPA2 private member functions>

      protected:
        @<ConicHPA2 protected member functions>

      public:
        @<ConicHPA2 public member functions>
    };
@end


@! ---------------------------------------------------------------------------
@subsubsection{{\tt CGAL\_ConicCPA2<PT,DA>}}
@! ---------------------------------------------------------------------------

@macro <ConicCPA2 interface and implementation> = @begin
    template < class _PT, class _DA>
    class ConicCPA2
    {
      public:   
        // types
        typedef           _PT      PT;
        typedef           _DA      DA;
        typedef  typename _DA::FT  FT;

      private:
        friend class Conic_2< CGAL::Cartesian<FT> >;
        friend class _Min_ellipse_2_adapterC2__Ellipse<PT,DA>;

        @<ConicCPA2 private data members>
        @<ConicCPA2 private member functions>

      protected:
        @<ConicCPA2 protected member functions>

      public:
        @<ConicCPA2 public member functions>

     };
@end

@! ---------------------------------------------------------------------------
@section{Implementation of Class templates {\tt CGAL\_ConicCPA2<PT,DA>} and 
{\tt CGAL\_ConicCPA2<PT,DA>}}
@! ---------------------------------------------------------------------------

We implement the classes in their interface to cope with insufficiencies
of the GNU compiler concerning scope operators in typenames. Because 
the implementations of most member functions are very similar for 
both representations, we always write them down in parallel.  

@! ---------------------------------------------------------------------------
@subsection{Private data members}
@! ---------------------------------------------------------------------------

An oriented conic $\C_{\r}$ is stored by its representation $\r$ and
certain @em{derived data}. These data are
\begin{itemize} 
\item the type of $\C_{\r}$, 
\item the orientation of $\C_{\r}$, 
\item degeneracy information, consisting of three flags indicating
whether $\C_{\r}$ is empty, trivial or degenerate.
\end{itemize}
Although type, orientation and degeneracy information can be retrieved from 
$\r$, it is more efficient to store them, because for example, repeated convex 
side tests on the same conic but with different points access these data
over and over again.  

@macro <ConicHPA2 private data members> = @begin
    DA                  dao;
    RT                  _r, _s, _t, _u, _v, _w;
    Conic_type          type;
    CGAL::Orientation   o;
    bool                empty, trivial, degenerate;

@end

@macro <ConicCPA2 private data members> = @begin
    DA                  dao;
    FT                  _r, _s, _t, _u, _v, _w;
    Conic_type          type;
    CGAL::Orientation   o;
    bool                empty, trivial, degenerate;

@end


@! ---------------------------------------------------------------------------
@subsection{Low-level private member functions}
@! ---------------------------------------------------------------------------
\label{private_methods}
Let's start with the low-level members that do not have a counterpart 
in the interface of the class @prg{Conic_2<R>}. 

@! ---------------------------------------------------------------------------
@subsubsection{Determinant}
@! ---------------------------------------------------------------------------

The function @prg{det} just 
computes $\det(\r)=4rs-t^2$, and as mentioned in subsection 
\ref{types_sec}, this value determines the type of the conic.

@macro <ConicHPA2 protected member functions> += @begin
    RT det () const
    {
        return RT(4)*s()*r() - t()*t();
    }

@end

@macro <ConicCPA2 protected member functions> += @begin 
    FT det () const
    {
        return FT(4)*s()*r() - t()*t();
    }

@end


@! ---------------------------------------------------------------------------
@subsubsection{Conic analysis}
@! ---------------------------------------------------------------------------

This method is the most important low-level method. It
initializes the derived data from the representation $\r$, by first
determining the conic's type and then handling the three possible types in
a @prg{case} statement. 

@macro <ConicHPA2 protected member functions> += @begin
    void analyse( ) 
    {
        RT d = det();
        type = (Conic_type)(CGAL_NTS sign(d));
        switch (type) {
        case HYPERBOLA: 
            {
                @<analyse hyperbola, homogeneous case>
            }
            break;
        case PARABOLA: 
            {
                @<analyse parabola, homogeneous case>
            }
            break;
        case ELLIPSE: 
            {
                @<analyse ellipse, homogeneous case>
            }
            break;
        }
    }

@end


@macro <ConicCPA2 protected member functions> += @begin
    void analyse( ) 
    {
        FT d = det();
        type = (Conic_type)(CGAL_NTS sign(d));
        switch (type) {
        case HYPERBOLA: 
            { 
                @<analyse hyperbola, Cartesian case>
            }
            break;
        case PARABOLA:
            {
                @<analyse parabola, Cartesian case>
            }
            break;
        case ELLIPSE:
            {
                @<analyse ellipse, Cartesian case>
            }
            break;
        }
    }

@end

 
Let us first deal with the case where
$\C_{\r}$ is a hyperbola or ellipse. Then we have seen in
subsection \ref{symmetry_sec} that $\C(\r)$ has a center of symmetry $c$ and
can be written in the form 
\[
\{p\mid (p-c)^T M (p-c) + 2w-c^TMc = 0\}, 
\]

Moreover, $\C(\r)$ is degenerate if the center lies on the conic.
This is the case if and only if $z := 2w-c^TMc=0$. To compute this 
value $z$, we go back to the formulas of subsection \ref{symmetry_sec}, where
we have seen that $$c = -M^{-1} \left(\begin{array}{c} u \\ v \end{array}
\right),$$
therefore
$$c^TMc = (u,v)M^{-1}\left(\begin{array}{c} u \\ v \end{array}
\right) = \frac{1}{\det(M)} (u,v) \left(
\begin{array}{cc} 2s & -t \\ -t & 2r \end{array}\right) 
\left(\begin{array}{c} u \\ v \end{array} \right) 
= \frac{1}{\det(M)}(2u^2 s + 2v^2 r - 2uvt).$$
This means, 
\begin{equation}
\label{z}
z = 2 \left(w - \frac{1}{\det(M)}(u^2 s + v^2 r - uvt)\right).
\end{equation}

To avoid divisions, we consider the value $z'=\det(M)z/2$ which satisfies 
$$
z' = \det(M) w - u^2 s - v^2 r + uvt.
$$
$z'$ has the same sign as $z$ in case of an ellipse, and opposite sign 
in case of a hyperbola. 

To proceed further, we note that only a parabola can be trivial. 
Moreover, a hyperbola cannot be empty (the matrix $M$ is indefinite,
therefore $x^TMx$ assumes arbitrary real values). In case of an ellipse,
$M$ is either positive or negative definite, meaning $r>0$ or $r<0$.
In the former case, the ellipse is empty if and only if $z>0$, in the
latter case, $z<0$ leads to an empty ellipse. Summarizing, the ellipse
is empty iff $rz>0$, equivalently $rz'>0$ 
  
Now consider orientation. A hyperbola is in positive 
orientation if and only if its center is on the negative
side, equivalently $z<0$, or $z'>0$. Similarly, the orientation 
is negative in case of $z'<0$. For $z=z'=0$, 
the hyperbola is degenerate, and so its orientation is zero. 

A non-degenerate ellipse, on the other hand, has positive orientation 
if and only if its center is in the positive side, equivalently, if
$z > 0$, or $z' > 0$. This is equivalent to $M$ being negative
definite, or $r<0$. In case of the degenerate ellipse $\E=\{c\}$, we have
$z=z'=0$, and the orientation is positive if and only if the negative side 
is nonempty (the positive side must agree with the empty convex side). As
before, this is equivalent to $M$ being negative definite, or $r<0$.   

In contrast to this, the empty ellipse has positive orientation if and 
only if the negative side is empty (the positive side must agree
with the convex side, which is the whole plane in this case). For
this, we get the equivalent condition $r>0$. 

@macro <analyse hyperbola, homogeneous case> = @begin
    trivial = empty = false;
    RT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();
    o = (CGAL::Orientation)(CGAL_NTS sign (z_prime));
    degenerate = (o == CGAL::ZERO);

@end


@macro <analyse hyperbola, Cartesian case> = @begin
    trivial = empty = false;
    FT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();
    o = (CGAL::Orientation)(CGAL_NTS sign (z_prime));
    degenerate = (o == CGAL::ZERO);

@end


@macro <analyse ellipse, homogeneous case> = @begin
    trivial = false;
    RT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();  
    if (CGAL_NTS is_positive (r())) {
        empty = CGAL_NTS is_positive(CGAL_NTS sign (z_prime));
        empty ? o = CGAL::POSITIVE : o = CGAL::NEGATIVE;
    } else {
        empty = CGAL_NTS is_negative(CGAL_NTS sign (z_prime));
        empty ? o = CGAL::NEGATIVE : o = CGAL::POSITIVE;
    }
    degenerate = empty || CGAL_NTS is_zero (z_prime);

@end

@macro <analyse ellipse, Cartesian case> = @begin
    trivial = false;
    FT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();
    if (CGAL_NTS is_positive (r())) {
        empty = CGAL_NTS is_positive(CGAL_NTS sign (z_prime));
        empty ? o = CGAL::POSITIVE : o = CGAL::NEGATIVE;
    } else {
        empty = CGAL_NTS is_negative(CGAL_NTS sign (z_prime));
        empty ? o = CGAL::NEGATIVE : o = CGAL::POSITIVE;
    }
    degenerate = empty || CGAL_NTS is_zero (z_prime);

@end

In the parabola case, we proceed as follows, first observing that
$\det(M)=0$ implies $r,s\geq 0$ or $r,s\leq 0$.  
Assume that $r\neq 0$. Then the conic equation (\ref{conic_def_hom}) can
be solved for $x$, obtaining
$$x = \frac{-ty - u \pm \sqrt{2(tu - 2rv)y +u^2 - 4rw}}{2r}.$$

The parabola is non-degenerate exactly if the factor $(tu - 2rv)$ is 
nonzero, where we obtain a curved object. For $tu=2rv$, the parabola
is either empty (this happens for $u^2<4rw$ when the radicant becomes
negative), a single line (if $u^2=4rw$), or a pair of lines (for
$u^2>4rw$). Let us treat the degenerate case first. 

In the empty case, the parabola has only one nonempty side, and the 
orientation is determined by $w$. $w>0$ means positive orientation
(because then the point $(0,0)$ is on the positive side, which therefore
equals the convex side in this case), $w<0$ means negative
orientation. Because of $u^2<4rw$, the case $w=0$ cannot occur. 

In case of $u^2 = 4rw$, the conic is given by 
$$\r(p) = r(x + \frac{ty + u}{2r})^2,$$
and the orientation is positive if and only if $r>0$ (in which case
every point $p=(x,y)$ is on the positive side, equivalently on the convex 
side).  

For $u^2 > 4rw$, we get a pair of parallel lines, of zero orientation
(because both positive and negative sides are nonempty, while the
non-convex side is empty). 

We can argue completely similar in the case $s\neq 0$. Here we get a
degeneracy exactly if $tv = 2su$, and the discriminant $v^2-4sw$ determines
whether the parabola is empty, equal to one line or to a pair of lines. 

If $r=s=0$ (implying $t=0$), we have the trivial conic 
if all other parameters are also zero. The trivial conic has
always positive orientation (because both positive and convex side
are empty). If $u=v=0$ but $w\neq 0$, we get the empty conic, where
the orientation is given by $w$. In any other case, we obtain a
single line, this time with zero orientation, because
it has nonempty positive and negative side but empty non-convex side. 

Now consider the non-degenerate case. We claim that the orientation is
positive if and only if $r,s\leq 0$. To see this, note that in this case,
the parabola can be written in the form
$$\r(p) = -(\sqrt{-r} x - \sqrt{-s} y)^2 + ux + vy + w.$$
$\r(p)$ is a concave function which implies that if $\r(p_1),\r(p_2)>0$,
then also $\r(p)>0$, $p$ a convex combination of $p_1,p_2$. This means
that the positive side is a convex set, thus equal to the convex side.
               
@macro <analyse parabola, homogeneous case> = @begin
    if (!CGAL_NTS is_zero (r())) {
        trivial         = false;
        degenerate      = (t()*u() == RT(2)*r()*v());
        if (degenerate) {
            CGAL::Sign discr = (CGAL::Sign)
                            CGAL_NTS sign(u()*u()-RT(4)*r()*w());
            switch (discr) {
                case CGAL::NEGATIVE:
                    empty = true;
                    o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                    break;
                case CGAL::ZERO:
                    empty = false;
                    o = (CGAL::Orientation)(CGAL_NTS sign (r()));
                    break;
                case CGAL::POSITIVE:
                    empty = false;
                    o = CGAL::ZERO;
                    break;
            }
        } else {
            empty = false;
            o = (CGAL::Orientation)(-CGAL_NTS sign (r()));
        }
    } else if (!CGAL_NTS is_zero (s())) {
        trivial         = false;
        degenerate      = (t()*v() == RT(2)*s()*u());
        if (degenerate) {
            CGAL::Sign discr = (CGAL::Sign)
                            CGAL_NTS sign(v()*v()-RT(4)*s()*w());
            switch (discr) {
                case CGAL::NEGATIVE:
                    empty = true;
                    o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                    break;
                case CGAL::ZERO:
                    empty = false;
                    o = (CGAL::Orientation)(CGAL_NTS sign (s()));
                    break;
                case CGAL::POSITIVE:
                    empty = false;
                    o = CGAL::ZERO;
                    break;
            }
        } else {
            empty = false;
            o = (CGAL::Orientation)(-CGAL_NTS sign (s()));
        }
    } else { // r=0, s=0 
        degenerate      = true;
        bool uv_zero    =    CGAL_NTS is_zero (u())
                          && CGAL_NTS is_zero (v());
        trivial         = uv_zero && CGAL_NTS is_zero (w());
        empty           = uv_zero && !trivial;
        if (empty)
            o = (CGAL::Orientation)(CGAL_NTS sign (w()));
        else if (trivial)
            o = CGAL::POSITIVE;
        else
            o = CGAL::ZERO;
    }

@end

@macro <analyse parabola, Cartesian case> = @begin
    if (!CGAL_NTS is_zero (r())) {
        trivial         = false;
        degenerate      = (t()*u() == FT(2)*r()*v());
        if (degenerate) {
            CGAL::Sign discr = (CGAL::Sign)
                            CGAL_NTS sign(u()*u()-FT(4)*r()*w());
            switch (discr) {
                case CGAL::NEGATIVE:
                    empty = true;
                    o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                    break;
                case CGAL::ZERO:
                    empty = false;
                    o = (CGAL::Orientation)(CGAL_NTS sign (r()));
                    break;
                case CGAL::POSITIVE:
                    empty = false;
                    o = CGAL::ZERO;
                    break;
            }
        } else {
            empty = false;
            o = (CGAL::Orientation)(-CGAL_NTS sign (r()));
        }
    } else if (!CGAL_NTS is_zero (s())) {
        trivial         = false;
        degenerate      = (t()*v() == FT(2)*s()*u());
        if (degenerate) {
            CGAL::Sign discr = (CGAL::Sign)
                            CGAL_NTS sign(v()*v()-FT(4)*s()*w());
            switch (discr) {
                case CGAL::NEGATIVE:
                    empty = true;
                    o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                    break;
                case CGAL::ZERO:
                    empty = false;
                    o = (CGAL::Orientation)(CGAL_NTS sign (s()));
                    break;
                case CGAL::POSITIVE:
                    empty = false;
                    o = CGAL::ZERO;
                    break;
            }
        } else {
            empty = false;
            o = (CGAL::Orientation)(-CGAL_NTS sign (s()));
        }
    } else { // r=0, s=0 
        degenerate      = true;
        bool uv_zero    =    CGAL_NTS is_zero (u())
                          && CGAL_NTS is_zero (v());
        trivial         = uv_zero && CGAL_NTS is_zero (w());
        empty           = uv_zero && !trivial;
        if (empty)
            o = (CGAL::Orientation)(CGAL_NTS sign (w()));
        else if (trivial)
            o = CGAL::POSITIVE;
        else
            o = CGAL::ZERO;
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Conic evaluation}
@! ---------------------------------------------------------------------------

For $\C_{\r}$ given by $R=(r,s,t,u,v,w)$ and a point $p=(x,y,h)$, the
homogeneous evaluation returns the value 
$\r(p) = rx^2+sy^2+txy+uxh+vyh+wh^2$.

@macro <ConicHPA2 protected member functions> += @begin
    RT evaluate (const PT& p) const
    {
        RT x, y, h;
        dao.get (p, x, y, h);
        return  r()*x*x + s()*y*y + t()*x*y + u()*x*h + v()*y*h + w()*h*h;
    }

@end

The Cartesian version is obtained for $h=1$, i.e. it computes the value
$\r(p) = rx^2+sy^2+txy+ux+vy+w$.

@macro <ConicCPA2 protected member functions> += @begin
    FT evaluate (const PT& p) const
    {
        FT x, y;
        dao.get (p, x, y);
        return r()*x*x + s()*y*y + t()*x*y + u()*x + v()*y + w();
    }

@end


@! ---------------------------------------------------------------------------
@subsection{High-level member functions}
@! ---------------------------------------------------------------------------


@! ---------------------------------------------------------------------------
@subsubsection{Construction}
@! ---------------------------------------------------------------------------

The construction from the representation proceeds by first setting
the vector components and then analysing the conic to obtain the derived
data.

@macro <ConicHPA2 public member functions> += @begin
    ConicHPA2 ( const DA& da = DA()) : dao( da) { }

    ConicHPA2 (RT r, RT s, RT t, RT u, RT v, RT w, const DA& da = DA())
        : dao( da), _r(r), _s(s), _t(t), _u(u), _v(v), _w(w)
    {
        analyse(); 
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    ConicCPA2 ( const DA& da = DA()) : dao( da) { }

    ConicCPA2 (FT r, FT s, FT t, FT u, FT v, FT w, const DA& da = DA())
        : dao( da), _r(r), _s(s), _t(t), _u(u), _v(v), _w(w)
    {
        analyse();
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Data Accessor}
@! ---------------------------------------------------------------------------

@macro <ConicHPA2 public member functions> += @begin
    const DA&  da() const
    {
        return dao;
    }

@end

@macro <ConicCPA2 public member functions> += @begin 
    const DA&  da() const
    {
        return dao;
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{General access}
@! ---------------------------------------------------------------------------

The coordinate access is straightforward.

@macro <ConicHPA2 public member functions> += @begin
    RT r() const { return _r;}
    RT s() const { return _s;}
    RT t() const { return _t;}
    RT u() const { return _u;}
    RT v() const { return _v;}
    RT w() const { return _w;}

@end


@macro <ConicCPA2 public member functions> += @begin
    FT r() const { return _r;}
    FT s() const { return _s;}
    FT t() const { return _t;}
    FT u() const { return _u;}
    FT v() const { return _v;}
    FT w() const { return _w;}

@end

To obtain the center, recall from subsection \ref{symmetry_sec} that it
is given by the point
\[
c = -M^{-1}
\left(\begin{array}{c} u \\ v \end{array}\right) 
= -\frac{1}{\det(M)}\left(\begin{array}{cc} 2s & -t \\ -t & 2r 
\end{array}\right) \left(\begin{array}{c} u \\ v \end{array}\right)
= -\frac{1}{\det(M)}\left(\begin{array}{cc} 2s\cdot u -t\cdot v \\ 
2r\cdot v - t\cdot u\end{array}\right).
\] 

In the homogeneous representation, the value $-\det(M)$ serves as the
$h$-component of the center, in the Cartesian representation, we divide
by it.

@macro <ConicHPA2 public member functions> += @begin
    PT center () const
    {
        CGAL_optimisation_precondition (type != PARABOLA);
        PT p;
        RT two = RT(2);
        dao.set( p, two*s()*u() - t()*v(), two*r()*v() - t()*u(), -det());
        return p;
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    PT center () const
    {
        CGAL_optimisation_precondition (type != PARABOLA);
        PT p;
        FT two = FT(2);
        FT div = -det();
        dao.set( p, (two*s()*u() - t()*v()) / div,
                    (two*r()*v() - t()*u()) / div);
        return p;
    }

@end


@! ---------------------------------------------------------------------------
@subsubsection{Type Related Access}
@! ---------------------------------------------------------------------------

Because the conic stores its type and degeneracy information, this is
straightforward.

@macro <ConicHPA2 public member functions> += @begin
    Conic_type conic_type () const
    {
        return type;
    }
            
    bool is_hyperbola () const
    {
        return (type == HYPERBOLA);
    }

    bool is_parabola () const
    {
        return (type == PARABOLA);
    }

    bool is_ellipse () const
    {
        return (type == ELLIPSE);
    }

    bool is_empty () const
    {
        return empty;
    }

    bool is_trivial () const
    {
        return trivial;
    }

    bool is_degenerate () const
    {
        return degenerate;
    }

@end


@macro <ConicCPA2 public member functions> += @begin
    Conic_type conic_type () const
    {
        return type;
    }
            
    bool is_hyperbola () const
    {
        return (type == HYPERBOLA);
    }

    bool is_parabola () const
    {
        return (type == PARABOLA);
    }

    bool is_ellipse () const
    {
        return (type == ELLIPSE);
    }

    bool is_empty () const
    {
        return empty;
    }

    bool is_trivial () const
    {
        return trivial;
    }

    bool is_degenerate () const
    {
        return degenerate;
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Orientation Related Access}
@! ---------------------------------------------------------------------------

@macro <ConicHPA2 public member functions> += @begin
    CGAL::Orientation orientation () const
    {
        return o;
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    CGAL::Orientation orientation () const
    {
        return o;
    }

@end


The orientation queries just evaluate the conic at the given point,
using the private method @prg{evaluate}. Recall that $p$ is in the
positive resp. negative side iff $\r(p)>0$ resp. $\r(p)<0$.  

@macro <ConicHPA2 public member functions> += @begin
    CGAL::Oriented_side oriented_side (const PT& p) const
    {
        return (CGAL::Oriented_side)(CGAL_NTS sign (evaluate (p)));
    }

    bool has_on_positive_side (const PT& p) const
    {    
        return (CGAL_NTS is_positive (evaluate(p)));
    }

    bool has_on_negative_side (const PT& p) const
    {
        return (CGAL_NTS is_negative (evaluate(p)));
    }

    bool has_on_boundary (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }

    bool has_on (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    CGAL::Oriented_side oriented_side (const PT& p) const
    {
        return (CGAL::Oriented_side)(CGAL_NTS sign (evaluate (p)));
    }

    bool has_on_positive_side (const PT& p) const
    {    
        return (CGAL_NTS is_positive (evaluate(p)));
    }

    bool has_on_negative_side (const PT& p) const
    {
        return (CGAL_NTS is_negative (evaluate(p)));
    }

    bool has_on_boundary (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }

    bool has_on (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }


@end

Then we have the convex side queries. Under nonzero orientation,
the side is determined by a conic evaluation. If the orientation is
zero, we know that the non-convex side is empty, see subsection
\ref{orientation_sec}.

@macro <ConicHPA2 public member functions> += @begin
    Convex_side convex_side (const PT& p) const 
    {
        switch (o) {
        case CGAL::POSITIVE:
            return (Convex_side)( CGAL_NTS sign (evaluate (p)));
        case CGAL::NEGATIVE:
            return (Convex_side)(-CGAL_NTS sign (evaluate (p)));
        case CGAL::ZERO:
            return (Convex_side)( CGAL_NTS sign (CGAL_NTS abs (evaluate(p))));
        }
        // keeps g++ happy
        return( Convex_side( 0));
    }

    bool has_on_convex_side (const PT& p) const
    {
        return (convex_side (p) == ON_CONVEX_SIDE);
    } 

    bool has_on_nonconvex_side (const PT& p) const
    {
        return (convex_side (p) == ON_NONCONVEX_SIDE);
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    Convex_side convex_side (const PT& p) const 
    {
        switch (o) {
        case CGAL::POSITIVE:
            return (Convex_side)( CGAL_NTS sign (evaluate (p)));
        case CGAL::NEGATIVE:
            return (Convex_side)(-CGAL_NTS sign (evaluate (p)));
        case CGAL::ZERO:
            return (Convex_side)( CGAL_NTS sign (CGAL_NTS abs (evaluate(p))));
        }
        // keeps g++ happy
        return( Convex_side( 0));
    }

    bool has_on_convex_side (const PT& p) const
    {
        return (convex_side (p) == ON_CONVEX_SIDE);
    } 

    bool has_on_nonconvex_side (const PT& p) const
    {
        return (convex_side (p) == ON_NONCONVEX_SIDE);
    }

@end

@! ---------------------------------------------------------------------------
@subsection{Comparison methods}
@! ---------------------------------------------------------------------------

We provide tests for equality and inequality of two conics.

@macro <ConicHPA2 public member functions> += @begin
    bool operator == ( const ConicHPA2<_PT,_DA>& c) const
    {
        // find coefficient != 0
        RT  factor1;
        if ( ! CGAL_NTS is_zero( r())) factor1 = r(); else
        if ( ! CGAL_NTS is_zero( s())) factor1 = s(); else
        if ( ! CGAL_NTS is_zero( t())) factor1 = t(); else
        if ( ! CGAL_NTS is_zero( u())) factor1 = u(); else
        if ( ! CGAL_NTS is_zero( v())) factor1 = v(); else
        if ( ! CGAL_NTS is_zero( w())) factor1 = w(); else
        CGAL_optimisation_assertion_msg( false, "all coefficients zero");

        // find coefficient != 0
        RT  factor2;
        if ( ! CGAL_NTS is_zero( c.r())) factor2 = c.r(); else
        if ( ! CGAL_NTS is_zero( c.s())) factor2 = c.s(); else
        if ( ! CGAL_NTS is_zero( c.t())) factor2 = c.t(); else
        if ( ! CGAL_NTS is_zero( c.u())) factor2 = c.u(); else
        if ( ! CGAL_NTS is_zero( c.v())) factor2 = c.v(); else
        if ( ! CGAL_NTS is_zero( c.w())) factor2 = c.w(); else
        CGAL_optimisation_assertion_msg( false, "all coefficients zero");

        return(    ( r()*factor2 == c.r()*factor1)
                && ( s()*factor2 == c.s()*factor1)
                && ( t()*factor2 == c.t()*factor1)
                && ( u()*factor2 == c.u()*factor1)
                && ( v()*factor2 == c.v()*factor1)
                && ( w()*factor2 == c.w()*factor1));
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    bool operator == ( const ConicCPA2<_PT,_DA>& c) const
    {
        // find coefficient != 0
        FT  factor1;
        if ( ! CGAL_NTS is_zero( r())) factor1 = r(); else
        if ( ! CGAL_NTS is_zero( s())) factor1 = s(); else
        if ( ! CGAL_NTS is_zero( t())) factor1 = t(); else
        if ( ! CGAL_NTS is_zero( u())) factor1 = u(); else
        if ( ! CGAL_NTS is_zero( v())) factor1 = v(); else
        if ( ! CGAL_NTS is_zero( w())) factor1 = w(); else
        CGAL_optimisation_assertion_msg( false, "all coefficients zero");

        // find coefficient != 0
        FT  factor2;
        if ( ! CGAL_NTS is_zero( c.r())) factor2 = c.r(); else
        if ( ! CGAL_NTS is_zero( c.s())) factor2 = c.s(); else
        if ( ! CGAL_NTS is_zero( c.t())) factor2 = c.t(); else
        if ( ! CGAL_NTS is_zero( c.u())) factor2 = c.u(); else
        if ( ! CGAL_NTS is_zero( c.v())) factor2 = c.v(); else
        if ( ! CGAL_NTS is_zero( c.w())) factor2 = c.w(); else
        CGAL_optimisation_assertion_msg( false, "all coefficients zero");

        return(    ( r()*factor2 == c.r()*factor1)
                && ( s()*factor2 == c.s()*factor1)
                && ( t()*factor2 == c.t()*factor1)
                && ( u()*factor2 == c.u()*factor1)
                && ( v()*factor2 == c.v()*factor1)
                && ( w()*factor2 == c.w()*factor1));
    }

@end

@! ---------------------------------------------------------------------------
@subsubsection{Private Methods} 
@! ---------------------------------------------------------------------------

A main difference between the public and the private set methods 
is that the private ones do not analyse the conic. If a conic
is constructed in a sequence of calls to set functions, it is more 
efficient to do an analysis only once for the final result, rather than
analysing any intermediate result. Therefore, the private set functions
also do not allow to specify an orientation, because in order to orient the
conic, an analysis would be necessary. Under this scheme, caution is in place,
of course: a conic constructed from a private set method is not fully
initialized, and calls to type and orientation related access functions
return undefined results. @prg{friend} classes need to care for this, 
whenever they call a private set method. 

@! ---------------------------------------------------------------------------
\paragraph{Linear combination of conics} ~
@! ---------------------------------------------------------------------------

@macro<ConicHPA2 private member functions> += @begin
    void 
    set_linear_combination (const RT& a1, const ConicHPA2<PT,DA>& c1,
                            const RT& a2, const ConicHPA2<PT,DA>& c2)
    {
        _r = a1 * c1.r() + a2 * c2.r();
        _s = a1 * c1.s() + a2 * c2.s();
        _t = a1 * c1.t() + a2 * c2.t();
        _u = a1 * c1.u() + a2 * c2.u();
        _v = a1 * c1.v() + a2 * c2.v();
        _w = a1 * c1.w() + a2 * c2.w();
    }

@end

@macro<ConicCPA2 private member functions> += @begin
    void 
    set_linear_combination (const FT& a1, const ConicCPA2<PT,DA>& c1,
                            const FT& a2, const ConicCPA2<PT,DA>& c2)
    {
        _r = a1 * c1.r() + a2 * c2.r();
        _s = a1 * c1.s() + a2 * c2.s();
        _t = a1 * c1.t() + a2 * c2.t();
        _u = a1 * c1.u() + a2 * c2.u();
        _v = a1 * c1.v() + a2 * c2.v();
        _w = a1 * c1.w() + a2 * c2.w();
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Two pairs of lines through four points} ~
@!----------------------------------------------------------------------------

Here is the method to get two line-pairs from four given points. It just
brings the points into (counter)clockwise order and then calls the
@prg{set_linepair} method for their two conic arguments, supplying the
points in the right order. The reordering needs access to the 
orientations of certain point triples. For points $p_i=(x_i,y_i,h_i)$,
$i=1\ldots3$, this orientation is given by the sign of the determinant
$$\det \left(
\begin{array}{ccc}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
h_1 & h_2 & h_3 
\end{array}\right),$$
see Lemma \ref{det_lemma} below. 
The following macros define this determinant (the Cartesian version is
obtained by setting $h_i=1,i=1\ldots 3$). 

@macro <h_orientation>(3) many = @begin
    (CGAL::Orientation)(CGAL_NTS sign
      (-h@1*x@3*y@2+h@3*x@1*y@2
       +h@1*x@2*y@3-h@2*x@1*y@3
       +h@2*x@3*y@1-h@3*x@2*y@1))
@end

@macro <c_orientation>(3) many = @begin
    (CGAL::Orientation)(CGAL_NTS sign
      (-x@3*y@2+x@1*y@2
       +x@2*y@3-x@1*y@3
       +x@3*y@1-x@2*y@1))
@end

@macro<ConicHPA2 private member functions> += @begin   
    static void set_two_linepairs (const PT& p1, 
                                   const PT& p2,
                                   const PT& p3,
                                   const PT& p4,
                                   ConicHPA2<PT,DA>& pair1, 
                                   ConicHPA2<PT,DA>& pair2)
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3, x4, y4, h4;
        const DA& da = pair1.da();
        da.get (p1, x1, y1, h1);
        da.get (p2, x2, y2, h2);
        da.get (p3, x3, y3, h3);
        da.get (p4, x4, y4, h4);

        CGAL::Orientation side1_24 = @<h_orientation>("2", "4", "1"),
                         side3_24 = @<h_orientation>("2", "4", "3");
        if (side1_24 != side3_24) {
            // (counter)clockwise order
            pair1.set_linepair (p1, p2, p3, p4);
            pair2.set_linepair (p2, p3, p4, p1);
        } else {
            CGAL::Orientation side1_32 = @<h_orientation>("3", "2", "1");
            if (side1_32 != side3_24) { 
                // p1, p2 need to be swapped
                pair1.set_linepair (p2, p1, p3, p4);
                pair2.set_linepair (p1, p3, p4, p2); 
            } else {
                // p2, p3 need to be swapped
                pair1.set_linepair (p1, p3, p2, p4);
                pair2.set_linepair (p3, p2, p4, p1);
            }
        }
    }

@end


@macro<ConicCPA2 private member functions> += @begin
    static void set_two_linepairs (const PT& p1, 
                                   const PT& p2,
                                   const PT& p3,
                                   const PT& p4,
                                   ConicCPA2<PT,DA>& pair1, 
                                   ConicCPA2<PT,DA>& pair2)
    {
        FT x1, y1, x2, y2, x3, y3, x4, y4;      
        const DA& da = pair1.da();
        da.get (p1, x1, y1);
        da.get (p2, x2, y2);
        da.get (p3, x3, y3);
        da.get (p4, x4, y4);

        CGAL::Orientation side1_24 = @<c_orientation>("2", "4", "1"),
                         side3_24 = @<c_orientation>("2", "4", "3");
        if (side1_24 != side3_24) {
            // (counter)clockwise order
            pair1.set_linepair (p1, p2, p3, p4);
            pair2.set_linepair (p2, p3, p4, p1);
        } else {
            CGAL::Orientation side1_32 = @<c_orientation>("3", "2", "1");
            if (side1_32 != side3_24) { 
                // p1, p2 need to be swapped
                pair1.set_linepair (p2, p1, p3, p4);
                pair2.set_linepair (p1, p3, p4, p2); 
            } else {
                // p2, p3 need to be swapped
                pair1.set_linepair (p1, p3, p2, p4);
                pair2.set_linepair (p3, p2, p4, p1);
            }
        }
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Some ellipse from two pairs of lines} ~
@!----------------------------------------------------------------------------

Assuming that we have constructed two line-pairs $\C_{\r_1}, \C_{\r_2}$ 
using the method @prg{set_two_linepairs} with points $p_1,p_2,p_3,p_4$
in convex position, an ellipse through the points can be obtained as a 
linear combination $\C_{\r}$, $\r = \lambda \r_1 + \mu \r_2$, where 
\begin{eqnarray*}
\lambda &=& \det({\r_2}) - 2(r_1s_2+r_2s_1) + t_1t_2, \\
\mu &=& \det({\r_1}) - 2(r_1s_2+r_2s_1) + t_1t_2,
\end{eqnarray*} 
$r_i,\ldots,w_i$ the components of $\r_i,i=1\ldots 2$.

@macro<ConicHPA2 private member functions> += @begin
    void set_ellipse (const ConicHPA2<PT,DA>& pair1, 
                      const ConicHPA2<PT,DA>& pair2)
    {
        RT b = RT(2) * (pair1.r() * pair2.s() + pair1.s() * pair2.r()) -
               pair1.t() * pair2.t();
        set_linear_combination (pair2.det()-b, pair1,
                                pair1.det()-b, pair2);  
    }

@end

@macro<ConicCPA2 private member functions> += @begin
    void set_ellipse (const ConicCPA2<PT,DA>& pair1, 
                      const ConicCPA2<PT,DA>& pair2)
    {
        FT b = FT(2) * (pair1.r() * pair2.s() + pair1.s() * pair2.r()) -
               pair1.t() * pair2.t();
        set_linear_combination (pair2.det()-b, pair1,
                                pair1.det()-b, pair2);
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Conic from two conics and a point} ~
@!----------------------------------------------------------------------------

If $\C_{\r_1}, \C_{\r_2}$ are two conics, $p$ some point, it is easy to see 
that the conic $\C_{\r}$ with $$\r = \r_2(p) \r_1 - \r_1(p) \r_2$$ is
a conic containing the set $(\C_{\r_1} \cap \C_{\r_2}) \cup \{p\}.$ Exactly
this conic is constructed here. 

@macro<ConicHPA2 private member functions> += @begin
    void set (const ConicHPA2<PT,DA>& c1,
              const ConicHPA2<PT,DA>& c2,
              const PT& p)
    {
        set_linear_combination (c2.evaluate(p), c1, -c1.evaluate(p), c2);
    }

@end

@macro<ConicCPA2 private member functions> += @begin
    void set (const ConicCPA2<PT,DA>& c1,
              const ConicCPA2<PT,DA>& c2,
              const PT& p)
    {
        set_linear_combination (c2.evaluate(p), c1, -c1.evaluate(p), c2);
    }

@end

@!----------------------------------------------------------------------------
\paragraph{Volume derivative of an ellipse} ~
@!----------------------------------------------------------------------------

Let $r(\tau),s(\tau), t(\tau), u(\tau), v(\tau), w(\tau)$ denote the
parameters of $\E(\tau)$. Omitting the parameter $\tau$ for the sake of
readability, we have
\begin{eqnarray*}
r &=& r_0 + \tau \d r, \\
s &=& s_0 + \tau \d s, \\
t &=& t_0 + \tau \d t, \\
u &=& u_0 + \tau \d u, \\
v &=& v_0 + \tau \d v, \\
w &=& w_0 + \tau \d w,
\end{eqnarray*}
$r_0,\ldots,w_0$ the representation of $\E$, $\d r,\ldots, \d w$ the 
given values $\d\r$.

Recall that $$\Vol(\E(\tau)) = \pi / \sqrt{\det(M / (2w - c^TMc))}.$$ This
implies
\[
\sgn\left( \frac{\d}{\d \tau} \left. \Vol(\E(\tau))\right|_{\tau=0}\right) =
-\sgn \left(\frac{\d}{\d \tau} \left. \det(M / (2w - c^TMc))
\right|_{\tau=0}\right).
\]
We have $\det(M / (2w - c^TMc)) = d / z^2$,
where
\begin{eqnarray*}
d &=& \det(M) = 4rs -t^2, \\
z &=& 2w - c^TMc = 2 \left(w - \frac{1}{d} (u^2 s - uvt + v^2r)\right),
\end{eqnarray*}
see also (\ref{z}). 
Hence
$$\frac{d}{z^2} = \frac{d^3}{4q^2}, \quad q = dw - u^2 s + uvt - v^2r.$$

It follows that 
\[
   \sgn\left(\frac{\d}{\d \tau} \Vol(\E(\tau))\right)
   = -\sgn\left(\frac{d^2(3d'q-2dq')}{4q^3}\right)
   = -\sgn \left( 3d'q-2dq'\right) \sgn (q).
\]

Consider the term $3d'q-2dq'=0$. We have

\begin{eqnarray*}
d &=& a_2\tau^2 + a_1\tau + a_0, \\
q &=& b_3\tau^3 + b_2\tau^2 + b_1\tau + b,
\end{eqnarray*}
where
\begin{eqnarray*}
a_2 &=& 4\d r\d s - \d t^2, \\
a_1 &=& 4r_0\d s + 4s_0\d r - 2t_0\d t, \\
a_0 &=& 4r_0s_0-t_0^2, \\
b_3 &=& (4\d r\d s - \d t^2)\d w - \d u^2 \d s - \d v^2 \d r+\d t \d u \d v, \\
b_2 &=& (4r_0\d s+4\d rs_0-2t_0\d t)\d w+(4\d r\d s-\d t^2)w_0-2u_0\d u\d s \\
& & -\d u^2s_0-2 v_0\d v\d r-\d v^2r_0+(u_0\d v+\d uv_0)\d t+\d u\d vt_0, \\
b_1 &=& (4r_0s_0-t_0^2)\d w+(4r_0\d s+4\d rs_0-2t_0\d t)w_0-u_0^2\d s \\
& & -2u_0\d u s_0-v_0^2\d r-2v_0\d vr_0+u_0v_0\d t+(u_0\d v+\d uv_0)t_0, \\
b_0 &=& (4r_0s_0-t_0^2)w_0-u_0^2s_0-v_0^2r_0+u_0v_0t_0
\end {eqnarray*}

Furthermore, 
$$3d'q-2dq' = c_3\tau^3 + c_2\tau^2 + c_1\tau + c_0,$$
with
\begin{eqnarray*}
c_3 &=& -3a_1b_3 + 2a_2b_2, \\
c_2 &=& -6a_0b_3 - a_1b_2 + 4a_2b_1, \\
c_1 &=& -4a_0b_2 + a_1b_1 + 6a_2b_0, \\
c_0 &=& -2a_0b_1 + 3a_1b_0.
\end{eqnarray*}

The desired sign of the volume derivative is obtained by
evaluating the expression $$-\sgn(3d'q-2dq')\sgn(q)$$ at $\tau=0$, 
where we get $3d'q-2dq' = c_0.$
Moreover, since $q = \det(M)(2w-c^TMc)/2$ with $\det(M)>0$, the sign of
$q$ is positive if and only if $2w-c^TMc$ is positive, equivalently
if the ellipse $\E(0)$ is in positive orientation (see also the description
of the method @prg{analyse}). This means, in case of positive orientation,
the sign of $-\sgn(3d'q-2dq')\sgn(q)$ equals the sign of $-c_0$, otherwise,
we return the sign of $c_0$. To compute $c_0$, we only need the values 
$a_1,a_0,b_1,b_0$ from above.

@macro<ConicHPA2 private member functions> += @begin
    CGAL::Sign vol_derivative (RT dr, RT ds, RT dt,
                              RT du, RT dv, RT dw) const
    {
        RT a1 = RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt,
           a0 = RT(4)*r()*s()-t()*t(),
           b1 = (RT(4)*r()*s()-t()*t())*dw+(RT(4)*r()*ds+RT(4)*dr*s()-
                RT(2)*t()*dt)*w()-u()*u()*ds -
                RT(2)*u()*du*s()-v()*v()*dr-RT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (RT(4)*r()*s()-t()*t())*w()
                -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c0 = -RT(2)*a0*b1 + RT(3)*a1*b0;

        return CGAL::Sign (-CGAL_NTS sign (c0)*o);
    }

@end

@macro<ConicCPA2 private member functions> += @begin
    CGAL::Sign vol_derivative (FT dr, FT ds, FT dt,
                              FT du, FT dv, FT dw) const
    {
        FT a1 = FT(4)*r()*ds+FT(4)*dr*s()-FT(2)*t()*dt,
           a0 = FT(4)*r()*s()-t()*t(),
           b1 = (FT(4)*r()*s()-t()*t())*dw+(FT(4)*r()*ds+FT(4)*dr*s()-
                FT(2)*t()*dt)*w()-u()*u()*ds -
                FT(2)*u()*du*s()-v()*v()*dr-FT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (FT(4)*r()*s()-t()*t())*w()
          -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c0 = -FT(2)*a0*b1 + FT(3)*a1*b0;

        return CGAL::Sign (-CGAL_NTS sign (c0)*o);
    }

@end

To find the value $\tau^*$ such that $E(\tau^*)$ is the ellipse of smallest
volume, we apply the Cardano formula to find the roots of the polynomial 
$p(\tau)=c_3\tau^3 + c_2\tau^2 + c_1\tau + c_0 = 0$. This is done 
by the function @prg{solve_cubic} which returns all (at most three) 
real roots. We then select the one which leads to largest (positive)
volume.

Here is an outline of the Cardano formula, for the polynomial 
$p(\tau)$, assuming that $c_3\neq 0$.
\begin{enumerate}
\item Divide by $c_3$ to normalize the equation, leading to
$$p(\tau) = \tau^3 + \gamma_2 \tau^2 + \gamma_1 \tau + \gamma_0,\quad 
\gamma_i := c_i/c_3, ~i=0,\ldots,2.$$ 
\item Eliminate the quadratic term by substituting $\tau := x - \gamma_2/3$.
This leads to
$$p(x) = x^3 + a x + b, \quad a = \gamma_1 - \frac{\gamma_2^2}{3},~ 
b = \frac{2}{27}\gamma_2^3 - \frac{1}{3}\gamma_1\gamma_2 + \gamma_0.$$ 
If $a=0$, $p$ has only one real root, namely $x_1 = \sqrt[3]{-b}$, so
let's assume $a\neq 0$.  
\item Define $$D := (a/3)^3 + (b/2)^2$$ and let $u$ be any number such
that $$u^3 = -b/2 + \sqrt{D}.$$ 
Note that $u^3$ is a solution of the quadratic equation
$$x^2 + bx - \left(\frac{a}{3}\right)^3 = 0.$$
This means, if $a\neq 0$, then also $u\neq 0$. 
If the discriminant $D$ is negative, $u$ is a complex number,
otherwise it must be chosen as the real number 
$$u := \sqrt[3]{-b/2 + \sqrt{D}}.$$
\item Let $u$ be of the form $u = u_R - u_I \I$ (we possibly have $u_I=0$).
$p$ has always a real root, given by 
$$x_1 = u_R \left(1 - \frac{a}{3\|u\|^2}\right).$$ 
\item If $D > 0$, the two other roots are complex. If 
$D \leq 0$, $p$ has two more real roots, 
given by 
\begin{eqnarray*}
x_2 & = & -\frac{1}{2}\left(1 - \frac{a}{3\|u\|^2}\right)(u_R-u_I\sqrt{3}),\\
x_3 & = & -\frac{1}{2}\left(1 - \frac{a}{3\|u\|^2}\right)(u_R+u_I\sqrt{3}).
\end{eqnarray*}
If $D=0$, we have $u_I=0$ and these two roots coincide. 
\end{enumerate} 
It follows that if $D\geq 0$, the roots of $p$ can be found by using only
$\sqrt{~}$ and $\sqrt[3]{~}$ operations over the real numbers. If $D<0$, 
it can be shown that this is not possible, and we need to approximate 
the values $u_R,u_I$ in some other way. For this, we need to solve the 
equation $$u^3 = C := -b/2 + i\sqrt{-D}$$
for $u$. Expressing $C$ in polar coordinates $(r,\phi)$ we get
\begin{eqnarray*}
r         &=& \|C\| = \sqrt{b^2/4 -D} = \sqrt{-(a/3)^3}, \\
\cos \phi &=& -\frac{b/2}{\|C\|}, 
\end{eqnarray*}
where $0\leq \phi < \pi$ because of $\sqrt{-D}>0$.
Therefore, a possible choice for $u$ is $u=(r', \phi')$ with
\begin{eqnarray*}
r'         &=& \sqrt{-a/3}, \\
\phi' &=&  \acos \left(-\frac{b/2}{\|C\|}\right)/3.
\end{eqnarray*}
This gives
\begin{eqnarray*}
u_R &=& r' \cos \phi', \\
u_I &=& r' \sin \phi',
\end{eqnarray*} 

Note that this implies $\|u\|^2 = -a/3$, showing that the roots of $p(x)$ 
assume the form 
\begin{eqnarray*}
x_1 &=& 2u_R, \\
x_2 &=& -(u_R -u_I\sqrt{3}),\\
x_3 &=& -(u_R +u_I\sqrt{3}).
\end{eqnarray*}
If $u^3$ is real, however, this is not true, and we really need to evaluate
the factor $$\alpha := \left(1 - \frac{a}{3\|u\|^2}\right)
= \left(1 - \frac{a}{3u^2}\right)$$ in this case to obtain
the roots of $p(x)$. In any case, we originally wanted to have the roots
of $p(\tau)$. Using the substitution formula $\tau=x-\gamma_2/3$, these 
roots are given by $\tau_i = x_i -\gamma_2/3, i=1,\ldots,3$.

The function @prg{solve_cubic} returns the number of distinct 
real roots of $p(\tau) = c_3\tau^3+c_2\tau^2+c_1\tau+c_0=0$ and stores 
them consecutively in @prg{r1}, @prg{r2} and @prg{r3}. Precondition is that
$p(\tau)$ is not a constant function.

@macro <function solve_cubic> zero = @begin
    template < class NT >
    int solve_cubic (NT c3, NT c2, NT c1, NT c0, 
                     NT& r1, NT& r2, NT& r3) 
    {
        if (c3 == 0.0) {
            // quadratic equation
            if (c2 == 0) {
                // linear equation
                CGAL_optimisation_precondition (c1 != 0);
                r1 = -c0/c1;
                return 1;
            }
            NT D = c1*c1-4*c2*c0;
            if (D < 0.0) 
                // only complex roots
                return 0;
            if (D == 0.0) {
                // one real root
                r1 = -c1/(2.0*c2);
                return 1;
            }
            // two real roots
            r1 = (-c1 + CGAL::sqrt(D))/(2.0*c2);
            r2 = (-c1 - CGAL::sqrt(D))/(2.0*c2);
            return 2;
        }

        // cubic equation
        // define the gamma_i
        NT g2 = c2/c3, 
           g1 = c1/c3, 
           g0 = c0/c3;

        // define a, b
        NT a = g1 - g2*g2/3.0,
           b = 2.0*g2*g2*g2/27.0 - g1*g2/3.0 + g0;

        if (a == 0) {
            // one real root
            /***** r1 = cbrt(-b) - g2/3.0; *****/
            r1 = exp(log(-b)/3.0) - g2/3.0;
            return 1;
        }

        // define D
        NT D  = a*a*a/27.0 + b*b/4.0;
        if (D >= 0.0) {
            // real case
            /***** NT u = cbrt(-b/2.0 + CGAL::sqrt(D)), *****/
            NT u = exp(log(-b/2.0 + CGAL::sqrt(D))),
                   alpha = 1.0 - a/(3.0*u*u);
            if (D == 0) {
                // two distinct real roots
                r1 =  u*alpha - g2/3.0;
                r2 =  -0.5*alpha*u - g2/3.0;
                return 2;
            }
            // one real root
            r1 = u*alpha - g2/3.0;
            return 1;
        }
        // complex case
        NT r_prime   = CGAL::sqrt(-a/3),
           phi_prime = acos (-b/(2.0*r_prime*r_prime*r_prime))/3.0,
           u_R       = r_prime * cos (phi_prime),
           u_I       = r_prime * sin (phi_prime);
        // three distinct real roots
        r1 = 2.0*u_R - g2/3.0;
        r2 = -u_R + u_I*CGAL::sqrt(3.0) - g2/3.0;
        r3 = -u_R - u_I*CGAL::sqrt(3.0) - g2/3.0;
        return 3;
    }
                    
@end

Here comes the actual computation of the volume minimum. To this end,
we compute the coefficients $c_3,c_2,c_1,c_0$ and find the roots
of $p(\tau)$. Among the roots we then select the one which leads 
to the smallest volume,
equivalently to the largest value of $\det(M/(2w-c^TMc))=d^3/4q^2$.
The coefficients of $d$ and $q$ have been computed before in order 
to find the roots, so we can directly evaluate the determinant, using
@prg{double} approximations of the coefficients.


@macro<ConicHPA2 private member functions> += @begin
    double vol_minimum (RT dr, RT ds, RT dt, RT du, RT dv, RT dw) const
    {
        RT a2 = RT(4)*dr*ds-dt*dt,
           a1 = RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt,
           a0 = RT(4)*r()*s()-t()*t(),
           b3 = (RT(4)*dr*ds-dt*dt)*dw-du*du*ds-dv*dv*dr+du*dv*dt,
           b2 = (RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt)*dw+
                (RT(4)*dr*ds-dt*dt)*w()-RT(2)*u()*du*ds-du*du*s()-
                RT(2)*v()*dv*dr-dv*dv*r()+(u()*dv+du*v())*dt+du*dv*t(),
           b1 = (RT(4)*r()*s()-t()*t())*dw+(RT(4)*r()*ds+RT(4)*dr*s()-
                RT(2)*t()*dt)*w()-u()*u()*ds -
                RT(2)*u()*du*s()-v()*v()*dr-RT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (RT(4)*r()*s()-t()*t())*w()
                -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c3 = -RT(3)*a1*b3 + RT(2)*a2*b2,
           c2 = -RT(6)*a0*b3 - a1*b2 + RT(4)*a2*b1,
           c1 = -RT(4)*a0*b2 + a1*b1 + RT(6)*a2*b0,
           c0 = -RT(2)*a0*b1 + RT(3)*a1*b0;

           if ( CGAL_NTS is_zero( c0)) return 0;// E(0) is the smallest ellipse

           double roots[3];
           int nr_roots = solve_cubic 
                                (CGAL::to_double(c3), CGAL::to_double(c2), 
                                 CGAL::to_double(c1), CGAL::to_double(c0), 
                                 roots[0], roots[1], roots[2]);
           CGAL_optimisation_precondition (nr_roots > 0); // minimum exists
           return best_value (roots, nr_roots,
                                 CGAL::to_double(a2), CGAL::to_double(a1),
                                 CGAL::to_double(a0), CGAL::to_double(b3),
                                 CGAL::to_double(b2), CGAL::to_double(b1),
                                 CGAL::to_double(b0));  
    }

@end


@macro<ConicCPA2 private member functions> += @begin
    double vol_minimum (FT dr, FT ds, FT dt, FT du, FT dv, FT dw) const
    {
        FT a2 = FT(4)*dr*ds-dt*dt,
           a1 = FT(4)*r()*ds+FT(4)*dr*s()-FT(2)*t()*dt,
           a0 = FT(4)*r()*s()-t()*t(),
           b3 = (FT(4)*dr*ds-dt*dt)*dw-du*du*ds-dv*dv*dr+du*dv*dt,
           b2 = (FT(4)*r()*ds+FT(4)*dr*s()-FT(2)*t()*dt)*dw+
                (FT(4)*dr*ds-dt*dt)*w()-FT(2)*u()*du*ds-du*du*s()-
                FT(2)*v()*dv*dr-dv*dv*r()+(u()*dv+du*v())*dt+du*dv*t(),
           b1 = (FT(4)*r()*s()-t()*t())*dw+(FT(4)*r()*ds+FT(4)*dr*s()-
                FT(2)*t()*dt)*w()-u()*u()*ds -
                FT(2)*u()*du*s()-v()*v()*dr-FT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (FT(4)*r()*s()-t()*t())*w()
                -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c3 = -FT(3)*a1*b3 + FT(2)*a2*b2,
           c2 = -FT(6)*a0*b3 - a1*b2 + FT(4)*a2*b1,
           c1 = -FT(4)*a0*b2 + a1*b1 + FT(6)*a2*b0,
           c0 = -FT(2)*a0*b1 + FT(3)*a1*b0;

           if ( CGAL_NTS is_zero( c0)) return 0;// E(0) is the smallest ellipse

           double roots[3];
           int nr_roots = solve_cubic
                                (CGAL::to_double(c3), CGAL::to_double(c2), 
                                 CGAL::to_double(c1), CGAL::to_double(c0), 
                                 roots[0], roots[1], roots[2]);
           CGAL_optimisation_precondition (nr_roots > 0); // minimum exists
           return best_value (roots, nr_roots,
                                 CGAL::to_double(a2), CGAL::to_double(a1),
                                 CGAL::to_double(a0), CGAL::to_double(b3),
                                 CGAL::to_double(b2), CGAL::to_double(b1),
                                 CGAL::to_double(b0));           
    }

@end

The function @prg{best_root} returns the value in its argument array
which leads to the largest determinant $d^3/q^2$. A precondition was that
an ellipse of smallest volume exists, so the largest determinant must 
be positive. 

@macro <function best_value> = @begin
    template < class NT >
    NT best_value (NT *values, int nr_values, 
                   NT a2, NT a1, NT a0,
                   NT b3, NT b2, NT b1, NT b0)
    {
        bool det_positive = false;
        NT d, q, max_det = 0, det, best = -1;
        for (int i=0; i<nr_values; ++i) {
            NT x = values[i];
            d = (a2*x+a1)*x+a0;
            q = ((b3*x+b2)*x+b1)*x+b0;
            det = d*d*d/(q*q);
            if (det > 0.0)
                if (!det_positive || (det > max_det)) {
                    max_det = det;
                    best = x;
                    det_positive = true;
                }
        }
        CGAL_optimisation_precondition (det_positive);
        return best;
    }
   
@end


@! ---------------------------------------------------------------------------
@subsubsection{Public Set Methods} 
@! ---------------------------------------------------------------------------

Here is the set method at coordinate level.

@macro <ConicHPA2 public member functions> += @begin
    void set (RT r_, RT s_, RT t_, RT u_, RT v_, RT w_)
    {
        _r = r_; _s = s_; _t = t_; _u = u_; _v = v_; _w = w_;
        analyse();
     }

@end

@macro <ConicCPA2 public member functions> += @begin
    void set (FT r_, FT s_, FT t_, FT u_, FT v_, FT w_)
    {
        _r = r_; _s = s_; _t = t_; _u = u_; _v = v_; _w = w_;
        analyse();
     }

@end


@! ---------------------------------------------------------------------------
\paragraph{Opposite Conic} ~
@! ---------------------------------------------------------------------------

The method @prg{set_opposite} just flips the representation $\r$ and the
orientation, all other derived data are taken over.  

@macro <ConicHPA2 public member functions> += @begin
    void set_opposite ()
    {
        _r = -r(); _s = -s(); _t = -t(); _u = -u(); _v = -v(); _w = -w();
        o = CGAL::opposite(orientation());
    }

@end

@macro <ConicCPA2 public member functions> += @begin
    void set_opposite ()
    {
        _r = -r(); _s = -s(); _t = -t(); _u = -u(); _v = -v(); _w = -w();
        o = CGAL::opposite(orientation());
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Pair of lines through four points} ~
@! ---------------------------------------------------------------------------
Given $p_1 = (x_1, y_1, h_1)$, 
$p_2 = (x_2, y_2, h_2)$, $p_3 = (x_3, y_3, h_3)$, 
$p_4 = (x_4, y_4, h_4)$, we develop a formula for representing the
pair of lines $\overline{p_1p_2}$, $\overline{p_3p_4}$ (forming a 
degenerate hyperbola) in the form of (\ref{conic_def_hom}). To this end, 
let $p=(x,y,h)$ be any point and define
\begin{equation}
\label{orientation_det}
[p_i, p_j, p] := \det 
\left(
\begin{array}{ccc}
x_i & x_j & x \\
y_1 & y_2 & y \\
h_1 & h_2 & h 
\end{array}\right).
\end{equation}
It is well known that $[p_i,p_j,p]$ records the orientation of the point
triple: let $\ell$ be the oriented line through $p_i$ and $p_j$. Then 
the following holds.
\begin{Lemma}
\label{det_lemma}
\[
p \mbox{~lies} 
\left\{ \begin{array}{c}
\mbox{to the left of} \\
\mbox{on} \\
\mbox{to the right of}
\end{array} \right\} \ell \Leftrightarrow
[p_i,p_j,p] 
\left\{\begin{array}{c}
>0 \\ = 0 \\ <0
\end{array}\right\}.
\]
\end{Lemma}

In particular, $p$ lies on $\overline{p_1p_2} \cup \overline{p_3p_4}$
iff $[p_1,p_2,p][p_3,p_4,p]=0$, and this expression turns out to be
of the form (\ref{conic_def_hom}), where @prg{Maple} gives us the concrete
values of $r,s,t,u,v,w$. Note that we must have $p_1\neq p_2$ and
$p_3\neq p_4$ to obtain reasonable results (and this was a precondition).

\begin{eqnarray*}
r &=& (y_1h_2-h_1y_2)(y_3h_4-h_3y_4), \\
s &=& (h_1x_2-x_1h_2)(h_3x_4-x_3h_4), \\
t &=& (h_1x_2-x_1h_2)(y_3h_4-h_3y_4)+(y_1h_2-h_1y_2)(h_3x_4-x_3h_4), \\
u &=& (-y_1x_2+x_1y_2)(y_3h_4-h_3y_4)+(y_1h_2-h_1y_2)(-y_3x_4+x_3y_4), \\
v &=& (-y_1x_2+x_1y_2)(h_3x_4-x_3h_4)+(h_1x_2-x_1h_2)(-y_3x_4+x_3y_4), \\
w &=& (-y_1x_2+x_1y_2)(-y_3x_4+x_3y_4).
\end{eqnarray*}

@macro <ConicHPA2 public member functions> += @begin
    void set_linepair (const PT& p1, const PT& p2, const PT& p3,
                       const PT& p4, const DA& da = DA()) 
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3, x4, y4, h4;
        da.get (p1, x1, y1, h1);
        da.get (p2, x2, y2, h2);
        da.get (p3, x3, y3, h3);
        da.get (p4, x4, y4, h4);
        
        // precondition: p1 != p2, p3 != p4
        CGAL_optimisation_precondition 
            ( ((x1*h2 != x2*h1) || (y1*h2 != y2*h1)) &&
              ((x3*h4 != x4*h3) || (y3*h4 != y4*h3)) );

        RT h1x2_x1h2 = h1*x2-x1*h2;
        RT h3x4_x3h4 = h3*x4-x3*h4;
        RT y1h2_h1y2 = y1*h2-h1*y2;
        RT y3h4_h3y4 = y3*h4-h3*y4;
        RT x1y2_y1x2 = x1*y2-y1*x2;
        RT x3y4_y3x4 = x3*y4-y3*x4;

        _r = y1h2_h1y2 * y3h4_h3y4;
        _s = h1x2_x1h2 * h3x4_x3h4; 
        _t = h1x2_x1h2 * y3h4_h3y4 + y1h2_h1y2 * h3x4_x3h4;
        _u = x1y2_y1x2 * y3h4_h3y4 + y1h2_h1y2 * x3y4_y3x4;
        _v = x1y2_y1x2 * h3x4_x3h4 + h1x2_x1h2 * x3y4_y3x4;
        _w = x1y2_y1x2 * x3y4_y3x4;

        analyse();
    }

@end

For the Cartesian representation we proceed completely similar,
replacing values $h_1,\ldots,h_4$ by 1.

@macro <ConicCPA2 public member functions> += @begin
    void set_linepair (const PT& p1, const PT& p2, const PT& p3, const PT& p4)
    {
        FT x1, y1, x2, y2, x3, y3, x4, y4;
        dao.get (p1, x1, y1);
        dao.get (p2, x2, y2);
        dao.get (p3, x3, y3);
        dao.get (p4, x4, y4);

        // precondition: p1 != p2, p3 != p4
        CGAL_optimisation_precondition 
            ( ((x1 != x2) || (y1 != y2)) &&
              ((x3 != x4) || (y3 != y4)) );

        FT x2_x1 = x2-x1;
        FT x4_x3 = x4-x3;
        FT y1_y2 = y1-y2;
        FT y3_y4 = y3-y4;
        FT x1y2_y1x2 = x1*y2-y1*x2;
        FT x3y4_y3x4 = x3*y4-y3*x4;

        _r = y1_y2 * y3_y4;
        _s = x2_x1 * x4_x3; 
        _t = x2_x1 * y3_y4 + y1_y2 * x4_x3;
        _u = x1y2_y1x2 * y3_y4 + y1_y2 * x3y4_y3x4;
        _v = x1y2_y1x2 * x4_x3 + x2_x1 * x3y4_y3x4;
        _w = x1y2_y1x2 * x3y4_y3x4;

        analyse();
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Smallest ellipse through three points} ~
@! ---------------------------------------------------------------------------

Given $p_1 = (x_1, y_1, h_1)$, 
$p_2 = (x_2, y_2, h_2)$, $p_3 = (x_3, y_3, h_3)$, we give a formula for
representing the ellipse of smallest volume containing $p_1,p_2,p_3$
in the form of (\ref{conic_def_hom}). For this, we use the following 
well-known formula for this ellipse in case of Cartesian points.
\begin{Lemma}
Let $q_1,q_2,q_3$ be non-collinear Cartesian points. Then the smallest
ellipse containing $q_1,q_2,q_3$ can be written as the set of points
$q=(x,y)$ satisfying
\begin{equation}
\label{ellipse_center_form}
(q-c)^T M (q-c) = 1,
\end{equation}
where $$c = \frac{1}{3}\sum_{i=1}^3q_i, \quad 
M^{-1} = \frac{2}{3} \sum_{i=1}^3 (q_i-c)(q_i-c)^T.$$
\end{Lemma}

To apply this Lemma to points $(p_1,p_2,p_3)$, we define
$q_i = (x_i/h_i, y_i/h_i)$, $i=1,\ldots,3$ and observe that
(\ref{ellipse_center_form}) can be written as
\[
q^T M q - 2q^TMc + c^TMc - 1 = 0,
\]
from which a representation in form of (\ref{conic_cart}) is obtained via
\[
\left( \begin{array}{cc} 2r & t \\ t & 2s \end{array}\right) := M,
\quad \left(\begin{array}{c} u \\ v \end{array}\right) := -Mc, 
\quad 2w := c^TMc - 1.
\]
Using @prg{Maple}, we get the following values for $\r=(r,s,t,u,v,w)$.

\begin{eqnarray*}
r &=& 3(y_1^2h_2^2h_3^2-y_1h_2h_3^2y_2h_1-y_1h_2^2h_3y_3h_1+y_2^2
h_1^2h_3^2-y_2h_1^2h_3y_3h_2+y_3^2h_1^2h_2^2) / d, \\
s &=& 3(x_1^2h_2^2h_3^2-x_1h_2h_3^2x_2h_1-x_1h_2^2h_3x_3h_1+x_2^2
h_1^2h_3^2-x_2h_1^2h_3x_3h_2+x_3^2h_1^2h_2^2) / d, \\
t &=& 3(-2x_1h_2^2h_3^2y_1+x_1h_2h_3^2y_2h_1+x_1h_2^2h_3y_3h_1+
x_2h_1h_3^2y_1h_2\\
&~&-2x_2h_1^2h_3^2y_2+x_2h_1^2h_3y_3h_2+x_3h_1h_2^2
y_1h_3+x_3h_1^2h_2y_2h_3-2x_3h_1^2h_2^2y_3) / d, \\
u &=& -3(y_2^2h_1^2h_3x_3-y_2h_1^2h_3y_3x_2-y_2h_1^2y_3h_2x_3+
y_3^2h_1h_2^2x_1+y_3^2h_1^2h_2x_2+y_1^2h_2h_3^2x_2\\
&~& +y_1^2h_2^2h_3x_3-y_1
h_2h_3^2y_2x_1-y_1h_3^2y_2h_1x_2-y_1h_2^2h_3y_3x_1-y_1h_2^2y_3
h_1x_3+y_2^2h_1h_3^2x_1) / d, \\
v &=& -3(-x_1h_2h_3^2y_1x_2-x_1h_2^2h_3y_1x_3+x_1^2h_2h_3^2y_2-x_1
h_3^2y_2h_1x_2+x_1^2h_2^2h_3y_3-x_1h_2^2y_3h_1x_3\\
&~&+x_2^2h_1
h_3^2y_1-x_2h_1^2h_3y_2x_3+x_2^2h_1^2h_3y_3-x_2h_1^2y_3h_2x_3+x_3^2h_1
h_2^2y_1+x_3^2h_1^2h_2y_2) / d, \\
w &=& 3(x_2h_3x_3h_2y_1^2-x_1h_2h_3x_3y_1y_2-x_2h_1h_3x_3y_1y_2+
x_3^2h_1h_2y_1y_2-x_1h_2h_3x_2y_1y_3+x_2^2h_1h_3y_1y_3\\
&~&-x_2h_1x_3
h_2y_1y_3+x_1h_3x_3h_1y_2^2+x_1^2h_2h_3y_2y_3-x_1h_3x_2h_1y_2
y_3-x_1h_2x_3h_1y_2y_3+x_1h_2x_2h_1y_3^2) / d,
\end{eqnarray*}
where 
$\delta = \left(\det \left(
\begin{array}{ccc}
x_1 & x_2 & x_3 \\
y_1 & y_2 & y_3 \\
h_1 & h_2 & h_3 
\end{array}\right)\right)^2 =
(-h_1x_3y_2+x_1y_2h_3+h_1x_2y_3-x_1y_3h_2+h_2x_3y_1-x_2y_1h_3)^2.$

After pre-computing the values
\[
x_i^2, y_i^2, x_ih_i, y_ih_i, h_i^2,
\] 
for $i=1,\ldots,3$, the components of the vector $\delta\r/3$ are easy
to obtain. Note that this vector is a legal representation of the ellipse
if $\delta\neq 0$. This is the case if and only if $p_1,p_2,p_3$ are
non-collinear, see Lemma \ref{det_lemma}. Moreover, one can show that
the formulas above determine an ellipse of negative orientation, regardless
of the point triple orientation. This means, if the orientation
was positive, we still need to flip the representation.  

@macro<ConicHPA2 public member functions> += @begin
    void set_ellipse (const PT& p1, const PT& p2, const PT& p3)
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3;
        dao.get (p1, x1, y1, h1);
        dao.get (p2, x2, y2, h2);
        dao.get (p3, x3, y3, h3);

        // precondition: p1, p2, p3 not collinear
        RT det = -h1*x3*y2+h3*x1*y2+h1*x2*y3-h2*x1*y3+h2*x3*y1-h3*x2*y1;
        CGAL_optimisation_precondition (!CGAL_NTS is_zero (det));

        RT x1x1 = x1*x1, y1y1 = y1*y1,
           x2x2 = x2*x2, y2y2 = y2*y2,
           x3x3 = x3*x3, y3y3 = y3*y3,  // x_i^2, y_i^2
           x1h1 = x1*h1, y1h1 = y1*h1, 
           x2h2 = x2*h2, y2h2 = y2*h2, 
           x3h3 = x3*h3, y3h3 = y3*h3,  // x_i h_i, y_i h_i
           h1h1 = h1*h1, 
           h2h2 = h2*h2, 
           h3h3 = h3*h3,                // h_i^2
           two = RT(2);                 // 2

        _r = y1y1*h2h2*h3h3 - y1h1*y2h2*h3h3 - y1h1*h2h2*y3h3 + 
             h1h1*y2y2*h3h3 - h1h1*y2h2*y3h3 + h1h1*h2h2*y3y3;

        _s = x1x1*h2h2*h3h3 - x1h1*x2h2*h3h3 - x1h1*h2h2*x3h3 +
             h1h1*x2x2*h3h3 - h1h1*x2h2*x3h3 + h1h1*h2h2*x3x3;

        _t = -two*x1*y1*h2h2*h3h3 + x1h1*y2h2*h3h3 + x1h1*h2h2*y3h3 +
                 y1h1*x2h2*h3h3 -two*h1h1*x2*y2*h3h3 + h1h1*x2h2*y3h3 +
                 y1h1*h2h2*x3h3 + h1h1*y2h2*x3h3 -two*h1h1*h2h2*x3*y3;

        _u = -(h1h1*y2y2*x3h3 - h1h1*x2*y2*y3h3 - h1h1*y2h2*x3*y3 +
                   x1h1*h2h2*y3y3 + h1h1*x2h2*y3y3 +y1y1*x2h2*h3h3 +
                   y1y1*h2h2*x3h3 - x1*y1*y2h2*h3h3 - y1h1*x2*y2*h3h3 -
                   x1*y1*h2h2*y3h3 - y1h1*h2h2*x3*y3 + x1h1*y2y2*h3h3);

        _v = -(h1h1*x2x2*y3h3 - h1h1*x2*y2*x3h3 - h1h1*x2h2*x3*y3 +
                   y1h1*h2h2*x3x3 + h1h1*y2h2*x3x3 +x1x1*y2h2*h3h3 +
                   x1x1*h2h2*y3h3 - x1*y1*x2h2*h3h3 - x1h1*x2*y2*h3h3 -
                   x1*y1*h2h2*x3h3 - x1h1*h2h2*x3*y3 + y1h1*x2x2*h3h3);

        _w = y1y1*x2h2*x3h3 - x1*y1*y2h2*x3h3 - y1h1*x2*y2*x3h3 +
             y1h1*y2h2*x3x3 - x1*y1*x2h2*y3h3 + y1h1*x2x2*y3h3 -
             y1h1*x2h2*x3*y3 + x1h1*y2y2*x3h3 + x1x1*y2h2*y3h3 -
             x1h1*x2*y2*y3h3 - x1h1*y2h2*x3*y3 + x1h1*x2h2*y3y3;

        type = ELLIPSE;
        degenerate = trivial = empty = false;
        o = CGAL::NEGATIVE;
        if (CGAL_NTS is_positive (det)) set_opposite (); 
        
    }

@end

As before, the Cartesian version is obtained by setting 
$h_1,h_2,h_3$ to 1.

@macro<ConicCPA2 public member functions> += @begin
    void set_ellipse (const PT& p1, const PT& p2, const PT& p3)
    {
        FT x1, y1, x2, y2, x3, y3;
        dao.get (p1, x1, y1);
        dao.get (p2, x2, y2);
        dao.get (p3, x3, y3);

        // precondition: p1, p2, p3 not collinear
        FT det = -x3*y2+x1*y2+x2*y3-x1*y3+x3*y1-x2*y1;
        CGAL_optimisation_precondition (!CGAL_NTS is_zero (det));

        FT x1x1 = x1*x1, y1y1 = y1*y1,
           x2x2 = x2*x2, y2y2 = y2*y2,
           x3x3 = x3*x3, y3y3 = y3*y3,  // x_i^2, y_i^2
           two = FT(2);

        _r = y1y1 - y1*y2 - y1*y3 + 
             y2y2 - y2*y3 + y3y3;

        _s = x1x1 - x1*x2 - x1*x3 +
             x2x2 - x2*x3 + x3x3;

        _t = -two*x1*y1 + x1*y2 + x1*y3 +
                 y1*x2 -two*x2*y2 + x2*y3 +
                 y1*x3 + y2*x3 -two*x3*y3;

        _u = -(y2y2*x3 - x2*y2*y3 - y2*x3*y3 +
                   x1*y3y3 + x2*y3y3 + y1y1*x2 +
                   y1y1*x3 - x1*y1*y2 - y1*x2*y2 -
                   x1*y1*y3 - y1*x3*y3 + x1*y2y2);

        _v = -(x2x2*y3 - x2*y2*x3 - x2*x3*y3 +
                   y1*x3x3 + y2*x3x3 + x1x1*y2 +
                   x1x1*y3 - x1*y1*x2 - x1*x2*y2 -
                   x1*y1*x3 - x1*x3*y3 + y1*x2x2);

        _w = y1y1*x2*x3 - x1*y1*y2*x3 - y1*x2*y2*x3 +
             y1*y2*x3x3 - x1*y1*x2*y3 + y1*x2x2*y3 -
             y1*x2*x3*y3 + x1*y2y2*x3 + x1x1*y2*y3 -
             x1*x2*y2*y3 - x1*y2*x3*y3 + x1*x2*y3y3;

        type = ELLIPSE;
        degenerate = trivial = empty = false;
        o = CGAL::NEGATIVE;
        if (CGAL_NTS is_positive (det)) set_opposite();
    }

@end


@! ---------------------------------------------------------------------------
\paragraph{Some ellipse through four points in convex position} ~
@! ---------------------------------------------------------------------------

This method builds on the private method to obtain an ellipse from two
pairs of lines through the four points. For constructing this pair, we
also have a method available.  

@macro<ConicHPA2 public member functions> += @begin
    void set_ellipse (const PT& p1, const PT& p2,
                      const PT& p3, const PT& p4,
                      CGAL::Orientation _o = POSITIVE) 
    {
        ConicHPA2<PT,DA> pair1, pair2;
        set_two_linepairs (p1, p2, p3, p4, pair1, pair2);
        set_ellipse (pair1, pair2);
        analyse();
        if (o != _o) set_opposite();
    }

@end


@macro<ConicCPA2 public member functions> += @begin
    void set_ellipse (const PT& p1, const PT& p2,
                      const PT& p3, const PT& p4,
                      CGAL::Orientation _o = POSITIVE) 
    {
        ConicCPA2<PT,DA> pair1, pair2;
        set_two_linepairs (p1, p2, p3, p4, pair1, pair2);
        set_ellipse (pair1, pair2);
        analyse();
        if (o != _o) set_opposite();
    }

@end

@! ---------------------------------------------------------------------------
\paragraph{Unique conic through five points} ~
@! ---------------------------------------------------------------------------

Using the previously defined methods, we implement the method to compute
the unique nontrivial conic through five given points 
$p_1,p_2,p_3,p_4,p_5$. For this,
we first compute the two conics 
$\C_1 = \overline{p_1p_2} \cup \overline{p_3p_4}$ and
$\C_2 = \overline{p_1p_4} \cup \overline{p_2p_3}$, using the 
@prg{set_linepair} method. This gives two conics having the points
$p_1,p_2,p_3,p_4$ in common. It follows that any linear combination
of them goes through $p_1,p_2,p_3,p_4$ as well. A particular linear 
combination is given by 
\begin{equation}
\label{unique_conic}
\C := \C_2(p_5) \C_1 - \C_1(p_5) \C_2,
\end{equation}
and it has the property that $\C(p_5)=0$, i.e. $\C$ goes through
$p_1,\ldots,p_5$. In case all points are distinct, this is the unique
nontrivial conic through the points. 

@macro<ConicHPA2 public member functions> += @begin
    void set (const PT& p1, const PT& p2, const PT& p3, const PT& p4,
              const PT& p5, CGAL::Orientation _o = POSITIVE)
    {
        ConicHPA2<PT,DA> c1; c1.set_linepair (p1, p2, p3, p4);
        ConicHPA2<PT,DA> c2; c2.set_linepair (p1, p4, p2, p3);
        set_linear_combination (c2.evaluate (p5), c1, 
                               -c1.evaluate (p5), c2);
        analyse();
        // precondition: all points distinct <=> conic nontrivial
        CGAL_optimisation_precondition (!is_trivial());
        if (o != _o) set_opposite();
    }

@end

@macro<ConicCPA2 public member functions> += @begin
    void set (const PT& p1, const PT& p2, const PT& p3, const PT& p4,
              const PT& p5, CGAL::Orientation _o = POSITIVE)
    {
        ConicCPA2<PT,DA> c1; c1.set_linepair (p1, p2, p3, p4);
        ConicCPA2<PT,DA> c2; c2.set_linepair (p1, p4, p2, p3);
        set_linear_combination (c2.evaluate (p5), c1, 
                               -c1.evaluate (p5), c2);
        analyse();
        // precondition: all points distinct <=> conic nontrivial
        CGAL_optimisation_precondition (!is_trivial());
        if (o != _o) set_opposite();
    }

@end


@! ---------------------------------------------------------------------------
@subsection{IO routines}
@! ---------------------------------------------------------------------------

@macro <ConicHPA2 I/O routines> = @begin
    template< class _PT, class _DA>
    std::ostream& operator << ( std::ostream& os, const ConicHPA2<_PT,_DA>& c)
    {
        return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
                   << c.u() << ' ' << c.v() << ' ' << c.w());
    }

    template< class _PT, class _DA>
    std::istream& operator >> ( std::istream& is, ConicHPA2<_PT,_DA>& c)
    {
        typedef           ConicHPA2<_PT,_DA>  Conic;
        typedef  typename _DA::RT                  RT;

        RT  r, s, t, u, v, w;
        is >> r >> s >> t >> u >> v >> w;
        c.set( r, s, t, u, v, w);

        return( is);
    }
@end

@macro <ConicCPA2 I/O routines> = @begin
    template< class _PT, class _DA>
    std::ostream& operator << ( std::ostream& os, const ConicCPA2<_PT,_DA>& c)
    {
        return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
                   << c.u() << ' ' << c.v() << ' ' << c.w());
    }

    template< class _PT, class _DA>
    std::istream& operator >> ( std::istream& is, ConicCPA2<_PT,_DA>& c)
    {
        typedef           ConicCPA2<_PT,_DA>  Conic;
        typedef  typename _DA::FT                  FT;

        FT  r, s, t, u, v, w;
        is >> r >> s >> t >> u >> v >> w;
        c.set( r, s, t, u, v, w);

        return( is);
    }
@end

@!----------------------------------------------------------------------------
@section{Files}
@!----------------------------------------------------------------------------

@i share/namespace.awi

@!----------------------------------------------------------------------------
@subsection{{\tt Conic\_misc.h}}
@!----------------------------------------------------------------------------

Here we collect all types and functions that are independent of the
representation type @prg{R}. These are the declarations of the enumeration
types @prg{Conic_type} and @prg{Convex_side}, the conic output
function and the functions in connection with the solution of cubic equations.

@file <include/CGAL/Conic_misc.h> = @begin
    @<file header>("include/CGAL/Conic_misc.h","2D Conic")

    #ifndef CGAL_CONIC_MISC_H
    #define CGAL_CONIC_MISC_H

    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")
    
    @<Conic_2 declaration>
    @<Conic_type declaration>
    @<Convex_side declaration>

    @!<Conic_2 Window_stream output>

    @<function best_value>
    @<function solve_cubic>

    @<namespace end>("CGAL")
    
    #endif // CGAL_CONIC_MISC_H

    @<end of file line>
@end

@!----------------------------------------------------------------------------
@subsection{{\tt Conic\_2.h}}
@!----------------------------------------------------------------------------

This file contains the implementation of the class @prg{Conic_2<R>}.
Depending on the loaded representation classes, the representation specific
classes @prg{ConicHPA2<PT,DA>} and/or  @prg{ConicCPA2<PT,DA>} are
included before that. 

@file <include/CGAL/Conic_2.h> = @begin
    @<file header>("include/CGAL/Conic_2.h","2D Conic")

    #ifndef CGAL_CONIC_2_H
    #define CGAL_CONIC_2_H
 
    #ifndef CGAL_REP_CLASS_DEFINED
    #  error  no representation class defined
    #endif // CGAL_REP_CLASS_DEFINED
 
    #ifdef CGAL_HOMOGENEOUS_H
    #  include <CGAL/ConicHPA2.h>
    #endif
 
    #ifdef CGAL_CARTESIAN_H
    #  include <CGAL/ConicCPA2.h>
    #endif

    #ifndef CGAL_POINT_2_H
    #  include <CGAL/Point_2.h>
    #endif

    #ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
    #include <CGAL/IO/forward_decl_window_stream.h>
    #endif

    @<namespace begin>("CGAL")
    
    @<Optimisation_ellipse_2 declaration>

    @<namespace end>("CGAL")
    
    @<Optimisation_ellipse_2 I/O operator declaration>

    @<namespace begin>("CGAL")

    @<Conic_2 interface and implementation>

    #ifndef CGAL_NO_OSTREAM_INSERT_CONIC_2
    @<Conic_2 I/O routines>
    #endif // CGAL_NO_OSTREAM_INSERT_CONIC_2

    @<namespace end>("CGAL")

    #endif // CGAL_CONIC_2_H

    @<end of file line>
@end

@!----------------------------------------------------------------------------
@subsection{{\tt ConicHPA2.h}}
@!----------------------------------------------------------------------------

Here is the class @prg{ConicHPA2<PT,DA>}\ldots

@file <include/CGAL/ConicHPA2.h> = @begin
    @<file header>("include/CGAL/ConicHPA2.h","2D Conic")

    #ifndef CGAL_CONICHPA2_H
    #define CGAL_CONICHPA2_H

    // includes
    #ifndef CGAL_CONIC_MISC_H
    #  include <CGAL/Conic_misc.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")

    @<ConicHPA2 declaration>
    @<ConicHPA2 interface and implementation>

    #ifndef CGAL_NO_OSTREAM_INSERT_CONICHPA2
    @<ConicHPA2 I/O routines>
    #endif // CGAL_NO_OSTREAM_INSERT_CONICHPA2

    @<namespace end>("CGAL")

    #endif // CGAL_CONICHPA2_H

    @<end of file line>
@end

@!----------------------------------------------------------------------------
@subsection{{\tt ConicCPA2.h}}
@!----------------------------------------------------------------------------

\ldots and the class @prg{ConicCPA2<PT,DA>}.

@file <include/CGAL/ConicCPA2.h> = @begin
    @<file header>("include/CGAL/ConicCPA2.h","2D Conic")

    #ifndef CGAL_CONICCPA2_H
    #define CGAL_CONICCPA2_H

    // includes
    #ifndef CGAL_CONIC_MISC_H
    #  include <CGAL/Conic_misc.h>
    #endif
    #ifndef CGAL_OPTIMISATION_ASSERTIONS_H
    #  include <CGAL/Optimisation/assertions.h>
    #endif

    @<namespace begin>("CGAL")
    
    @<ConicCPA2 declaration>
    @<ConicCPA2 interface and implementation>

    #ifndef CGAL_NO_OSTREAM_INSERT_CONICCPA2
    @<ConicCPA2 I/O routines>
    #endif // CGAL_NO_OSTREAM_INSERT_CONICCPA2

    @<namespace end>("CGAL")

    #endif // CGAL_CONICCPA2_H

    @<end of file line>
@end

@! ----------------------------------------------------------------------------
@! Conic_2_Window_stream.h
@! ----------------------------------------------------------------------------

\subsection{Conic\_2\_Window\_stream.h}
@file <include/CGAL/IO/Conic_2_Window_stream.h> = @begin
    @<file header>(
        "include/CGAL/IO/Conic_2_Window_stream.h",
        "graphical output to `leda_window' for Conic_2 algo.")

    // Each of the following operators is individually 
    // protected against multiple inclusion.

    // Window_stream I/O operators
    // ===========================

    // Conic_2
    // -------
    @<Conic_2 graphical output operator>

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
        "Min_ellipse_2",
        "Geometric Optimisation",
        "Conic_2",
        "$Revision$","$Date$",
        "Bernd Gärtner, Sven Schönherr <sven@@inf.ethz.ch>",
        "ETH Zürich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
        "@2")
@end

@! ===== EOF ==================================================================
