@! ============================================================================
@! The CGAL Project
@! Implementation: dD Smallest Enclosing Sphere
@! ----------------------------------------------------------------------------
@! file  : web/Optimisation/Min_sphere_d.aw
@! author: Bernd Gärtner (gaertner@inf.fu-berlin.de)
@! ----------------------------------------------------------------------------
@! $Revision$
@! $Date$
@! ============================================================================
 
@p maximum_input_line_length = 180
@p maximum_output_line_length = 180

@documentclass[twoside]{article}
@usepackage[latin1]{inputenc}
@usepackage{a4wide2}
@usepackage{amssymb}
@usepackage{cc_manual}
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
\newcommand{\WHILE}{\keyword{WHILE} }
\newcommand{\REPEAT}{\keyword{REPEAT} \+ \\ }
\newcommand{\UNTIL}{\< \keyword{UNTIL} \- }

\newcommand{\ct}[1]{\hfill (* #1 *)}
\newcommand{\q}{\mathbf{q}}
\newcommand{\C}{\mathbf{c}}

\newtheorem{Definition}{Definition}[section]
\newtheorem{Lemma}[Definition]{Lemma}
\newtheorem{Proposition}[Definition]{Proposition}
\newtheorem{Algorithm}[Definition]{Algorithm}

\newcommand{\ms}{\textrm{ms}}
\newcommand{\msb}{\textrm{msb}}
\newcommand{\welzlms}{\texttt{welzl\_ms}}
\newcommand{\mtfms}{\texttt{mtf\_ms}}
\newcommand{\pivotms}{\texttt{pivot\_ms}}
\newcommand{\msbasis}{\texttt{ms\_basis}}

\newcommand{\linebreakByHand}{\ccTexHtml{\\}{}}
\newcommand{\SaveSpaceByHand}{}  %%%%% [2]{\ccTexHtml{#1}{#2}}

@! ============================================================================
@! Title
@! ============================================================================
\RCSdef{\rcsrevision}{$Revision$}
\RCSdefDate{\rcsdate}{$Date$}

@t vskip 5 mm
@t title titlefont centre "CGAL -- Smallest Enclosing Sphere*"
@t vskip 1 mm
@t title smalltitlefont centre "Implementation Documentation"
@t vskip 5 mm
@t title smalltitlefont centre "Bernd Gärtner and Sven Schönherr"
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

We provide an implementation of an optimisation algorithm for computing
the smallest (w.r.t.\ volume) enclosing sphere of a finite point set $P$
in $d$-dimensional Euclidean space. The class template 
\ccc{Min_sphere_d} is implemented
as a semi-dynamic data structure, thus allowing to insert points while
maintaining the smallest enclosing sphere. It is parameterized with a
traits class that defines the interface between the
optimisation algorithm and the point class representing members of $P$. We
provide ready-to-use traits classes for computing smallest spheres 
of two-, three- and $d$-dimensional points of the \cgal\ kernel.

This document is organized as follows. Section \ref{sec:Min_sphere_d_specs}
contains the specifications and traits class requirements as they appear 
in the CGAL Reference Manual. Section \ref{sec:Min_sphere_d_impl} describes 
the algorithms we use and gives their implementation. 
Sections \ref{sec:opt_sphere_d}, \ref{sec:opt_sphere_dC} and 
\ref{sec:opt_sphere_dH} implements the sphere types that encapsulate 
the basic primitive operations of the \ccc{Min_sphere_d} algorithms. 
Sections \ref{sec:Min_sphere_d_traits_2}, \ref{sec:Min_sphere_d_traits_3} and
\ref{sec:Min_sphere_d_traits_d} provide the traits implementations for the
points of the \cgal\ kernel. Section \ref{sec:misc} contains miscellanea.
 
In Section \ref{sec:test} we provide a test program which performs 
some correctness checks. Finally, the product files are created in 
Section \ref{sec:files}.

\tableofcontents

@! ============================================================================
@! Specifications
@! ============================================================================

\clearpage
\section{Specifications and Requirements} 
\label{sec:Min_sphere_d_specs}

\emph{Note:} Below some references are undefined, they refer to sections
in the \cgal\ Reference Manual.

\renewcommand{\ccSection}{\ccSubsection}
\renewcommand{\ccFont}{\tt}
\renewcommand{\ccEndFont}{}
\newcommand{\cgalSetMinSphereLayout}{\ccTexHtml{%
\ccSetThreeColumns{Support_point_iterator}{}{creates a variable
    \ccc{min_sphere} of type \ccc{Min_sphere_d<Traits>}.}
\ccPropagateThreeToTwoColumns}{}}
\newcommand{\cgalColumnLayout}{\ccTexHtml{%
  \ccSetThreeColumns{Oriented_side}{}{\hspace*{8.5cm}}
  \ccPropagateThreeToTwoColumns}{}}

\input{Min_sphere_d.tex}

@! ============================================================================
@! Implementations
@! ============================================================================

\section{Class Template \texttt{CGAL\_Min\_sphere\_d<Traits>}} 
\label{sec:Min_sphere_d_impl}

\subsection{Introduction and Algorithms}

Given an $n$-point set $P=\{p_0,\ldots,p_{n-1}\}$ in $d$-space, we are looking
for the sphere of smallest radius that contains $P$. In other words, 
we are looking for the minimum value of $r$ such a point $c\in\R^d$ exists 
with $\|c-p_i\|\leq r$, for $i=0,\ldots,n-1$.

For sets $Q,B\subseteq P$, $Q\cap B=\emptyset$, let $\ms(Q,B)$ denote 
the smallest sphere that
contains $Q$ and has $B$ on its boundary (if it exists). The solution to
the problem is then given by $\ms(P) := \ms(P,\emptyset)$. 
Finally, let $\msb(B)$ denote $\ms(\emptyset,B)$, the smallest sphere
with the points of $B$ on its boundary. 

In this section, we describe algorithms to compute
$\ms(Q,B)$ for given sets $Q,B$.

\paragraph{Welzl's Method.}
The following randomized recursive algorithm, due to Welzl 
\cite{w-sedbe-91a}, computes
$\ms(Q,B)$ (if it exists). The expected runtime is $O(n)$ if $d$ is constant. 
The algorithm assumes the existence of a primitive operation
to compute $\msb(B)$, with $B$ a set of at most $d+1$ points. Section
\ref{sec:primitives} describes how this primitive can be realized. 

\begin{Algorithm} ~

\begin{pseudocode}{$\welzlms(Q,B)$:}
\IF $Q=\emptyset$ \OR $|B|=d+1$ \THEN 
        \RETURN $\msb(B)$ \\
\ELSE
        choose $p\in Q$ uniformly at random \\
        $\ms :=\welzlms(Q-\{p\},B)$ \\
        \IF $p\not\in \ms$ \THEN 
                \RETURN $\welzlms(Q-\{p\},B\cup\{p\})$ \\
        \ELSE
                \RETURN $\ms$ \\
        \END
\END
\end{pseudocode}
\end{Algorithm}

Note that if $\ms(Q,B)$ exists, then also $\msb(B)$ and $\ms(Q-\{p\},B)$
exist. Moreover, in \cite{w-sedbe-91a} it is shown that in case 
of $p\not\in \ms$, $\ms(Q,B)=\ms(Q-\{p\},B\cup\{p\})$ holds, and so 
the algorithm is always called with pairs $(Q,B)$ for which $\ms(Q,B)$
exists. The correctness also follows, once the termination criterion 
`$|B|=d+1$' is proven to be correct. For this, we also refer 
to \cite{w-sedbe-91a}. 

The call to $\msb(B)$ is termed the @em{basis computation}, 
$B$ itself is the @em{basis}. The query `$p\not\in \ms$' is the
@em{violation test}. It is clear that these
are the only operations in the algorithm where actual 
computations take place; moreover, given center and squared 
radius of the sphere, the violation test is straightforward.

The algorithm computes even more than just $\ms(Q,B)$. It
is also able to return a @em{support set} $S$ which is an inclusion-minimal
subset $S\subseteq Q$ such that $\ms(S,B)=\ms(Q,B)$ 
holds. This implies $\ms(S,B)=\msb(S\cup B)$. It is not hard to see that 
the last basis $B'$ for which the call to $\msb(B')$ was executed during 
the algorithm, defines a support set $S=B'-B$.

\paragraph{Welzl's Move-to-Front Method.}
To increase its efficiency in practice, Welzl's algorithm 
is usually enhanced with the
@em{move-to-front} heuristic. This heuristic keeps the points of $Q$ as an
ordered sequence $Q_k = (q_0,\ldots,q_{k-1})$, $k=|Q|$ and builds 
$\ms(Q,B)$ incrementally,
by adding one point after another from the sequence. In case $p=q_i$ is 
found to be outside the current ball, it is moved to the front of the 
sequence after the recursive call to compute $\ms(Q_{i},B\cup\{p\})$ with
$Q_i := (q_0,\ldots,q_{i-1})$ has been completed. The algorithm then 
looks as follows ($Q_k$ is initially set to a random 
permutation of the points in $Q$).  

\begin{Algorithm} ~

\begin{pseudocode}{$\mtfms(Q_k,B)$:}
$\ms := \msb(B)$ \\
\IF $|B|=d+1$ \THEN 
        \RETURN $\ms$ \\
\END
\FOR $i=0$ \TO $k-1$ \DO
        \IF $q_i\not\in \ms$ \THEN
                $\ms := \mtfms(Q_{i},B\cup\{q_i\})$ \\
                update $Q_k$ by moving $q_i$ to the front \\
        \END
\END
\RETURN $\ms$
\end{pseudocode}
\end{Algorithm}

Here is an alternative
view of the algorithm in which the pair $(\ms,B)$ is
maintained globally. This is the version we implement 
below.

\begin{Algorithm} ~

\label{alg:mtf_ms}
\begin{pseudocode}{$\mtfms(Q_k)$:}
\ct{Precondition: $\ms = \msb(B)$} \\
\ct{Postcondition: $\ms = \ms(Q_k)$} \\
\IF $|B|=d+1$ \THEN
        \RETURN \\
\END 
\FOR $i=0$ \TO $k-1$ \DO
        \IF $q_i\not\in \ms$ \THEN
                $B := B\cup\{q_i\}$ \\
                $\ms := \ms(B)$ \\
                $\mtfms(Q_{i})$ \\
                $B := B\setminus\{q_i\}$ \\
                update $Q_k$ by moving $q_i$ to the front \\
        \END
\END
\end{pseudocode}
\end{Algorithm}

By induction, one can prove that upon termination of $\mtfms(Q_k)$,
a support set $S$ appears as a prefix of the list $Q_k$. 

\paragraph{A Pivoting Variant.}
\label{sec:pivot_ms}
For larger dimension ($5\leq d\leq 20$), the following
variant is a dramatic improvement in practice. It calls \mtfms\ as a 
subroutine. As before, we assume that the pair $(\ms, B)$ is maintained
globally. In addition, the algorithm maintains the size $s$ of the current
support set after each call to \mtfms, where we assume that the procdure 
\mtfms\ itself sets this size. By the remark just made, the actual support 
set is then the length-$s$-prefix of the current list. 

Because we will use the variant only for $Q_n=P$, $P$ the whole
point set, we omit the parameter $Q_k$. 

\begin{Algorithm} ~

\label{alg:pivot_ms}
\begin{pseudocode}{$\pivotms()$:}
\ct{Precondition: $\ms = \msb(B)$} \\
\ct{Postcondition: $\ms = \ms(P)=\msb(B\cup S), S=\{q_0,\ldots,q_{s-1}\}$} \\
$t := \min (d+1,n)$ \\
$\mtfms (Q_{t})$  \\
\REPEAT
        \ct{$\ms = \msb(B\cup S) = (c, r^2), S = \{q_0,\ldots, q_{s-1}\}$} \\
        $excess := 0$ \\
        \FOR $i=t$ \TO $n-1$ \DO
                $e := \|q_i-c\|^2 -r^2$ \\
                \IF $e > excess$ \THEN
                        $excess = e$ \\
                        $pivot = i$ \\
                \END
        \END
        \IF $excess > 0$ \THEN
                $t := s+1$ \\
                $B:=B\cup\{q_i\}$ \\
                $\ms := \msb(B)$ \\
                $\mtfms (Q_s)$ \\
                $B := B\setminus\{q_i\}$ \\
                update $Q_k$ by moving $q_{pivot}$ to the front \\
        \END
\UNTIL $excess = 0$
\end{pseudocode}
\end{Algorithm} 
Note that starting the \FOR-loop with $i=t$ ensures that only points
are checked which have not already been handled by the previous call
to \mtfms (these points are contained in the current sphere anyway). 

Because $\mtfms$ is only called for small sets, the randomization is no
longer very useful and will be omitted. 

\subsection{Implementation}

\subsubsection{Types and Data Members}

@macro <Min_sphere_d types>  many = @begin


    public:     
        typedef typename Traits::Rep_tag        Rep_tag;
        typedef typename Traits::RT             RT;
        typedef typename Traits::FT             FT;
        typedef typename Traits::Point_d        Point; // Point type
          
        typedef  typename Traits::Access_dimension_d
          Access_dimension_d;
        typedef  typename Traits::Access_coordinates_begin_d
          Access_coordinates_begin_d;
        typedef  typename Traits::Construct_point_d
          Construct_point_d;
        
        typedef typename std::list<Point>::const_iterator       
                Point_iterator;
        typedef typename std::list<Point>::const_iterator
                Support_point_iterator;

    private:
        typedef typename std::list<Point>::iterator             It;

@end

The point set $P=Q_n$ is stored as an STL list, the pair $(\ms,B)$
is maintained as a sphere object of type 
@prg{Optimisation_sphere_d<R,Point,DA>}. 
To manipulate the pair $(\ms,B)$ as in the
algorithms above, the sphere offers functions \texttt{push} and 
\texttt{pop}, motivated by the fact that $B$ behaves like a stack
throughout the algorithm.

The support set $S$ is implicitly stored by keeping an iterator 
delimiting its defining prefix in $Q_n$.

The parameters of the sphere ($c$ and $r^2$) are only needed to compute 
the `excess' of $q_i$ w.r.t. $\ms$ -- the class 
@prg{Optimisation_sphere_d<R,Point,DA>} can do that more efficiently 
by itself and offers a corresponding member function. 

Finally, we store traits and data accessor objects and the ambient 
dimension, i.e. the dimension of the points in $P$.  

@macro <Min_sphere_d data members>  many = @begin
    private:
        int                                     d;            // ambient dim 
        std::list<Point>                        points;       // keeps P = Q_n
        Traits                                  tco;          // traits object
        Optimisation_sphere_d<Rep_tag, FT, RT, Point,Traits>        
                                                ms_basis; // keeps  miniball
        It                                      support_end;  // delimites S

#ifdef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION
        #define ms_basis(X) ms_basis(typename Traits::Rep_tag(), X)
#endif

@end

\subsubsection{Creation}
The default constructor creates $\ms(\emptyset)$, of ambient dimension
$-1$.

@macro<Min_sphere_d constructors>  += @begin
    Min_sphere_d ()
    : d(-1), tco( Traits()), ms_basis (tco),
      support_end(points.begin())
        {}

    
    Min_sphere_d (const Traits& traits)
    : d(-1), tco( traits), ms_basis (tco),
      support_end(points.begin())
        {}

    
        
@end

Here is the constructor to create $\ms(P)$ from a range of points. After
copying the points into the internal representation, it deduces the
ambient dimension from the first point, checks whether all points have
the same dimension, sets the size of \msbasis\ to the appropriate dimension
and calls the private method
\pivotms. 

@macro<Min_sphere_d constructors> += @begin
    // STL-like constructor (member template) 
    template <class InputIterator>
    Min_sphere_d( InputIterator first,
                       InputIterator last)
        : d(-1), points( first, last), tco( Traits()), ms_basis (tco),
          support_end(points.begin())
    {
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim(d));
            ms_basis.get_sphere(Rep_tag()).set_size (d); 
            pivot_ms();
        } 
    }

    template <class InputIterator>
    Min_sphere_d( InputIterator first,
                       InputIterator last,
                       const Traits& traits)
        : d(-1), points( first, last), tco( traits), ms_basis (tco),
          support_end(points.begin())
    {
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim(d));
            ms_basis.get_sphere(Rep_tag()).set_size (d); 
            pivot_ms();
        } 
    }
    
@end

The copy constructor just copies the points and reconstructs the sphere.
This is a linear-time operation, because (as argued above) the support
points come first in the point list, leading to a fast reconstruction
algorithm. 

@macro <Min_sphere_d constructors> += @begin
    Min_sphere_d (const Min_sphere_d& msph) 
    : d(msph.ambient_dimension()),
      points (msph.points_begin(), msph.points_end()), tco (msph.tco),
          ms_basis (tco), support_end (points.begin())
          
    {
        if (d != -1) {
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        }
    }

@end

The assignment operator removes the old point set before reconstructing. 

@macro <Min_sphere_d constructors> += @begin
    Min_sphere_d& operator=(const Min_sphere_d& msph)
    {
        if (this != &msph) {
            points.erase (points.begin(), points.end());
            d = msph.ambient_dimension();
            points.insert
              (points.begin(), msph.points_begin(), msph.points_end());
              ms_basis.get_sphere(Rep_tag()).set_tco(msph.tco);
            support_end = points.begin();
            tco = msph.tco;
            if (d != -1) {
                ms_basis.get_sphere(Rep_tag()).set_size (d);
                pivot_ms();
            }
        } 
        return *this;   
    }

@end


\subsubsection{Access Functions}
These are either straightforward or rely on corresponding
access functions of the type @prg{Optimisation_sphere_d<R>}.

@macro <Min_sphere_d access functions>  many = @begin
    int number_of_points() const
    {
        return points.size();
    }

    int number_of_support_points() const
    {
        return ms_basis.get_sphere(Rep_tag()).number_of_support_points();
    }

    Point_iterator points_begin () const
    {
        return points.begin();
    }
 
    Point_iterator points_end () const
    {
        return points.end();
    }

    Support_point_iterator support_points_begin () const
    {
        return points.begin();
    }
 
    Support_point_iterator support_points_end () const
    {
        return support_end; 
    }

    int ambient_dimension () const
    {
        return d;
    }

    Point center () const
    {
        CGAL_optimisation_precondition (!is_empty());
        return ms_basis.get_sphere(Rep_tag()).center();
    }

    FT squared_radius () const
    {
        CGAL_optimisation_precondition (!is_empty());
        return ms_basis.get_sphere(Rep_tag()).squared_radius();
    }

@end

\subsubsection{Predicates}
These are as easy as the access functions. All sidedness tests rely on
the @prg{excess} predicate of the sphere, which measures $\|p-c\|^2-r^2$
for given point $p$. For the sidedness-tests, we check the required
preconditions.

@macro <Min_sphere_d predicates>  many = @begin
    Bounded_side bounded_side (const Point& p) const
    {
        if (d == -1)
           return ON_UNBOUNDED_SIDE;
        else {
          CGAL_optimisation_precondition
           (d == tco.access_dimension_d_object()(p));
           return (Bounded_side
               (-CGAL::sign (ms_basis.get_sphere(Rep_tag()).excess (p))));
        }
    }

    bool has_on_bounded_side (const Point& p) const
    {
        if (d == -1)
           return false;
        else {
          CGAL_optimisation_precondition
           (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_negative (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool has_on_unbounded_side (const Point& p) const
    {
        if (d == -1)
           return true;
        else {
          CGAL_optimisation_precondition
          (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_positive (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool has_on_boundary (const Point& p) const
    {
        if (d == -1)
           return false;
        else {
          CGAL_optimisation_precondition
          (d == tco.access_dimension_d_object()(p));
           return (CGAL_NTS is_zero (ms_basis.get_sphere(Rep_tag()).excess (p)));
        }
    }

    bool is_empty () const
    {
        return (d == -1);
    }

    bool is_degenerate () const
    {
        return (ms_basis.get_sphere(Rep_tag()).number_of_support_points() < 2);
    }

@end

\subsubsection{Modifiers}

@macro <Min_sphere_d modifiers>  += @begin
   void clear ()
   {
        d = -1;
        points.erase (points.begin(), points.end());
        ms_basis.get_sphere(Rep_tag()).set_size (-1);
        support_end = points.begin();
   }

@end 

The set function is similar to the constructor above, after erasing 
the old point set. 

@macro <Min_sphere_d modifiers> += @begin
    // STL-like set(member template) 
    template <class InputIterator>
    void set ( InputIterator first,
               InputIterator last)
    {
        points.erase (points.begin(), points.end());
        points.insert (points.begin(), first, last);
        support_end = points.begin();
        if (points.size()>0) {
            d = tco.access_dimension_d_object() (points.front());
            CGAL_optimisation_precondition ((d>=0) && all_points_have_dim (d));
            ms_basis.get_sphere(Rep_tag()).set_size (d);
            pivot_ms();
        } else {
            d = -1;
            ms_basis.get_sphere(Rep_tag()).set_size (-1);
        }
    }

@end

To insert a single point $p$ into $\ms(P)$, we distinguish two cases. 
If $p$ is contained in the sphere \ms\ currently maintained, we just need 
to append $p$ to the list $P$. Otherwise, we know that $\ms(P\cup\{p\})$
equals $\ms(P,\{p\})$, which we compute using the method \pivotms. 
Postcondition is as in @prg{mtf_ms} and @prg{pivot_ms} that a support 
set appears as a prefix of the point list.

@macro <Min_sphere_d modifiers> += @begin
    void insert (const Point& p)
    {
        if (has_on_unbounded_side (p)) {
            if (is_empty()) {
                d = tco.access_dimension_d_object() (p);
                CGAL_optimisation_precondition (d>=0);
                ms_basis.get_sphere(Rep_tag()).set_size (d);
            }
            // ensure precondition of pivot_ms
            ms_basis.get_sphere(Rep_tag()).push (p);
            pivot_ms ();
            ms_basis.get_sphere(Rep_tag()).pop ();
            points.push_front (p);  // ensure postcondition of insert
        } else
            points.push_back (p);   // just append p
            if (support_end == points.end()) --support_end;
    }

@end

To insert a range, we call @prg{insert} for all points in the range.

@macro <Min_sphere_d modifiers> += @begin
    template <class InputIterator>
    void insert (InputIterator first, InputIterator last)
    {
        for (InputIterator i=first; i!=last; ++i)
            insert (*i);
    }

@end

\subsubsection{Validity check}
The validity check performs two tests. First, it checks whether all points
$p\in P$ are contained in the current sphere, and second, it calls the 
validity check of the sphere to verify that the current sphere is indeed 
equal to $\msb(S)$, $S$ the current support set, and that $S$ is minimal 
with this property.

To perform the first check, the natural thing to do would be to get
the current center and squared radius and then compute the excess of 
every point $p\in P$ by hand. However, this cannot be done in a
representation-independent way - the formulas look different for
cartesain and homogeneous coordinates. This could be fixed by
requiring the data accessor to be able to return a cartesian (or a
homogeneous) representation in any case. To put this burden on the
implementor just for the purpose of checking, however, seems not 
appropriate. Instead we assume the correctness of the
@prg{has_on_unbounded_side} predicate, which is realized via the @prg{excess}
method of the
representation-dependent type @prg{Sphere}. Inspect the corresponding
code if you have doubts about this assumption.
 
@macro <Min_sphere_d validity check>  many = @begin
    bool is_valid (bool verbose = false, int level = 0) const
    {
        Verbose_ostream verr (verbose);

        // sphere verification
        verr << "  (a) sphere verification..." << std::flush;
        if (ms_basis.get_sphere(Rep_tag()).is_valid (verbose))
            verr << "passed." << std::endl;
        else 
            return false;

        // containment check
        verr << "  (b) containment check..." << std::flush;

        // non-support-points
        Point_iterator i;
        for (i=support_end; i!=points.end(); ++i)
            if (has_on_unbounded_side (*i))
                return (_optimisation_is_valid_fail (verr, 
                  "sphere does not contain all points"));

        // support points
        for (i=points.begin(); i!=support_end; ++i)
            if (!has_on_boundary (*i))
                return (_optimisation_is_valid_fail (verr, 
                  "sphere does not have all support points on boundary"));
  
        verr << "passed." << std::endl;
        verr << "object is valid!" << std::endl;
        return true;
    }

@end

\subsubsection{Miscellaneous}

@macro <Min_sphere_d miscellaneous>  many = @begin
   const Traits& traits() const
   {
        return tco;
   }

@end

\subsubsection{I/O}
          
@macro <Min_sphere_d I/O operators declaration>  many = @begin
    template < class Traits >
    std::ostream&
    operator << ( std::ostream& os, const Min_sphere_d<Traits>& ms);

    template < class Traits >
    std::istream&
    operator >> ( std::istream& is, Min_sphere_d<Traits> & ms);
@end

@macro <Min_sphere_d I/O operators>  many = @begin
    template < class Traits >
    std::ostream&
    operator << ( std::ostream& os, const Min_sphere_d<Traits>& min_sphere)
    {
        typedef typename Min_sphere_d<Traits>::Point  Point;

        switch ( get_mode( os)) {

          case IO::PRETTY:
            os << std::endl;
            os << "Min_sphere_d( |P| = " << min_sphere.number_of_points()
               << ", |S| = " << min_sphere.number_of_support_points()
               << std::endl;
            os << "  P = {" << std::endl;
            os << "    ";
            std::copy( min_sphere.points_begin(), min_sphere.points_end(),
                  std::ostream_iterator<Point>( os, ",\n    "));
            os << "}" << std::endl;
            os << "  S = {" << std::endl;
            os << "    ";
            std::copy( min_sphere.support_points_begin(),
                  min_sphere.support_points_end(),
                  std::ostream_iterator<Point>( os, ",\n    "));
            os << "}" << std::endl;
            os << "  center = " << min_sphere.center() << std::endl;
            os << "  squared radius = " << min_sphere.squared_radius()
               << std::endl;
            os << ")" << std::endl;
            break;

          case IO::ASCII: 
            os << min_sphere.number_of_points() << std::endl;
            std::copy( min_sphere.points_begin(), min_sphere.points_end(),
                  std::ostream_iterator<Point>( os, "\n"));
            break;

          case IO::BINARY:
            os << min_sphere.number_of_points() << " ";
            std::copy( min_sphere.points_begin(), min_sphere.points_end(),
                  std::ostream_iterator<Point>( os, " "));
            break;

          default:
            CGAL_optimisation_assertion_msg
                ( false, "get_mode( os) invalid!");
            break; }

        return( os);
    }

    template < class Traits >
    std::istream&
    operator >> ( std::istream& is, Min_sphere_d<Traits>& min_sphere)
    {
        switch ( get_mode( is)) {

          case IO::PRETTY:
            std::cerr << std::endl;
            std::cerr << "Stream must be in ascii or binary mode" << std::endl;
            break;

          case IO::ASCII:
          case IO::BINARY:
          {
            min_sphere.clear();
            int n; is >> n;
            typename Min_sphere_d<Traits>::Point p;
            for (int i=0; i<n; ++i) {
                is >> p;
                min_sphere.insert (p);
            }
          } break;

          default:
            CGAL_optimisation_assertion_msg( false, "IO::mode invalid!");
            break; 
     }

        return( is);
    }

@end

\subsubsection{Private member function \texttt{mtf\_ms}}
This method realizes the pseudocode of Algorithm \ref{alg:mtf_ms}. As a
postcondition, @prg{support_end} is a past-the-end iterator of the 
support points in the list.

@macro <Min_sphere_d private method mtf_ms>  many = @begin
    void mtf_ms (It k)
    {
        support_end = points.begin();
        if (ms_basis.get_sphere(Rep_tag()).size_of_basis()==d+1) return;
        for (It i = points.begin(); i!=k;) {
            It j = i++; 
            if (CGAL_NTS is_positive (ms_basis.get_sphere(Rep_tag()).excess(*j))) {
                ms_basis.get_sphere(Rep_tag()).push (*j);
                mtf_ms (j);
                ms_basis.get_sphere(Rep_tag()).pop();
                move_to_front (j);
            } 
        }
    }

@end
         
\subsubsection{Private member function \texttt{pivot\_ms}}
This method realizes the pseudocode of Algorithm \ref{alg:pivot_ms}. As a
postcondition, @prg{support_end} is a past-the-end iterator of the 
support points
in the list. You might wonder why the @prg{excess} method returns a ring
type element; here the hack is that @prg{excess} returns only a positive
multiple of the real excess, where the multiplicative factor is the same
for all points. This functionality is completely sufficient for our 
purposes and more efficient to realize.

@macro <Min_sphere_d private method pivot_ms>  many = @begin
    void pivot_ms ()
    {
        It t = points.begin();
        std::advance (t, std::min (d+1, (int)points.size()));
        mtf_ms (t);

        RT excess, e;
        do {
            excess = RT(0);
            It pivot;
            for (It i=t; i!=points.end(); ++i) {
                e = ms_basis.get_sphere(Rep_tag()).excess(*i); 
                if (e > excess) {
                   excess = e;
                   pivot = i;
                }
            }
            if (CGAL_NTS is_positive (excess)) {
                t = support_end;
                if (t==pivot) ++t; //  inserted from the esa code       
                ms_basis.get_sphere(Rep_tag()).push (*pivot);
                mtf_ms (support_end); 
                ms_basis.get_sphere(Rep_tag()).pop();
                move_to_front (pivot);
            }
        } while (CGAL_NTS is_positive (excess));
    }

@end

\subsubsection{Private member function \texttt{move\_to\_front}}
Here we call the @prg{splice} method of the point list which does
exactly the move-to-front operation. However, if the current
@prg{support_end} iterator happens to point to the item to be
moved, this iterator has to be advanced first. 

@macro <Min_sphere_d private method move_to_front> many = @begin
    void move_to_front (It j)
    {
        if (support_end == j) 
           ++support_end;
        points.splice (points.begin(), points, j); 
    }

@end
        
\subsubsection{Private member function \texttt{all\_points\_have\_dim}}
This function checks whether all points currently stored have the
given dimension $d$. This is more efficient than checking the dimension
of a point each time it is handled during the algorithm (which can happen
several times).

@macro <Min_sphere_d private method all_points_have_dim> many = @begin
    bool all_points_have_dim (int dim) const
    {
        for (Point_iterator i=points.begin(); i!=points.end(); ++i)
            if (tco.access_dimension_d_object()(*i) != dim)
                return false;
        return true;
    }

@end

\section{Class Template
\texttt{CGAL\_Optimisation\_sphere\_d<Rep\_tag,FT,RT,PT,DA>}}
\label{sec:opt_sphere_d}

\subsection{Introduction and the Main Primitive}
\label{sec:primitives}
An object of the type @prg{Optimisation_sphere_d<Rep_tag,FT,RT,PT,DA>} is a 
pair $(\ms, B)$, where $\ms = \msb(B\cup S)$ for some $S$ disjoint from $B$.
The set $B$ behaves like a stack throughout the algorithm. 

The sphere has an ambient dimension associated with it, this is the
dimension of the points that form the sets $B$. This dimension is 
$-1$ after construction, meaning that it is not known.  

The class has representation-dependent
specializations (over the @prg{Rep_tag}) which realize the actual
functionality. 

The main functionality provided is the computation of $\msb(B)$, the
smallest sphere with a point set $B$ on the boundary. The algorithm
guarantees that $B$ is always a set of affinely independent points, 
from which $|B|\leq d+1$ follows. In that case, 
$\msb(B)$ is the unique circumsphere of the points in $B$ with center
restricted to the affine hull of $B$. This means, the center
$\C$ and squared radius $r^2$ can be computed from the following system 
of equations, where $B=\{\q_0,\ldots,\q_{m-1}\}, m\leq d+1$. 
\begin{eqnarray*}
\|\q_i-\C\|^2 &=& r^2, \quad i=0,\ldots m-1, \\
\sum_{i=0}^{m-1} \lambda_i \q_i &=& \C, \\
\sum_{i=0}^{m-1} \lambda_i &=& 1.
\end{eqnarray*}
Assume $\q_i$ and $\C$ have cartesian coordinates $q_i, c\in \R^d$. 
Defining $\alpha := r^2-c^Tc$ and substituting $\C$ in the first equation
with $\sum_{i=0}^{m-1} \lambda_i q_i$, we obtain a linear system in 
the variables $\alpha, \lambda_0,\ldots, \lambda_{m-1}$ which we can 
write in matrix form as follows.

\begin{equation}
\label{smb_system}
\left(\begin{array}{cccc}
0 & 1 & \cdots & 1 \\
1 & 2q_0^Tq_0 & \cdots &2q_0^Tq_{m-1} \\
\vdots & \vdots & & \vdots \\
1 & 2q_0^Tq_{m-1} & \cdots &2q_{m-1}^Tq_{m-1}
\end{array}\right) 
\left(\begin{array}{c} \alpha \\ \lambda_0 \\ \vdots \\ \lambda_{m-1}
\end{array}\right) =
\left(\begin{array}{c} 1 \\ q_0^Tq_0 \\ \vdots \\ 
q_{m-1}^Tq_{m-1} \end{array}\right).
\end{equation}

Computing the values of 
$\alpha, \lambda_0,\ldots, \lambda_{m-1}$ 
amounts to solving the linear system (\ref{smb_system}). $\C$ and 
$r^2$ are then easily obtained via
\begin{eqnarray}
\label{eq:center}
\C &=&  \sum_{i=0}^{m-1} \lambda_i \q_i, \\
\label{eq:sqr_radius}
r^2 &=& \alpha + \|\C\|^2.
\end{eqnarray}

In Section \ref{sec:opt_sphere_dC} we discuss 
how to solve system (\ref{smb_system}) in the cartesian
case; Section \ref{sec:opt_sphere_dH} presents a similar system for points
with homogeneous coordinates and a 
corresponding solution method. In both case, it is crucial for the
efficiency that the stack-like behavior of $B$ is exploited. If
a point $\q_m$ is added to $B$, we do not solve the resulting
system (\ref{smb_system}) from scratch but update our solution 
device from the previous $B$.  

\section{Specialized Class Template
\texttt{CGAL\_Optimisation\_sphere\_d<Cartesian\_tag,FT,RT,PT,DA>}}
\label{sec:opt_sphere_dC}

\subsection{Realizing the Main Primitive}

The device for solving system (\ref{smb_system}) will be the explicit
inverse $A^{-1}_B$ of the matrix 
\[ A_B := 
\left(\begin{array}{cccc}
0 & 1 & \cdots & 1 \\
1 & 2q_0^Tq_0 & \cdots &2q_0^Tq_{m-1} \\
\vdots & \vdots & & \vdots \\
1 & 2q_0^Tq_{m-1} & \cdots &2q_{m-1}^Tq_{m-1}
\end{array}\right) 
\]
along with the vector
\[ v_B := 
\left(\begin{array}{c} 1 \\ q_0^Tq_0 \\ \vdots \\ 
q_{m-1}^Tq_{m-1} \end{array}\right).
\]

Given $A^{-1}_B$ and $v_B$, it takes merely a matrix-vector 
multiplication to obtain the solution to (\ref{smb_system}).

\subsection{Updating the Inverse}
\label{cartesian_update}
Assume, point $\q_m$ is pushed onto $B$, resulting in the set $B'$.
Computing $A^{-1}_{B'}$ from scratch would take time $O(dm + m^3)$ using
Gauss elimination ($O(dm)$ time is needed to compute the new values
$q_j^Tq_m, j=0,\ldots,m$). We are going to get $A^{-1}_{B'}$ faster and 
without any Gauss elimination, by updating the previous inverse
$A^{-1}_{B}$.

It is not hard to check that
\[A_{B'} =  \left(\begin{array}{ccc|c}
1 &        &   & 0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & 0 \\ \hline
a_0 & \cdots & a_m & 1
\end{array}\right) 
\left(\begin{array}{ccc|c}
  &        &   & 0 \\
  &   A_B  &   & \vdots \\
  &        &   & 0 \\ \hline
0 & \cdots & 0 & z
\end{array}\right)
\left(\begin{array}{ccc|c}
1 &        &   & a_0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & a_m \\ \hline
0 & \cdots & 0 & 1
\end{array}\right),
\]
where 
\[
\left(\begin{array}{c} a_0 \\ \vdots \\ a_m \end{array}\right) =
A^{-1}_B \left(\begin{array}{c} 1 \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right)\]
and
\begin{eqnarray*}
z &=& 2q_m^Tq_m - \left(1, 2q_0^Tq_m, \ldots, 2q_{m-1}^Tq_m\right ) 
A^{-1}_B \left(\begin{array}{c} 1 \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right) \\
&=& 2q_m^Tq_m - a_0 - 2a_1q_0^Tq_m\cdots - 2a_mq_{m-1}^Tq_m .
\end{eqnarray*}

This implies
\[A^{-1}_{B'} = 
\left(\begin{array}{ccc|c}
1 &        &   & -a_0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & -a_j \\ \hline
0 & \cdots & 0 & 1
\end{array}\right)
\left(\begin{array}{ccc|c}
  &        &   & 0 \\
  &   A^{-1}_B  &   & \vdots \\
  &        &   & 0 \\ \hline
0 & \cdots & 0 & 1/z
\end{array}\right)
\left(\begin{array}{ccc|c}
1 &        &   & 0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & 0 \\ \hline
-a_0 & \cdots & -a_j & 1
\end{array}\right),
\]
equivalently
\begin{eqnarray}
\nonumber
A^{-1}_{B'} &=&
\frac{1}{z}
\left(\begin{array}{c|c}
z A^{-1}_B + 
\left(\begin{array}{cccc}
a_0^2 & a_0a_1 & \cdots & a_0a_m \\
a_0a_1 & a_1^2 & \cdots & a_1a_m \\
\vdots & \vdots    &        & \vdots \\
a_0a_m & a_1a_m & \cdots & a_m^2
\end{array}\right) &
\begin{array}{c}
-a_0 \\ -a_1 \\ \vdots \\ -a_m
\end{array} \\ \hline
\begin{array}{cccc}
-a_0 & -a_1 & \cdots & -a_m
\end{array} & 1
\end{array}
\right) \\
\label{eq:cart_update}
&=& 
\left(\begin{array}{c|c}
A^{-1}_B + 
\left(\begin{array}{cccc}
a_0^2 & a_0a_1 & \cdots & a_0a_m \\
a_0a_1 & a_1^2 & \cdots & a_1a_m \\
\vdots & \vdots    &        & \vdots \\
a_0a_m & a_1a_m & \cdots & a_m^2
\end{array}\right) / z &
\begin{array}{c}
-a_0/z \\ -a_1/z \\ \vdots \\ -a_m/z
\end{array} \\ \hline
\begin{array}{cccc}
-a_0/z & -a_1/z & \cdots & -a_m/z
\end{array} & 1/z
\end{array}
\right).
\end{eqnarray}

To perform this update, we need to compute $q_j^Tq_m, j=0,\ldots,m$ as
before in time $O(dm)$. The values $a_0,\ldots,a_m$ are obtained in time 
$O(m^2)$ by a matrix multiplication, $z$ can then be derived in time
$O(m)$. The final update according to
(\ref{eq:cart_update}) takes another $O(m^2)$ time, for a total of 
$O(dm + m^2)$ time. This means,
we save an $O(m)$-factor compared to the complete reinversion of
$A_{B'}$; moreover, the update is easy to code. 

\subsection{Implementation}

@macro <Optimisation_sphereCd declaration> many = @begin
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphere_d<Cartesian_tag, FT, RT, PT, Traits>;

@end

\subsubsection{Data Members}
In the algorithm, $B=\{\q_0,\ldots,\q_{m-1}\}$ is an ordered set (elements are 
ordered by time of addition to $B$). Define $B^j := \{\q_0,\ldots,\q_{j}\}$.
For any $j$ in the range $0\leq j\leq m-1$ we store
\begin{itemize}
\item $q_j$, the cartesian representation of $\q_j$,  
\item the matrix $A^{-1}_{B^j}$ (because this matrix is symmetric, 
the entries below and including the main diagonal suffice) 
\end{itemize}
In addition, we keep the vector $v_B=(1,q_0^Tq_0,\ldots,q_{m-1}^Tq_{m-1})^T$
(and therefore have all vectors $v_{B^j}$ implicitly available).

As indicated above, from this it takes merely a matrix-vector multiplication 
to obtain the solution to (\ref{smb_system}). For this, we provide
a (private) multiplication method and reserve storage for the 
solution vector $x=(\alpha_0,\lambda_0,\ldots,\lambda_{m-1})^T$.

For doing the matrix update, we reserve space for 
the vector $v=(1, 2q_0^Tq_m,\ldots,2q_{m-1}^Tq_m)^T$. To store the $a_i$,
we also use $x$. 

Let's check the required dimensions of the data: $m$ can be at most
$d+1$, so $j$ possibly assumes $d+1$ values between $0$ and $d$. 
$A^{-1}_{B^j}$ has size $(j+2)\times (j+2)$, so the vectors $v_B$, $x$ and $v$
have dimension at most $d+2$. Because $d$ is not known at compile time,
we store only pointers which we allocate in the constructor below that
receives a dimension argument. 

The set $B\cup S$ will be of the form $\{q_0,\ldots,q_{s-1}\},s\geq m$,
and we keep its size. In addition, we store center and squared radius 
of the current ball $\msb(B\cup S)$.

Finally, we have a data accessor object.

@macro <Optimisation_sphereCd types and data members>  many = @begin
    typedef
      typename Traits::Access_coordinates_begin_d::Coordinate_iterator  It;

    // hack needed for egcs, see also test programs
    typedef             FT                              FT_;

    int                 d;                      // dimension
    int                 m;                      // |B|
    int                 s;                      // |B\cup S|
    
    FT**                q;                      // the q_j's
    FT***               inv;                    // the A^{-1}_{B^j}'s
    FT*                 v_basis;                // the vector v_B
     
    FT*                 x;                      // solution vector
    FT*                 v;                      // auxiliary vector

    FT*                 c;                      // center, for internal use
    PT                  ctr;                    // center, for external use
    FT                  sqr_r;                  // squared_radius

    Traits              tco;

    

@end

\subsubsection{Construction}

The default constructor generates an empty sphere of ambient dimension 
$-1$. Only combinatorial information (number of support points, size of
basis, validity) can be retrieved from this sphere. Excess, center and
squared radius are undefined.   

@macro <Optimisation_sphereCd constructor>  many = @begin
    Optimisation_sphere_d (const Traits& t = Traits())
        : d(-1), m(0), s(0), tco (t)
    {}

    void set_tco (const Traits& tcobj)
    {
        tco = tcobj;
    }

@end

The method @prg{init} allocates space for the arrays declared above,
according to the dimensions given above. Moreover, it sets $c=0$ and
$r^2=-1$ to create an empty sphere and puts the leading $1$ into the
vector $v_B$. init will only be called for @prg{ambient_dimension} at least
0.
 
@macro <Optimisation_sphereCd init method>  many = @begin
    void init (int ambient_dimension)
    {
        d = ambient_dimension;
        m = 0;
        s = 0;
        sqr_r = -FT(1);

        q =             new FT*[d+1];
        inv =           new FT**[d+1];
        v_basis =       new FT[d+2];
        x =             new FT[d+2];
        v =             new FT[d+2];
        c =             new FT[d];

        for (int j=0; j<d+1; ++j) {
            q[j] =      new FT[d];
            inv[j] =    new FT*[j+2];
            for (int row=0; row<j+2; ++row)
                inv[j][row] = new FT[row+1];
        }
        
        for (int i=0; i<d; ++i) 
            c[i] = FT(0);
        v_basis[0] = FT(1);
    }

@end    

\subsubsection{Destruction}
Here we just get rid of the arrays, by calling the method @prg{destroy}.

@macro <Optimisation_sphereCd destructor> many = @begin
    ~Optimisation_sphere_d ()
    {
        if (d != -1)
           destroy();
    }

@end

@macro <Optimisation_sphereCd destroy method>  many = @begin
    void destroy ()
    {
        for (int j=0; j<d+1; ++j) {
            for (int row=0; row<j+2; ++row)
                delete[] inv[j][row];
            delete[] inv[j];
            delete[] q[j];
        }
        delete[] c;
        delete[] v;
        delete[] x;
        delete[] v_basis;
        delete[] inv;
        delete[] q;
    }

@end

\subsubsection{Method \texttt{set\_size}}

Resets the sphere to an empty sphere of the given dimension (which
can also be $-1$). 

@macro <Optimisation_sphereCd set_size method> many = @begin
    void set_size (int ambient_dimension)
    {
        if (d != -1)
            destroy();
        if (ambient_dimension != -1)
            init(ambient_dimension);
        else {
            d = -1;
            m = 0;
            s = 0;
        }
    }

@end
            

\subsubsection{Method \texttt{push}}
This method implements the update from Section \ref{cartesian_update} and
recomputes center and squared radius. Note that if currently $m=0$ holds, 
we have no matrix to update. In this case, we generate the matrix
\[ A^{-1}_{B^0} = \left( \begin{array}{cc}
                        0 & 1 \\
                        1 & 2q_0^Tq_0
                  \end{array}\right)^{-1} = 
                  \left( \begin{array}{cc}
                        - 2q_0^Tq_0& 1 \\
                        1 & 0
                  \end{array}\right)
\]
`by hand'.

@macro <Optimisation_sphereCd push method>  many = @begin
    void push (const PT& p)
    {   
        // store q_m = p by copying its cartesian coordinates into q[m]
        It i(tco.access_coordinates_begin_d_object()(p)); FT *o;
        for (o=q[m]; o<q[m]+d; *(o++)=*(i++));

        // update v_basis by appending q_m^Tq_m
        v_basis[m+1] = prod(q[m],q[m],d);

        if (m==0) 
        {
            // set up A^{-1}_{B^0} directly
            FT** M = inv[0];
            M[0][0] = -FT_(2)*v_basis[1];
            M[1][0] = FT_(1);
            M[1][1] = FT_(0);
        } else {
            // set up vector v by computing 2q_j^T q_m, j=0,...,m-1
            v[0] = FT_(1);
            for (int j=0; j<m; ++j)
                v[j+1] = FT_(2)*prod(q[j],q[m],d);
                 
            // compute a_0,...,a_m
            multiply (m-1, v, x);               // x[j]=a_j, j=0,...,m

            // compute z
            FT z = FT_(2)*v_basis[m+1] - prod(v,x,m+1);
            CGAL_optimisation_assertion (!CGAL_NTS is_zero (z));
            FT inv_z = FT_(1)/z;

            // set up A^{-1}_{B^m}
            FT** M = inv[m-1];          // A^{-1}_B, old matrix
            FT** M_new = inv[m];        // A^{-1}_{B'}, new matrix

            // first m rows
            int row, col;
            for (row=0; row<m+1; ++row)
                for (col=0; col<row+1; ++col)
                    M_new [row][col] = M[row][col] + x[row]*x[col]*inv_z;

            // last row
            for (col=0; col<m+1; ++col)
                M_new [m+1][col] = -x[col]*inv_z; 
            M_new [m+1][m+1] = inv_z;
        } 
        s = ++m; 
        compute_c_and_sqr_r();  // side effect: sets x
    }   

@end

\subsubsection{Method \texttt{pop}}

We decrease $m$; note that the sphere $\ms(B\cup S)$ 
(represented by $s$, the center and squared radius is unaffected 
by the @prg{pop} operation by definition, see Section \ref{sec:pivot_ms}.

@macro <Optimisation_sphereCd pop method>  many = @begin
    void pop ()
    {
        --m;
    }

@end

\subsubsection{Method \texttt{excess}}

@macro <Optimisation_sphereCd excess method>  many = @begin
    FT excess (const PT& p) const
    {
        // compute (c-p)^2
        FT sqr_dist (FT(0));
        It i(tco.access_coordinates_begin_d_object()(p));
        FT *j;
        for (j=c; j<c+d; ++i, ++j)
            sqr_dist += CGAL_NTS square(*i-*j);
        return sqr_dist - sqr_r;
     }

@end

\subsubsection{Access Methods}

@macro <Optimisation_sphereCd access methods>  many = @begin

   PT center () const
   {
        return tco.construct_point_d_object()(d, c, c+d);
   }

   FT squared_radius () const
   {
        return sqr_r;
   }

   int number_of_support_points () const
   {
        return s;
   }

   int size_of_basis () const
   {
        return m;
   }

@end

\subsubsection{Method \texttt{is\_valid}}

The sphere is valid if the center is strictly in the convex hull of the 
support points and no support point is redundant. The latter holds if the
matrix $A_B$ was never singular. We implicitly checked this already
with the assertion @prg{!is_zero(z)} during the matrix update formula.
The convex hull property is easy to check when going back to the formula 
for the center
given by (\ref{eq:center}). The vector @prg{x} set up in computing the
center still stores the multipliers $\lambda_0,\ldots,\lambda_{m-1}$ of the
support points that define the center, so we just need to check whether 
they are positive. Note that $\lambda_0$ has index 1 in $x$. 

@macro <Optimisation_sphereCd is_valid method>  many = @begin
    bool is_valid (bool verbose = false, int level = true) const
    {
        Verbose_ostream verr (verbose);
        for (int j=1; j<m+1; ++j)
            if (!CGAL_NTS is_positive (x[j]))
                return (_optimisation_is_valid_fail
                    (verr, "center not in convex hull of support points"));
        return (true);
    }
        
@end        

\subsubsection{Private Method \texttt{multiply}}
This method computes the product \[A^{-1}_{B^j} v\]
and stores it in \texttt{x[0],...,x[j+1]}. Recall that we 
only keep the entries below the diagonal of $A^{-1}_{B^j}$.

@macro <Optimisation_sphereCd multiply method>  many = @begin
    void multiply (int j, const FT* vec, FT* res)
    {
        FT** M = inv[j];
        for (int row=0; row<j+2; ++row) {
            res[row] = prod(M[row],vec,row+1);
            for (int col = row+1; col<j+2; ++col)
                res[row] += M[col][row]*vec[col];
        }
    }

@end

\subsubsection{Private Method \texttt{compute\_c\_and\_sqr\_r}}
This is easy, using formulas (\ref{eq:center}) and (\ref{eq:sqr_radius}).
The multipliers $\lambda_j$ are stored in @prg{x[1],...,x[m]}, $\alpha$
appears in @prg{x[0]}.

@macro <Optimisation_sphereCd compute_c_and_sqr_r method>  many = @begin
    void compute_c_and_sqr_r ()
    {
        multiply (m-1, v_basis, x);

        for (int i=0; i<d; ++i) c[i] = FT(0);
        for (int j=0; j<m; ++j) {
            FT l = x[j+1], *q_j = q[j];
            for (int i=0; i<d; ++i) 
                c[i] += l*q_j[i];
        }
        sqr_r = x[0] + prod(c,c,d);
    }

@end
            

\subsubsection{Private Method \texttt{prod}}
Computes the dot product of the first $k$ entries of $v$ and $w$

@macro <Optimisation_sphereCd prod method>  many = @begin
    FT prod (const FT* v1, const FT* v2, int k) const
    {
        FT res(FT(0));
        for (const FT *i=v1, *j=v2; i<v1+k; res += (*(i++))*(*(j++)));
        return res;
    }

@end 

\section{Specialized Class Template
\texttt{CGAL\_Optimisation\_sphere\_d<Homogeneous\_tag,FT,RT,PT,DA>}}
\label{sec:opt_sphere_dH}

\subsection{Realizing the Main Primitive}

Let's reconsider the equation system that led to (\ref{smb_system}). 
If $\q_i$ is in homogeneous representation, then its cartesian
representation is of the form $q_i/h_i\in \R^d, q_i\in \R^d, h_i\in \R$.
We say that $q_i$ is the @em{cartesian part} of $\q_i$. In that case, system
(\ref{smb_system}) can be written as

\begin{equation}
\label{smb_system_hom}
\left(\begin{array}{cccc}
0 & h_0 & \cdots & h_{m-1} \\
h_0 & 2q_0^Tq_0 & \cdots &2q_0^Tq_{m-1} \\
\vdots & \vdots & & \vdots \\
h_{m-1} & 2q_0^Tq_{m-1} & \cdots &2q_{m-1}^Tq_{m-1}
\end{array}\right) 
\left(\begin{array}{c} \alpha \\ \lambda_0/h_0 \\ \vdots \\
\lambda_{m-1}/h_{m-1}
\end{array}\right) =
\left(\begin{array}{c} 1 \\ q_0^Tq_0/h_0 \\ \vdots \\ 
q_{m-1}^Tq_{m-1}/h_{m-1} \end{array}\right).
\end{equation}

Again, our device for solving this will be the explicit inverse $A^{-1}_B$
of the system's matrix, represented in rational form
\[A^{-1}_B = \frac{\tilde{A}^{-1}_B}{D},\]
where the entries of $\tilde{A}^{-1}_B$ are elements of the ring type
underlying the homogeneous representation of the points, and $D$ is the
determinant $D=\det(A_B)$. (Such a representation exists as easily follows
from Cramer's rule). Let 
\[ h = \prod_{i=0}^{m-1}h_i\]
and
\[\left(\begin{array}{c} \tilde{\alpha} \\
\tilde{\lambda}_0 \\ \vdots \\ \tilde{\lambda}_{m-1}\end{array}\right)
= \tilde{A}^{-1}_B\cdot h\cdot \left(\begin{array}{c} 1 \\ q_0^Tq_0/h_0 
\\ \vdots \\ q_{m-1}^Tq_{m-1}/h_{m-1} \end{array}\right).\]
Then the vector 
$(\tilde{\alpha},\tilde{\lambda}_0, \ldots, \tilde{\lambda}_{m-1})^T$
is a vector of ring elements, and the center $\C$ of the sphere can 
be obtained as a homogeneous point with cartesian part

\begin{equation}
\label{eq:centerh}
c =  \sum_{i=0}^{m-1} \tilde{\lambda}_i q_i
\end{equation} and homogenizing coordinate $hD$.

The squared radius is the field type element
\begin{equation}
\label{eq:sqr_radiush}
r^2 = \alpha + \|\C\|^2 = \frac{\tilde{\alpha}}{hD} + \frac{c^Tc}{(hD)^2}
= \frac{\tilde{\alpha}hD + c^Tc}{(hD)^2}.
\end{equation}
        
\subsection{Updating the Inverse}
\label{homogeneous_update}
As before in the cartesian case, we are concerned with computing the 
inverse $A^{-1}_{B'}$ from the `old' inverse $A^{-1}_{B}$, 
$B' = B\cup\{\q_m\}.$

Recall that
\[
A_{B'} = \left(\begin{array}{cccc|c}
&  &        &   & h_m \\
&  &        &   & 2q_0^Tq_{m} \\
&  &   A_B  &   & \vdots \\
&  &        &   & 2q_{m-1}^Tq_{m} \\ \hline
h_m & 2q_0^Tq_{m} & \cdots & 2q_{m-1}^Tq_{m} & 2q_{m}^Tq_{m}
\end{array}\right).\]

This implies
\begin{equation}
\label{eq:update_formula} 
A_{B'} =  \left(\begin{array}{ccc|c}
1 &        &   & 0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & 0 \\ \hline
a_0 & \cdots & a_m & 1
\end{array}\right) 
\left(\begin{array}{ccc|c}
  &        &   & 0 \\
  &   A_B  &   & \vdots \\
  &        &   & 0 \\ \hline
0 & \cdots & 0 & z
\end{array}\right)
\left(\begin{array}{ccc|c}
1 &        &   & a_0 \\
  &  \ddots     &   & \vdots \\
  &        & 1  & a_m \\ \hline
0 & \cdots & 0 & 1
\end{array}\right),
\end{equation}
where 
\[
\left(\begin{array}{c} a_0 \\ \vdots \\ a_m \end{array}\right) =
A^{-1}_B \left(\begin{array}{c} h_m \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right)\]
and
\begin{eqnarray*}
z &=& 2q_m^Tq_m - \left(h_m, 2q_0^Tq_m, \ldots, 2q_{m-1}^Tq_m\right ) 
A^{-1}_B \left(\begin{array}{c} h_m \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right) \\
&=& 2q_m^Tq_m - a_0h_m - 2a_1q_0^Tq_m\cdots - 2a_mq_{m-1}^Tq_m .
\end{eqnarray*}

As in the cartesian case, we get
\[A^{-1}_{B'} =
\frac{1}{z}
\left(\begin{array}{c|c}
z A^{-1}_B + 
\left(\begin{array}{cccc}
a_0^2 & a_0a_1 & \cdots & a_0a_m \\
a_0a_1 & a_1^2 & \cdots & a_1a_m \\
\vdots & \vdots    &        & \vdots \\
a_0a_m & a_1a_m & \cdots & a_m^2
\end{array}\right) &
\begin{array}{c}
-a_0 \\ -a_1 \\ \vdots \\ -a_m
\end{array} \\ \hline
\begin{array}{cccc}
-a_0 & -a_1 & \cdots & -a_m
\end{array} & 1
\end{array}
\right).\]

Now recall that we store $A^{-1}_B$ in the form
\[A^{-1}_B = \frac{\tilde{A}^{-1}_B}{D}.\]

Define 
\[
\left(\begin{array}{c} \tilde{a}_0 \\ \vdots \\ \tilde{a}_m \end{array}\right) 
= \tilde{A}^{-1}_B \left(\begin{array}{c} h_m \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right)\]

and 
\begin{eqnarray*}
\tilde{z} &=& D\cdot 2q_m^Tq_m - 
\left(h_m, 2q_0^Tq_m, \ldots, 2q_{m-1}^Tq_m\right ) 
\tilde{A}^{-1}_B \left(\begin{array}{c} h_m \\ 2q_0^Tq_{m} \\ 
\vdots \\2q_{m-1}^Tq_{m} \end{array}\right) \\
&=& D\cdot 2q_m^Tq_m - \tilde{a}_0h_m - 2\tilde{a}_1q_0^Tq_m\cdots -
2\tilde{a}_mq_{m-1}^Tq_m.
\end{eqnarray*}

Then the $\tilde{a}_i$ and $\tilde{z}$ are ring type values with
$\tilde{a}_i = Da_i$ for all $i$ and $\tilde{z}=Dz$, and we can write
\begin{eqnarray*}
A^{-1}_{B'} &=& 
\frac{D}{\tilde{z}}
\left(\begin{array}{c|c}
\frac{\tilde{z} \tilde{A}^{-1}_B}{D^2} + 
\left(\begin{array}{cccc}
\tilde{a}_0^2 & \tilde{a}_0\tilde{a}_1 & \cdots & \tilde{a}_0\tilde{a}_m \\
\tilde{a}_0\tilde{a}_1 & \tilde{a}_1^2 & \cdots & \tilde{a}_1\tilde{a}_m \\
\vdots & \vdots    &        & \vdots \\
\tilde{a}_0\tilde{a}_m & \tilde{a}_1\tilde{a}_m & \cdots & \tilde{a}_m^2
\end{array}\right) / D^2 &
\begin{array}{c}
-\tilde{a}_0/D \\ -\tilde{a}_1/D \\ \vdots \\ -\tilde{a}_m/D
\end{array} \\ \hline
\begin{array}{cccc}
-\tilde{a}_0/D & -\tilde{a}_1/D & \cdots & -\tilde{a}_m/D
\end{array} & 1
\end{array}
\right) \\
&=&
\frac{1}{D\tilde{z}}
\left(\begin{array}{c|c}
\tilde{z} \tilde{A}^{-1}_B + 
\left(\begin{array}{cccc}
\tilde{a}_0^2 & \tilde{a}_0\tilde{a}_1 & \cdots & \tilde{a}_0\tilde{a}_m \\
\tilde{a}_0\tilde{a}_1 & \tilde{a}_1^2 & \cdots & \tilde{a}_1\tilde{a}_m \\
\vdots & \vdots    &        & \vdots \\
\tilde{a}_0\tilde{a}_m & \tilde{a}_1\tilde{a}_m & \cdots & \tilde{a}_m^2
\end{array}\right) &
\begin{array}{c}
-\tilde{a}_0 D \\ -\tilde{a}_1 D \\ \vdots \\ -\tilde{a}_m D
\end{array} \\ \hline
\begin{array}{cccc}
-\tilde{a}_0 D & -\tilde{a}_1 D & \cdots & -\tilde{a}_m D
\end{array} & D^2
\end{array}
\right).
\end{eqnarray*}
Because we have $\tilde{z}=\det(A_{B'})$, as easily follows from 
(\ref{eq:update_formula}), all ring elements in the matrix must be 
divisible by $D$, and performing this division (which the ring type
supports exactly), we can store $A^{-1}_{B'}$ as a matrix of ring 
elements, along with a common denominator $\tilde{z}$.

\subsection{Implementation}

@macro <Optimisation_sphereHd declaration> many = @begin 
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphere_d<Homogeneous_tag, FT, RT, PT, Traits>;

@end


\subsubsection{Data Members}

Considering the solution- and update formulas above, it becomes clear
that we need some additional data members, compared to the cartesian
case. As before, we store the cartesian part $q_j$ of the point $q_j$
and the matrix ${A}^{-1}_{B^j}$, for every $j\in\{0,\ldots,m-1\}$, in
the rational format $\tilde{A}^{-1}_{B^j}/D_j$, this means we need one
more array to hold the denominators. 

Storing the vector 
\[\tilde{v}_{B^{m-1}} = h\cdot \left(\begin{array}{c} 1 \\ q_0^Tq_0/h_0 \\ 
\vdots \\ q_{m-1}^Tq_{m-1}/h_{m-1} \end{array}\right)\]
that directly enables us to compute the values $\tilde{\alpha}$ and
$\tilde{\lambda}_0,\ldots,\tilde{\lambda}_{m-1}$ by matrix multiplication
is substantially more expensive now than in the cartesian case, where we 
needed only one array which was updated upon insertion of a new point 
$\q_m$. Here we need one vector 
\[\tilde{v}_{B^j} = \prod_{k=0}^j h_j\cdot 
\left(\begin{array}{c} 1 \\ q_0^Tq_0/h_0 \\ 
\vdots \\ q_{j}^Tq_{j}/h_{j} \end{array}\right)\]
for every value of $j$ to allow 
an efficient update. The center is stored as a homogeneous point now, with
homogenizing coordinate $hD$, as suggested by the formula (\ref{eq:centerh})
above. The squared radius is kept by its ring type numerator (its denominator 
is $(hD)^2$).  

@macro <Optimisation_sphereHd types and data members>  many = @begin
    typedef
      typename Traits::Access_coordinates_begin_d::Coordinate_iterator  It;

    // hack needed for egcs, see also test programs
    typedef             RT                              RT_;

    int                 d;                      // dimension
    int                 m;                      // |B|
    int                 s;                      // |B\cup S|
    
    RT**                q;                      // the q_j's
    RT***               inv;                    // the \tilde{A}^{-1}_{B^j}'s
    RT*                 denom;                  // corresponding denominators
    RT**                v_basis;                // the \tilde{v}_B^j
     
    RT*                 x;                      // solution vector
    mutable RT*         v;                      // auxiliary vector

    RT*                 c;                      // center, for internal use
    PT                  ctr;                    // center, for external use
    RT                  sqr_r;                  // numerator of squared_radius

    Traits              tco; 

@end

\subsubsection{Construction}

The constructor allocates space for the arrays declared above,
according to the dimensions as given in the cartesian case, with the exception
of @prg{v_basis} which is an array of vectors now and @prg{c} which gets
$d+1$ coordinates. Plus, there is an additional array of denominators.
As before, the constructor sets $c=0$ and $r^2=-1$ to create an empty 
sphere. 
 
@macro <Optimisation_sphereHd constructor>  many = @begin
    Optimisation_sphere_d (const Traits& t = Traits())
        : d(-1), m(0), s(0), tco(t)
    {}

    void set_tco (const Traits& tcobj)
    {
        tco = tcobj;
    }

@end

@macro <Optimisation_sphereHd init method>  many = @begin
    void init (int ambient_dimension)
    {
        d = ambient_dimension;
        m = 0;
        s = 0;
        sqr_r = -RT(1);

        q =             new RT*[d+1];
        inv =           new RT**[d+1];
        denom =         new RT[d+1];
        v_basis =       new RT*[d+1];
        x =             new RT[d+2];
        v =             new RT[d+2];
        c =             new RT[d+1];

        for (int j=0; j<d+1; ++j) {
            q[j] =      new RT[d];
            inv[j] =    new RT*[j+2];
            v_basis[j] =new RT[j+2];
            for (int row=0; row<j+2; ++row)
                inv[j][row] = new RT[row+1];
        }
        
        for (int i=0; i<d; ++i) 
            c[i] = RT(0);
        c[d] = RT(1);
    }

@end    

\subsubsection{Destruction}
Here we just get rid of the arrays, by calling @prg{destroy}

@macro <Optimisation_sphereHd destructor>  many = @begin
    ~Optimisation_sphere_d ()
    {
        if (d != -1) 
            destroy();
    }

@end

@macro <Optimisation_sphereHd destroy method>  many = @begin
    void destroy ()
    {
        for (int j=0; j<d+1; ++j) {
            for (int row=0; row<j+2; ++row)
                delete[] inv[j][row];
            delete[] v_basis[j];
            delete[] inv[j];
            delete[] q[j];
        }
        delete[] c;
        delete[] v;
        delete[] x;
        delete[] v_basis;
        delete[] denom;
        delete[] inv;
        delete[] q;
    }

@end

\subsubsection{Method \texttt{set\_size}}

Resets the sphere to an empty sphere of the given dimension (which
can also be $-1$). 

@macro <Optimisation_sphereHd set_size method> many = @begin
    void set_size (int ambient_dimension)
    {
        if (d != -1)
            destroy();
        if (ambient_dimension != -1)
            init(ambient_dimension);
        else {
            d = -1;
            m = 0;
            s = 0;
        } 
    }

@end

\subsubsection{Method \texttt{push}}
This method implements the update from Section \ref{homogeneous_update} and
recomputes center and squared radius. 
We do not only update the matrix here but also set the vector
$\tilde{v}_{B^m}$. Note that if currently $m=0$ holds, we have no matrix 
and vector to update. In this case, we generate the matrix
\[ A^{-1}_{B^0} = \left( \begin{array}{cc}
                        0 & h_0 \\
                        h_0 & 2q_0^Tq_0
                  \end{array}\right)^{-1} = 
                  \frac{1}{-h_0^2} 
                  \left( \begin{array}{cc}
                        2q_0^Tq_0& -h_0 \\
                        -h_0 & 0
                  \end{array}\right)
\]
and the vector
\[v_{B^0} = \left(\begin{array}{c} h_0 \\ q_0^Tq_0\end{array}\right)\]
`by hand'.

@macro <Optimisation_sphereHd push method>  many = @begin
    void push (const PT& p)
    {
        // store q_m = p by copying its cartesian part into q[m]
        It i(tco.access_coordinates_begin_d_object()(p)); RT *o;
        for (o=q[m]; o<q[m]+d; *(o++)=*(i++));

        // get homogenizing coordinate
        RT hom = *(i++);

        if (m==0) 
        {
            // set up v_{B^0} directly
            v_basis[0][0] = hom;
            v_basis[0][1] = prod(q[0],q[0],d);
            
            // set up \tilde{A}^{-1}_{B^0} directly
            RT** M = inv[0];
            M[0][0] = RT_(2)*v_basis[0][1];
            M[1][0] = -hom;
            M[1][1] = RT_(0);
            denom[0] = -CGAL_NTS square(hom);  // det(\tilde{A}_{B^0})
            
        } else {
            // set up v_{B^m}
            int j;
            RT sqr_q_m = prod(q[m],q[m],d);
            v_basis[m][m+1] = v_basis[m-1][0]*sqr_q_m;
            for (j=0; j<m+1; ++j)
                v_basis[m][j] = hom*v_basis[m-1][j];

            
            // set up vector v by computing 2q_j^T q_m, j=0,...,m-1
            v[0] = hom;
            for (j=0; j<m; ++j)
                v[j+1] = RT_(2)*prod(q[j],q[m],d);
                 
            // compute \tilde{a}_0,...,\tilde{a}_m
            multiply (m-1, v, x);               // x[j]=\tilde{a}_j, j=0,...,m

            // compute \tilde{z}
            RT old_denom = denom[m-1];
            RT z = old_denom*RT_(2)*sqr_q_m - prod(v,x,m+1);
            CGAL_optimisation_assertion (!CGAL_NTS is_zero (z));

            // set up \tilde{A}^{-1}_{B^m}
            RT** M = inv[m-1];          // \tilde{A}^{-1}_B, old matrix
            RT** M_new = inv[m];        // \tilde{A}^{-1}_{B'}, new matrix

            // first m rows
            int row, col;
            for (row=0; row<m+1; ++row)
                for (col=0; col<row+1; ++col)
                    M_new [row][col] 
                        = (z*M[row][col] + x[row]*x[col])/old_denom;

            // last row
            for (col=0; col<m+1; ++col)
                M_new [m+1][col] = -x[col]; 
            M_new [m+1][m+1] = old_denom;

            // new denominator
            denom[m] = z;
        } 
        s = ++m; 
        compute_c_and_sqr_r();
    }   

@end

\subsubsection{Method \texttt{pop}}

@macro <Optimisation_sphereHd pop method>  many = @begin
    void pop ()
    {
        --m;
    }

@end

\subsubsection{Method \texttt{excess}}

The value computed is the $(h_phD)^2$-fold multiple of the true excess, where
$h_p$ is the homogeneous coordinate of $p$. This way we avoid computing 
with field type numbers. Recall that $hD$ is stored in the homogenizing 
coordinate of the center @prg{c}.

@macro <Optimisation_sphereHd excess method>  many = @begin
    RT excess (const PT& p) const
    {
        // store hD times the cartesian part of p in v
        RT hD = c[d];
        It i(tco.access_coordinates_begin_d_object()(p)); RT *o;

#ifndef CGAL_CFG_NO_MUTABLE

        for ( o=v; o<v+d; *(o++)=hD*(*(i++)));
#else
        for (o=(RT*)v; o<(RT*)v+d; *(o++)=hD*(*(i++)));

#endif // CGAL_CFG_NO_MUTABLE

        // get h_p
        RT h_p = *(i++);
        CGAL_optimisation_precondition (!CGAL_NTS is_zero (h_p));

        // compute (h_p h D)^2 (c-p)^2
        RT sqr_dist(RT(0));
        for (int k=0; k<d; ++k)
            sqr_dist += CGAL_NTS square(h_p*c[k]-v[k]);

        // compute excess
        return sqr_dist - CGAL_NTS square(h_p)*sqr_r;
     }

@end

\subsubsection{Access Methods}

@macro <Optimisation_sphereHd access methods>  many = @begin

   PT center () const
   {
        return tco.construct_point_d_object()(d,c,c+d+1);
   }

   FT squared_radius () const
   {
        return FT(sqr_r)/FT(CGAL_NTS square(c[d]));
   }

   int number_of_support_points () const
   {
        return s;
   }

   int size_of_basis () const
   {
        return m;
   }

@end

\subsubsection{Method \texttt{is\_valid}}
As above, we still need to check whether the $\lambda_j$ that define
the center are strictly positive. However, we do not have them explicitly.
The vector @prg{x} stores the values $\tilde{\lambda}_j$, where
\[\tilde{\lambda}_j = \frac{hD\lambda_j}{h_j}, \]
equivalently
\[\lambda_j = \frac{\tilde{\lambda}_jh_j}{hD}.\]
We can immediately access the signs of $\tilde{\lambda}_j$ and $hD$ (which
is the homogenizing coordinate of the center @prg{c}). To get the sign of
$h_j$, we observe that $h_j$ is the quotient of the first entries of
$v_{B^j}$ and $v_{B^{j-1}}$, and these entries are explicitly stored.

@macro <Optimisation_sphereHd is_valid method>  many = @begin
    bool is_valid (bool verbose = false, int level = true) const
    {
        if (d==-1) return true;
        Verbose_ostream verr (verbose);
        int sign_hD = CGAL::sign(c[d]), s_old = 1, 
            s_new = CGAL::sign(v_basis[0][0]), signum;
        for (int j=1; j<m+1; ++j) {
            signum = sign_hD * s_old * s_new * CGAL::sign(x[j]);
            if (!CGAL_NTS is_positive (signum))
                return (_optimisation_is_valid_fail
                    (verr, "center not in convex hull of support points"));
            s_old = s_new; s_new = CGAL::sign(v_basis[j][0]);
        }
        return true;    
    }
        
@end        

\subsubsection{Private Method \texttt{multiply}}
This method computes the product \[\tilde{A}^{-1}_{B^j} v\]
and stores it in \texttt{x[0],...,x[j+1]}. Recall that we 
only keep the entries below the diagonal of $\tilde{A}^{-1}_{B^j}$.

@macro <Optimisation_sphereHd multiply method>  many = @begin
    void multiply (int j, const RT* vec, RT* res)
    {
        RT** M = inv[j];
        for (int row=0; row<j+2; ++row) {
            res[row] = prod(M[row],vec,row+1);
            for (int col = row+1; col<j+2; ++col)
                res[row] += M[col][row]*vec[col];
        }
    }

@end

\subsubsection{Private Method \texttt{compute\_c\_and\_sqr\_r}}
This is easy, using formulas (\ref{eq:centerh}) and (\ref{eq:sqr_radiush}).
The multipliers $\tilde{\lambda}_j$ are stored in @prg{x[1],...,x[m]}, 
$\tilde{\alpha}$ appears in @prg{x[0]}.

@macro <Optimisation_sphereHd compute_c_and_sqr_r method>  many = @begin
    void compute_c_and_sqr_r ()
    {
        // solve
        multiply (m-1, v_basis[m-1], x);

        // set cartesian part
        for (int i=0; i<d; ++i) c[i] = RT(0);
            for (int j=0; j<m; ++j) {
                RT l = x[j+1], *q_j = q[j];
                for (int i=0; i<d; ++i)
                    c[i] += l*q_j[i];
        }
        c[d] = v_basis[m-1][0]*denom[m-1];                // hD
        sqr_r = x[0]*c[d] + prod(c,c,d); // \tilde{\alpha}hD+c^Tc
    }

@end

\subsubsection{Private Method \texttt{prod}}
Computes the dot product of the first $k$ entries of $v$ and $w$

@macro <Optimisation_sphereHd prod method>  many = @begin
    RT prod (const RT* v1, const RT* v2, int k) const
    {
        RT res = RT(0);
        for (const RT *i=v1, *j=v2; i<v1+k; res += (*(i++))*(*(j++)));
        return res;
    }

@end 


            
@! ============================================================================
@! Test Programs
@! ============================================================================
\section{Test Programs}
\label{sec:test}

Here we run test programs, with the aim of code coverage and consistency
checking. However, even if the programs goes through without errors, the
implementation might still contain hidden bugs -- not everything can be
checked here.

We generate a set of random input points, in cartesian and homogeneous 
representation (with homogenizing coordinates not equal to one). Then
we compute their smallest enclosing spheres using several methods,
compare the results and apply the remaining functions. 

\subsection{Testing the $d$-dimensional default traits}

\subsubsection{Typedefs and Constants}

@macro <Min_sphere_d test program - definitions> += @begin
    typedef leda_rational                       NT;
    CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_rational);
    typedef NT                                  FT;
    typedef NT                                  RT;
    typedef Cartesian<FT>                       C;
    typedef Homogeneous<RT>                     H;
    typedef Optimisation_d_traits_d<C>          Cartesian_traits;
    typedef Optimisation_d_traits_d<H>          Homogeneous_traits;
    typedef Min_sphere_d<Cartesian_traits>      Min_sphereC;
    typedef Min_sphere_d<Homogeneous_traits>    Min_sphereH;
    typedef Point_d<C>                          PointC;
    typedef Point_d<H>                          PointH;


    // checked traits
    Cartesian_traits                                    tC;
    Homogeneous_traits                                  tH;

    const int n = 10;                           // number of points
    const int D = 2;                            // dimension of points
    const int r = 128;                          // coordinate range
    bool verbose = true;

    // NT == Quotient<NT> ?
    bool equals (NT x, Quotient<NT> y)
    {
        return (Quotient<NT>(x) == y);
    }

    // PointC == Point H ?
    bool equals (PointC p, PointH q) 
    {
        int dim =  p.dimension();
        if (q.dimension() != dim) return false;
        for (int j=0; j<dim; ++j)
           if (!equals(p.cartesian(j),q.cartesian(j)))
                return false;
        return true;
    }

    
@end

\subsubsection{Generation of Points}

@macro <Min_sphere_d test program - generate points> many = @begin
    Point_d<C>  PC[n];          // cartesian point set
    Point_d<H>  PH[n];          // homogeneous point set
    int         coord[D];       // arrays for random coordinates
    FT          coordC[D];      // ... cartesian version
    RT          coordH[D+1];    // ... homogeneous version
    Random      my_random;      // random number generator

   
    for (int i=0; i<n; ++i) {
        int j;
        // random coordinates
        for (j=0; j<D; ++j)
            coord[j] = my_random (r);

        // cartesian point
        for (j=0; j<D; ++j)
            coordC[j] = FT (coord[j]);
        PC[i] = Point_d<C>(D, coordC, coordC+D);

        // homogeneous point
        for (j=0; j<D; ++j)
            coordH[j] = RT(2*coord[j]);
        coordH[D] = RT(2);
        PH[i] = Point_d<H>(D, coordH, coordH+D+1);
    }

@end

\subsubsection{Constructors}

@macro <Min_sphere_d test program - test constructors> many = @begin
    // test constructors

    // default
    #ifdef __BORLANDC__
    // problems with default traits argument in constructor calls
    Min_sphereC         msC_empty(tC); assert (msC_empty.is_valid(verbose));
    Min_sphereH         msH_empty(tH); assert (msH_empty.is_valid(verbose));
    #else
    Min_sphereC         msC_empty; assert (msC_empty.is_valid(verbose));
    Min_sphereH         msH_empty; assert (msH_empty.is_valid(verbose));
    #endif
    
    // from range
    Min_sphereC         msC (PC, PC+n, tC); assert (msC.is_valid(verbose));
    Min_sphereH         msH (PH, PH+n, tH); assert (msH.is_valid(verbose));

    // copy
    Min_sphereC         msC1 (msC); assert (msC1.is_valid(verbose));
    Min_sphereH         msH1 (msH); assert (msH1.is_valid(verbose));
    
    PointC              centerC (msC.center()), centerC1 (msC1.center());
    PointH              centerH (msH.center()), centerH1 (msH1.center());

    FT                  radiusC (msC.squared_radius()),
                        radiusC1 (msC1.squared_radius());
    Quotient<RT>        radiusH (msH.squared_radius()),
                        radiusH1 (msH1.squared_radius());

    assert (equals (centerC, centerH));
    assert (centerC == centerC1); assert(centerH == centerH1);
    assert (equals (radiusC, radiusH));
    assert (radiusC == radiusC1); assert(radiusH == radiusH1);

    // assignment 
    msC1 = msC; msH1 = msH;
    assert (centerC == centerC1); assert(centerH == centerH1);
    assert (radiusC == radiusC1); assert(radiusH == radiusH1);

@end

\subsubsection{Modifiers}

@macro <Min_sphere_d test program - test modifiers> many = @begin
   // test set method
   msC.set (PC, PC+n); assert (msC.is_valid(verbose));
   msH.set (PH, PH+n); assert (msH.is_valid(verbose));
   assert (centerC == msC.center());
   assert (centerH == msH.center());
   assert (radiusC == msC.squared_radius());
   assert (radiusH == msH.squared_radius());

   // test clear and insert methods
   msC.clear(); assert (msC.is_valid(verbose));
   msH.clear(); assert (msH.is_valid(verbose));
   msC.insert (PC, PC+n); assert (msC.is_valid(verbose));
   msH.insert (PH, PH+n); assert (msH.is_valid(verbose));
   assert (centerC == msC.center());
   assert (centerH == msH.center());
   assert (radiusC == msC.squared_radius());
   assert (radiusH == msH.squared_radius());

   // combined set and insert
   msC.set (PC, PC+n/2); msC.insert (PC+n/2, PC+n); 
   assert (msC.is_valid(verbose));
   msH.set (PH, PH+n/2); msH.insert (PH+n/2, PH+n); 
   assert (msH.is_valid(verbose));
   assert (centerC == msC.center());
   assert (centerH == msH.center());
   assert (radiusC == msC.squared_radius());
   assert (radiusH == msH.squared_radius());

@end

\subsubsection{Access Functions}

@macro <Min_sphere_d test program - test access functions> many = @begin
   // test access functions

   // number_of_points
   assert (msC.number_of_points() == n);
   assert (msH.number_of_points() == n);
 
   // number_of_support_points
   assert (msC.number_of_support_points() == msH.number_of_support_points());

   // points_begin, points_end
   int m; 
   m = std::distance (msC.points_begin(), msC.points_end());
   assert (m == n);
   m = std::distance (msH.points_begin(), msH.points_end());
   assert (m == n);

   // support_points_begin, support_points_end
   m = 
   std::distance (msC.support_points_begin(), msC.support_points_end());
   assert (m == msC.number_of_support_points());
   m =
   std::distance (msH.support_points_begin(), msH.support_points_end());
   assert (m == msH.number_of_support_points());

   // ambient dim
   assert (msC.ambient_dimension() == D);
   assert (msH.ambient_dimension() == D);

   // center and squared radius already tested

@end

\subsubsection{Predicates}

@macro <Min_sphere_d test program - test predicates> many = @begin
   // test predicates

   // bounded_side
   assert (msC.bounded_side (centerC) == ON_BOUNDED_SIDE);
   assert (msH.bounded_side (centerH) == ON_BOUNDED_SIDE);

   // has_on_bounded_side
   assert (msC.has_on_bounded_side (centerC));
   assert (msH.has_on_bounded_side (centerH));

   // has_on_boundary already tested in is_valid method
   // has_on_unbounded_side already tested in is_valid method

   // is_empty
   assert (!msC.is_empty());
   assert (!msH.is_empty());

   // is_degenerate
   assert (!msC.is_degenerate());
   assert (!msH.is_degenerate());

@end
  
\subsubsection{Miscellaneous}

\subsubsection{I/O}
We write the spheres into a string and read them back in to compare
it with the written one.

@macro <Min_sphere_d test program - test I/O> many = @begin
   int size = 100*n;    
   char* buffer = new char[size];  // suffices to store the sphere

   std::ostrstream ost (buffer, size);  // output string
   set_ascii_mode (ost);
   ost << msC << msH << std::endl;      // write spheres

   std::istrstream ist (buffer, size);  // input string
   set_ascii_mode (ist);           
   ist >> msC >> msH;              // read spheres

   assert (centerC == msC.center());
   assert (centerH == msH.center());
   assert (radiusC == msC.squared_radius());
   assert (radiusH == msH.squared_radius());

   delete[] buffer;

@end

\subsubsection{The Complete Test Program}
 
Here is the complete test program.

@macro <Min_sphere_d test program> many = @begin
    // hack to overcome internal compiler error of egcs
    #define FT_(x) FT(x,1)
    #define RT_(x) RT(x,1)

    // hack to overcome external linkage conflict warning in MIPS
    #define __SGI_STL_INTERNAL_RELOPS

    #include<CGAL/basic.h>

    // only do something if LEDA available
    #ifdef CGAL_USE_LEDA
    
    #include<CGAL/Random.h>
    #include<CGAL/Cartesian.h>
    #include<CGAL/Homogeneous.h>
    #include<strstream>
    #include <cassert>
    #include<CGAL/Optimisation_d_traits_d.h>
    #include<CGAL/Min_sphere_d.h> 
    #include<CGAL/leda_rational.h>

    using namespace CGAL;

    @<Min_sphere_d test program - definitions>
    
    int main ()
    {
        @<Min_sphere_d test program - generate points>
        @<Min_sphere_d test program - test constructors>
        @<Min_sphere_d test program - test modifiers>
        @<Min_sphere_d test program - test access functions>
        @<Min_sphere_d test program - test predicates>
        @<Min_sphere_d test program - test I/O>

        return 0;
    }

    #else
    int main ()
    {
      return 0;
    }
    #endif

@end 

\subsection{The two-dimensional traits}
Here we basically aim at code coverage for the traits class -- the 
@prg{Min_sphere_d} functionality has extensively been tested above.

@macro <Min_sphere_d_traits_2 test program> many = @begin
    // hack to overcome internal compiler error of egcs
    #define FT_(x) FT(x,1)
    #define RT_(x) RT(x,1)

    // hack to overcome external linkage conflict warning in MIPS
    #define __SGI_STL_INTERNAL_RELOPS

    #include <CGAL/basic.h>

    // only do something if LEDA available
    #ifdef CGAL_USE_LEDA
    
    #include <cassert>
    #include<CGAL/Random.h>
    #include<CGAL/Cartesian.h>
    #include<CGAL/Homogeneous.h>
    #include<CGAL/Optimisation_d_traits_2.h>
    #include<CGAL/Min_sphere_d.h>
    #include<CGAL/Min_circle_2.h>
    #include<CGAL/Min_circle_2_traits_2.h>  //why is this not in Min_circle_2.h
    #include<CGAL/leda_rational.h>
    CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_rational);
    using namespace CGAL;

    typedef leda_rational                               NT;
    typedef NT                                          FT;
    typedef NT                                          RT;
    typedef Cartesian<FT>                               C;
    typedef Homogeneous<RT>                     H;
    typedef Optimisation_d_traits_2<C>          Cartesian_traits;
    typedef Optimisation_d_traits_2<H>          Homogeneous_traits;
    typedef Min_sphere_d<Cartesian_traits>      Min_sphereC;
    typedef Min_sphere_d<Homogeneous_traits>    Min_sphereH;
    typedef C::Point_2                          PointC;
    typedef H::Point_2                          PointH;


    // Min_circle_2 stuff
    typedef Min_circle_2_traits_2<C>            Circle_traitsC;
    typedef Min_circle_2_traits_2<H>            Circle_traitsH;
    typedef Min_circle_2<Circle_traitsC>                Min_circleC;
    typedef Min_circle_2<Circle_traitsH>                Min_circleH;

    // checked traits
    Cartesian_traits                                    tC;
    Homogeneous_traits                                  tH;

    const int n = 10;                           // number of points
    const int D = 2;
    const int r = 128;                          // coordinate range
    bool verbose = true;

    // NT == Quotient<NT> ?
    bool equals (NT x, Quotient<NT> y)
    {
        return (Quotient<NT>(x) == y);
    }

    // PointC == Point H ?
    bool equals (PointC p, PointH q) 
    {
        int dim =  p.dimension();
        if (q.dimension() != dim) return false;
        for (int j=0; j<dim; ++j)
           if (!equals (p.cartesian(j),q.cartesian(j)))
                return false;
        return true;
    }

    int main ()
    {
        Point_2<C>      PC[n];          // cartesian point set
        Point_2<H>      PH[n];          // homogeneous point set
        int             coord[D];       // arrays for random coordinates
        FT              coordC[D];      // ... cartesian version
        RT              coordH[D+1];    // ... homogeneous version
        Random          my_random;              // random number generator

        

        for (int i=0; i<n; ++i) {
            int j;
            // random coordinates
            for (j=0; j<D; ++j)
                coord[j] = my_random (r);

            // cartesian point
            for (j=0; j<D; ++j)
                coordC[j] = FT (coord[j]);
            PC[i] = Point_2<C>(coordC[0], coordC[1]);

            // homogeneous point
            for (j=0; j<D; ++j)
                coordH[j] = RT(2*coord[j]);
            coordH[D] = RT(2);
            PH[i] = Point_2<H>(coordH[0], coordH[1], coordH[2]);
        }

        Min_sphereC     msC (PC, PC+n, tC); assert (msC.is_valid(verbose));
        Min_sphereH     msH (PH, PH+n, tH); assert (msH.is_valid(verbose));

        PointC                  centerC (msC.center());
        PointH                  centerH (msH.center());

        FT                      radiusC (msC.squared_radius());
        Quotient<RT>    radiusH (msH.squared_radius());

        assert (equals (centerC, centerH));
        assert (equals (radiusC, radiusH));

        // test for equality with Min_circle, no randomization
        Min_circleC     mcC (PC, PC+n, false); 
        Min_circleH     mcH (PH, PH+n, false);

        assert (centerC == mcC.circle().center());
        assert (centerH == mcH.circle().center());
        assert (radiusC == mcC.circle().squared_radius());
        assert (radiusH == mcH.circle().squared_radius());

        return 0;
    }

    #else

    int main()
    {
      return 0;
    }

    #endif

@end

\subsection{The three-dimensional traits}
Here we basically aim at code coverage for the traits class -- the 
@prg{Min_sphere_d} functionality has extensively been tested above.

@macro <Min_sphere_d_traits_3 test program> many = @begin
    // hack to overcome internal compiler error of egcs
    #define FT_(x) FT(x,1)
    #define RT_(x) RT(x,1)

    // hack to overcome external linkage conflict warning in MIPS
    #define __SGI_STL_INTERNAL_RELOPS

    #include <CGAL/basic.h>

    // only do something if LEDA available
    #ifdef CGAL_USE_LEDA
    
    #include <cassert>
    #include<CGAL/Random.h>
    #include<CGAL/Cartesian.h>
    #include<CGAL/Homogeneous.h>
    #include<CGAL/Optimisation_d_traits_3.h>
    #include<CGAL/Min_sphere_d.h>
    #include<CGAL/leda_rational.h>
    CGAL_DEFINE_ITERATOR_TRAITS_POINTER_SPEC( leda_rational);
    using namespace CGAL;

    typedef leda_rational                               NT;
    typedef NT                                          FT;
    typedef NT                                          RT;
    typedef Cartesian<FT>                               C;
    typedef Homogeneous<RT>                     H;
    typedef Optimisation_d_traits_3<C>          Cartesian_traits;
    typedef Optimisation_d_traits_3<H>          Homogeneous_traits;
    typedef Min_sphere_d<Cartesian_traits>      Min_sphereC;
    typedef Min_sphere_d<Homogeneous_traits>    Min_sphereH;
    typedef C::Point_3                          PointC;
    typedef H::Point_3                          PointH;
    

    // checked traits
    Cartesian_traits                                    tC;
    Homogeneous_traits                                  tH;

    const int n = 10;                           // number of points
    const int D = 3;
    const int r = 128;                          // coordinate range
    bool verbose = true;

    // NT == Quotient<NT> ?
    bool equals (NT x, Quotient<NT> y)
    {
        return (Quotient<NT>(x) == y);
    }

    // PointC == Point H ?
    bool equals (PointC p, PointH q) 
    {
        int dim =  p.dimension();
        if (q.dimension() != dim) return false;
        for (int j=0; j<dim; ++j)
           if (!equals(p.cartesian(j),q.cartesian(j)))
                return false;
        return true;
    }

    int main ()
    {
        Point_3<C>      PC[n];          // cartesian point set
        Point_3<H>      PH[n];          // homogeneous point set
        int             coord[D];       // arrays for random coordinates
        FT              coordC[D];      // ... cartesian version
        RT              coordH[D+1];    // ... homogeneous version
        Random          my_random;              // random number generator

       
        for (int i=0; i<n; ++i) {
            int j;
            // random coordinates
            for (j=0; j<D; ++j)
                coord[j] = my_random (r);

            // cartesian point
            for (j=0; j<D; ++j)
                coordC[j] = FT (coord[j]);
            PC[i] = Point_3<C>(coordC[0], coordC[1], coordC[2]);

            // homogeneous point
            for (j=0; j<D; ++j)
                coordH[j] = RT(2*coord[j]);
            coordH[D] = RT(2);
            PH[i] = Point_3<H>(coordH[0], coordH[1], coordH[2], coordH[3]);
        }       

        Min_sphereC     msC (PC, PC+n, tC); assert (msC.is_valid(verbose));
        Min_sphereH     msH (PH, PH+n, tH); assert (msH.is_valid(verbose));

        PointC                  centerC (msC.center()); 
        PointH                  centerH (msH.center()); 

        FT                      radiusC (msC.squared_radius());
        Quotient<RT>    radiusH (msH.squared_radius());

        assert (equals (centerC, centerH));
        assert (equals (radiusC, radiusH));

        return 0;
    }

    #else

    int main()
    {
      return 0;
    }

    #endif

@end
        
@! ============================================================================
@! Files
@! ============================================================================
\section{File organisation}
\label{sec:files}


@! ----------------------------------------------------------------------------
@! main_Min_sphere_d.tex
@! ----------------------------------------------------------------------------
\subsection{main_Min\_sphere\_d.tex}

@file <main_Min_sphere_d.tex>  = @begin

% =============================================================================
% The CGAL Reference Manual
% Chapter: Geometric Optimisation
% Section: Smallest Enclosing Sphere
% -----------------------------------------------------------------------------
% file  : doc_tex/basic/Optimisation/Optimisation_ref/main_Min_sphere_d.tex
% author: Bernd Gärtner (gaertner@@inf.ethz.ch)
% -----------------------------------------------------------------------------
% $Revision$
% $Date$
% $CGAL_Package: Min_sphere_d WIP $
% =============================================================================
    
\newcommand{\cgalSetMinSphereLayout}{\ccTexHtml{%
    \ccSetThreeColumns{Support_point_iterator}{}{creates a variable
      \ccc{min_sphere} of type \ccc{CGAL_Min_sphere_d<Traits>}.}
      \ccPropagateThreeToTwoColumns}{}}
      
    \input{Optimisation_ref/Min_sphere_d}
@end

@! ----------------------------------------------------------------------------
@! Min_sphere_d.tex
@! ----------------------------------------------------------------------------
\subsection{Min\_sphere\_d.tex}

@file <Min_sphere_d.tex>  = @begin

% =============================================================================
% The CGAL Reference Manual
% Chapter: Geometric Optimisation
% Section: Smallest Enclosing Sphere
% -----------------------------------------------------------------------------
% file  : doc_tex/basic/Optimisation/Optimisation_ref/Min_sphere_d.tex
% author: Bernd Gärtner (gaertner@@inf.ethz.ch)
% -----------------------------------------------------------------------------
% $Revision$
% $Date$
% $CGAL_Package: Min_sphere_d WIP $
% =============================================================================

    \begin{ccRefClass}{Min_sphere_d<Traits>}
    \ccIndexSubitem[t]{sphere}{smallest enclosing}
    \ccIndexSubitem[t]{smallest enclosing}{sphere}
    \ccIndexSubitem[t]{bounding volumes}{smallest enclosing sphere}
    \ccIndexSubitemSeeAlso[t]{annulus}{smallest enclosing sphere}
    \ccIndexSubitemSeeAlso[t]{circle}{smallest enclosing sphere}
\cgalSetMinSphereLayout

% -----------------------------------------------------------------------------
\ccDefinition

An object of the class \ccRefName\ is the unique sphere of
smallest volume enclosing a finite (multi)set of points in $d$-dimensional
Euclidean space $\E_d$. For a set $P$ we denote by $ms(P)$ the
smallest sphere that contains all points of $P$. $ms(P)$ can
be degenerate, i.e.\ $ms(P)=\mbox{\ccTexHtml{$\;\emptyset$}{&Oslash;}}$
if $P=\mbox{\ccTexHtml{$\;\emptyset$}{&Oslash;}}$ and $ms(P)=\{p\}$ if
$P=\{p\}$.

An inclusion-minimal subset $S$ of $P$ with $ms(S)=ms(P)$ is called a
\emph{support set}, the points in $S$ are the \emph{support points}.
A support set has size at most $d+1$, and all its points lie on the
boundary of $ms(P)$. In general, neither the support set nor its size 
are unique.

The algorithm
computes a support set $S$ which remains fixed until the next insert
or clear operation.

    \ccInclude{CGAL/Min_sphere_d.h}

    \ccRequirements
    \ccIndexRequirements
    
    The class \ccRefName\ expects a model of the concept
    \ccc{OptimisationDTraits} as its template argument.
      We provide the models \ccc{Optimisation_d_traits_2},
      \ccc{Optimisation_d_traits_3} and \ccc{Optimisation_d_traits_d}
    for two-, three-, and $d$-dimensional points respectively. 

\ccTypes
\ccIndexClassTypes

\ccSetThreeColumns{typedef Traits::Point}{Sphere;}{}

\ccNestedType{Traits}{}
\ccNestedType{FT}{typedef to \ccc{Traits::FT}.}
\ccNestedType{Point}{typedef to \ccc{Traits::Point}.}


\ccSetTwoColumns{Min_sphere_d<Traits>:: Support_point_iterator}{}

\ccGlueBegin
\ccUnchecked
    \ccNestedType{Point_iterator}{non-mutable model of the STL
    concept \ccc{BidirectionalIterator} with value type \ccc{Point}. Used
    to access the points used to build the smallest enclosing sphere.}
    \ccNestedType{Support_point_iterator}{non-mutable model of the STL
    concept \ccc{BidirectionalIterator} with value type \ccc{Point}. Used
    to access the support points defining the smallest enclosing sphere.}
\ccGlueEnd

\cgalSetMinSphereLayout

% -----------------------------------------------------------------------------
\ccCreation
\ccIndexClassCreation
\ccCreationVariable{min_sphere}

\ccUnchecked
\ccConstructor{ Min_sphere_d (const Traits& traits = Traits());}{
        creates a variable of type \ccRefName\ and 
        initializes it to $ms(\mbox{\ccTexHtml{$\;\emptyset$}{&Oslash;}})$.
          If the traits parameter is not supplied, the class \ccc{Traits}
          must provide a default constructor.}

\ccUnchecked
\ccConstructor{ template < class InputIterator >
                Min_sphere_d( InputIterator  first,
                                   InputIterator  last,
                                   const Traits&  traits = Traits());}{
        creates a variable \ccVar\ of type \ccRefName.
        It is initialized to $ms(P)$ with $P$ being the set of points
        in the range [\ccc{first},\ccc{last}).
        \ccRequire The value type of \ccc{first} and \ccc{last}
        is \ccc{Point}. If the traits parameter is not supplied,
        the class \ccc{Traits} must provide a default constructor.
        \ccPrecond All points have the same dimension.}

\ccUnchecked
\ccHidden
\ccConstructor{ Min_sphere_d( const Min_sphere_d<Traits>&);}{
        copy constructor.}

\ccHidden
\ccConstructor{ ~Min_sphere_d( );}{
        destructor.}

\ccUnchecked
\ccHidden
\ccMemberFunction{ Min_sphere_d<Traits>&
                   operator = ( const Min_sphere_d<Traits>&);}{
        assignment operator.}

% -----------------------------------------------------------------------------\ccAccessFunctions
\begin{ccIndexMemberFunctions}
\ccIndexMemberFunctionGroup{access}
    
\ccIndexSubitem[t]{support set}{\ccFont Min_sphere_d}
\ccMemberFunction{ int  number_of_points( ) const;}{
        returns the number of points of \ccVar, i.e.\ $|P|$.}

\ccMemberFunction{ int  number_of_support_points( ) const;}{
        returns the number of support points of \ccVar, i.e.\ $|S|$.}

\ccGlueBegin
\ccMemberFunction{ Point_iterator  points_begin() const;}{
        returns an iterator referring to the first point of \ccVar.}
%
\ccMemberFunction{ Point_iterator  points_end() const;}{
        returns the corresponding past-the-end iterator.}
\ccGlueEnd

\ccGlueBegin
\ccMemberFunction{ Support_point_iterator  support_points_begin() const;}{
        returns an iterator referring to the first support point of \ccVar.}
%
\ccMemberFunction{ Support_point_iterator  support_points_end() const;}{
        returns the corresponding past-the-end iterator.}
\ccGlueEnd

\ccMemberFunction{ int ambient_dimension() const;}{
        returns the dimension of the points in $P$. If \ccVar\
        is empty, the ambient dimension is $-1$.}

\ccMemberFunction{ const Point&  center( ) const;}{
        returns the center of \ccVar.
        \ccPrecond \ccVar\ is not empty.}

        \ccMemberFunction{FT squared_radius( ) const;}{
        returns the squared radius of \ccVar.
        \ccPrecond \ccVar\ is not empty.}

% -----------------------------------------------------------------------------
\ccPredicates
\ccIndexMemberFunctionGroup{predicates}

By definition, an empty \ccRefName\ has no boundary and no
bounded side, i.e.\ its unbounded side equals the whole space $\E_d$.

\ccMemberFunction{ Bounded_side
                   bounded_side( const Point& p) const;}{
        returns \ccc{CGAL::ON_BOUNDED_SIDE}, \ccc{CGAL::ON_BOUNDARY}, or
        \ccc{CGAL::ON_UNBOUNDED_SIDE} iff \ccc{p} lies properly inside,
        on the boundary, or properly outside of \ccVar, resp.
        \ccPrecond if \ccVar\ is not empty, the dimension of $p$ 
        equals \ccc{ambient_dimension()}.}

\ccMemberFunction{ bool  has_on_bounded_side( const Point& p) const;}{
        returns \ccc{true}, iff \ccc{p} lies properly inside \ccVar.
        \ccPrecond if \ccVar\ is not empty, the dimension of $p$ 
        equals \ccc{ambient_dimension()}.}

\ccMemberFunction{ bool  has_on_boundary( const Point& p) const;}{
        returns \ccc{true}, iff \ccc{p} lies on the boundary
        of \ccVar.
        \ccPrecond if \ccVar\ is not empty, the dimension of $p$ 
        equals \ccc{ambient_dimension()}.}

\ccMemberFunction{ bool  has_on_unbounded_side( const Point& p) const;}{
        returns \ccc{true}, iff \ccc{p} lies properly outside of \ccVar.
        \ccPrecond if \ccVar\ is not empty, the dimension of $p$ 
        equals \ccc{ambient_dimension()}.}

\ccMemberFunction{ bool  is_empty( ) const;}{
        returns \ccc{true}, iff \ccVar\ is empty (this implies
        degeneracy).}

\ccIndexSubitem[t]{degeneracies}{\ccFont Min_sphere_d}
\ccMemberFunction{ bool  is_degenerate( ) const;}{
        returns \ccc{true}, iff \ccVar\ is degenerate, i.e.\ if
        \ccVar\ is empty or equal to a single point, equivalently if
        the number of support points is less than 2.}

% -----------------------------------------------------------------------------
\ccModifiers
\ccIndexMemberFunctionGroup{modifiers}
    
\ccMemberFunction{ void clear ();}{
        resets \ccVar\ to $ms(\mbox{\ccTexHtml{$\;\emptyset$}{&Oslash;}})$.}

\ccUnchecked
\ccMemberFunction{ template < class InputIterator >
                   void set( InputIterator first,
                             InputIterator last );}{
        sets \ccVar\ to the $ms(P)$, where $P$ is the set of points
        in the range [\ccc{first},\ccc{last}).
        \ccRequire The value type of \ccc{first} and \ccc{last} is
        \ccc{Point}.
        \ccPrecond All points have the same dimension.}

\ccMemberFunction{ void  insert( const Point& p);}{
        inserts \ccc{p} into \ccVar. If \ccc{p} lies inside the
        current sphere, this is a constant-time operation, otherwise
        it might take longer, but usually substantially less than  
        recomputing the smallest enclosing sphere from scratch.
        \ccPrecond The dimension of \ccc{p} equals \ccc{ambient_dimension()}
        if \ccVar\ is not empty.}

\ccUnchecked
\ccMemberFunction{ template < class InputIterator >
                   void  insert( InputIterator  first,
                                 InputIterator  last );}{
        inserts the points in the range [\ccc{first},\ccc{last})
        into \ccVar\ and recomputes the smallest enclosing sphere, by
        calling \ccc{insert} for all points in the range.
        \ccRequire  The value type of \ccc{first} and \ccc{last} is
      \ccc{Point}.
        \ccPrecond All points have the same dimension. If
        \ccVar\ is not empty, this dimension must be equal to
        \ccc{ambient_dimension()}.}

% -----------------------------------------------------------------------------
    \ccHeading{Validity Check}
    \ccIndexMemberFunctionGroup{validity check}
  \ccIndexSubitem[t]{validity check}{\ccFont Min_sphere_d}
  An object \ccVar\ is valid, iff
  \begin{itemize}
    \item \ccVar\ contains all points of its defining set $P$,
    \item \ccVar\ is the smallest sphere containing its support set $S$, and
    \item $S$ is minimal, i.e.\ no support point is redundant.
  \end{itemize}

  \emph{Note:} Under inexact arithmetic, the result of the
  validation is not realiable, because the checker itself can suffer
  from numerical problems.

  \ccMemberFunction{ bool is_valid( bool verbose = false,
                                    int  level   = 0    ) const;}{
        returns \ccc{true}, iff \ccVar\ is valid. If \ccc{verbose}
        is \ccc{true}, some messages concerning the performed checks
        are written to standard error stream. The second parameter
        \ccc{level} is not used, we provide it only for consistency
        with interfaces of other classes.}

% -----------------------------------------------------------------------------
\ccHeading{Miscellaneous}
\ccIndexMemberFunctionGroup{miscellaneous}

\def\ccTagRmConstRefPair{\ccFalse}
\ccMemberFunction{ const Traits&  traits( ) const;}{
        returns a const reference to the traits class object.}
        \def\ccTagRmConstRefPair{\ccTrue}

 \end{ccIndexMemberFunctions}
       
% -----------------------------------------------------------------------------
\ccHeading{I/O}
\begin{ccIndexGlobalFunctions}
    
\ccIndexGlobalFunctionGroup{output}
\ccFunction{ std::ostream& operator << ( std::ostream& os,
                                    const Min_sphere_d<Traits>&
                                        min_sphere);}{
        writes \ccVar\ to output stream \ccc{os}.
        \ccRequire The output operator is defined for \ccc{Point}.}

\ccIndexGlobalFunctionGroup{input}
        
\ccFunction{ std::istream& operator >> ( std::istream& is,
                                    Min_sphere_d<Traits> min_sphere&);}{
        reads \ccVar\ from input stream \ccc{is}.
          \ccRequire The input operator is defined for \ccc{Point}.}
    
\end{ccIndexGlobalFunctions}

% -----------------------------------------------------------------------------
\ccSeeAlso

    \ccRefIdfierPage{CGAL::Optimisation_d_traits_2<R,ET,NT>}\\
    \ccRefIdfierPage{CGAL::Optimisation_d_traits_3<R,ET,NT>}\\
    \ccRefIdfierPage{CGAL::Optimisation_d_traits_d<R,ET,NT>}\\
    \ccRefConceptPage{OptimisationDTraits}\\
    \ccRefIdfierPage{CGAL::Min_circle_2<Traits>}\\
    \ccRefIdfierPage{CGAL::Min_annulus_d<Traits>}

% -----------------------------------------------------------------------------
\ccImplementation
\ccIndexImplementation

    \ccIndexSubitem[t]{incremental algorithm}{\ccFont Min_sphere_d}
    \ccIndexSubitem[t]{move-to-front heuristic}{\ccFont Min_sphere_d}
We implement the algorithm of Welzl with move-to-front
heuristic~\cite{w-sedbe-91a} for small point sets, combined with a new
efficient method for large sets, which is particularly tuned for
moderately large dimension ($d \leq 20$) \cite{g-frseb-99}.
The creation time is almost
always linear in the number of points. Access functions and predicates
take constant time, inserting a point might take up to linear time,
but substantially less than computing the new smallest enclosing
sphere from scratch. The clear operation and the check for validity
    each take linear time.

    \ccExample
    \ccIncludeVerbatim{Optimisation_ref/min_sphere_d_example.C}
    
\end{ccRefClass}

% ===== EOF ===================================================================

@end

@! ----------------------------------------------------------------------------
@! min_sphere_d_example.C
@! ----------------------------------------------------------------------------

@file <min_sphere_d_example.C> = @begin

#include<CGAL/Cartesian.h>
#include<iostream>
#include<cstdlib>
#include<CGAL/Random.h>
#include<CGAL/Optimisation_d_traits_d.h>    
#include<CGAL/Min_sphere_d.h>

typedef CGAL::Cartesian<double>                R;
typedef CGAL::Optimisation_d_traits_d<R>       Traits;              
typedef CGAL::Min_sphere_d<Traits>             Min_sphere;
typedef R::Point_d                             Point;

const int n = 10;                        // number of points
const int d = 5;                         // dimension of points    

int main ()
{
    Point         P[n];                  // n points
    double        coord[d];              // d coordinates
    CGAL::Random  r;                     // random number generator

    for (int i=0; i<n; ++i) {
        for (int j=0; j<d; ++j)
            coord[j] = r.get_double();       
        P[i] = Point(d, coord, coord+d); // random point
    } 

    Min_sphere  ms (P, P+n);             // smallest enclosing sphere

    CGAL::set_pretty_mode (std::cout);
    std::cout << ms;                     // output the sphere

    return 0;
}

@end 



@! ----------------------------------------------------------------------------
@! Min_sphere_d.h
@! ----------------------------------------------------------------------------
\subsection{Min\_sphere\_d.h}

@file <include/CGAL/Min_sphere_d.h> = @begin
    @<Min_sphere_d header>("include/CGAL/Min_sphere_d.h")

    #ifndef CGAL_MIN_SPHERE_D_H
    #define CGAL_MIN_SPHERE_D_H

    // Class declarations
    // ==================

    // Class interface and implementation
    // ==================================
    // includes

    #  include <CGAL/basic.h>
   
    #  include <CGAL/Optimisation/assertions.h>
    
    #  include <CGAL/Optimisation/basic.h>
    
    #  include <CGAL/Optimisation_sphere_d.h>
    

    #ifndef CGAL_PROTECT_LIST
    #  include <list>
    #  define CGAL_PROTECT_LIST
    #endif

    #ifndef CGAL_PROTECT_IOSTREAM
    #  include <iostream>
    #  define CGAL_PROTECT_IOSTREAM
    #endif

    #ifdef _MSC_VER
    #ifdef _VIRTUAL
    #pragma message("Min_sphere_d.h: this code may not compile due to a problem")
    #pragma message("in xlocnum. Please consult the installation guide for a fix.")
    #endif
    #endif

    CGAL_BEGIN_NAMESPACE


    template <class Traits>
    class Min_sphere_d
    {
        @<Min_sphere_d types>
        @<Min_sphere_d data members>
    public:
        @<Min_sphere_d constructors>
        @<Min_sphere_d access functions>
        @<Min_sphere_d predicates>
        @<Min_sphere_d modifiers>
        @<Min_sphere_d validity check>
        @<Min_sphere_d miscellaneous>
    private:
        @<Min_sphere_d private method mtf_ms>
        @<Min_sphere_d private method pivot_ms>
        @<Min_sphere_d private method move_to_front>
        @<Min_sphere_d private method all_points_have_dim>
    };

    // Function declarations
    // =====================
    // I/O
    // ---

    @<Min_sphere_d I/O operators declaration>

    CGAL_END_NAMESPACE

    #ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
    #  include <CGAL/Min_sphere_d.C>
    #endif      

    #endif // CGAL_MIN_SPHERE_D_H

    @<end of file line>

@end

@! ----------------------------------------------------------------------------
@! Min_sphere_d.C
@! ----------------------------------------------------------------------------
\subsection{Min\_sphere\_d.C}

@file <include/CGAL/Min_sphere_d.C> = @begin
    @<Min_sphere_d header>("include/CGAL/Min_sphere_d.C")

    CGAL_BEGIN_NAMESPACE

    // Class implementation (continued)
    // ================================
    // I/O
    // ---
    @<Min_sphere_d I/O operators>

    CGAL_END_NAMESPACE

    @<end of file line>

@end

@! ----------------------------------------------------------------------------
@! Optimisation_sphere_d.h
@! ----------------------------------------------------------------------------
\subsection{Optimisation\_sphere\_d.h}

@file <include/CGAL/Optimisation_sphere_d.h> = @begin
    @<Min_sphere_d header>("include/CGAL/Optimisation_sphere_d.h")

    #ifndef CGAL_OPTIMISATION_SPHERE_D_H
    #define CGAL_OPTIMISATION_SPHERE_D_H

    CGAL_BEGIN_NAMESPACE

    // Class declarations
    // ==================
    // general template
    template <class Rep_tag, class FT, class RT, class PT, class Traits>
    class Optimisation_sphere_d;

#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

    @<Optimisation_sphereHd declaration>
    @<Optimisation_sphereCd declaration>

#else

   template <class FT, class RT, class PT, class Traits>
   class Optimisation_sphereCd;

   template < class FT, class RT, class PT, class Traits>
   class Optimisation_sphereHd;

#endif

    CGAL_END_NAMESPACE

    // Class interfaces and implementation
    // ==================================
    // includes

    
    #include <CGAL/Optimisation/basic.h>
    
    #include <CGAL/Optimisation/assertions.h>
   

    CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_PARTIAL_CLASS_TEMPLATE_SPECIALISATION

    // Cartesian version
    // -----------------
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphere_d<Cartesian_tag, FT, RT, PT, Traits>
    {
    private:
        @<Optimisation_sphereCd types and data members>

    public:

        Optimisation_sphere_d& get_sphere (Cartesian_tag t)
        {return *this;}

        const Optimisation_sphere_d&
          get_sphere (Cartesian_tag t) const
        {return *this;}

        @<Optimisation_sphereCd constructor>
        @<Optimisation_sphereCd init method>
        @<Optimisation_sphereCd destructor>
        @<Optimisation_sphereCd destroy method>
        @<Optimisation_sphereCd set_size method>
        @<Optimisation_sphereCd push method>
        @<Optimisation_sphereCd pop method>
        @<Optimisation_sphereCd excess method>
        @<Optimisation_sphereCd access methods>
        @<Optimisation_sphereCd is_valid method>
    private:
        @<Optimisation_sphereCd multiply method>
        @<Optimisation_sphereCd compute_c_and_sqr_r method>
        @<Optimisation_sphereCd prod method>
    };


    // Homogeneous version
    // -----------------
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphere_d<Homogeneous_tag, FT, RT, PT, Traits>
    {
    private:
        @<Optimisation_sphereHd types and data members>

    public:

        Optimisation_sphere_d& get_sphere (Homogeneous_tag t)
        {return *this;}

        const Optimisation_sphere_d&
          get_sphere (Homogeneous_tag t) const
        {return *this;}

        @<Optimisation_sphereHd constructor>
        @<Optimisation_sphereHd init method>
        @<Optimisation_sphereHd destructor>
        @<Optimisation_sphereHd destroy method>
        @<Optimisation_sphereHd set_size method>
        @<Optimisation_sphereHd push method>
        @<Optimisation_sphereHd pop method>
        @<Optimisation_sphereHd excess method>
        @<Optimisation_sphereHd access methods>
        @<Optimisation_sphereHd is_valid method>
    private:
        @<Optimisation_sphereHd multiply method>
        @<Optimisation_sphereHd compute_c_and_sqr_r method>
        @<Optimisation_sphereHd prod method>
    };

#else

    // general template 
    template <class Rep_tag, class FT, class RT, class PT,class Traits>
    class Optimisation_sphere_d
    {
        Optimisation_sphereCd<FT,RT,PT,Traits> ms_cart;
        Optimisation_sphereHd<FT,RT,PT,Traits> ms_hom;

    public:

        Optimisation_sphere_d (Cartesian_tag t, const Traits& tr = Traits())
        : ms_cart (tr) {}

        Optimisation_sphere_d (Homogeneous_tag t, const Traits& tr = Traits())
        : ms_hom (tr) {}

        Optimisation_sphereCd<FT,RT,PT,Traits>& get_sphere (Cartesian_tag t)
        {return ms_cart;}

        Optimisation_sphereHd<FT,RT,PT,Traits>& get_sphere (Homogeneous_tag t)
        {return ms_hom;}

        const Optimisation_sphereCd<FT,RT,PT,Traits>& get_sphere
        (Cartesian_tag t) const
        {return ms_cart;}

        const Optimisation_sphereHd<FT,RT,PT,Traits>& get_sphere
        (Homogeneous_tag t) const
        {return ms_hom;}

    };

    // Cartesian version
    // -----------------
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphereCd

    #define Optimisation_sphere_d Optimisation_sphereCd
    {
    private:
        @<Optimisation_sphereCd types and data members>

    public:
        @<Optimisation_sphereCd constructor>
        @<Optimisation_sphereCd init method>
        @<Optimisation_sphereCd destructor>
        @<Optimisation_sphereCd destroy method>
        @<Optimisation_sphereCd set_size method>
        @<Optimisation_sphereCd push method>
        @<Optimisation_sphereCd pop method>
        @<Optimisation_sphereCd excess method>
        @<Optimisation_sphereCd access methods>
        @<Optimisation_sphereCd is_valid method>
    private:
        @<Optimisation_sphereCd multiply method>
        @<Optimisation_sphereCd compute_c_and_sqr_r method>
        @<Optimisation_sphereCd prod method>
    };

    #undef Optimisation_sphere_d

   
    // Homogeneous version
    // -----------------
    template <class FT, class RT, class PT, class Traits>
    class Optimisation_sphereHd
    
    #define Optimisation_sphere_d Optimisation_sphereHd
    {
    private:
        @<Optimisation_sphereHd types and data members>

    public:
        @<Optimisation_sphereHd constructor>
        @<Optimisation_sphereHd init method>
        @<Optimisation_sphereHd destructor>
        @<Optimisation_sphereHd destroy method>
        @<Optimisation_sphereHd set_size method>
        @<Optimisation_sphereHd push method>
        @<Optimisation_sphereHd pop method>
        @<Optimisation_sphereHd excess method>
        @<Optimisation_sphereHd access methods>
        @<Optimisation_sphereHd is_valid method>
    private:
        @<Optimisation_sphereHd multiply method>
        @<Optimisation_sphereHd compute_c_and_sqr_r method>
        @<Optimisation_sphereHd prod method>
    };

    #undef Optimisation_sphere_d

#endif 

     CGAL_END_NAMESPACE

    #endif // CGAL_OPTIMISATION_SPHERE_D_H

   

    @<end of file line>

@end

@! ----------------------------------------------------------------------------
@! min_sphere_test.C
@! ----------------------------------------------------------------------------
\subsection{min\_sphere\_test.C}

@file <min_sphere_test.C>  = @begin
    @<Min_sphere_d header>("min_sphere_test.C")
    @<Min_sphere_d test program>

    @<end of file line>

@end

@! ----------------------------------------------------------------------------
@! min_sphere_traits_2_test.C
@! ----------------------------------------------------------------------------
\subsection{min\_sphere\_traits\_2\_test.C}

@file <min_sphere_traits_2_test.C> = @begin
    @<Min_sphere_d header>("min_sphere_traits_2_test.C")
    @<Min_sphere_d_traits_2 test program>

    @<end of file line>

@end

@! ----------------------------------------------------------------------------
@! min_sphere_traits_3_test.C
@! ----------------------------------------------------------------------------
\subsection{min\_sphere\_traits\_3\_test.C}

@file <min_sphere_traits_3_test.C> = @begin
    @<Min_sphere_d header>("min_sphere_traits_3_test.C")
    @<Min_sphere_d_traits_3 test program>

    @<end of file line>

@end
 
\subsection*{File Header}

A formatted file header allows easy identification of the files. It is
parameterized with the title of the implementation, the product file
name, the source file name, the author name, and the RCS variables
@em{Revision} and @em{Date} of the source file.

@macro <file header>(8) zero many = @begin
// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-wip $
// release_date  : $CGAL_Date$
//
// chapter       : $CGAL_Chapter: Optimisation $
// package       : $CGAL_Package: MinSphere $
// file          : @2
// source        : web/@3.aw
// revision      : @7
// revision_date : @8
// author(s)     : @4
//                 @5
//
// coordinator   : @6
//
// implementation: @1
// ============================================================================
@end

@macro <dividing line> zero many = @begin
// ============================================================================
@end
 
@macro <end of file line> zero many = @begin
// ===== EOF ==================================================================
@end
 
@macro <Min_sphere_d header>(1) many = @begin
    @<file header>("dD Smallest Enclosing Sphere",@1,
                   "Optimisation/Min_sphere_d",
                   "Sven Schönherr <sven@@inf.fu-berlin.de>",
                   "Bernd Gärtner",
                   "ETH Zurich (Bernd Gärtner <gaertner@@inf.ethz.ch>)",
                   "$Revision$","$Date$")
@end



@! ============================================================================
@! Bibliography
@! ============================================================================

\clearpage
\bibliographystyle{plain}
\bibliography{geom,cgal}

@! ===== EOF ==================================================================
