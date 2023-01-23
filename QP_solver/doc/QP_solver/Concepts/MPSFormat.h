
/*!
\ingroup PkgQPSolverConcepts
\cgalConcept

MPS is a commonly used file format for storing linear and quadratic
programs according to the concepts `QuadraticProgram`,
`LinearProgram`,
`NonnegativeQuadraticProgram`, and
`NonnegativeLinearProgram`, see also
<A HREF="https://en.wikipedia.org/wiki/MPS_(format)"><TT>https://en.wikipedia.org/wiki/MPS_(format)</TT></A>.

\cgal supports a large subset of
this format, but there are MPS files around that we cannot read (for
example, files that encode integrality constraints on the variables).
Also, there might be some other MPS-based solvers that will not be able
to read the MPS files written by \cgal, since we do not strictly
adhere to the very rigid layout requirements of the original MPS
format.

Let's look at an example first. The quadratic program

\f[
\begin{array}{lrcl}
\mbox{minimize} & x^2 + 4(y-4)^2 &(=& x^2 + 4y^2 - 32y + 64) \\
\mbox{subject to} & x + y &\leq& 7 \\
& -x + 2y &\leq& 4 \\
& x &\geq& 0 \\
& y &\geq& 0 \\
& y &\leq& 4
\end{array}
\f]

has the following description in MPS format.

\verbatim
NAME first_qp
ROWS
N obj
L c0
L c1
COLUMNS
x0 c0 1
x0 c1 -1
x1 obj -32
x1 c0 1
x1 c1 2
RHS
rhs obj -64
rhs c0 7
rhs c1 4
BOUNDS
UP BND x1 4
QMATRIX
x0 x0 2
x1 x1 8
\endverbatim

Here comes a semiformal description of the format in general.

\cgalHeading{NAME Section}

This (mandatory) section consists of a single line
starting with <TT>NAME</TT>. Everything starting from the
first non-whitespace after that until the end of the line
constitutes the name of the problem.

\cgalHeading{ROWS Section}

In the (mandatory) <TT>ROW</TT> section, you find one line for every
constraint, where the letter <TT>L</TT> indicates relation \f$ \leq\f$,
letter <TT>G</TT> stands for \f$ \geq\f$, and <TT>E</TT> for \f$ =\f$. In
addition, there is a row for the linear objective function (indicated
by letter <TT>N</TT>). In that section, names are assigned to the
constraints (here: <TT>c0, c1</TT>) and the objective function (here:
<TT>obj</TT>). An MPS file may encode several linear objective
functions by using several rows starting with <TT>N</TT>, but we ignore
all but the first.

\cgalHeading{COLUMNS Section}

The (mandatory) <TT>COLUMNS</TT> section encodes the constraint matrix
\f$ A\f$ and the linear objective function vector \f$ c\f$. Every line consists
of one or two sequences of three tokens \f$ j i val\f$, where \f$ j\f$ is the
name of a variable (here, we have variables <TT>x0,x1</TT>), \f$ i\f$ is
the name of a constraint or the objective function, and \f$ val\f$ is the
value \f$ A_{ij}\f$ (if \f$ i\f$ names a constraint), or \f$ c_j\f$ (if \f$ i\f$ names the
linear objective function). Values for pairs \f$ (i,j)\f$ that are not
specified in this section default to \f$ 0\f$. Otherwise, for every pair
\f$ (i,j)\f$, the <I>last</I> specified value determines \f$ A_{ij}\f$ or \f$ c_j\f$.

\cgalHeading{RHS Section}

This (mandatory) section encodes the right-hand side vector \f$ b\f$ and
the constant term \f$ c_0\f$ in the objective function. The first token in
every line is an identifier (here: <TT>rhs</TT>). An MPS file may
encode several right-hand sides \f$ b\f$ by using several such identifiers,
but we ignore all lines having an identifier different from that of
the first line.

The right-hand side identifier is succeeded by one or two sequences
of tokens \f$ i val\f$, where \f$ i\f$ names a constraint or the linear
objective function, and \f$ val\f$ specifies the value \f$ b_i\f$ (if
\f$ i\f$ names a constraint), or \f$ -c_0\f$ (if \f$ i\f$ names the linear objective
function). Values that are not specified in this section default to \f$ 0\f$.
Otherwise, for every \f$ i\f$, the <I>last</I> specified value determines
\f$ b_{i}\f$ or \f$ -c_0\f$.

\cgalHeading{BOUNDS Section}

This (optional) section encodes the lower and upper bound vectors \f$ l\f$
and \f$ u\f$ for the variables. The default bounds for any variable \f$ x_j\f$ are
\f$ 0\leq x_j\leq \infty\f$; the
<TT>BOUNDS</TT> section is used to override these defaults. In particular,
if there is no <TT>BOUNDS</TT> section, the program is nonnegative and
actually a model of the concept `NonnegativeQuadraticProgram`
or `NonnegativeLinearProgram`.

The first token in every line is succeeded by an (optional) identifier
(here: <TT>BND</TT>). An MPS file may encode several bound vectors \f$ l\f$
and \f$ u\f$ by using several such identifiers, but we ignore all lines
having an identifier different from that of the first line. The first
token \f$ t\f$ itself determines the type of the bound, and the token \f$ j\f$
after the bound identifier names the variable to which the bound applies
In case of bound types <TT>FX</TT>, <TT>LO</TT>, and
<TT>UP</TT>, there is another token \f$ val\f$ that specifies the bound
value. Here is how bound type and value determine a bound for variable
\f$ x_j\f$. There may be several bound specifications for a single variable, and
they are processed in order of appearance.

bound type | resulting bound
---------- | ----------------
FX         | \f$x_j = val\f$ (\f$x_j\f$ becomes a fixed variable)
LO         | \f$x_j \geq val\f$ (upper bound remains unchanged)
UP         | \f$x_j \leq val\f$ (lower bound remains unchanged, except if \f$val<0\f$; then, a zero lower bound is reset to \f$-\infty\f$)
FR         | \f$-\infty \leq x_j\leq\infty\f$ (previous bounds are discarded)
MI         | \f$x_j\geq -\infty\f$ (upper bound remains unchanged)
PL         | \f$x_j\leq \infty\f$ (lower bound remains unchanged)

\cgalHeading{QMATRIX / QUADOBJ / DMATRIX Section}

This (optional) section encodes the quadratic objective
function matrix \f$ D\f$. Every line is a sequence \f$ i j val\f$ of
three tokens, where both \f$ i\f$ and \f$ j\f$ name variables, and
\f$ val\f$ is the value \f$ 2D_{i,j}\f$ (in case of <TT>QMATRIX</TT>
or <TT>QUADOBJ</TT>), or \f$ D_{ij}\f$ (in case of <TT>DMATRIX</TT>).

In case of <TT>QMATRIX</TT> and <TT>DMATRIX</TT>, <I>all</I> nonzero
entries must be specified: if there is a line \f$ i j val\f$, then there
must also be a line \f$ j i val\f$, since \f$ D\f$ is required to be symmetric.
In case of <TT>QUADOBJ</TT>, only the entries of \f$ 2D\f$ on or below the
diagonal must be specified, entries above the diagonal are deduced from
symmetry. It is not allowed to specify two or more <I>different</I>
nonzero values for an unordered pair \f$ \{i,j\}\f$.

If this section is missing or does not contain nonzero values, the
program is a model of the concept `LinearProgram`.

\cgalHeading{Miscellaneous}

Our MPS format also supports an (optional) <TT>RANGES</TT> section,
but we don't explain this here.

\sa `CGAL::Quadratic_program_from_mps<NT>`

*/

class MPSFormat {
public:

}; /* end MPSFormat */

