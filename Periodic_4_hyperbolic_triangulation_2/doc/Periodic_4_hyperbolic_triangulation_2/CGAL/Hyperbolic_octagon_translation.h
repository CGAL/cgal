// Copyright (c) 1999-2018   INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2HyperbolicTranslationsClasses

The class `Hyperbolic_octagon_translation` defines an object to represent a hyperbolic
translation of the fundamental group of the Bolza surface \f$\mathcal M\f$. It accepts
one template parameter:

\tparam FT         %Field number type. Must provide exact computations with algebraic numbers,
                        notably with nested square roots. The default value for this parameter is
                        `CORE::Expr`.

A translation \f$g\f$ in \f$\mathcal G\f$ is a mapping acting on the hyperbolic plane
\f$\mathbb H^2\f$. It has the form
\f[ g(z) = \frac{ \alpha\cdot z + \beta }{ \overline{\beta}\cdot z + \overline{\alpha} }, \qquad
        \alpha,\beta \in \mathbb C, \qquad z \in \mathbb H^2, \qquad |\alpha|^2 - |\beta|^2 = 1, \f]
where \f$\overline{\alpha}\f$ ane \f$\overline{\beta}\f$ are the complex conjugates of
\f$\alpha\f$ and \f$\beta\f$ respectively. In this implementation, the translation \f$g\f$
is uniquely defined by its coefficients \f$\alpha\f$ and \f$\beta\f$.

Considering the set of generators \f$\mathcal A = [a, \overline{b}, c, \overline{d}, \overline{a},
b, \overline{c}, d]\f$ as an alphabet, a translation \f$g\f$ in \f$\mathcal G\f$ can be seen as a
word on the alphabet \f$\mathcal A\f$. Each letter of this alphabet is represented as an unsigned
integer from from 0 to 7, and each word (i.e., translation) is a sequence of letters.
*/

template <typename FT = CORE::Expr>
class Hyperbolic_octagon_translation {

public:


        /// \name Types

        /*!
                %Field number type of the coefficients of \f$\alpha\f$ and \f$\beta\f$.
                 This number type must be the same as the field number type `FT` of the
                 class `Periodic_4_hyperbolic_Delaunay_triangulation_2_traits`.
                 \sa `alpha()`
                 \sa `beta()`
        */
        typedef unspecified_type                         FT;

        /*!
                Represents a single letter of the alphabet \f$\mathcal A\f$.
                By extension, represents a generator of the group \f$\mathcal G\f$.
        */
        typedef unsigned short int                         Word_letter;

        /*!
                Represents a word on the alphabet \f$\mathcal A\f$. By extension, represents
                a hyperbolic translation in the group \f$\mathcal A\f$.
        */
        typedef std::vector<Word_letter>         Word;

        /*!
                Enumeration type for the alphabet \f$\mathcal A\f$. This enumeration can
                be used to recover the generators of the group \f$\mathcal G\f$.
                \sa `generator()`
        */
        enum Generator: Word_letter {
                /*! translation \f$a\f$ */
                A                 = 0,
                /*! translation \f$\overline{b}\f$ */
                B_BAR         = 1,
                /*! translation \f$c\f$ */
                C                 = 2,
                /*! translation \f$\overline{d}\f$ */
                D_BAR         = 3,
                /*! translation \f$\overline{a}\f$ */
                A_BAR         = 4,
                /*! translation \f$b\f$ */
                B                  = 5,
                /*! translation \f$\overline{c}\f$ */
                C_BAR         = 6,
                /*! translation \f$d\f$ */
                D                 = 7
        };

        /// \name Creation

        /// Default constructor. Creates the identity translation of the group \f$\mathcal G\f$.
        Hyperbolic_octagon_translation();

        /// Creates the translation described by the word `w`.
        Hyperbolic_octagon_translation(Word w);

        /// Creates the translation described by the one-letter word `l`.
        Hyperbolic_octagon_translation(Word_letter l);

        /// Creates the translation described by the two-letter word `lm`.
        Hyperbolic_octagon_translation(Word_letter l, Word_letter m);

        /// Creates the translation described by the three-letter word `lmn`.
        Hyperbolic_octagon_translation(Word_letter l, Word_letter m, Word_letter n);

        /// Creates the translation described by the four-letter word `lmno`.
        Hyperbolic_octagon_translation(Word_letter l, Word_letter m, Word_letter n, Word_letter o);


        /// \name Operations

        /// Multiplication operator; composes the current translation with `other`.
        Hyperbolic_octagon_translation operator*(const Hyperbolic_octagon_translation& other) const;

        ///        Difference operator; the difference of two translations \f$v\f$ and \f$w\f$ is defined as \f$v * w^{-1}\f$.
          Hyperbolic_octagon_translation operator-(const Hyperbolic_octagon_translation& other) const;

          /// Assignment operator; modifying the translation after the assignment leaves `other` unaffected.
        Hyperbolic_octagon_translation& operator=(const Hyperbolic_octagon_translation& other);

        /// Equality operator.
        bool operator==(const Hyperbolic_octagon_translation<FT>& other) const;

        /// Inequality operator.
        bool operator!=(const Hyperbolic_octagon_translation<FT>& other) const;

        /*!
                Comparison operator.
                 Each translation \f$g\f$ of \f$\mathcal G\f$, when applied to the octagon \f$\mathcal D_O\f$,
                 produces a copy of \f$\mathcal D_O\f$ labeled by the translation \f$g\f$. The copies of
                 \f$\mathcal D_O\f$ incident to \f$\mathcal D_O\f$ are naturally ordered counter-clockwise around
                 \f$\mathcal D_O\f$. The comparison operator compares two translations based on the ordering
                 of the copies of \f$\mathcal D_O\f$ that they produce. For more details, see Section
                 \ref P4HT2_representation of the User manual.
        */
        bool operator<(const Hyperbolic_octagon_translation<FT>& other) const;

        /// Returns the inverse of the current translation.
        Hyperbolic_octagon_translation inverse() const;


        /// \name Access Functions

        /// Returns the coefficient \f$\alpha\f$ of the translation. The first element of the returned `pair`
        /// contains the real part of \f$\alpha\f$. The second element contains its imaginary part.
        std::pair<FT,FT> alpha() const;

        /// Returns the coefficient \f$\beta\f$ of the translation. The first element of the returned `pair`
        /// contains the real part of \f$\beta\f$. The second element contains its imaginary part.
        std::pair<FT,FT> beta() const;

        /// Returns `true` if the current translation represents the identity element of the group \f$\mathcal G\f$.
        bool is_identity();


        /// \name Static Access Functions

        /*!
                Return the generator `wl` of the group \f$\mathcal G\f$.
                Note that `wl` can be an element of the enumeration set Generator. The calls
                `generator(0)` and `generator(A)` will both return the translation \f$a\f$.
        */
        static Self generator(const Word_letter wl) {
                return Self(wl);
        }

        /*!
                Returns the set of generators of \f$\mathcal G\f$ and their inverses in a `std::vector`.
                The generators are given in the order:
                \f$ [a, \overline{b}, c, \overline{d}, \overline{a}, b, \overline{c}, d].\f$
         */
        static void generators(std::vector< Hyperbolic_octagon_translation<FT> >& gens);


        /// \name Utility Functions
        /// @{

                /*!
                        Returns a string representation of the translation, containing its `Word_letter`s.
                        This function is given as utility for printing/debugging purposes. For example,
                        for the translation \f$abcd\f$, this function returns `0527`, and for the identity
                        it returns `_`.
                */
                std::string to_string() const;
        /// @}
};

} // namespace CGAL
