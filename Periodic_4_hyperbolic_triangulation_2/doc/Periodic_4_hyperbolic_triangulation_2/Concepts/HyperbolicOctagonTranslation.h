// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept
\cgalModifBegin
The concept `HyperbolicOctagonTranslation` describes the requirements for an object that represents 
a hyperbolic translation in the fundamental group \f$\mathcal G\f$ of the the Bolza surface. 

Generally, a hyperbolic translation \f$g\f$ is an orientation-preserving isometry acting 
on the hyperbolic plane, and is specified by two complex numbers \f$\alpha\f$ and \f$\beta\f$:
\f[
	g(z) = \frac{\alpha \cdot z + \beta}{\overline{\beta} \cdot z + \overline{\alpha}},
	\qquad |\alpha|^2 - |\beta|^2 = 1,
\f]
where \f$\overline{a}\f$ and \f$\overline{b}\f$ denote the complex conjugates of \f$a\f$
and \f$b\f$ respectively. A hyperbolic translation is also a word on the generators of the 
fundamental group of a hyperbolic surface.

Specifically for the group \f$\mathcal G\f$, let \f$[g_0, g_1, \ldots, g_7]\f$ be the set of its 
generators, ordered counter-clockwise around the origin, with \f$g_k\f$ translating the origin 
\f$O\f$ along the real axis in the positive direction. The coefficients \f$\alpha_k\f$ and \f$\beta_k\f$ 
of each generator of this ordered set are given as 
\f[
	\alpha_k = 1 + \sqrt{2}, \qquad \beta_k = \exp\left(\tfrac{ik\pi}{4}\right)\sqrt{2}\sqrt{1+\sqrt{2}}.
\f]
Note that with the notation given in Section \ref P4HT2_thespace of the User manual, the ordered 
set of generators can be written also as
\f[
	[g_k, k = 0, \ldots, 7] = [ a, \overline{b}, c, \overline{d}, \overline{a}, b, \overline{c}, d ].
\f] 
\cgalModifEnd

\cgalHasModel CGAL::Hyperbolic_octagon_translation
*/


class HyperbolicOctagonTranslation {

public:
	
	
	/// \name Types
	
	/// Number type of the coefficients of \f$\alpha\f$ and \f$\beta\f$.
	/// This number type must be compatible with the field type `FT` of the concept 
	/// `Periodic_4HyperbolicDelaunayTriangulationTraits_2`.
	typedef unspecified_type 			NT;

	/// \cgalModifBegin Represents a single letter of the word corresponding to the translation. \cgalModifEnd
	typedef unsigned short int 			Word_letter;

	/// Word representation of the translation, a sequence of elements of type `Word_letter`.
	typedef std::vector<Word_letter> 	Word;

	/// \cgalModifBegin Represents a complex number with coefficients of type `NT`.\cgalModifEnd
	typedef std::pair<NT,NT> 			Complex;

	/// \cgalModifBegin Represents the coefficients \f$\alpha\f$ and \f$\beta\f$ of a hyperbolic translation.\cgalModifEnd
	typedef std::pair<Complex,Complex> 	Coefficients;

	/// \name Creation

	/// Default constructor. Creates the identity translation.
	HyperbolicOctagonTranslation();

	/// Creates the translation described by the word `w`.
	HyperbolicOctagonTranslation(Word w);

	/// Creates the translation described by the one-letter word `l`.
	HyperbolicOctagonTranslation(Word_letter l);

	/// Creates the translation described by the two-letter word `lm`.
	HyperbolicOctagonTranslation(Word_letter l, Word_letter m);

	/// Creates the translation described by the three-letter word `lmn`.
	HyperbolicOctagonTranslation(Word_letter l, Word_letter m, Word_letter n);

	/// Creates the translation described by the four-letter word `lmno`.
	HyperbolicOctagonTranslation(Word_letter l, Word_letter m, Word_letter n, Word_letter o);


	/// \name Operations
	
	/// Multiplication operator; composes the current translation with `other`. 
	HyperbolicOctagonTranslation operator*(const HyperbolicOctagonTranslation& other) const;

	///	Difference operator; the difference of two translations \f$v\f$ and \f$w\f$ is defined as \f$v * w^{-1}\f$. 
  	HyperbolicOctagonTranslation operator-(const HyperbolicOctagonTranslation& other) const;

  	/// Assignment operator; modifying the translation after the assignment leaves `other` unaffected.
	HyperbolicOctagonTranslation& operator=(const HyperbolicOctagonTranslation& other);

	/// Equality operator.
	bool operator==(const Hyperbolic_octagon_translation<NT>& other) const;

	/// Inequality operator.
	bool operator!=(const Hyperbolic_octagon_translation<NT>& other) const;

	/// Comparison operator. 
	/// The comparison is done on the ordering of the translations around the original domain.
	bool operator<(const Hyperbolic_octagon_translation<NT>& other) const;

	/// Returns the inverse of the current translation.
	HyperbolicOctagonTranslation inverse() const;


	/// \name Access Functions

	/// Returns the coefficient \f$\alpha\f$ of the translation's matrix. The first element of the returned `pair` 
	/// contains the real part of \f$\alpha\f$. The second element contains its imaginary part.
	Complex alpha() const;

	/// Returns the coefficient \f$\beta\f$ of the translation's matrix. The first element of the returned `pair` 
	/// contains the real part of \f$\beta\f$. The second element contains its imaginary part.
	Complex beta() const;

	/// Returns `true` if the current translation represents the identity element of the group.
	bool is_identity();


	/// \name Static Access Functions

	/*!
		\cgalModifBegin
		Returns a vector containing the coefficients of the generators of \f$\mathcal G\f$ and of their inverses.
		The translations must be given in the order:
		\f[ [g_0, g_1, \ldots, g_7]. \f]
		\cgalModifEnd
 	*/
	static void get_generator_coefficients(std::vector< Coefficients >& gens);


	/*!
		\cgalModifBegin
		Returns a vector containing the generators of \f$\mathcal G\f$ and their inverses.
		The translations must be given in the order:
		\f[ [g_0, g_1, \ldots, g_7]. \f]
		\cgalModifEnd
 	*/
	static void get_generators(std::vector<HyperbolicOctagonTranslation>& gens);

}; // end of class HyperbolicOctagonTranslation

