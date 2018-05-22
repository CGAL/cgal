// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

The concept `HyperbolicOctagonTranslation` describes the requirements for the class of 
hyperbolic translations specific to the Bolza surface. 

Generally, a hyperbolic translation \f$g\f$ is an orientation-preserving isometry acting 
on the hyperbolic plane, and is specified by two complex numbers \f$\alpha\f$ and \f$\beta\f$:
\f[
	g(z) = \frac{\alpha \cdot z + \beta}{\overline{\beta} \cdot z + \overline{\alpha}},
	\qquad |\alpha|^2 - |\beta|^2 = 1.
\f]
A hyperbolic translation is also a word on the generators of the fundamental group 
of a hyperbolic surface.

In our particular case, we consider an application specific to the Bolza surface. For our 
application, we only need to deal with a finite set of hyperbolic translations from the 
fundamental group of the Bolza surface. A notable property of these translations is that
their word length is at most four; the composition of two translation always gives a
translation that can be simplified so that its length is at most four. For more details, 
see \cgalCite{cgal:btv-dtosl-16} and \cgalCite{cgal:it-idtbs-17}.

\cgalHasModel CGAL::Hyperbolic_octagon_translation
*/


class HyperbolicOctagonTranslation {

public:
	
	
	/// \name Types
	
	/// Number type of the coefficients of \f$\alpha\f$ and \f$\beta\f$.
	typedef unspecified_type 			NT;

	/// Represents a single letter of the word corresponding to the translation.
	typedef unspecified_type 			Word_letter;

	/// Word representation of the translation, a sequence of elements of type `Word_letter`.
	typedef std::vector<Word_letter> 	Word;

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
	std::pair<NT,NT> alpha() const;

	/// Returns the coefficient \f$\beta\f$ of the translation's matrix. The first element of the returned `pair` 
	/// contains the real part of \f$\beta\f$. The second element contains its imaginary part.
	std::pair<NT,NT> beta() const;

	/// Returns `true` if the current translation represents the identity element of the group.
	bool is_identity();


}; // end of class HyperbolicOctagonTranslation

