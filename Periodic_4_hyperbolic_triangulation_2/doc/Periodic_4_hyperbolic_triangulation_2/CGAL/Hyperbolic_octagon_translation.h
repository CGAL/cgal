// Copyright (c) 1999-2018   INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2HyperbolicTranslationsClasses

The class `Hyperbolic_octagon_translation` is the default model for the 
concept `HyperbolicOctagonTranslation`. It accepts one template parameter:

\tparam NT 	Number Type. Must provide exact computations with algebraic numbers, 
			notably with nested square roots. The default value for this parameter 
			is `CORE::Expr`.

\cgalModels HyperbolicOctagonTranslation
*/

template <typename NT>
class Hyperbolic_octagon_translation {

public:

	/// \name Functions
	/// @{

		/*!
			Returns a string representation of the word that the translation represents.
		*/
		std::string to_string() const;
	/// @}
};

} // namespace CGAL
