// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.

/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2Concepts
\cgalConcept

\cgalModifBegin
\cgalRefines TriangulationDSFaceBase_2

A refinement of the concept `TriangulationDSFaceBase_2` that adds an interface for hyperbolic translations.
\cgalModifEnd

At the base level, a face stores handles to its incident vertices and to its neighboring faces. 
Compare with Section \ref Section_2D_Triangulations_Software_Design of the 2D Triangulations 
package. The vertices and neighbors are indexed counter-clockwise 0, 1, and 2. Neighbor `i` lies
opposite to vertex `i`. 

For periodic hyperbolic triangulations, the face base class needs to store three hyperbolic translations, 
one for each vertex. Applying each translation to the point stored in the corresponding vertex produces 
the canonical representative of the face in the hyperbolic plane.

\cgalHasModel CGAL::Periodic_4_hyperbolic_triangulation_ds_face_base_2

\sa `TriangulationDataStructure_2`
\sa `Periodic_4HyperbolicTriangulationDSVertexBase_2`

*/

class Periodic_4HyperbolicTriangulationDSFaceBase_2 
{
public:

	/// \name Types
	/// @{

		/*!
			This must be a model of the concept HyperbolicOctagonTranslation.
		*/
		typedef undefined_type 						Hyperbolic_translation;
	///@}


	/// \name Creation
	/// @{

		/*!
		Default constructor.
		*/
		Periodic_4HyperbolicTriangulationDSFaceBase_2();

		/*!
		Creates a face with vertices `v0, v1` and `v1`.
		*/
		Periodic_4HyperbolicTriangulationDSFaceBase_2(
			const Vertex_handle& v0, const Vertex_handle& v1,
			const Vertex_handle& v2); 

		/*!
		Creates a face with vertices `v0, v1` and `v1`, setting 
		also neighborhood relations with `n0, n1` and `n2`.
		*/
		Periodic_4HyperbolicTriangulationDSFaceBase_2(
			const Vertex_handle& v0, const Vertex_handle& v1,
			const Vertex_handle& v2, const Face_handle&   n0,
			const Face_handle&   n1, const Face_handle&   n2); 
	/// @}


	/// \name Access functions
	/// @{	

		/*!
			Returns the translation corresponding to the vertex `i`.
			\pre \f$ 0 \leq i \leq 2\f$.
		*/
		Hyperbolic_translation translation(int i) const;


		/*!
			Sets the `i`-th translation to `new_tr`.
			\pre \f$ 0 \leq i \leq 2\f$.
		*/
		void set_translation(const int& i, const Hyperbolic_translation& new_tr);


		/*!
			Changes the orientation of the face by exchanging `vertex(0)` with `vertex(1)`, 
			`neighbor(0)` with `neighbor(1)`, and `translation(0)` with `translation(1)`.

			\sa TriangulationDSFaceBase_2::reorient() 
		*/
		void reorient();

	/// @}

};