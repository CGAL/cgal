// Copyright (c) 2016-2017 INRIA Nancy - Grand Est (France).
// All rights reserved.


namespace CGAL {


/*!
\ingroup PkgPeriodic4HyperbolicTriangulation2MainClasses

The class `Periodic_4_hyperbolic_Delaunay_triangulation_2` allows the construction and
handling of Delaunay triangulations of the Bolza surface. 

The class expects two template parameters.

\tparam GT 	Geometric Traits parameter. Must be a model of the concept 
			`Periodic_4HyperbolicDelaunayTriangulationTraits_2`. This parameter has no
			default value.
\tparam	TDS %Triangulation Data Structure parameter. Must be a model of the concept
		`TriangulationDataStructure_2` with some additional functionality in the vertices 
		and faces, following the concepts `Periodic_4HyperbolicTriangulationVertexBase_2` and
		`Periodic_4HyperbolicTriangulationFaceBase_2`, respectively. The default value for 
		this parameter is
		`Triangulation_data_structure_2< Periodic_4_hyperbolic_triangulation_vertex_base_2<GT>, Periodic_4_hyperbolic_triangulation_face_base_2<GT> >`

*/

  template <  class GT, class TDS >
  class Periodic_4_hyperbolic_Delaunay_triangulation_2: public Periodic_4_hyperbolic_triangulation_2<GT, TDS> {

  public:

  	/// \name Creation
  	/// @{
		/*!
			Default constructor with an optional parameter for the geometric traits object.
		*/	
		Periodic_4_hyperbolic_Delaunay_triangulation_2(const Geom_traits &gt = Geom_traits());

		/*!
			Copy constructor.
		*/
		Periodic_4_hyperbolic_Delaunay_triangulation_2(const Periodic_4_hyperbolic_Delaunay_triangulation_2& tr);

		/*!
		\cgalModifBegin
			Point range constructor.
			Initializes a triangulation and then inserts all points in the `InputIterator`.
		\cgalModifEnd
		*/
		template < class InputIterator >
		Periodic_4_hyperbolic_Delaunay_triangulation_2( InputIterator first, InputIterator last, const Geom_traits& gt = Geom_traits() );
	/// @}


	/// \name Clearing the triangulation
	/// @{

		/*!
			Deletes all faces and vertices of the triangulation, and re-initializes 
			the triangulation data structure with the set of dummy points.
		*/
		void clear();
	/// @}


	/// \name Queries
	/// @{
	 	
	 	/*!
	 		\cgalModifBegin
			Computes the conflict zone induced by `p`.
			The output iterator `it` contains all faces in conflict with `p`.
			The optional parameters `start` and `ltr`, if given, must be such 
			that `start` translated by `ltr` is in conflict with `p`. 
			\cgalModifEnd
	 	*/
	 	template<class OutputFaceIterator>
		void
		find_conflicts(	const Point& p,
						OutputFaceIterator it,
						Face_handle start = Face_handle(),
						Hyperbolic_translation ltr = Hyperbolic_translation()) const;
 

		/*!
			Checks if the vertex `vh` is part of the set of vertices storing dummy points.
		*/
		bool is_dummy_vertex(Vertex_handle vh) const;

		/*!
			Returns the number of dummy points currently existing in the triangulation.
		*/
		int number_of_dummy_points() const;

		/*!
		\cgalModifBegin
			Checks the combinatorial validity of the triangulation, and also the validity of its geometric embedding. 
		\cgalModifEnd
		*/
		bool is_valid(bool verbose = false) const;
	/// @}


	/// \name Point insertion
	/// @{

		/*!
			\cgalModifBegin
			Inserts the point `p` in the triangulation.
			The face `start`, if given, is used as a starting place for the location of the point.
			Note that this function does not remove unnecessary dummy points. 
			The removal, if desired, should be done by manually calling the function 
			`clean_dummy_points()`.
			\cgalModifEnd
		*/
		Vertex_handle insert(const Point  &p, Face_handle start = Face_handle());

		/*!
			\cgalModifBegin
			Inserts all points in the input iterator into the triangulation.
			Note that this function by default tries to remove unnecessary dummy points 
			at the end of the insertion process. This behavior is controlled by the 
			optional boolean parameter `flag_clean_dummy_points`; if automatic removal
			of the dummy points is not desired, set the flag to `false`.
			\cgalModifEnd
		*/
		template < class InputIterator >
		std::ptrdiff_t
		insert(InputIterator first, InputIterator last, const bool flag_clean_dummy_points = true);
	/// @}

	/// \name Vertex removal
	/// @{

		/*!
			\cgalModifBegin
			Removes the vertex `v` from the triangulation.
			Note that `v` is removed from the triangulation only if the resulting
			triangulation remains a simplicial complex. The function returns `true`
			if the vertex `v` is removed from the triangulation, otherwise it returns 
			`false`.
			\cgalModifEnd
		*/
		bool remove(Vertex_handle v);
	/// @}

	/// \name Removal of dummy points
	/// @{

		/*!
			\cgalModifBegin
			Tries to remove the dummy points from the triangulation one by one.
			Returns the number of dummy points that could not be removed and are 
			still present in the triangulation.
			\cgalModifEnd
		*/
		int clean_dummy_points();
	/// @}
	

	/// \name Voronoi diagram
	/// @{
		
		/*!
			Returns the dual of face `f`, i.e., the hyperbolic center of the disk 
			defined by the vertices of `f`. If the optional parameter `nbtr` is given, 
			then this function returns the dual of the  periodic triangle `(f, nbtr)`.
		*/
		Voronoi_point dual (Face_handle f, Hyperbolic_translation nbtr = Hyperbolic_translation()) const;

		/*!
			Returns the hyperbolic segment that is dual to the edge `e`.
		*/
		Segment dual(const Edge &e) const;
	/// @}


};  // class Periodic_4_hyperbolic_Delaunay_triangulation_2



} // namespace CGAL





