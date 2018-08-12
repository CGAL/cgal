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
		and faces, following the models `Periodic_4HyperbolicTriangulationDSVertexBase_2` and
		`Periodic_4HyperbolicTriangulationDSFaceBase_2`, respectively. The default value for 
		this parameter is
		`Triangulation_data_structure_2< Triangulation_vertex_base_2<GT, Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT> >, Triangulation_face_base_2<GT, Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT> > >`

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
			Computes recursively the conflict zone induced by `p`.
			The output iterator `it` contains all faces in conflict with `p`.
	 	*/
	 	template<class OutputFaceIterator>
		void
		get_conflicts(	const Point& p,
						OutputFaceIterator it) const;
 

		/*!
			\cgalModifBegin
			Same as above. 
			This is the function called recursively by the overloaded version above.
			The parameter `cf` is the face to examine if in conflict with `p`.
			The hyperbolic translation `tr` is such that the periodic triangle
			`(cf, tr)` is adjacent to the triangle containing `p`. The set `visited` 
			contains all faces that have already been examined so far, to avoid 
			unnecessary controls. The faces in conflict with `p` are stored in the 
			output iterator `it`.
			\cgalModifEnd
		*/
	  	template<class OutputFaceIterator>
		void
		get_conflicts(	const Point& p,
						const Face_handle cf,
						Hyperbolic_translation tr,
						std::set<Face_handle>& visited,
						OutputFaceIterator it) const;

		/*!
			Checks if the vertex `vh` is part of the set of vertices storing dummy points.
		*/
		bool is_dummy_vertex(Vertex_handle vh) const;

		/*!
			Returns the number of dummy points currently existing in the triangulation.
		*/
		int number_of_dummy_points() const;

	/// @}


	/// \name Point insertion
	/// @{

		/*!
			\cgalModifBegin
			Inserts the point `p` in the triangulation.
			The face `start`, if given, is used as a starting place for the location of the point.
			If `batch_insertion` is set to `false` (the default value), then after the insertion 
			of `p` in the triangulation, all unnecessary dummy points in the triangulation (if any) 
			are removed. If `batch_insertion` is set to `false`, then this operation is omitted.
			\cgalModifEnd
		*/
		Vertex_handle insert(const Point  &p, Face_handle start = Face_handle(), bool batch_insertion = false);

		/*!
			\cgalModifBegin
			Inserts all points in the input iterator into the triangulation.
			At the end of the insertion, all unnecessary dummy points in the triangulation (if any)
			are removed.
			\cgalModifEnd
		*/
		template < class InputIterator >
		std::ptrdiff_t
		insert(InputIterator first, InputIterator last);
	/// @}


	/// \name Vertex removal
	/// @{

		/*!
			Removes the vertex `v` from the triangulation.
			The triangulation is fixed locally by triangulating the resulting hole. 
			Note that `v` is removed from the triangulation only if the resulting
			triangulation would be a simplicial complex. The function returns `true`
			if the vertex `v` is removed from the triangulation, otherwise it returns 
			`false`.
		*/
		bool remove(Vertex_handle v);
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





