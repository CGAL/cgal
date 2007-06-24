#ifndef CGAL_ARR_DO_INTERSECT_ZONE_VISITOR_H
#define CGAL_ARR_DO_INTERSECT_ZONE_VISITOR_H

/*! \file
 * Definition of the Arr_do_intersect_zone_visitor_2 class.
 */

CGAL_BEGIN_NAMESPACE

/*! \class
 * A visitor class for Arrangement_zone_2, which check whether
 * a given x-monotone curve intersects the arrangment.
 * The class shouldbe templated by an Arrangement_2 class.
 */
template <class Arrangement_>
class Arr_do_intersect_zone_visitor
{
public:

  typedef Arrangement_                                Arrangement_2;
  typedef typename Arrangement_2::Traits_2            Traits_2;

  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;

  typedef typename Arrangement_2::Point_2             Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef std::pair<Halfedge_handle, bool>            Result;

private:

  const Halfedge_handle         invalid_he;    // Invalid halfedge.
  const Vertex_handle           invalid_v;     // Invalid vertex.

  bool                          m_intersect;   // Boolean to hold the answer.

public:

  /*! Constructor. */
  Arr_do_intersect_zone_visitor () :
    invalid_he (),
    invalid_v (),
    m_intersect (false)
  {}

  /*! Initialize the visitor with an arrangement object. */
  void init (Arrangement_2 *)
  {
    m_intersect = false;
  }

  /*!
   * Handle the a subcurve located in the interior of a given face.
   * \param cv The subcurve.
   * \param face The face containing cv's interior.
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param left_he The halfedge that contains the left endpoint of cv
   *               (or an invalid handle if no such halfedge exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_he The halfedge that contains the right endpoint of cv
   *                 (or an invalid handle if no such halfedge exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         subcurve into the arrangement.
   */
  Result found_subcurve (const X_monotone_curve_2& cv,
                         Face_handle face,
                         Vertex_handle left_v, Halfedge_handle left_he,
                         Vertex_handle right_v, Halfedge_handle right_he)
  { 
    if ((left_v == invalid_v) && (right_v == invalid_v) &&
        (left_he == invalid_he) && (right_he == invalid_he))
    {
      // The current subcurve just lies inside the given face, and its
      // endpoints are not incident to any valid vertex or edge, so it does
      // not intersect the arrangement.
      return (Result (invalid_he, false));
    }

    // We found an intersection. Note we return a result indicating that the
    // zone-computation can stop here.
    m_intersect = true;
    return (Result (invalid_he, true));
  }

  /*!
   * Handle the a subcurve that overlaps a given edge.
   * \param cv The overlapping subcurve.
   * \param he The overlapped halfedge (directed from left to right).
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         overlapping subcurve into the arrangement.
   */
  Result found_overlap (const X_monotone_curve_2& cv,
                        Halfedge_handle he,
                        Vertex_handle left_v, Vertex_handle right_v)
  {
    // We found an overlap (hence an intersection). Note we return a result
    // indicating that the zone-computation can stop here.
    m_intersect = true;
    return (Result (invalid_he, true));
  }

  bool do_intersect () const
  {
    return (m_intersect);
  }
};


//-----------------------------------------------------------------------------
// Checks if the given x-monotone curve intersects the existing arrangement.
//
template <class Traits, class Dcel, class PointLocation>
bool do_intersect_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr, 
                                    const typename Traits::X_monotone_curve_2& c,
                                    const PointLocation& pl)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  // Define a zone-computation object an a visitor that performs the
  // intersection check.
  typedef Arr_do_intersect_zone_visitor<Arrangement_2>  Zone_visitor;
  
  Zone_visitor                                     visitor;
  Arrangement_zone_2<Arrangement_2, Zone_visitor>  arr_zone (arr, &visitor);

  arr_zone.init (c, pl);
  arr_zone.compute_zone();

  return (visitor.do_intersect());
}

//-----------------------------------------------------------------------------
// Checks if the given x-monotone curve intersects the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel>
bool do_intersect_x_monotone_curve (Arrangement_2<Traits,Dcel>& arr, 
                                    const typename Traits::X_monotone_curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  //insert the curve using the walk point location
  do_intersect_x_monotone_curve (arr, c, walk_pl);
  return;
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
//
template <class Traits, class Dcel, class PointLocation>
bool do_intersect_curve (Arrangement_2<Traits,Dcel>& arr, 
                         const typename Traits::Curve_2& c,
                         const PointLocation& pl)
{
  // Obtain an arrangement accessor.
  typedef Arrangement_2<Traits,Dcel>                     Arrangement_2;

  // Break the input curve into x-monotone subcurves and isolated points.
  typedef Arr_traits_adaptor_2<Traits>                   Traits_adaptor_2;

  Traits_adaptor_2   *traits =
                        static_cast<Traits_adaptor_2*> (arr.get_traits());

  std::list<CGAL::Object>                     x_objects;
  std::list<CGAL::Object>::const_iterator     obj_iter;
  const typename Traits::X_monotone_curve_2  *x_curve;
  const typename Traits::Point_2             *iso_p;

  traits->make_x_monotone_2_object() (c,
                                      std::back_inserter (x_objects));

  // Insert each x-monotone curve into the arrangement.
  for (obj_iter = x_objects.begin(); obj_iter != x_objects.end(); ++obj_iter)
  {
    // Act according to the type of the current object.
    x_curve = object_cast<typename Traits::X_monotone_curve_2> (&(*obj_iter));
    if (x_curve != NULL)
    {
      // Check if the x-monotone subcurve intersects the arrangement.
      if (do_intersect_x_monotone_curve(arr, *x_curve, pl) == true)
        return true;
    }
    else
    {
      iso_p = object_cast<typename Traits::Point_2> (&(*obj_iter));
      CGAL_assertion (iso_p != NULL);

      // Check whether the isolated point lies inside a face (otherwise,
      // it conincides with a vertex or an edge).
      CGAL::Object  obj = pl.locate (*iso_p);

      return (object_cast<typename Arrangement_2::Face_const_handle>(&obj) !=
              NULL);
    }
  }

  // If we reached here, the curve does not intersect the arrangement.
  return (false);
}

//-----------------------------------------------------------------------------
// Checks if the given curve intersects the existing arrangement.
// Overloaded version with no point location object - the walk point-location
// strategy is used as default.
//
template <class Traits, class Dcel>
bool do_intersect_curve (Arrangement_2<Traits, Dcel>& arr, 
                         const typename Traits::Curve_2& c)
{
  typedef Arrangement_2<Traits, Dcel>                          Arrangement_2;
  typedef Arr_walk_along_line_point_location<Arrangement_2>    Walk_pl;
  
  // create walk point location object
  Walk_pl    walk_pl(arr);

  // check if the curve intersects the arrangement using the walk point 
  // location.
  return do_intersect_curve (arr, c, walk_pl);
}

CGAL_END_NAMESPACE

#endif
