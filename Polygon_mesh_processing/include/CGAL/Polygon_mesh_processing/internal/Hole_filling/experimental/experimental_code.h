// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ilker O. Yaz

/************************************************************************
 * Currently not useful code pieces, in case there will be any need in the future.
 * Also they might not work when they plugged-in due to recent changes.
 ************************************************************************/

// It can produce a patch from both complete and incomplete lambda
// WARNING: Not working good for all cases
// For holes, this code first close them then erase them.
// However the algorithm might produce holes which are invalid to close (closing them breaks edge manifoldness, so erasing doesn't work)
template<class Polyhedron, class OutputIteratorPatch, class OutputIteratorHole>
struct Tracer_polyhedron_incomplete
{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle    Facet_handle;

  Tracer_polyhedron_incomplete(OutputIteratorPatch out,
                               OutputIteratorHole out_hole,
                               Polyhedron& polyhedron,
                               std::vector<Halfedge_handle>& P)
    : out(out), out_hole(out_hole), polyhedron(polyhedron), P(P)
  { }

  template <class LookupTable>
  void
  operator()(const LookupTable& lambda, int i, int k)
  {
    std::vector<Facet_handle> facets_to_delete;
    (*this)(lambda, i, k, facets_to_delete, true);
    for(typename std::vector<Facet_handle>::iterator it = facets_to_delete.begin();
        it != facets_to_delete.end(); ++it)
    {
      *out_hole++=(*it)->halfedge(); // each deleted facet corresponds to a new hole
      polyhedron.erase_facet((*it)->halfedge());
    }
  }

private:
  template <class LookupTable>
  Halfedge_handle
  operator()(const LookupTable& lambda,
             int i, int k,
             std::vector<Facet_handle>& facets_to_delete,
             bool last)
  {
    if(i + 1 == k) { return P[i+1]; }

    Halfedge_handle h, g;
    if(i+2 == k){

      if(last)
      { h = polyhedron.fill_hole(P[i+1]); }
      else
      { h = polyhedron.add_facet_to_border(P[i+1]->prev(), P[i+2/*k*/]); }

      CGAL_assertion(h->facet() != Facet_handle());

      int la = lambda.get(i, k);
      if(la == -1) {
        facets_to_delete.push_back(h->facet());
      }
      else {
        *out++ = h->facet();
      }
      return h->opposite();
    }
    else
    {
      int la = lambda.get(i, k);
      if(la == -1) {
        if(last)
        { h = polyhedron.fill_hole(P[i+1]); }
        else
        { h = polyhedron.add_facet_to_border(P[i+1]->prev(), P[i+2/*k*/]); }
        facets_to_delete.push_back(h->facet());
        return h->opposite();
      }
      else {
        h = operator()(lambda, i, la, facets_to_delete, false);
        g = operator()(lambda, la, k, facets_to_delete, false);

        if(last)
        { h = polyhedron.fill_hole(g); }
        else
        { h = polyhedron.add_facet_to_border(h->prev(), g); }

        CGAL_assertion(h->facet() != Facet_handle());
        *out++ = h->facet();
        return h->opposite();
      }
    }
  }

public:
  OutputIteratorPatch out;
  OutputIteratorHole out_hole;
  Polyhedron& polyhedron;
  std::vector<Halfedge_handle>& P;
};

// Try closing holes by gathering incomplete patches together (an external approach)
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline_incomplete(InputIterator pbegin, InputIterator pend,
                                     InputIterator qbegin, InputIterator qend,
                                     OutputIterator out)
{

  typedef typename std::iterator_traits<InputIterator>::value_type Point_3;
  typedef Weight_incomplete<Weight_min_max_dihedral_and_area>      Weight;
  typedef Weight_calculator<Weight, Is_valid_degenerate_triangle>  WC;

  typedef std::vector<boost::tuple<int, int, int> >                   Facet_vector; /* deliberately not OutputIteratorValueType*/
  typedef std::back_insert_iterator<Facet_vector>                     OutIt;
  typedef Tracer_polyline_incomplete<Facet_vector::value_type, OutIt> Tracer;
  typedef std::pair<int, int> Range;

  std::vector<Point_3> P(pbegin, pend);
  std::vector<Point_3> Q(qbegin, qend);

  if(P.front() != P.back()){
    P.push_back(P.front());
    if( !Q.empty() && P.size() > Q.size()) {
      Q.push_back(Q.front());
    }
  }

  std::vector<OutputIteratorValueType>   patch_facets;
  std::stack<Range>                      remaining_holes;

  remaining_holes.push(Range(0, P.size() -2));

  while(!remaining_holes.empty()) {
    Range h = remaining_holes.top();
    remaining_holes.pop();

    std::vector<Point_3> P_r(&P[h.first], (&P[h.second]) + 1);
    std::vector<Point_3> Q_r;
    if(!Q.empty()) { Q_r.insert(Q_r.begin(), &Q[h.first], (&Q[h.second]) + 1); };

    Facet_vector new_facets;
    OutIt new_facets_out(new_facets);
    std::vector<Range> new_holes;
    std::back_insert_iterator<std::vector<Range> > new_holes_out(new_holes);
    Tracer tracer(new_facets_out, new_holes_out);

    triangulate_hole_polyline(P_r, Q_r, tracer, WC(), true, true);

    if(new_facets.empty()) {
      new_holes.clear();
      //triangulate_hole_polyline(P_r, Q_r, tracer, WC(), false, true);
      if(new_facets.empty()) {
        // if no patch facets created and also we are using brute force approach, then there is nothing to do,
        // leave `out` intact and return
        CGAL_warning_msg(false, "Returning no output. Filling hole with incomplete patches is not successful!");
        return out;
      }
    }
    // put new borders to remaining_holes
    for(typename std::vector<Range>::iterator it = new_holes.begin(); it != new_holes.end(); ++it) {
      remaining_holes.push(std::make_pair(it->first + h.first, it->second + h.first));
    }
    tracer.remaining_holes.clear();
    for(Facet_vector::iterator it = new_facets.begin(); it != new_facets.end(); ++it) {
      patch_facets.push_back(
        OutputIteratorValueType(it->get<0>() + h.first, it->get<1>() + h.first, it->get<2>() + h.first));
    }
  }

  return std::copy(patch_facets.begin(), patch_facets.end(), out);
}

// (for Polyhedron_3) Try closing holes by gathering incomplete patches together (an external approach)
template<class Polyhedron, class OutputIterator>
std::pair<OutputIterator, Weight_min_max_dihedral_and_area>
triangulate_hole_Polyhedron_incomplete(Polyhedron& polyhedron,
                                       typename Polyhedron::Halfedge_handle border_halfedge,
                                       OutputIterator out)
{
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle    Facet_handle;

  Weight_min_max_dihedral_and_area weight_total = Weight_min_max_dihedral_and_area::DEFAULT();
  std::vector<Facet_handle>   patch_facets;
  std::stack<Halfedge_handle> remaining_holes;
  remaining_holes.push(border_halfedge);

  while(!remaining_holes.empty()) {
    Halfedge_handle h = remaining_holes.top();
    remaining_holes.pop();
    std::vector<Halfedge_handle> holes_h;

    std::size_t patch_facets_before = patch_facets.size();
    Weight_min_max_dihedral_and_area w =
      triangulate_hole_Polyhedron(polyhedron, h, back_inserter(patch_facets), back_inserter(holes_h), true);

    if(patch_facets_before == patch_facets.size()) {
      holes_h.clear();
      patch_facets_before = patch_facets.size();
      w = triangulate_hole_Polyhedron(polyhedron, h, back_inserter(patch_facets), back_inserter(holes_h), false);
      if(patch_facets_before == patch_facets.size()) {
        // if no patch facets created and also we are using brute force approach, then there is nothing to do,
        // leave `out` intact and return
        CGAL_warning_msg(false, "Returning no output. Filling hole with incomplete patches is not successful!");
        return std::make_pair(out, Weight_min_max_dihedral_and_area::NOT_VALID());
      }
    }
    // put new borders to remaining_holes and update weight
    for(typename std::vector<Halfedge_handle>::iterator it = holes_h.begin(); it != holes_h.end(); ++it) {
      //remaining_holes.push(*it);
    }
    weight_total = weight_total + w;
  }

  out = std::copy(patch_facets.begin(), patch_facets.end(), out);
  return std::make_pair(out, weight_total);
}
