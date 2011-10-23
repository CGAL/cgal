// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <iostream>
#include <fstream>

typedef CGAL::Combinatorial_map_with_points<3>     LCC_3;
typedef LCC_3::Dart_handle                         Dart_handle;
typedef LCC_3::Vertex_attribute                    Vertex;
typedef LCC_3::Point                               Point;
typedef LCC_3::Vector                              Vector;
typedef LCC_3::FT                                  FT;

#define PI 3.1415926535897932

// Smoth a vertex depending on the vertices of its incident facet.
class Smooth_old_vertex
{
public:
  /**  Contructor.
   * @param amap is the map to smooth
   * @param amark is a mark designing old darts (i.e. darts not created during
   *        the triangulation step)
   */
  Smooth_old_vertex(LCC_3 &amap, unsigned int amark) : mmap(amap)
  {}
  
  Vertex operator()( Vertex& v) const
  {
    Dart_handle d = v.dart();
    CGAL_assertion(d!=NULL);

    int degree=0;
    bool open = false;
    
    LCC_3::One_dart_per_incident_cell_range<1,0>::iterator 
      it=mmap.one_dart_per_incident_cell<1,0>(d).begin(),
      itend(mmap.one_dart_per_incident_cell<1,0>(d).end());
    for ( ; it!=itend; ++it )
      {
	++degree;
	if ( it->is_free(2) ) open = true;
      }        

    if ( open ) return v;

    LCC_3::FT alpha = (4.0f - 2.0f *
                     (LCC_3::FT)cos( 2.0f * PI / (LCC_3::FT)degree)) / 9.0f;
    LCC_3::Vector vec = (v.point() - CGAL::ORIGIN) * ( 1.0f - alpha);

    for (it.rewind(); it!=mmap.one_dart_per_incident_cell<1,0>(d).end(); ++it)
      {
	CGAL_assertion(!it->is_free(2));

        vec = vec + (mmap.point(it->other_extremity()) - CGAL::ORIGIN)
	  * alpha / degree;
      }
    
    Vertex res(CGAL::ORIGIN + vec);
    res.set_dart(d);
    return res;
  }
private:
  LCC_3& mmap;
};

// Flip an edge, work in 2D and in 3D.
Dart_handle flip_edge(LCC_3 &m, Dart_handle d)
{
  CGAL_assertion( d!=NULL && !d->is_free(2) );

  if ( !CGAL::is_removable<LCC_3,1>(m,d) ) return NULL;

  Dart_handle d2 = d->beta(1)->beta(1);
  CGAL::remove_cell<LCC_3,1>(m, d);
  CGAL::insert_cell_1_in_cell_2(m, d2, d2->beta(1)->beta(1));

  return d2->beta(0);
}

// Subdivide each facet of the map by using sqrt(3)-subdivision.
void subdivide_map_3(LCC_3& m)
{
  if (m.is_empty() )
    return;

  unsigned int mark    = m.get_new_mark();
  unsigned int treated = m.get_new_mark();
  m.negate_mark(mark); // All the old darts are marked in O(1).

  // 1) We smoth the old vertices.
  std::vector<Vertex> vertices; // smooth the old vertices
  vertices.reserve(m.number_of_vertex_attributes()); // get intermediate space
  std::transform( m.vertex_attributes().begin(), m.vertex_attributes().end(),
                  std::back_inserter(vertices), Smooth_old_vertex(m,mark) );

  // 2) We subdivide each facet.
  m.negate_mark(treated); // All the darts are marked in O(1).
  unsigned int nb=0;
  for ( LCC_3::Dart_range::iterator it=m.darts().begin(); 
	m.number_of_marked_darts(treated)>0; ++it )
    {
      ++nb;
      if ( m.is_marked(it, treated) )
	{
	  // We unmark the darts of the facet to process only once dart/facet.
	  CGAL::unmark_cell<LCC_3,2>(m, it, treated);
	  // We triangulate the facet.
    m.insert_barycenter_in_cell<2>(it);
	}
    }

  CGAL_assertion( m.is_whole_map_unmarked(treated) );
  CGAL_assertion( m.is_valid() );
  m.free_mark(treated);

  // 3) We update the coordinates of old vertices.
  for(std::vector<Vertex>::iterator vit=vertices.begin();
      vit!=vertices.end();++vit)
    {
      CGAL_assertion(vit->dart()!=NULL);
      m.point(vit->dart())=vit->point();
    }

  // 4) We flip all the old edges.
  m.negate_mark(mark); // Now only new darts are marked.
  Dart_handle d2 = NULL;
  for (LCC_3::Dart_range::iterator it=m.darts().begin(); it != m.darts().end(); )
    {
      d2 = it++;
      CGAL_assertion(d2!=NULL);
      if (!m.is_marked(d2, mark))   // This is an old dart.
	{
	  // We process only the last dart of a same edge.
	  if (!d2->is_free(2) && (d2->beta(2)->beta(3)==d2->beta(3)->beta(2)))
	    {
	      if ( m.is_marked(d2->beta(2), mark) &&
		   (d2->is_free(3) ||
		    (m.is_marked(d2->beta(3), mark) &&
		     m.is_marked(d2->beta(2)->beta(3), mark))) )
		{
		  m.negate_mark(mark); // thus new darts will be marked
		  flip_edge(m, d2);
		  m.negate_mark(mark);
		}
	      else
		m.mark(d2, mark);
	    }
	  else
	    m.mark(d2, mark);
	}
    }
  CGAL_assertion( m.is_whole_map_marked(mark) );
  m.free_mark(mark);

  CGAL_postcondition(m.is_valid());
}

Dart_handle make_iso_cuboid(LCC_3& lcc, const Point& basepoint, FT lg)
{
	return lcc.make_hexahedron(basepoint,
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(lg,0,0)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(lg,lg,0)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(0,lg,0)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(0,lg,lg)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(0,0,lg)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(lg,0,lg)),
                             LCC_3::Construct_translated_point()(basepoint,
                                                                 LCC_3::Vector(lg,lg,lg)));
}

int main(int narg, char** argv)
{
  if ( narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
    {
      std::cout<<"Usage : a.out [-h -?] [filename1 ... filenamek]"<<std::endl;
      std::cout<<"  Without parameter, this program creates two cubes, 3-sew them "
	"and subdivides them by using sqrt(3)-subdivision."<<std::endl;
      std::cout<<"  With parameters, load all the filename into a same map "
	"each file must be an off file) and subdivides the map by using "
	"sqrt(3)-subdivision."<<std::endl;
      exit(EXIT_SUCCESS);
    }

  LCC_3 cm;

  if ( narg==1 )
    {       
      // Create 2 cubes.
      Dart_handle d1 = make_iso_cuboid(cm, Point(-2, 0, 0), 1);
      Dart_handle d2 = make_iso_cuboid(cm, Point(0, 0, 0), 1);
       
      // 3-Sew the 2 cubes along one facet
      cm.sew<3>(d1->beta(1)->beta(1)->beta(2), d2->beta(2));
    }
  else
    {
      for (int i=1; i<narg; ++i)
	{
	  std::ifstream is(argv[i]);
	  CGAL::import_from_polyhedron_flux<LCC_3>(cm,is);
	  is.close();
	}
    }

  std::cout << "****** Initial object ******" << std::endl;
  // Display the vertices.
  CGAL::set_ascii_mode(std::cout);
  std::cout << "Vertices: ";
  for (LCC_3::Vertex_attribute_range::iterator
	 v(cm.vertex_attributes().begin()), vend(cm.vertex_attributes().end());
       v!=vend; ++v)
    std::cout << v->point() << "; ";
  std::cout << std::endl;

  // Display the m characteristics.
  cm.display_characteristics(std::cout)
    << ", valid=" << cm.is_valid() << std::endl;
      
  // Subdivide the volumes.
  subdivide_map_3(cm);

  std::cout << "****** Subdivided object ******" << std::endl;
  // Display the vertices.
  CGAL::set_ascii_mode(std::cout);
  std::cout << "Vertices: ";
  for (LCC_3::Vertex_attribute_range::iterator
	 v(cm.vertex_attributes().begin()), vend(cm.vertex_attributes().end());
       v!=vend; ++v)
    std::cout << v->point() << "; ";
  std::cout << std::endl;

  // Display the m characteristics.
  cm.display_characteristics(std::cout)
    << ", valid=" << cm.is_valid() << std::endl;

  return 1;
}

