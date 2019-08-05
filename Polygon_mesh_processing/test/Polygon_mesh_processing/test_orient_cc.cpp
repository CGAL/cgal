#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> SMesh;

template<class TriangleMesh, class NamedParameters>
bool test_orientation(TriangleMesh& tm, bool is_positive, const NamedParameters& np)
{
  typedef boost::graph_traits<TriangleMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::face_descriptor face_descriptor;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<TriangleMesh,
      NamedParameters>::const_type Vpm;
  typedef typename CGAL::Polygon_mesh_processing::GetFaceIndexMap<TriangleMesh,
      NamedParameters>::const_type Fid_map;

  Vpm vpm = boost::choose_param(get_param(np, CGAL::internal_np::vertex_point),
                                CGAL::get_const_property_map(boost::vertex_point, tm));

  Fid_map fid_map = boost::choose_param(get_param(np, CGAL::internal_np::face_index),
                                        CGAL::get_const_property_map(boost::face_index, tm));

  std::vector<std::size_t> face_cc(num_faces(tm), std::size_t(-1));

  // set the connected component id of each face
  std::size_t nb_cc = PMP::connected_components(tm,
                                           CGAL::bind_property_maps(fid_map,CGAL::make_property_map(face_cc)),
                                           PMP::parameters::face_index_map(fid_map));

  // extract a vertex with max z coordinate for each connected component
  std::vector<vertex_descriptor> xtrm_vertices(nb_cc, Graph_traits::null_vertex());
  for(vertex_descriptor vd : vertices(tm))
  {
    face_descriptor test_face = face(halfedge(vd, tm), tm);
    if(test_face == Graph_traits::null_face())
      test_face = face(opposite(halfedge(vd, tm), tm), tm);
    std::size_t cc_id = face_cc[get(fid_map,test_face )];
    if (xtrm_vertices[cc_id]==Graph_traits::null_vertex())
      xtrm_vertices[cc_id]=vd;
    else
      if (get(vpm, vd).z()>get(vpm,xtrm_vertices[cc_id]).z())
        xtrm_vertices[cc_id]=vd;
  }
  std::vector<std::vector<face_descriptor> > ccs(nb_cc);
  for(face_descriptor fd : faces(tm))
  {
    ccs[face_cc[get(fid_map,fd)]].push_back(fd);
  }

  //test ccs orientation
  for(std::size_t id=0; id<nb_cc; ++id)
  {
    if((!PMP::internal::is_outward_oriented(xtrm_vertices[id], tm, np)
         &&  is_positive)
       || (PMP::internal::is_outward_oriented(xtrm_vertices[id], tm, np)
           &&  !is_positive))
    {
      std::cerr<<" the orientation failed"<<std::endl;
      return false;
    }
  }
  return true;
}


int main()
{

  std::ifstream input("data-coref/nested_cubes_invalid_volume.off");
  assert(input);
  SMesh sm1, sm2, sm3, sm4, volume;
  input >> sm1;
  sm2 = sm1;
  sm3 = sm1;
  sm4 = sm1;
  volume = sm1;
  PMP::orient(sm1);
  if(!test_orientation(sm1, true, PMP::parameters::all_default()))
    return 1;
  typedef boost::property_map<SMesh, CGAL::vertex_point_t>::type Ppmap;
  typedef boost::property_map<SMesh, CGAL::face_index_t>::type Fidmap;
  Ppmap vpmap2 = get(CGAL::vertex_point, sm2);
  Fidmap fidmap2 = get(CGAL::face_index, sm2);

  PMP::orient(sm2, PMP::parameters::vertex_point_map(vpmap2)
                                   .face_index_map(fidmap2));
  if(!test_orientation(sm2, true, PMP::parameters::vertex_point_map(vpmap2)
                       .face_index_map(fidmap2)))
  {
    std::cerr << "ERROR for test1\n";
    return 1;
  }

  PMP::orient(sm3, PMP::parameters::outward_orientation(false));
  if(!test_orientation(sm3, false, PMP::parameters::all_default()))
  {
    std::cerr << "ERROR for test2\n";
    return 1;
  }

  Ppmap vpmap4 = get(CGAL::vertex_point, sm4);
  Fidmap fidmap4 = get(CGAL::face_index, sm4);

  PMP::orient(sm4, PMP::parameters::vertex_point_map(vpmap4)
                                   .face_index_map(fidmap4)
                                   .outward_orientation(false));
  if(!test_orientation(sm4, false, PMP::parameters::vertex_point_map(vpmap4)
                       .face_index_map(fidmap4)))
  {
    std::cerr << "ERROR for test3\n";
    return 1;
  }

  PMP::orient_to_bound_a_volume(volume);
  if( !PMP::does_bound_a_volume(volume))
  {
    std::cerr << "ERROR for test4\n";
    return 1;
  }

  return 0;
}
