#ifndef _TRANSLATE_
#define _TRANSLATE_

// translate each vertex or a polyhedron by a specified vector

template <class Polyhedron, class Kernel>
void translate(Polyhedron& polyhedron,
               const typename Kernel::Vector_3& translate)
{
  typename Polyhedron::Vertex_iterator v;
  for(v = polyhedron.vertices_begin();
      v != polyhedron.vertices_end();
      v++)
    v->point() = v->point() + translate;
}

#endif // _TRANSLATE_
