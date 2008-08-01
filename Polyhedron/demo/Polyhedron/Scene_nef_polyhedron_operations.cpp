#include "Scene.h"
#include "Nef_type.h"

Nef_polyhedron* Scene::new_nef_polyhedron() 
{
  return new Nef_polyhedron;
}

Nef_polyhedron* Scene::copy_nef_polyhedron(Nef_polyhedron* poly)
{
  return new Nef_polyhedron(*poly);
}

void Scene::destroy_nef_polyhedron(Nef_polyhedron* poly)
{
  delete poly;
}
