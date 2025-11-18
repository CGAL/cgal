namespace CGAL {
/*!
  \ingroup PkgHyperbolicSurfaceTriangulation2DirichletDomain

  \return a vector with the vertices of a Dirichlet domain whose base point is translated to the origin of the Poincaré disk.

  `domain` must be a fundamental domain whose vertices represent a same point on the corresponding hyperbolic surface.
*/
std::vector<typename Traits::Hyperbolic_Voronoi_point_2> Dirichlet_vertices(Hyperbolic_fundamental_domain_2<Traits> & domain);

} // namespace CGAL