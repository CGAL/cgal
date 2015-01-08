struct Distance {
  typedef Point Query_item;
  typedef double FT;
  typedef CGAL::Dimension_tag<3> D;

  double transformed_distance(const Point& p1, const Point& p2) const {
    double distx= p1.x()-p2.x();
    double disty= p1.y()-p2.y();
    double distz= p1.z()-p2.z();
    return distx*distx+disty*disty+distz*distz;
  }

  double min_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<FT,D>& b) const {
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }

  double min_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<double>& dists){   
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)){
      dists[0] = (b.min_coord(0)-h);
      distance += dists[0]*dists[0];
    }
    if (h > b.max_coord(0)){
      dists[0] = (h-b.max_coord(0));
      distance += dists[0]*dists[0];
    }
    h=p.y();
    if (h < b.min_coord(1)){
      dists[1] = (b.min_coord(1)-h);
      distance += dists[1]*dists[1];
    }
    if (h > b.max_coord(1)){
      dists[1] = (h-b.max_coord(1));
      distance += dists[1]*dists[1];
    }
    h=p.z();
    if (h < b.min_coord(2)){
      dists[2] = (b.min_coord(2)-h);
      distance += dists[2]*dists[2];
    }
    if (h > b.max_coord(2)){
      dists[2] = (h-b.max_coord(2));
      distance += dists[2]*dists[2];
    }
    return distance;
  }


  double max_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<FT,D>& b) const {
    double h = p.x();

    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.y();
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }

  double max_distance_to_rectangle(const Point& p,
				   const CGAL::Kd_tree_rectangle<FT,D>& b,std::vector<double>& dists){   
    double h = p.x();

    dists[0] = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0)) : (b.max_coord(0)-h);
    
    h=p.y();
    dists[1] = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1)) : (b.max_coord(1)-h);
    h=p.z();
    dists[2] = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2)) : (b.max_coord(2)-h);
    return dists[0] * dists[0] + dists[1] * dists[1] + dists[2] * dists[2];
  }

  double new_distance(double& dist, double old_off, double new_off,
		      int /* cutting_dimension */)  const {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return std::sqrt(d); }

}; // end of struct Distance
