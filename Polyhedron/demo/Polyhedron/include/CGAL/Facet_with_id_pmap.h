template<class FaceGraph, class ValueType>
struct Facet_with_id_pmap
    : public boost::put_get_helper<ValueType&,
                                   Facet_with_id_pmap<FaceGraph,ValueType> >
{
  typedef typename boost::graph_traits<FaceGraph>::face_descriptor key_type;
  typedef ValueType value_type;
  typedef value_type& reference;
  typedef boost::lvalue_property_map_tag category;
  
  Facet_with_id_pmap(const FaceGraph& fg,
                     std::vector<ValueType>& internal_vector)
    : internal_vector(internal_vector), fg(fg) { }
  
  reference operator[](key_type key) const
  { return internal_vector[get(CGAL::face_index, fg, key)]; }
private:
  std::vector<ValueType>& internal_vector;
  const FaceGraph& fg; 
};
