#ifndef MULTI_INFO_H
#define MULTI_INFO_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <list>

template<class Info_t>
class Multi_info
{
public:
  typedef Info_t                              Info;
  typedef std::list<Info>                     Info_list;
  typedef typename Info_list::const_iterator  Info_iterator;
  typedef typename Info_list::size_type       size_type;


  Multi_info() {}

  Multi_info(Info info) {
    info_list_.push_back(info);
  }

  Multi_info(Info_list info_list) : info_list_(info_list) {}

  template<typename InputIterator>
  Multi_info(InputIterator first, InputIterator beyond)
    : info_list_(first, beyond) {}

  size_type size() const { return info_list_.size(); }
  bool is_empty() const { return info_list_.empty(); }

  const Info_list& info_list() const { return info_list_; }

  Info_iterator info_begin() const {
    return info_list_.begin();
  }

  Info_iterator info_end() const {
    return info_list_.end();
  }

private:
  Info_list info_list_;
};

template<class Info>
std::ostream&
operator<<(std::ostream& os, const Multi_info<Info>& info)
{
  if ( info.is_empty() ) {
    return os << "{}";
  }

  typedef typename Multi_info<Info>::Info_iterator iterator;
  iterator last = --info.info_end();
  os << "{";
  for (iterator it = info.info_begin(); it != last; ++it) {
    os << *it << ", ";
  }
  os << *last << "}";

  return os;
}


template<class Info_t>
struct Multi_info_convert_info
{
  typedef Info_t             Info;
  typedef Multi_info<Info>   result_type;

  inline
  Multi_info<Info> operator()(const Multi_info<Info>& minfo0, bool) const
  {
    return minfo0;
  }

  inline
  Multi_info<Info> operator()(const Multi_info<Info>& minfo0,
                              const Multi_info<Info>& minfo1, bool) const
  {
    return minfo0;
  }
};

template<class Info_t>
struct Multi_info_merge_info
{
  typedef Info_t            Info;
  typedef Multi_info<Info>  result_type;

  inline
  Multi_info<Info> operator()(const Multi_info<Info>& minfo0,
                              const Multi_info<Info>& minfo1) const
  {
    typedef typename Multi_info<Info>::Info_list       Info_list;

    Info_list merged_info = minfo0.info_list();
    Info_list copy = minfo1.info_list();

    merged_info.splice(merged_info.end(), copy);
    return merged_info;
  }
};


#endif // MULTI_INFO_H
