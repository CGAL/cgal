#ifndef TRANSLATION_INFO
#define TRANSLATION_INFO

#include <CGAL/Hyperbolic_isometry_2.h>

template<typename String>
class TranslationInfo
{
public:
  
  TranslationInfo() : color(0)
  {
  }
  
  TranslationInfo(const String& name, int new_color = 0)
    : name_of_translation(name), color(new_color)
  {
  }
  
  void setString(const String& name)
  {
    name_of_translation = name;
  }
  
  const String& toString() const
  {
    return name_of_translation;
  }
  
  void setColor(int new_color)
  {
    color = new_color;
  }
  
  int getColor() const
  {
    return color;
  }
   
private:
  String name_of_translation;
  
  int color; 
};

template<typename String>
TranslationInfo<String> operator + (const TranslationInfo<String>& l, const TranslationInfo<String>& r)
{
  TranslationInfo<String> result;
  result.setString(l.toString() + r.toString());
  
  return result;
}

template<typename Gt, typename Info>
class IsometryWithInfo : public CGAL::Hyperbolic_isometry_2<Gt>
{
public:
  typedef CGAL::Hyperbolic_isometry_2<Gt> Base;
  //typedef Base::Point Point;
  
  IsometryWithInfo(const Base& base, const Info& new_info = Info()) : Base(base), _info(new_info)
  {
  }
  
  IsometryWithInfo operator * (const IsometryWithInfo& other) const
  {
    Base base = this->Base::operator *(other);
    Info new_info = info() + other.info();
    
    return IsometryWithInfo(base, new_info);
  }
  
  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
  
private:
  Info _info;
};

#endif // TRANSLATION_INFO
