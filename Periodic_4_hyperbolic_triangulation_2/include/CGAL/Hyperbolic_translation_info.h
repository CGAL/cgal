#ifndef HYPERBOLIC_TRANSLATION_INFO
#define HYPERBOLIC_TRANSLATION_INFO

template<typename String>
class Hyperbolic_translation_info
{
public:
  
  Hyperbolic_translation_info() : color(0)
  {
  }
  
  Hyperbolic_translation_info(const String& name, int new_color = 0)
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
Hyperbolic_translation_info<String> operator + (const Hyperbolic_translation_info<String>& l, const Hyperbolic_translation_info<String>& r)
{
  Hyperbolic_translation_info<String> result;
  result.setString(l.toString() + r.toString());
  
  return result;
}



#endif // HYPERBOLIC_TRANSLATION_INFO
