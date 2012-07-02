template<typename String>
class TranslationInfo
{
public:
  
  TranslationInfo() : color(0)
  {
  }
  
  TranslationInfo(const String& name)
    : name_of_translation(name), color(0)
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
