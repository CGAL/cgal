#ifndef CGAL_FILTRED_CIRCULATOR_H
#define CGAL_FILTRED_CIRCULATOR_H

CGAL_BEGIN_NAMESPACE

template <class Circ, class Pred>
class Filtred_circulator : public Circ
{
  bool is_null;
  Pred test;

public:
  typedef Filtred_circulator<Circ,Pred> Self;


  Filtred_circulator(const Pred p=Pred()): is_null(true), test(p) {};

  Filtred_circulator(const Self& c): Circ(c), is_null(c.is_null),
    test(c.test) {};

  Self& operator=(const Self& c)
    {
      //this->Circ::operator=(c); // This does not work with bcc
      //*this = c;  // This does not work with VC++6
      static_cast<Circ&>(*this) = static_cast<const Circ&>(c);
      is_null=c.is_null;
      test=c.test;
      return *this;
    }

  Filtred_circulator(const Circ& c, const Pred& p=Pred())
    : Circ(c), is_null(false), test(p)
    {
      if(test(**this))
	is_null=false;
      else
	{
	  Self end(*this);
	  do { 
	    this->Circ::operator++();
	  } while( !test(**this) && (*this)!=end );
	  if((*this)==end)
	    is_null=true;
	}
    };

  bool operator==( CGAL_NULL_TYPE ) const {
    return is_null;
  }

  bool operator!=( CGAL_NULL_TYPE ) const {
    return !is_null;
  }

  bool operator==(const Self& c) const 
    {
      return is_null==c.is_null && this->Circ::operator==(c);
    }

  bool operator!=( const Self& c) const { return !(*this == c); }

  Self& operator++() {
    CGAL_assertion(!is_null);
    do {
      this->Circ::operator++();
    } while( !test(**this) );
    return *this;
  }

  Self  operator++(int) {
    Self tmp= *this;
    ++*this;
    return tmp;
  }
};

CGAL_END_NAMESPACE

#endif
