#ifndef CGAL_HYPERBOLIC_TRANSLATIONS_2_H
#define CGAL_HYPERBOLIC_TRANSLATIONS_2_H

#include <CGAL/Hyperbolic_isometry_2.h>

template<typename Translation>
struct Element
{
    typedef typename std::list<Element> List;
    typedef CGAL::Circulator_from_container<List> Circulator;
    
    Translation g;
    
    // circulator iterator to an inverse translation in the list
    Circulator inverse;
};

template<typename Gt>
class Diametric_translations
{
public:
    typedef typename Gt::FT FT;
    typedef typename std::complex<FT> complex;
    typedef CGAL::Hyperbolic_isometry_2<Gt> Hyperbolic_isometry;
    
    typedef Element<Hyperbolic_isometry> Element_t;
    typedef std::list<Element_t> List;
    typedef typename List::iterator List_iterator;
    typedef CGAL::Circulator_from_container<List> Circulator;
    
    typedef std::pair<Hyperbolic_isometry, int> Node;
    typedef std::vector<Node> Vector;
    typedef typename Vector::iterator Vector_iterator;
    
    
    static Hyperbolic_isometry& a()
    {
        compute();
        return g[0].first;
    }
    
    static Hyperbolic_isometry& b()
    {
        compute();
        return g[5].first;
    } 
    
    static Hyperbolic_isometry& c()
    {
        compute();
        return g[2].first;
    }
    
    static Hyperbolic_isometry& d()
    {
        compute();
        return g[7].first;
    }
    
    static Hyperbolic_isometry& inva()
    {
        compute();
        return g[4].first;
    }
    
    static Hyperbolic_isometry& invb()
    {
        compute();
        return g[1].first;
    }
    
    static Hyperbolic_isometry& invc()
    {
        compute();
        return g[6].first;
    }
    
    static Hyperbolic_isometry& invd()
    {
        compute();
        return g[3].first;
    }
    
    static const Vector& get_vector_of_translations()
    {
        compute();
        return g;
    }
    
    static Vector_iterator vector_begin()
    {
        compute();
        return g.begin();
    }
    
    static Vector_iterator vector_end()
    {
        compute();
        return g.end();
    }
    
    static List_iterator list_begin()
    {
        compute();
        return l.begin();
    }
    
    static List_iterator list_end()
    {
        compute();
        return l.end();
    }
    
    static List& list()
    {
        compute();
        return l;
    }
    
private:
    
    static void compute_g()
    {
        const FT a1 = FT(1) + CGAL::sqrt(2.);
        const FT a2 = FT(0);
        
        std::complex<FT> a(a1, a2);
        
        // This is the multiplicative factor for all b's
        const FT bb = CGAL::sqrt(2.) * CGAL::sqrt( a1 );
        
        std::vector< std::complex<FT> > b;
        
        // Euler's identity: exp(i k \theta) = cos(k \theta) + i sin(k \theta)
        for (int k = -2; k < 6; k++) {
            b.push_back( std::complex<FT>( bb * cos( FT(k) * CGAL_PI / FT(4)), bb * sin( FT(k) * CGAL_PI / FT(4) ) ) );
        }
        
        g.resize(8);
        
        
         
         // a
        int index = 0;
        int invindex;
        
        g[index].first = Hyperbolic_isometry(a, b[index]);
        g[index].second = 4;        // Note that second = (index + 4) % 8
        invindex = g[index].second;
        g[invindex].first = Hyperbolic_isometry(a, b[invindex]);
        g[invindex].second = index;
         
        // b
        index = 5;
        g[index].first = Hyperbolic_isometry(a, b[index]);
        g[index].second = 1;
        invindex = g[index].second;
        g[invindex].first = Hyperbolic_isometry(a, b[invindex]);
        g[invindex].second = index;
         
        // c
        index = 2;
        g[index].first = Hyperbolic_isometry(a, b[index]);
        g[index].second = 6;
        invindex = g[index].second;
        g[invindex].first = Hyperbolic_isometry(a, b[invindex]);
        g[invindex].second = index;
         
        // d
        index = 7;
        g[index].first = Hyperbolic_isometry(a, b[index]);
        g[index].second = 3;
        invindex = g[index].second;
        g[invindex].first = Hyperbolic_isometry(a, b[invindex]);
        g[invindex].second = index;
        
    }
    
    
    static void compute_l()
    {
        l.resize(g.size());
        
        std::vector<Circulator> aux_list;
        aux_list.reserve(8);
        
        for(List_iterator li = l.begin(); li != l.end(); li++) {
            aux_list.push_back( Circulator(&l, li) );
        }
        
        for(typename List::size_type i = 0; i < aux_list.size(); i++) {
            aux_list[i]->g = g[i].first;
            aux_list[i]->inverse = aux_list[g[i].second];
        }
    }
    
    static void compute()
    {
        static bool computed = false;
        if(!computed) {
            compute_g();
            compute_l();
            computed = true;
        }
    }
    
    static Vector g;
    
    static List l;
};

// default initialization

template<typename Gt>
typename Diametric_translations<Gt>::Vector
Diametric_translations<Gt>::g;

template<typename Gt>
typename Diametric_translations<Gt>::List
Diametric_translations<Gt>::l;

#endif // CGAL_HYPERBOLIC_TRANSLATIONS_2_H
