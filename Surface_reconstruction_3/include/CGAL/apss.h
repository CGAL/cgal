
#ifndef _APSS_h_
#define _APSS_h_

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Surface_mesher/Implicit_surface_oracle_3.h>
#include <CGAL/Min_sphere_d.h>
#include <CGAL/Optimisation_d_traits_3.h>

#include "apss_oracle.h"

namespace CGAL
{

template<
    typename GT,
    typename PointList,
    typename NormalList
    >
class APSS 
{
public:
    typedef GT Geom_traits;
    typedef typename Geom_traits::Sphere_3 Sphere_3;
    typedef typename Geom_traits::FT Scalar;
    typedef typename Geom_traits::Point_3 Point;
    typedef APSS<Geom_traits,PointList,NormalList> Self;
    
//     typedef Surface_mesher::Implicit_surface_oracle_3< Geom_traits, Self> Surface_mesher_traits_3;
    typedef ApssOracle<Geom_traits, Self> Surface_mesher_traits_3;
    friend class ApssOracle<Geom_traits, Self> ;

private:

    class KdTreeElement : public Geom_traits::Point_3
    {
    public:
        unsigned int index;
        KdTreeElement() {}
        
        KdTreeElement(const typename Geom_traits::Point_3& p, unsigned int id=0)
            : Geom_traits::Point_3(p), index(id)
        {}
    };

    class KdTreeGT : public Geom_traits
    {
    public:
        typedef KdTreeElement Point_3;
    };

    typedef Search_traits_3<KdTreeGT> TreeTraits;
    typedef Orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
    typedef typename Neighbor_search::Tree Tree;
    

public:

    APSS( const PointList& apoints, const NormalList& anormals)
    {
        m = new Private(apoints,anormals);
        
        unsigned int i=0;
        m->treeElements.reserve(m->points.size());
        for (typename PointList::const_iterator it=m->points.begin() ; it<m->points.end() ; ++it,++i)
        {
            m->treeElements.push_back(KdTreeElement(*it,i));
        }
        m->tree = new Tree(m->treeElements.begin(), m->treeElements.end());
        
        // Compute the radius of each points.
        // The union of these ball defines the surface definition domain.
        i = 0;
        for (typename PointList::const_iterator it=m->points.begin() ; it<m->points.end() ; ++it,++i)
        {
            KdTreeElement query(*it);
            Neighbor_search search(*(m->tree), query, 16);
            Scalar maxdist2 = 0.;
            for (typename Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
            {
                if (maxdist2<it->second)
                    maxdist2 = it->second;
            }
            m->radii.push_back(2.*sqrt(maxdist2/16.));
        }
        
        Min_sphere_d< CGAL::Optimisation_d_traits_3<GT> > ms3(m->points.begin(), m->points.end());
        
        m->boundingSphere = Sphere_3(ms3.center(), 1.2*sqrt(ms3.squared_radius()));
//         m->boundingSphere = Sphere_3(CGAL::ORIGIN, 2);
        
        m->sqError = 0.0000001 * GT().compute_squared_radius_3_object()(m->boundingSphere);
    }
    
    APSS(const APSS& other) {
        m = other.m;
        m->count++;
    }
    
    APSS& operator = (const APSS& other) {
        m = other.m;
        m->count++;
    }
    
    ~APSS() {
        if (--(m->count)==0)
            delete m;
    }
    
    void setNofNeighbors(unsigned int k) {m->nofNeighbors = k;}
    
    /** Fit an algebraic sphere on a set of neigbors in a Moving Least Square sense.
        The weight function is scaled such that the weight of the furthest neighbor is 0.
    */
    void fit(Neighbor_search& neighbors) const
    {
        typename Neighbor_search::iterator last = neighbors.end();
        last--;
        Scalar r2 = last->second;
        
        Point sumP(0,0,0);
        Point sumN(0,0,0);
        Scalar sumDotPP = 0.;
        Scalar sumDotPN = 0.;
        Scalar sumW = 0.;
        
        r2 *= 1.001;
        Scalar invr2 = 1./r2;
        for (typename Neighbor_search::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            const Point& p = m->points[it->first.index];
            const Point& n = m->normals[it->first.index];
            Scalar w = 1. - it->second*invr2;
            w = w*w; w = w*w;
            
            sumP = add(sumP,mul(w,p));
            sumN = add(sumN,mul(w,n));
            sumDotPP += w * dot(p,p);
            sumDotPN += w * dot(p,n);
            sumW += w;
        }
        
        Scalar invSumW = 1./sumW;
        m->as.u4 = 0.5 * (sumDotPN - invSumW*dot(sumP,sumN))/(sumDotPP - invSumW*dot(sumP,sumP));
        m->as.u13 = mul(invSumW,add(sumN,mul(-2.*m->as.u4,sumP)));
        m->as.u0 = -invSumW*(dot(m->as.u13,sumP)+m->as.u4*sumDotPP);
        m->as.finalize();
    }
    
    /** Check whether the point p is close to the input points or not.
    */
    inline bool isValid(const Neighbor_search& neighbors, const Point& p) const
    {
        for (typename Neighbor_search::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
        {
            Scalar r = 2. * m->radii[it->first.index];
            if (r*r > it->second)
                return true;
        }
        return false;
    }
    
    /** Compute the potentiel value.
    */
    Scalar operator()(const Point& p) {
        KdTreeElement query(p);
        Neighbor_search search(*(m->tree), query, m->nofNeighbors);
        unsigned int closestId = search.begin()->first.index;
        
        if (!isValid(search,p))
        {
            Point n;
            Point pp = p;
            project(pp,n,1);
            pp = sub(p,pp);
            return length(pp) * ( dot(n,pp)>0. ? 1. : -1.);
        }
        
        fit(search);
        return m->as.euclideanDistance(p);
    }
    
    /** Projects the point p onto the MLS surface.
    */
    void project(Point& p, unsigned int maxNofIterations = 20) const
    {
        Point n;
        project(p,n,maxNofIterations);
    }
    
    /** Projects the point p onto the MLS surface, and returns an approximate normal.
    */
    void project(Point& p, Point& n, unsigned int maxNofIterations = 20) const
    {
        Point source = p;
        
        Scalar delta2 = 0.;
        unsigned int countIter = 0;
        do {
        
            Neighbor_search search(*(m->tree), p, m->nofNeighbors);
            
            // if p is far away the input point cloud, start with the closest point.
            if (!isValid(search,p))
            {
                p = search.begin()->first;
                n = m->normals[search.begin()->first.index];
                delta2 = search.begin()->second;
            }
            else
            {
                fit(search);
                
                Point oldP = p;
                if (m->as.state==AlgebraicSphere::SPHERE)
                {
                    // projection onto a sphere.
                    Point dir = normalize(sub(source,m->as.center));
                    p = add(m->as.center,mul(m->as.radius,dir));
                    Scalar flipN = m->as.u4<0. ? -1. : 1.;
                    if (!isValid(search,p))
                    {
                        // if the closest intersection is far away the input points,
                        // then we take the other one.
                        p = add(m->as.center,mul(-m->as.radius,dir));
                        flipN = -flipN;
                    }
                    n = mul(flipN,dir);
                    
                    if (!isValid(search,p))
                    {
                        std::cout << "Invalid projection\n";
                    }
                }
                else if (m->as.state==AlgebraicSphere::PLANE)
                {
                    // projection onto a plane.
                    p = sub(source, mul(dot(m->as.normal,source)+m->as.d,m->as.normal));
                }
                else
                {
                    // iterative projection onto the algebraic sphere.
                    p = m->as.iProject(source);
                }
                
                oldP = sub(oldP,p);
                delta2 = dot(oldP,oldP);
            }
            
        } while ( ((++countIter)<maxNofIterations) && (delta2<m->sqError) );
    }

    const Scalar& squared_error_bound() const
    {
        return m->sqError;
    }

    const Sphere_3& bounding_sphere() const
    {
        return m->boundingSphere;
    }

public:

    inline static Scalar dot(const Point& a, const Point& b) {
        return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
    }
    inline static Scalar length(const Point& a) {
        return sqrt(dot(a,a));
    }
    inline static Point mul(Scalar s, const Point& p) {
        return Point(p.x()*s, p.y()*s, p.z()*s);
    }
    inline static Point add(Scalar s, const Point& p) {
        return Point(p.x()+s, p.y()+s, p.z()+s);
    }
    inline static Point add(const Point& a, const Point& b) {
        return Point(a.x()+b.x(), a.y()+b.y(), a.z()+b.z());
    }
    inline static Point sub(const Point& a, const Point& b) {
        return Point(a.x()-b.x(), a.y()-b.y(), a.z()-b.z());
    }
    inline static Point normalize(const Point& p) {
        Scalar s = 1. / length(p);
        return mul(s,p);
    }
    inline static Point cross(const Point& a, const Point& b) {
        return Point(a.y()*b.z() - a.z()*b.y(),
            a.z()*b.x() - a.x()*b.z(),
            a.x()*b.y() - a.y()*b.x());
    }

private:

    struct AlgebraicSphere {
        Scalar u0, u4;
        Point u13;
        enum State {UNDETERMINED=0,PLANE=1,SPHERE=2};
        State state;
        struct {Point center; Scalar radius;};
        struct {Point normal; Scalar d;};
        
        AlgebraicSphere() : state(UNDETERMINED) {}
        
        /** Converts the algebraix sphere to an explicit sphere or plane.
        */
        void finalize(void) {
            if (fabs(u4)>1e-9)
            {
                state = SPHERE;
                Scalar b = 1./u4;
                center = mul(-0.5*b,u13);
                radius = sqrt(dot(center,center) - b*u0);
            }
            else if (u4==0.)
            {
                state = PLANE;
                Scalar s = 1./length(u13);
                normal = mul(s,u13);
                d = u0*s;
            }
            else
            {
                state = UNDETERMINED;
            }
        }
        
        /** Compute the Euclidean distance between the algebraic surface and a point.
        */
        Scalar euclideanDistance(const Point& p) {
            if (state==SPHERE)
            {
                Scalar aux = length(sub(center,p)) - radius;
                if (u4<0.)
                    aux = -aux;
                return aux;
            }
            
            if (state==PLANE)
                return dot(p,normal) + d;
            
            // else, tedious case, fall back to an iterative method:
            return iEuclideanDistance(p);
        }
        
        /** Euclidean distance via an iterative projection procedure.
            This is an optimized version of distance(x,iProject(x)).
        */
        inline Scalar iEuclideanDistance(const Point& x) const
        {
            Scalar d = 0.;
            Point grad;
            Point dir = add(mul(2.*u4,x),u13);
            Scalar ilg = 1./length(dir);
            dir = mul(ilg,dir);
            Scalar ad = u0 + dot(u13,x) + u4 * dot(x,x);
            Scalar delta = -ad*(ilg<1.?ilg:1.);
            Point p = add(x, mul(delta,dir));
            d += delta;
            for (int i=0 ; i<5 ; ++i)
            {
                grad = add(mul(2.*u4,p),u13);
                ilg = 1./length(grad);
                delta = -(u0 + dot(u13,p) + u4 * dot(p,p))*(ilg<1.?ilg:1.);
                p = add(p,mul(delta,dir));
                d += delta;
            }
            return -d;
        }
        
        /** Iterative projection.
        */
        inline Point iProject(const Point& x) const
        {
            Point grad;
            Point dir = add(mul(2.*u4,x),u13);
            Scalar ilg = 1./length(dir);
            dir = mul(ilg,dir);
            Scalar ad = u0 + dot(u13,x) + u4 * dot(x,x);
            Scalar delta = -ad*(ilg<1.?ilg:1.);
            Point p = add(x, mul(delta,dir));
            for (int i=0 ; i<5 ; ++i)
            {
                grad = add(mul(2.*u4,p),u13);
                ilg = 1./length(grad);
                delta = -(u0 + dot(u13,p) + u4 * dot(p,p))*(ilg<1.?ilg:1.);
                p = add(p,mul(delta,dir));
            }
            return p;
        }
    };

    struct Private {
        Private( const PointList& apoints, const NormalList& anormals)
            : points(apoints), normals(anormals), nofNeighbors(12), count(1)
        {}
        Tree* tree;
        const PointList& points;
        const NormalList& normals;
        std::vector<KdTreeElement> treeElements;
        std::vector<Scalar> radii;
        Sphere_3 boundingSphere;
        Scalar sqError;
        unsigned int nofNeighbors;
        mutable AlgebraicSphere as;
        int count;
    };
    Private* m;
};

template <typename GT, typename PointList, typename NormalList>
APSS<GT, PointList, NormalList>
make_apss(GT, PointList points, NormalList normals)
{
    typedef APSS<GT, PointList, NormalList> Surface;
    return Surface(points, normals);
}


}

#endif
