CGAL_BEGIN_NAMESPACE

#ifndef CGAL_ALGEBRAIC_1_IMPL_H
#define CGAL_ALGEBRAIC_1_IMPL_H

template<class T>bool Algebraic_1::operator==(const T &n2)const{
	if(contains(n2))
		if(is_point())
			return true;
		else
			overlap();
	return false;
};

template<class T>bool Algebraic_1::operator!=(const T &n2)const{
	return !(operator==(n2));
};

template<class T>bool Algebraic_1::operator<=(const T &n2)const{
	return ((operator==(n2))||(operator<(n2)));
};

template<class T>bool Algebraic_1::operator>=(const T &n2)const{
	return ((operator==(n2))||(operator>(n2)));
};

template<class T>bool operator==(const T &n1,const Algebraic_1 &n2){
	return (n2==n1);
}

template<class T>bool operator!=(const T &n1,const Algebraic_1 &n2){
	return (n2!=n1);
}

template<class T>bool operator<(const T &n1,const Algebraic_1 &n2){
	return (n2>n1);
}

template<class T>bool operator>(const T &n1,const Algebraic_1 &n2){
	return (n2<n1);
}

template<class T>bool operator<=(const T &n1,const Algebraic_1 &n2){
	return (n2>=n1);
}

template<class T>bool operator>=(const T &n1,const Algebraic_1 &n2){
	return (n2<=n1);
}

#endif // CGAL_ALGEBRAIC_1_IMPL_H

CGAL_END_NAMESPACE
