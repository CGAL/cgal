inline
Oriented_side
power_testH2_SAF(
    const Static_filter_error &phx,
    const Static_filter_error &phy,
    const Static_filter_error &phw,
    const Static_filter_error &pwt,
    const Static_filter_error &qhx,
    const Static_filter_error &qhy,
    const Static_filter_error &qhw,
    const Static_filter_error &qwt,
    const Static_filter_error &rhx,
    const Static_filter_error &rhy,
    const Static_filter_error &rhw,
    const Static_filter_error &rwt,
    const Static_filter_error &thx,
    const Static_filter_error &thy,
    const Static_filter_error &thw,
    const Static_filter_error &twt,
    double & epsilon_0)
{
  typedef Static_filter_error RT;

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhw = square(rhw);
    RT drz = square(rhx) + square(rhy) - rwt*drhw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) - twt*dthw;

    return Oriented_side(sign_of_determinant4x4_SAF(dphx, dphy, dpz, dphw,
	                                        dqhx, dqhy, dqz, dqhw,
	                                        drhx, drhy, drz, drhw,
	                                        dthx, dthy, dtz, dthw,
		epsilon_0));
}

inline
Oriented_side
power_testH2_SAF(
    const Restricted_double &phx,
    const Restricted_double &phy,
    const Restricted_double &phw,
    const Restricted_double &pwt,
    const Restricted_double &qhx,
    const Restricted_double &qhy,
    const Restricted_double &qhw,
    const Restricted_double &qwt,
    const Restricted_double &rhx,
    const Restricted_double &rhy,
    const Restricted_double &rhw,
    const Restricted_double &rwt,
    const Restricted_double &thx,
    const Restricted_double &thy,
    const Restricted_double &thw,
    const Restricted_double &twt,
    const double & epsilon_0)
{
  typedef Restricted_double RT;

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhw = square(rhw);
    RT drz = square(rhx) + square(rhy) - rwt*drhw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) - twt*dthw;

    return Oriented_side(sign_of_determinant4x4_SAF(dphx, dphy, dpz, dphw,
	                                        dqhx, dqhy, dqz, dqhw,
	                                        drhx, drhy, drz, drhw,
	                                        dthx, dthy, dtz, dthw,
		epsilon_0));
}

inline
Oriented_side
power_testH2_SAF(
    const Static_filter_error &phx,
    const Static_filter_error &phy,
    const Static_filter_error &phw,
    const Static_filter_error &pwt,
    const Static_filter_error &qhx,
    const Static_filter_error &qhy,
    const Static_filter_error &qhw,
    const Static_filter_error &qwt,
    const Static_filter_error &thx,
    const Static_filter_error &thy,
    const Static_filter_error &thw,
    const Static_filter_error &twt,
    double & epsilon_0,
    double & epsilon_1)
{
  typedef Static_filter_error RT;

    
    
    RT pa, qa, ta;

    if (phx * qhw != qhx * phw )
    {
	pa = phx*phw;
	qa = qhx*qhw;
	ta = thx*thw;
    }
    else
    {   
	pa = phy*phw;
	qa = qhy*qhw;
	ta = thy*thw;
    }

    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) - pwt*dphw;

    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) - qwt*dqhw;

    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) - twt*dthw;

    return Oriented_side(CGAL::compare_SAF(pa, qa,
		epsilon_0) *
	                 sign_of_determinant3x3_SAF(pa, dpz, dphw,
				                qa, dqz, dqhw,
				                ta, dtz, dthw,
		epsilon_1));
}

inline
Oriented_side
power_testH2_SAF(
    const Restricted_double &phx,
    const Restricted_double &phy,
    const Restricted_double &phw,
    const Restricted_double &pwt,
    const Restricted_double &qhx,
    const Restricted_double &qhy,
    const Restricted_double &qhw,
    const Restricted_double &qwt,
    const Restricted_double &thx,
    const Restricted_double &thy,
    const Restricted_double &thw,
    const Restricted_double &twt,
    const double & epsilon_0,
    const double & epsilon_1)
{
  typedef Restricted_double RT;

    
    
    RT pa, qa, ta;

    if (phx * qhw != qhx * phw )
    {
	pa = phx*phw;
	qa = qhx*qhw;
	ta = thx*thw;
    }
    else
    {   
	pa = phy*phw;
	qa = qhy*qhw;
	ta = thy*thw;
    }

    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) - pwt*dphw;

    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) - qwt*dqhw;

    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) - twt*dthw;

    return Oriented_side(CGAL::compare_SAF(pa, qa,
		epsilon_0) *
	                 sign_of_determinant3x3_SAF(pa, dpz, dphw,
				                qa, dqz, dqhw,
				                ta, dtz, dthw,
		epsilon_1));
}

