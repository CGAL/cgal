inline
Oriented_side
power_testH3_SAF(
    const Static_filter_error &phx,
    const Static_filter_error &phy,
    const Static_filter_error &phz,
    const Static_filter_error &phw,
    const Static_filter_error &pwt,
    const Static_filter_error &qhx,
    const Static_filter_error &qhy,
    const Static_filter_error &qhz,
    const Static_filter_error &qhw,
    const Static_filter_error &qwt,
    const Static_filter_error &rhx,
    const Static_filter_error &rhy,
    const Static_filter_error &rhz,
    const Static_filter_error &rhw,
    const Static_filter_error &rwt,
    const Static_filter_error &shx,
    const Static_filter_error &shy,
    const Static_filter_error &shz,
    const Static_filter_error &shw,
    const Static_filter_error &swt,
    const Static_filter_error &thx,
    const Static_filter_error &thy,
    const Static_filter_error &thz,
    const Static_filter_error &thw,
    const Static_filter_error &twt,
    double & epsilon_0)
{
  typedef Static_filter_error RT;

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphz = phz*phw;
    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) + square(phz) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhz = qhz*qhw;
    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) + square(qhz) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhz = rhz*rhw;
    RT drhw = square(rhw);
    RT drz = square(rhx) + square(rhy) + square(rhz) - rwt*drhw;

    RT dshx = shx*shw;
    RT dshy = shy*shw;
    RT dshz = shz*shw;
    RT dshw = square(shw);
    RT dsz = square(shx) + square(shy) + square(shz) - swt*dshw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthz = thz*thw;
    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) + square(thz) - twt*dthw;

    return Oriented_side(- sign_of_determinant5x5_SAF(dphx, dphy, dphz, dpz, dphw,
	                                        dqhx, dqhy, dqhz, dqz, dqhw,
	                                        drhx, drhy, drhz, drz, drhw,
	                                        dshx, dshy, dshz, dsz, dshw,
	                                        dthx, dthy, dthz, dtz, dthw,
		epsilon_0));
}

inline
Oriented_side
power_testH3_SAF(
    const Restricted_double &phx,
    const Restricted_double &phy,
    const Restricted_double &phz,
    const Restricted_double &phw,
    const Restricted_double &pwt,
    const Restricted_double &qhx,
    const Restricted_double &qhy,
    const Restricted_double &qhz,
    const Restricted_double &qhw,
    const Restricted_double &qwt,
    const Restricted_double &rhx,
    const Restricted_double &rhy,
    const Restricted_double &rhz,
    const Restricted_double &rhw,
    const Restricted_double &rwt,
    const Restricted_double &shx,
    const Restricted_double &shy,
    const Restricted_double &shz,
    const Restricted_double &shw,
    const Restricted_double &swt,
    const Restricted_double &thx,
    const Restricted_double &thy,
    const Restricted_double &thz,
    const Restricted_double &thw,
    const Restricted_double &twt,
    const double & epsilon_0)
{
  typedef Restricted_double RT;

    RT dphx = phx*phw;
    RT dphy = phy*phw;
    RT dphz = phz*phw;
    RT dphw = square(phw);
    RT dpz = square(phx) + square(phy) + square(phz) - pwt*dphw;

    RT dqhx = qhx*qhw;
    RT dqhy = qhy*qhw;
    RT dqhz = qhz*qhw;
    RT dqhw = square(qhw);
    RT dqz = square(qhx) + square(qhy) + square(qhz) - qwt*dqhw;

    RT drhx = rhx*rhw;
    RT drhy = rhy*rhw;
    RT drhz = rhz*rhw;
    RT drhw = square(rhw);
    RT drz = square(rhx) + square(rhy) + square(rhz) - rwt*drhw;

    RT dshx = shx*shw;
    RT dshy = shy*shw;
    RT dshz = shz*shw;
    RT dshw = square(shw);
    RT dsz = square(shx) + square(shy) + square(shz) - swt*dshw;

    RT dthx = thx*thw;
    RT dthy = thy*thw;
    RT dthz = thz*thw;
    RT dthw = square(thw);
    RT dtz = square(thx) + square(thy) + square(thz) - twt*dthw;

    return Oriented_side(- sign_of_determinant5x5_SAF(dphx, dphy, dphz, dpz, dphw,
	                                        dqhx, dqhy, dqhz, dqz, dqhw,
	                                        drhx, drhy, drhz, drz, drhw,
	                                        dshx, dshy, dshz, dsz, dshw,
	                                        dthx, dthy, dthz, dtz, dthw,
		epsilon_0));
}

