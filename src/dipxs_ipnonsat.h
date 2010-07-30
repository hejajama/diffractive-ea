#ifndef DIPXS_IPNONSAT_H
#define DIPXS_IPNONSAT_H

/*
 * Dipole cross section in IP non sat model
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
#include "dipxs.h"
#include "nucleus.h"
 
class DipXS_IPNonSat : public DipXS
{
    public:
        DipXS_IPNonSat(Nucleus nucleus_);
        REAL DipXSection_b(REAL r, REAL xbj, Vec b ); // Impact parameter representation
        REAL DipXSection_b(REAL r, REAL xbj, REAL b ); // b=|b|
    
    private:

};


#endif
