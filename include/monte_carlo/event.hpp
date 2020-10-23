// Tiny container struct to hold four vectors into some accept/reject function
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#ifndef _PEVENT_
#define _PEVENT_

#include "TLorentzVector.h"

namespace jpacPhoto
{
    class event
    {
        public:
        event(){};

        event(TLorentzVector ep, TLorentzVector em, TLorentzVector prec)
        : p_ep(ep), p_em(em), p_prec(prec)
        {};

        TLorentzVector p_ep, p_em, p_prec;

        inline void update(TLorentzVector ep, TLorentzVector em, TLorentzVector prec)
        {
            p_ep = ep; p_em = em; p_prec = prec;
        };
    };
};

#endif
