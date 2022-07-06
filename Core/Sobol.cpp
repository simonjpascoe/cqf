#include "stdafx.h"
#include <time.h>

#include "Sobol.h"

namespace Common {

//
// Implementation translated from Numerical Receipes in C (Press et al)
//
Sobol::Sobol(unsigned long dimensions) : m_dimensions(dimensions)
{
    if (m_dimensions > MAXDIM) { throw "Requested dimensions for Sobol generator too large.";}

    m_ip[0] = 0; m_ip[1] = 1; m_ip[2] = 1; m_ip[3] = 2; m_ip[4] = 1; m_ip[5] = 4;

    m_mdeg[0] = 1; m_mdeg[1] = 2; m_mdeg[2] = 3; m_mdeg[3] = 3; m_mdeg[4] = 4; m_mdeg[5] = 4;

    m_iv[0]  = 1;  m_iv[1]  = 1;  m_iv[2]  = 1; m_iv[3]  = 1;  m_iv[4]  = 1;  m_iv[5] = 1;
    m_iv[6]  = 3;  m_iv[7]  = 1;  m_iv[8]  = 3; m_iv[9]  = 1;  m_iv[10] = 3;  m_iv[11] = 1;
    m_iv[12] = 5;  m_iv[13] = 7;  m_iv[14] = 7; m_iv[15] = 5;  m_iv[16] = 1;  m_iv[17] = 3;
    m_iv[18] = 15; m_iv[19] = 11; m_iv[20] = 5; m_iv[21] = 3;  m_iv[22] = 1;  m_iv[23] = 7;

    int j,k,l;
    unsigned i,ipp;
    for (k=0; k<MAXDIM; k++) { m_ix.push_back(0); }
    m_in=0;
    m_fac = 1.0 / (1 << MAXBIT);
    for (j=0,k=0; j<MAXBIT; j++,k+=MAXDIM) { m_iu[j] = &m_iv[k];}
    for (k=0; k<MAXDIM; k++) {
        for (j=0; j<m_mdeg[k]; j++) { m_iu[j][k] <<= (MAXBIT-1-j);}
        for (j=m_mdeg[k]; j<MAXBIT; j++) {
            ipp=m_ip[k];
            i=m_iu[j-m_mdeg[k]][k];
            i ^= (i >> m_mdeg[k]);
            for (l=m_mdeg[k]-1; l>=1; l--) {
                if (ipp & 1) { i ^= m_iu[j-l][k];}
                ipp >>= 1;
            }
            m_iu[j][k]=i;
        }
    }
}

Sobol::~Sobol(void)
{
}

std::vector<double> Sobol::Draw()
{
    std::vector<double> x(m_dimensions);

    auto im = m_in++;
    int j;
    for (j=0; j<MAXBIT; j++) {
        if (!(im & 1)) break;
        im >>= 1;
    }

    im=j*MAXDIM;
    for (unsigned k=0; k<m_dimensions; k++)
    {
        m_ix[k] ^= m_iv[im+k];
        x[k]=m_ix[k]*m_fac;
    }

    return x;
}

// skip the first n draws. this is the same as Draw, except doesn't return the 
// x vector, and doesn't even allocate.
void Sobol::Skip(unsigned long n)
{
    for (unsigned long c=0; c<n; c++) {
        auto im = m_in++;
        int j;
        for (j=0; j<MAXBIT; j++) {
            if (!(im & 1)) break;
            im >>= 1;
        }

        im=j*MAXDIM;
        for (unsigned k=0; k<m_dimensions; k++)
        {
            m_ix[k] ^= m_iv[im+k];
        }
    }
}

}