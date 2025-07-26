#include <cassert>
#include <algorithm>
#include <cstdio>
#include "GradGen.h"

bool g_bExGEnd_MAG = true;

// interpolation function
template<typename T>
static T intp(double dXEv, double dX0, T tY0, double dX1, T tY1)
{
    if (std::fabs(dXEv-dX0) < std::fabs(dXEv-dX1))
    {
        return tY0 + (tY1-tY0)/(dX1-dX0) * (dXEv-dX0);
    }
    else
    {
        return tY1 + (tY0-tY1)/(dX0-dX1) * (dXEv-dX1);
    }
}

GradGen::GradGen
    (
        const TrajFunc* ptTraj,
        double dSLim, double dGLim,
        double dDt, int64_t lOs, 
        double dG0Norm, double dG1Norm
    ):
    m_ptTraj(ptTraj),
    m_dSLim(dSLim), 
    m_dGLim(dGLim), 
    m_dDt(dDt), 
    m_lOs(lOs), 
    m_dG0Norm(dG0Norm), 
    m_dG1Norm(dG1Norm)
{

}

GradGen::~GradGen()
{

}

bool GradGen::sovQDE(double* pdSol0, double* pdSol1, double dA, double dB, double dC)
{
    double dDelta = dB*dB - 4*dA*dC;
    if (dDelta<0) dDelta = 0;
    if (pdSol0) *pdSol0 = (-dB-std::sqrt(dDelta))/(2*dA);
    if (pdSol1) *pdSol1 = (-dB+std::sqrt(dDelta))/(2*dA);
    return true;
}

double GradGen::getCurRad(double dP)
{
    v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dP);
    v3 v3D2kDp2; m_ptTraj->getD2kDp2(&v3D2kDp2, dP);
    double dNume = pow(v3::norm(v3DkDp), 3e0);
    double dDeno = v3::norm(v3::cross(v3DkDp, v3D2kDp2));
    return dNume/dDeno;
}

double GradGen::getDp(const v3& v3G, double dDt, double dP, double dSignDp)
{
    // solve `ΔP` by RK2
    double dDl = v3::norm(v3G)*dDt;
    // k1
    double dK1;
    {
        v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dP);
        double dDlDp = v3::norm(v3DkDp)*dSignDp;
        dK1 = 1e0/dDlDp;
    }
    // k2
    double dK2;
    {
        v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dP+dK1*dDl);
        double dDlDp = v3::norm(v3DkDp)*dSignDp;
        dK2 = 1e0/dDlDp;
    }
    double dDp = dDl*(dK1 + dK2)/2e0;
    return dDp;
}

#if 1

double GradGen::getDp(const v3& v3GPrev, const v3& v3GThis, double dDt, double dPPrev, double dPThis, double dSignDp)
{
    // solve `ΔP` by RK2
    double dDl = v3::norm(v3GThis)*dDt;
    // k1
    double dK1;
    {
        v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dPThis);
        double dDlDp = v3::norm(v3DkDp)*dSignDp;
        dK1 = 1e0/dDlDp;
    }
    // k2
    double dK2;
    {
        v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dPThis+dK1*dDl);
        double dDlDp = v3::norm(v3DkDp)*dSignDp;
        dK2 = 1e0/dDlDp;
    }
    double dDp = dDl*(0.5*dK1 + 0.5*dK2);
    return dDp;
}

#else // less accurate due to estimation of PNext

double GradGen::getDp(const v3& v3GPrev, const v3& v3GThis, double dDt, double dPPrev, double dPThis, double dSignDp)
{
    // solve `ΔP` by RK2
    double dDl = v3::norm(v3GThis)*dDt;
    v3 v3DkDp0; m_ptTraj->getDkDp(&v3DkDp0, dPThis);
    v3 v3DkDp1; m_ptTraj->getDkDp(&v3DkDp1, dPThis*2e0-dPPrev);
    double dDlDp0 = v3::norm(v3DkDp0)*dSignDp;
    double dDlDp1 = v3::norm(v3DkDp1)*dSignDp;
    return dDl*(1e0/dDlDp0 + 1e0/dDlDp1)/2e0;
}

#endif

bool GradGen::step(v3* pv3GUnit, double* pdGNormMin, double* pdGNormMax, double dP, double dSignDp, const v3& v3G, double dSLim, double dDt)
{
    // current gradient direction
    v3 v3DkDp; m_ptTraj->getDkDp(&v3DkDp, dP);
    double dDlDp = v3::norm(v3DkDp)*dSignDp;
    if (pv3GUnit) *pv3GUnit = v3DkDp/dDlDp;
    
    // current gradient magnitude
    return sovQDE
    (
        pdGNormMin, pdGNormMax,
        1e0,
        -2e0*v3::inner(v3G, *pv3GUnit),
        v3::inner(v3G, v3G) - std::pow(dSLim*dDt, 2e0)
    );
}

bool GradGen::compute(lv3* plv3G, ld* pldP)
{
    bool bRet = true;
    double dP0 = m_ptTraj->getP0();
    double dP1 = m_ptTraj->getP1();

    int64_t lNit = 0;

    // backward
    v3 v3G1Unit; bRet &= m_ptTraj->getDkDp(&v3G1Unit, dP1);
    v3G1Unit = v3G1Unit * (dP0>dP1?1e0:-1e0);
    v3G1Unit = v3G1Unit / v3::norm(v3G1Unit);
    double dG1Norm = m_dG1Norm;
    dG1Norm = std::min(dG1Norm, m_dGLim);
    dG1Norm = std::min(dG1Norm, std::sqrt(m_dSLim*getCurRad(dP1)));
    v3 v3G1 = v3G1Unit * dG1Norm;

    ld ldP_Bac; ldP_Bac.push_back(dP1);
    lv3 lv3G_Bac; lv3G_Bac.push_back(v3G1);
    ld ldGNorm_Bac; ldGNorm_Bac.push_back(v3::norm(v3G1));
    while (1)
    {
        double dP = *ldP_Bac.rbegin();
        v3 v3G = *lv3G_Bac.rbegin();
        // update grad
        v3 v3GUnit;
        double dGNorm;
        bRet &= step(&v3GUnit, NULL, &dGNorm, dP, (dP0-dP1)/std::fabs(dP0-dP1), v3G, m_dSLim, m_dDt/m_lOs);
        dGNorm = std::min(dGNorm, m_dGLim);
        dGNorm = std::min(dGNorm, std::sqrt(m_dSLim*getCurRad(dP)));
        v3G = v3GUnit*dGNorm;

        // update para
        dP += getDp(*lv3G_Bac.rbegin(), v3G, m_dDt/m_lOs, *ldP_Bac.rbegin(), dP, (dP0-dP1)/std::fabs(dP0-dP1));
        // dP += getDp(v3G, m_dDt/m_lOs, dP, (dP0-dP1)/std::fabs(dP0-dP1));

        // stop or append
        if (std::fabs(*ldP_Bac.rbegin() - dP1) >= (1-1e-6)*std::fabs(dP0 - dP1))
        {
            break;
        }
        else
        {
            // printf("bac: dP = %lf\n", dP); // test
            ldP_Bac.push_back(dP);
            lv3G_Bac.push_back(v3G);
            ldGNorm_Bac.push_back(v3::norm(v3G));
        }
    }
    vd vdP_Bac(ldP_Bac.rbegin(), ldP_Bac.rend());
    vd vdGNorm_Bac(ldGNorm_Bac.rbegin(), ldGNorm_Bac.rend());
    
    lNit += ldP_Bac.size();

    // forward
    v3 v3G0Unit; bRet &= m_ptTraj->getDkDp(&v3G0Unit, dP0);
    v3G0Unit = v3G0Unit * (dP1>dP0?1e0:-1e0);
    v3G0Unit = v3G0Unit / v3::norm(v3G0Unit);
    double dG0Norm = m_dG0Norm;
    dG0Norm = std::min(dG0Norm, m_dGLim);
    dG0Norm = std::min(dG0Norm, std::sqrt(m_dSLim*getCurRad(dP0)));
    dG0Norm = std::min(dG0Norm, intp(dP0, vdP_Bac[0], vdGNorm_Bac[0], vdP_Bac[1], vdGNorm_Bac[1]));
    v3 v3G0 = v3G0Unit * dG0Norm;

    ld ldP; if (pldP==NULL) pldP = &ldP;
    pldP->clear(); pldP->push_back(dP0);
    plv3G->clear(); plv3G->push_back(v3G0);
    ld::reverse_iterator ildP_Bac = std::next(ldP_Bac.rbegin());
    ld::reverse_iterator ildGNorm_Bac = std::next(ldGNorm_Bac.rbegin());
    while (1)
    {
        double dP = *pldP->rbegin();
        v3 v3G = *plv3G->rbegin();

        // update grad
        v3 v3GUnit;
        double dGNorm;
        bRet &= step(&v3GUnit, NULL, &dGNorm, dP, (dP1-dP0)/std::fabs(dP1-dP0), v3G, m_dSLim, m_dDt/m_lOs);
        dGNorm = std::min(dGNorm, m_dGLim);
        dGNorm = std::min(dGNorm, std::sqrt(m_dSLim*getCurRad(dP)));

        // find index for interpolation
        while (std::fabs(dP-dP0) > std::fabs(*ildP_Bac-dP0))
        {
            if (std::next(ildP_Bac)!=ldP_Bac.rend())
            {
                ++ildP_Bac;
                ++ildGNorm_Bac;
            }
            else break;
        }

        // interpolation
        dGNorm = std::min(dGNorm, intp(dP, *std::prev(ildP_Bac), *std::prev(ildGNorm_Bac), *ildP_Bac, *ildGNorm_Bac));
        v3G = v3GUnit*dGNorm;

        // update para
        dP += getDp(*plv3G->rbegin(), v3G, m_dDt/m_lOs, *pldP->rbegin(), dP, (dP1-dP0)/std::fabs(dP1-dP0));
        // dP += getDp(v3G, m_dDt/m_lOs, dP, (dP1-dP0)/std::fabs(dP1-dP0));

        // stop or append
        if (std::fabs(*pldP->rbegin() - dP0) >= (1-1e-6)*std::fabs(dP1 - dP0) || dGNorm <= 0)
        {
            break;
        }
        else
        {
            // printf("for: dP = %lf\n", dP); // test
            pldP->push_back(dP);
            plv3G->push_back(v3G);
        }
    }
    v3G1 = (v3::norm(v3G1)!=0 ? v3G1/v3::norm(v3G1) : v3(0,0,0)) * std::min(v3::norm(v3G1), v3::norm(*plv3G->rbegin()));
    
    lNit += pldP->size();
    printf("MAG Nit: %ld\n", (int64_t)lNit);

    // deoversamp the para. vec.
    {
        ld::iterator ildP = pldP->begin();
        int64_t n = pldP->size();
        for (int64_t i = 0; i < n; ++i)
        {
            if(i%m_lOs!=m_lOs/2) ildP = pldP->erase(ildP);
            else ++ildP;
        }
    }

    // derive gradient
    {
        plv3G->clear();
        ld::iterator ildP = std::next(pldP->begin());
        int64_t n = pldP->size();
        for (int64_t i = 1; i < n; ++i)
        {
            v3 v3K1; bRet &= m_ptTraj->getK(&v3K1, *ildP);
            v3 v3K0; bRet &= m_ptTraj->getK(&v3K0, *std::prev(ildP));
            plv3G->push_back((v3K1 - v3K0)/m_dDt);
            ++ildP;
        }
    }
    pldP->pop_front();

    if (g_bExGEnd_MAG)
    {
        // add ramp gradient to satisfy desired Gstart and Gfinal
        lv3 lv3GRampFront; bRet &= GradGen::ramp_front(&lv3GRampFront, *plv3G->begin(), v3G0, m_dSLim, m_dDt);
        lv3 lv3GRampBack; bRet &= GradGen::ramp_back(&lv3GRampBack, *plv3G->rbegin(), v3G1*-1, m_dSLim, m_dDt);
        
        // corresponding parameter sequence
        for (int64_t i = 0; i < (int64_t)lv3GRampFront.size(); ++i) pldP->push_front(dP0);
        for (int64_t i = 0; i < (int64_t)lv3GRampBack.size(); ++i) pldP->push_back(dP1);

        // concate ramp gradient
        plv3G->splice(plv3G->begin(), lv3GRampFront);
        plv3G->splice(plv3G->end(), lv3GRampBack);
    }

    return bRet;
}

bool GradGen::ramp_front(lv3* plv3GRamp, const v3& v3G0, const v3& v3G0Des, double dSLim, double dDt)
{
    v3 v3Dg = v3G0Des - v3G0;
    v3 v3DgUnit = v3::norm(v3Dg)!=0 ? v3Dg/v3::norm(v3Dg) : v3(0,0,0);
    int64_t lNSamp = (int64_t)std::ceil(v3::norm(v3Dg)/(dSLim*dDt));

    // derive ramp gradient
    plv3GRamp->clear();
    for (int64_t i = 1; i < lNSamp; ++i)
    {
        plv3GRamp->push_front(v3G0 + v3DgUnit * (dSLim*dDt) * i);
    }
    if (lNSamp>0) plv3GRamp->push_front(v3G0Des);
    
    return true;
}

double GradGen::ramp_front(lv3* plv3GRamp, const v3& v3G0, const v3& v3G0Des, int64_t lNSamp, double dDt)
{
    v3 v3Dg = v3G0Des - v3G0;
    v3 v3DgUnit = v3::norm(v3Dg)!=0 ? v3Dg/v3::norm(v3Dg) : v3(0,0,0);
    double dSLim = v3::norm(v3Dg)/(lNSamp*dDt);

    // derive ramp gradient
    plv3GRamp->clear();
    for (int64_t i = 1; i < lNSamp; ++i)
    {
        plv3GRamp->push_front(v3G0 + v3DgUnit * (dSLim*dDt) * i);
    }
    if (lNSamp>0) plv3GRamp->push_front(v3G0Des);

    return dSLim;
}

bool GradGen::ramp_back(lv3* plv3GRamp, const v3& v3G1, const v3& v3G1Des, double dSLim, double dDt)
{
    v3 v3Dg = v3G1Des - v3G1;
    v3 v3DgUnit = v3::norm(v3Dg)!=0 ? v3Dg/v3::norm(v3Dg) : v3(0,0,0);
    int64_t lNSamp = (int64_t)std::ceil(v3::norm(v3Dg)/(dSLim*dDt));

    // derive ramp gradient
    plv3GRamp->clear();
    for (int64_t i = 1; i < lNSamp; ++i)
    {
        plv3GRamp->push_back(v3G1 + v3DgUnit * (dSLim*dDt) * i);
    }
    if (lNSamp>0) plv3GRamp->push_back(v3G1Des);
    
    return true;
}

double GradGen::ramp_back(lv3* plv3GRamp, const v3& v3G1, const v3& v3G1Des, int64_t lNSamp, double dDt)
{
    v3 v3Dg = v3G1Des - v3G1;
    v3 v3DgUnit = v3::norm(v3Dg)!=0 ? v3Dg/v3::norm(v3Dg) : v3(0,0,0);
    double dSLim = v3::norm(v3Dg)/(lNSamp*dDt);

    // derive ramp gradient
    plv3GRamp->clear();
    for (int64_t i = 1; i < lNSamp; ++i)
    {
        plv3GRamp->push_back(v3G1 + v3DgUnit * (dSLim*dDt) * i);
    }
    if (lNSamp>0) plv3GRamp->push_back(v3G1Des);
    
    return dSLim;
}

bool GradGen::revGrad(v3* pv3M0Dst, lv3* plv3Dst, const v3& v3M0Src, const lv3& lv3Src, double dDt)
{
    bool bRet = true;

    if(lv3Src.size() <= 1) bRet = false;

    // derive Total M0
    *pv3M0Dst = v3M0Src;
    lv3::const_iterator ilv3Src = lv3Src.begin();
    for(int64_t i = 0; i < (int64_t)lv3Src.size()-1; ++i)
    {
        *pv3M0Dst += (*ilv3Src + *std::next(ilv3Src))*dDt/2e0;
        ++ilv3Src;
    }

    // reverse gradient
    *plv3Dst = lv3(lv3Src.rbegin(), lv3Src.rend());
    lv3::iterator ilv3Dst = plv3Dst->begin();
    while (ilv3Dst != plv3Dst->end())
    {
        *ilv3Dst *= -1;
        ++ilv3Dst;
    }

    return bRet;
}

bool GradGen::calM0(v3* pv3M0, const lv3& lv3Grad, double dDt, const v3& v3GBegin, const v3& v3GEnd)
{
    *pv3M0 = v3(0,0,0);
    const v3* pv3Grad = &v3GBegin;
    lv3::const_iterator ilv3Grad = lv3Grad.begin();
    while (ilv3Grad != lv3Grad.end())
    {
        *pv3M0 += (*pv3Grad + *ilv3Grad)*dDt/2e0;
        pv3Grad = &*ilv3Grad;
        ++ilv3Grad;
    }
    *pv3M0 += (*pv3Grad + v3GEnd)*dDt/2e0;

    return true;
}