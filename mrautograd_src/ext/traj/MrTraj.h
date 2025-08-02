#pragma once

#include "TrajFunc.h"
#include "../mag/GradGen.h"
#include <string>
#include <stdexcept>
#include <ctime>
#include "Spline.h"

#ifdef MTG_EXIST
    #include "../mtg/header.h"
#endif

bool g_bUseMtg_MrTraj = false; // use Lustig's MinTimeGrad solver
int64_t g_lOsOw_MrTraj = -1; // oversample ratio, overwrite set value

#define GOLDRAT ((1e0+std::sqrt(5e0))/2e0)
#define GOLDANG ((3e0-std::sqrt(5e0))*M_PI)


class Spline_TrajFunc: public TrajFunc
{
public:
    typedef std::vector<v3> vv3;
    typedef std::vector<double> vd;

    Spline_TrajFunc(const vv3& vv3K)
    {
        int64_t lNTrajSamp = vv3K.size();

        vd vdP(lNTrajSamp);
        for (int64_t i = 0; i < lNTrajSamp; ++i)
        {
            vdP[i] = i;
        }

        vd vdX(lNTrajSamp), vdY(lNTrajSamp), vdZ(lNTrajSamp);
        for (int64_t i = 0; i < lNTrajSamp; ++i)
        {
           vdX[i] =  vv3K[i].m_dX;
           vdY[i] =  vv3K[i].m_dY;
           vdZ[i] =  vv3K[i].m_dZ;
        }

        m_intpX.fit(vdP, vdX);
        m_intpY.fit(vdP, vdY);
        m_intpZ.fit(vdP, vdZ);

        m_dP0 = *vdP.begin();
        m_dP1 = *vdP.rbegin();
    }
    
    bool getK(v3* pv3K, double dP) const
    {
        pv3K->m_dX = m_intpX.eval(dP);
        pv3K->m_dY = m_intpY.eval(dP);
        pv3K->m_dZ = m_intpZ.eval(dP);

        return true;
    }

    bool getDkDp(v3* pv3K, double dP) const
    {
        pv3K->m_dX = m_intpX.eval(dP, 1);
        pv3K->m_dY = m_intpY.eval(dP, 1);
        pv3K->m_dZ = m_intpZ.eval(dP, 1);

        return true;
    }

    bool getD2kDp2(v3* pv3K, double dP) const
    {
        pv3K->m_dX = m_intpX.eval(dP, 2);
        pv3K->m_dY = m_intpY.eval(dP, 2);
        pv3K->m_dZ = m_intpZ.eval(dP, 2);

        return true;
    }
protected:
    Spline m_intpX, m_intpY, m_intpZ;
};

/* 
 * A set of trajectories sufficient to fully-sample the k-space
 * defined by:
 * 1. some Base trajectories with different shapes.
 * 2. a acquisition plan which decides which Base trajectory to use,
 *    and how to transform to the desired gradient of a particular acquisition.
 * 
 * notice:
 * 1. Different Base trajectoires share the same traj. func. getK(),
 *    but the behaviour of getK() may differ due to constant parameters.
 */

class MrTraj
{
public:
    typedef std::vector<int64_t> vl;
    typedef std::list<int64_t> ll;
    typedef std::vector<double> vd;
    typedef std::list<double> ld;
    typedef std::string str;
    typedef std::vector<v3> vv3;
    typedef std::list<v3> lv3;
    typedef std::vector<vv3> vvv3;
    typedef std::vector<lv3> vlv3;
    typedef std::list<vv3> lvv3;
    typedef std::vector<TrajFunc*> vptf;
    typedef struct
    {
        bool bIs3D;
        double dFov;
        int64_t lNPix;
    } GeoPara;
    typedef struct
    {
        double dSLim;
        double dGLim;
        double dDt;
        bool bMaxG0;
        bool bMaxG1;
    } GradPara;
    
    const double dGyoMagRat; // Hz/T
    
    MrTraj():
        dGyoMagRat(42.5756e6)
    {}
    
    virtual ~MrTraj()
    {}
    
    virtual bool getM0PE(v3* pv3M0PE, int64_t lIAcq) const = 0;
    
    virtual bool getGRO(lv3* plv3GRO, int64_t lIAcq) const = 0;
    
    virtual int64_t getNWait(int64_t lIAcq) const = 0;
    
    virtual int64_t getNSamp(int64_t lIAcq) const = 0;

    const GeoPara& getGeoPara() const
    { return m_sGeoPara; }

    const GradPara& getGradPara() const
    { return m_sGradPara; }
    
    int64_t getNAcq() const
    { return m_lNAcq; }

    static bool genRandIdx(vl* pvlIdx, int64_t lN)
    {
        ll llSeq;
        for(int64_t i = 0; i < lN; ++i)
        {
            llSeq.push_back(i);
        }
        ll::iterator illSeq = llSeq.begin();
        pvlIdx->resize(lN);
        
        int64_t lIntv = (int64_t)round(llSeq.size()*(llSeq.size()+1)*(1e0/(GOLDRAT+1e0))); lIntv %= std::max(llSeq.size(), (size_t)1);
        for(int64_t i = 0; i < lN; ++i)
        {
            for(int64_t j = 0; j < lIntv; ++j)
            {
                illSeq++;
                if(illSeq==llSeq.end()) illSeq=llSeq.begin();
            }
            pvlIdx->at(i) = *illSeq;
            illSeq = llSeq.erase(illSeq);
            if(illSeq==llSeq.end()) illSeq=llSeq.begin();
            lIntv = (int64_t)round(llSeq.size()*(llSeq.size()+1)*(1e0/(GOLDRAT+1e0))); lIntv %= std::max(llSeq.size(), (size_t)1);
        }

        return true;
    }

    static bool calM0SP(v3* pv3M0SP, const v3& v3M0PE, const lv3 lv3GRO)
    {
        bool bRet = true;
        *pv3M0SP = v3M0PE;
        lv3::const_iterator ilv3GRO = lv3GRO.begin();
        lv3::const_iterator ilv3NonZero = ilv3GRO;
        for(int64_t i = 0; i < (int64_t)lv3GRO.size()-1; ++i)
        {
            *pv3M0SP += (*ilv3GRO + *std::next(ilv3GRO))*1e0/2e0; // assume `dt` to be 1
            if (v3::norm(*ilv3GRO)!=0e0) ilv3NonZero = ilv3GRO;
            ++ilv3GRO;
        }

        if (v3::norm(*pv3M0SP) != 0e0)
        {
            *pv3M0SP /= v3::norm(*pv3M0SP);
        }
        else if (v3::norm(*ilv3NonZero) != 0e0)
        {
            *pv3M0SP = *ilv3NonZero/v3::norm(*ilv3NonZero);
        }
        else
        {
            *pv3M0SP = v3(1,0,0);
            bRet = false;
        }

        return bRet;
    }

protected:
    GeoPara m_sGeoPara;
    GradPara m_sGradPara;
    int64_t m_lNAcq;
    
    // calculate required num. of rot. to satisfy Nyquist sampling (for spiral only)
    static int64_t calNRot(double dRhoRotAng, int64_t lNPix)
    {
        return (int64_t)std::ceil(lNPix*2e0*M_PI*dRhoRotAng);
    }

    // calculate required num. of rot. to satisfy Nyquist sampling
    static int64_t calNRot(TrajFunc* ptraj, double dP0, double dP1, int64_t lNPix, int64_t lNSamp=1000)
    {
        /*
         * Note:
         * This method is base on the derivative of trajectory function,
         * only local sampling is considered, so there are limitations for
         * this method
         * 
         * Applicable Trajectories:
         * Spiral, Cones, Rosette (single petal)
         * 
         * Non-Applicable Trajectories:
         * Yarnball, Rosette (multi petal)
         */

        // calculate and find min. rot. ang.
        double dNyqIntv = 1e0/lNPix;
        double dMinRotAng = 2e0*M_PI;
        for (int64_t lIk = 1; lIk < lNSamp-1; ++lIk)
        {
            double dP = dP0 + ((double)lIk/lNSamp)*(dP1-dP0);
            v3 v3K; ptraj->getK(&v3K, dP);
            v3 v3DkDpara;
            {
                v3 v3K_Nx; ptraj->getK(&v3K_Nx, dP+1e-7);
                v3DkDpara = v3K_Nx - v3K;
            }
            if (v3::norm(v3DkDpara)==0) continue;
            v3 v3DkDphi;
            {
                v3 v3K_Nx; v3::rotate(&v3K_Nx, 2, 1e-7, v3K);
                v3DkDphi = v3K_Nx - v3K;
            }
            if (v3::norm(v3DkDphi)==0) continue;
            double dRho = std::sqrt(v3K.m_dX*v3K.m_dX + v3K.m_dY*v3K.m_dY);
            double dCos = v3::inner(v3DkDpara, v3DkDphi) / (v3::norm(v3DkDpara) * v3::norm(v3DkDphi));
            double dSin = std::sqrt(1e0 - std::min(dCos*dCos,1e0));
            double dRotAng = (dNyqIntv/dSin) / (dRho);
            dMinRotAng = std::min(dMinRotAng, dRotAng);
        }

        // ensure the rot. Num. is a integer
        return (int64_t)std::ceil(2e0*M_PI/dMinRotAng);
    }
    
    static double calRotAngInc(int64_t lNRot)
    {
        return 2e0*M_PI/lNRot;
    }

    static bool calGRO_MAG(lv3* plv3G, ld* pldP, const TrajFunc& tf, const GradPara& sGradPara, int64_t lOs=8)
    {
        bool bRet = true;
        const double& dSLim = sGradPara.dSLim;
        const double& dGLim = sGradPara.dGLim;
        const double& dDt = sGradPara.dDt;
        const bool& bMaxG0 = sGradPara.bMaxG0;
        const bool& bMaxG1 = sGradPara.bMaxG1;

        GradGen gg(&tf, dSLim, dGLim, dDt, lOs, bMaxG0?1e15:0e0, bMaxG1?1e15:0e0);
        bRet &= gg.compute(plv3G, pldP);

        return bRet;
    }

    static bool calGRO_MTG(lv3* plv3G, ld* pldP, const vd& vdC, const GradPara& sGradPara)
    {
#ifdef MTG_EXIST
        bool bRet = true;
        const double& dSLim = sGradPara.dSLim;
        const double& dGLim = sGradPara.dGLim;
        const double& dDt = sGradPara.dDt;
        const bool& bMaxG0 = sGradPara.bMaxG0;
        const bool& bMaxG1 = sGradPara.bMaxG1;

        // Prepare arg. for Lustig's function
        double g0 = bMaxG0?1e15:0e0, gfin = bMaxG1?1e15:0e0, gmax = dGLim, smax = dSLim, T = dDt, ds = -1;

        double *p_Cx = nullptr, *p_Cy = nullptr, *p_Cz = nullptr;
        double *p_gx = nullptr, *p_gy = nullptr, *p_gz = nullptr;
        double *p_p = nullptr;
        double *p_sx = nullptr, *p_sy = nullptr, *p_sz = nullptr;
        double *p_kx = nullptr, *p_ky = nullptr, *p_kz = nullptr;
        double *p_sdot = nullptr, *p_sta = nullptr, *p_stb = nullptr;
        
        double time = 0;
        int size_interpolated = 0, size_sdot = 0, size_st = 0;
        int gfin_empty = 0, ds_empty = 1;

        // Call Lustig's function (assume it is linked in or compiled as C)
        minTimeGradientRIV(
            vdC.data(), vdC.size()/3, 3, g0, gfin, gmax, smax, T, ds,
            &p_Cx, &p_Cy, &p_Cz, &p_gx, &p_gy, &p_gz, &p_p,
            &p_sx, &p_sy, &p_sz, &p_kx, &p_ky, &p_kz, &p_sdot, &p_sta, &p_stb, &time,
            &size_interpolated, &size_sdot, &size_st, gfin_empty, ds_empty);

        // Copy results to C++ outputs
        plv3G->clear();
        for (int i = 0; i < size_interpolated; ++i)
        {
            plv3G->push_back(v3(p_gx[i], p_gy[i], p_gz[i]));
        }
        if (pldP)
        {
            pldP->clear();
            for (int i = 0; i < size_interpolated; ++i)
            {
                pldP->push_back(p_p[i]);
            }
        }
        
        free(p_Cx);    free(p_Cy);    free(p_Cz);
        free(p_gx);    free(p_gy);    free(p_gz);
        free(p_p);
        free(p_sx);    free(p_sy);    free(p_sz);
        free(p_kx);    free(p_ky);    free(p_kz);
        free(p_sdot);  free(p_sta);   free(p_stb);

        return bRet;
#else
        char sErrMsg[] = "[ERROR] MTG removed due to copyright.";
        puts(sErrMsg);
        throw std::runtime_error(sErrMsg);
        return false;
#endif
    }

    static bool calGRO(lv3* plv3G, ld* pldP, const TrajFunc& tf, const GradPara& sGradPara, int64_t lOs=8)
    {
        bool bRet = true;

        // calculate gradient
        if(!g_bUseMtg_MrTraj)
        {
            bRet &= calGRO_MAG(plv3G, pldP, tf, sGradPara, lOs);
        }
        else
        {
            // Prepare trajectory sampling
            int64_t lNTrajSamp = 1000;
            vd vdC(lNTrajSamp * 3, 0.0);

            // Sample the trajectory at N points
            double dP0 = tf.getP0();
            double dP1 = tf.getP1();
            for (int i = 0; i < lNTrajSamp; ++i)
            {
                double dP = dP0 + (dP1-dP0)* (i)/double(lNTrajSamp-1);
                v3 v3K; tf.getK(&v3K, dP);
                v3K *= 4.257; // k is defined by k*4.257 in Lustig's method
                vdC[i] = v3K.m_dX;
                vdC[i + lNTrajSamp] = v3K.m_dY;
                vdC[i + 2*lNTrajSamp] = v3K.m_dZ;
            }

            bRet &= calGRO_MTG(plv3G, pldP, vdC, sGradPara);
        }

        return bRet;
    }

    static bool calGRO(lv3* plv3G, ld* pldP, const vv3& vv3TrajSamp, const GradPara& sGradPara, int64_t lOs=8)
    {
        bool bRet = true;
        int64_t lNTrajSamp = vv3TrajSamp.size();

        // calculate gradient
        if(!g_bUseMtg_MrTraj)
        {
            Spline_TrajFunc sptf(vv3TrajSamp);
            bRet &= calGRO_MAG(plv3G, pldP, sptf, sGradPara, lOs);
        }
        else
        {
            // Prepare trajectory sampling
            vd vdC(lNTrajSamp*3);

            // Sample the trajectory at N points
            for (int i = 0; i < lNTrajSamp; ++i)
            {
                v3 v3K = vv3TrajSamp[i]*4.257; // k is defined by k*4.257 in Lustig's method
                vdC[i] = v3K.m_dX;
                vdC[i + lNTrajSamp] = v3K.m_dY;
                vdC[i + 2*lNTrajSamp] = v3K.m_dZ;
            }

            bRet &= calGRO_MTG(plv3G, pldP, vdC, sGradPara);
        }

        return bRet;
    }
    
    static bool calGrad(v3* pv3M0PE, lv3* plv3GRO, ld* pldP, int64_t* plNWait, int64_t* plNSamp, const TrajFunc& tfBaseTraj, const GradPara& sGradPara, int64_t lOs=8)
    {
        bool bRet = true;
        const double& dGLim = sGradPara.dGLim;
        const double& dSLim = sGradPara.dSLim;
        const double& dDt = sGradPara.dDt;
        double dTRampFront = sGradPara.bMaxG0 ? dGLim/dSLim : 0e0;
        double dTRampBack = sGradPara.bMaxG1 ? dGLim/dSLim : 0e0;
        if (g_lOsOw_MrTraj>0) lOs = g_lOsOw_MrTraj;
        
        // calculate GRO with ramp-up and ramp-down
        calGRO(plv3GRO, pldP, tfBaseTraj, sGradPara, lOs);

        lv3 lv3GRampFront; bRet &= dSLim>=GradGen::ramp_front(&lv3GRampFront, *plv3GRO->begin(), v3(0,0,0), int64_t(dTRampFront/dDt), dDt);
        lv3 lv3GRampBack; bRet &= dSLim>=GradGen::ramp_back(&lv3GRampBack, *plv3GRO->rbegin(), v3(0,0,0), int64_t(dTRampBack/dDt), dDt);

        *plNWait = lv3GRampFront.size();
        *plNSamp = plv3GRO->size();

        // calculate M0 of PE
        bRet &= tfBaseTraj.getK0(pv3M0PE);
        v3 v3M0Ramp; GradGen::calM0(&v3M0Ramp, lv3GRampFront, dDt, v3(0,0,0), *plv3GRO->begin());
        *pv3M0PE -= v3M0Ramp;

        // concate ramp gradient
        plv3GRO->splice(plv3GRO->begin(), lv3GRampFront);
        plv3GRO->splice(plv3GRO->end(), lv3GRampBack);
        
        return bRet;
    }
    
    static bool calGrad(v3* pv3M0PE, lv3* plv3GRO, ld* pldP, int64_t* plNWait, int64_t* plNSamp, const vv3& vv3TrajSamp, const GradPara& sGradPara, int64_t lOs=8)
    {
        bool bRet = true;
        const double& dGLim = sGradPara.dGLim;
        const double& dSLim = sGradPara.dSLim;
        const double& dDt = sGradPara.dDt;
        double dTRampFront = sGradPara.bMaxG0 ? dGLim/dSLim : 0e0;
        double dTRampBack = sGradPara.bMaxG1 ? dGLim/dSLim : 0e0;
        if (g_lOsOw_MrTraj>0) lOs = g_lOsOw_MrTraj;
        
        // calculate GRO with ramp-up and ramp-down
        calGRO(plv3GRO, pldP, vv3TrajSamp, sGradPara, lOs);

        lv3 lv3GRampFront; bRet &= dSLim>=GradGen::ramp_front(&lv3GRampFront, *plv3GRO->begin(), v3(0,0,0), int64_t(dTRampFront/dDt), dDt);
        lv3 lv3GRampBack; bRet &= dSLim>=GradGen::ramp_back(&lv3GRampBack, *plv3GRO->rbegin(), v3(0,0,0), int64_t(dTRampBack/dDt), dDt);

        *plNWait = lv3GRampFront.size();
        *plNSamp = plv3GRO->size();

        // calculate M0 of PE
        *pv3M0PE = *vv3TrajSamp.begin();
        v3 v3M0Ramp; GradGen::calM0(&v3M0Ramp, lv3GRampFront, dDt, v3(0,0,0), *plv3GRO->begin());
        *pv3M0PE -= v3M0Ramp;

        // concate ramp gradient
        plv3GRO->splice(plv3GRO->begin(), lv3GRampFront);
        plv3GRO->splice(plv3GRO->end(), lv3GRampBack);
        
        return bRet;
    }
};
