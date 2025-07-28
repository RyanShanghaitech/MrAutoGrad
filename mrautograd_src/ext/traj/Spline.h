#pragma once

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>

class Spline
{
public:
    typedef std::vector<double> vd;

    Spline()
    {}

    Spline(const vd& vdX, const vd& vdY)
    {
        fit(vdX, vdY);
    }
        
    bool fit(const vd& vdX, const vd& vdY)
    {
        bool bRet = true;
        m_vdX = vd(vdX);
        m_vdY = vd(vdY);

        const int64_t lN = m_vdX.size();
        if (lN < 2 || m_vdY.size() != lN)
        {
            throw std::invalid_argument("Input vectors must have the same size and at least 2 points.");
            bRet &= false;
            return bRet;
        }

        vd vdH(lN - 1);
        for (int64_t i = 0; i < lN - 1; ++i)
            vdH[i] = m_vdX[i + 1] - m_vdX[i];

        // Step 1: Set up the tridiagonal system
        vd vdAlpha(lN, 0.0);
        for (int64_t i = 1; i < lN - 1; ++i)
            vdAlpha[i] = (3.0 / vdH[i]) * (m_vdY[i + 1] - m_vdY[i]) - (3.0 / vdH[i - 1]) * (m_vdY[i] - m_vdY[i - 1]);

        // Step 2: Solve tridiagonal system for c (second derivatives)
        vd vdL(lN, 1.0), vdMu(lN, 0.0), vdZ(lN, 0.0);
        m_vdC.resize(lN, 0.0);
        m_vdB.resize(lN - 1, 0.0);
        m_vdD.resize(lN - 1, 0.0);
        m_vdA = m_vdY;

        for (int64_t i = 1; i < lN - 1; ++i)
        {
            vdL[i] = 2.0 * (m_vdX[i + 1] - m_vdX[i - 1]) - vdH[i - 1] * vdMu[i - 1];
            vdMu[i] = vdH[i] / vdL[i];
            vdZ[i] = (vdAlpha[i] - vdH[i - 1] * vdZ[i - 1]) / vdL[i];
        }

        // Natural spline boundary conditions
        vdL[lN - 1] = 1.0;
        vdZ[lN - 1] = 0.0;
        m_vdC[lN - 1] = 0.0;

        // Back substitution
        for (int i = int(lN) - 2; i >= 0; --i)
        {
            m_vdC[i] = vdZ[i] - vdMu[i] * m_vdC[i + 1];
            m_vdB[i] = (m_vdA[i + 1] - m_vdA[i]) / vdH[i] - vdH[i] * (m_vdC[i + 1] + 2.0 * m_vdC[i]) / 3.0;
            m_vdD[i] = (m_vdC[i + 1] - m_vdC[i]) / (3.0 * vdH[i]);
        }
        
        return bRet;
    }

    double eval(double dXEval, int64_t lOrder=0) const // order: order of derivation, default is 0 (function value)
    {
        // Find the interval using bisection
        int64_t lN = m_vdX.size();
        if (lN < 2)
        {
            throw std::invalid_argument("Input vectors must have the same size and at least 2 points.");
        }

        int64_t low = 0;
        int64_t high = lN - 1;
        while (high - low > 1)
        {
            int64_t mid = (low + high) / 2;
            if (m_vdX[mid] > dXEval) high = mid;
            else             low = mid;
        }

        double dDx = dXEval - m_vdX[low];
        if (lOrder==0) return
        (
            m_vdA[low]
            + m_vdB[low] * dDx
            + m_vdC[low] * dDx*dDx
            + m_vdD[low] * dDx*dDx*dDx
        );
        if (lOrder==1) return
        (
            m_vdB[low]
            + m_vdC[low] * 2e0*dDx
            + m_vdD[low] * 3e0*dDx*dDx
        );
        if (lOrder==2) return
        (
            m_vdC[low] * 2e0
            + m_vdD[low] * 6e0*dDx
        );
        if (lOrder==3) return
        (
            m_vdD[low] * 6e0
        );
        return 0e0;
    }

private:
    vd m_vdX, m_vdY, m_vdA, m_vdB, m_vdC, m_vdD;
};