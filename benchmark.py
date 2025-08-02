import mrautograd as mag
from numpy import *
from matplotlib.pyplot import *
from numpy.linalg import norm

# select solver
mag.setSolverMtg(1) # 0 for MAG solver, 1 for MTG solver

# arguments
gamma = 42.5756e6
fov = 0.256
nPix = 256
nAx = -1
sLim = 100 * gamma * fov/nPix
gLim = 120e-3 * gamma * fov/nPix
dtGrad = 10e-6
dtADC = 2.5e-6
argCom = dict(dFov=fov, lNPix=nPix, dSLim=sLim, dGLim=gLim, dDt=dtGrad)

def eval():
    global nAx
    # nAx = 2; return mag.getG_Spiral(bIs3D=0, **argCom)
    # nAx = 2; return mag.getG_VarDenSpiral(bIs3D=0, **argCom)
    nAx = 2; return mag.getG_Rosette(bIs3D=0, **argCom)
    # nAx = 2; return mag.getG_Rosette_Trad(**argCom, dOm1=10*pi, dOm2=8*pi, dTmax=1, dTacq=2e-03)
    # nAx = 3; return mag.getG_Shell3d(**argCom)
    # nAx = 3; return mag.getG_Yarnball(**argCom)
    # nAx = 3; return mag.getG_Seiffert(**argCom)
    # nAx = 3; return mag.getG_Cones(**argCom)

mag.setTrajRev(0) # don't reverse trajectory
mag.setGoldAng(1) # use golden angle interleaving if possible
mag.setShuf(0) # disable shuffle
mag.setMaxG0(0) # don't maximize G0
mag.setMaxG1(0) # don't maximize G1
mag.setExGEnd(0) # don't ensure exact G0 and G1 for MAG, for fair comparison
mag.setMagOs(8) # set temporal oversampling of MAG to 8

for i in range(10):
    lstArrK0, lstArrGrad = eval()

# derive shape parameter
if nAx==2:
    lstArrK0 = [arrK0[:2] for arrK0 in lstArrK0]
    lstArrGrad = [arrG[:,:2] for arrG in lstArrGrad]
nRO, nAx = lstArrGrad[0].shape

# derive slewrate
lstArrSlew = [diff(arrG, axis=0)/dtGrad for arrG in lstArrGrad]
sMax = max(norm(concatenate(lstArrSlew)/gamma*nPix/fov,axis=-1))
gMax = max(norm(concatenate(lstArrGrad)/gamma*nPix/fov,axis=-1))

# derive trajectory
lstArrK = []
for arrK0, arrGrad in zip(lstArrK0, lstArrGrad):
    arrK, _ = mag.cvtGrad2Traj(arrGrad, dtGrad, dtADC)
    arrK += arrK0
    lstArrK.append(arrK)

# k-space and g-space
iArrK = len(lstArrK)*2//3

figure(figsize=(18,9), dpi=120)

subplot(221, projection="3d" if nAx==3 else None)
plot(*lstArrK[iArrK].T, ".-")
axis("equal")
grid("on")
title(f"kspace {1}/{len(lstArrGrad)}")

subplot(223, projection="3d" if nAx==3 else None)
plot(*lstArrGrad[iArrK].T, ".-")
axis("equal")
grid("on")
title("gspace")

# gradient and slewrate
subplot(222)
for iAx in range(nAx):
    plot(lstArrGrad[iArrK][:,iAx]/gamma*nPix/fov, ".-")
grid("on")
title(f"Gradient")

subplot(224)
plot(norm(lstArrSlew[iArrK],axis=-1)/gamma*nPix/fov, ".-", c="tab:blue")
ylim(sLim/gamma*nPix/fov*0.9, sLim/gamma*nPix/fov*1.1)
grid("on")

twinx()
plot(norm(lstArrGrad[iArrK],axis=-1)/gamma*nPix/fov*1e3, ".-", c="tab:orange")
ylim(gLim/gamma*nPix/fov*0.9*1e3, gLim/gamma*nPix/fov*1.1*1e3)
grid("on")
title(f"Grad & Slew amp., max grad:{gMax*1e3:.3f}, max slew:{sMax:.3f}")

subplots_adjust(0.05,0.1,0.95,0.9, 0.2, 0.2)

savefig("figure.pdf")
show()
