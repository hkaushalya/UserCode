import FWCore.ParameterSet.Config as cms
from UserCode.MetSigForQcd.metsig_cfi import *

usePrescaleWeight = 0
applyLumiweighing = 0
applyEventweighing = 0


njetmin = 3
njetmax = 1000
verbose = False
njetCut = 1
dPhiCut = 1
dMinMetSig = 0

print 'metsigall_cfi: usePrescaleWeight  is set to = ', usePrescaleWeight
print 'metsigall_cfi: applyLumiweighing  is set to = ', applyLumiweighing
print 'metsigall_cfi: applyEventweighing is set to = ', applyEventweighing
print 'metsigall_cfi: njetmin = ', njetmin
print 'metsigall_cfi: njetmax = ', njetmax
print 'metsigall_cfi: njetcut = ', njetCut
print 'metsigall_cfi: dphicut = ', dPhiCut


metsignomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsigmht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)

metsigmht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)

metsigmht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing),
                      verbose = cms.untracked.int32(verbose)

)

#metsig variations

dMinMetSig25 = 25.0

metsig25nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsig25mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsig25mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig25mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig25),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


dMinMetSig50 = 50.0

metsig50nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig50mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig50mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig50mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig50mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig50),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

dMinMetSig75 = 75.0

metsig75nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsig75mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig75mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig75mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsig75mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig75),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)



dMinMetSig100 = 100.0

metsig100nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsig100mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig100mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig100mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig100),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


dMinMetSig200 = 200

metsig200nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig200mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig200),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


dMinMetSig300 = 300

metsig300nomht = metsig.clone(
                      dMinMHT = cms.untracked.double(0.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig300mht25 = metsig.clone(
                      dMinMHT = cms.untracked.double(25.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht50 = metsig.clone(
                      dMinMHT = cms.untracked.double(50.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht75 = metsig.clone(
                      dMinMHT = cms.untracked.double(75.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht100 = metsig.clone(
                      dMinMHT = cms.untracked.double(100.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht125 = metsig.clone(
                      dMinMHT = cms.untracked.double(125.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht150 = metsig.clone(
                      dMinMHT = cms.untracked.double(150.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)
metsig300mht175 = metsig.clone(
                      dMinMHT = cms.untracked.double(175.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig300mht200 = metsig.clone(
                      dMinMHT = cms.untracked.double(200.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig300mht350 = metsig.clone(
                      dMinMHT = cms.untracked.double(350.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)

metsig300mht500 = metsig.clone(
                      dMinMHT = cms.untracked.double(500.0),   
                      dMinMetSig = cms.untracked.double(dMinMetSig300),   
							 dMinNjet = cms.untracked.double(njetmin),
							 dMaxNjet = cms.untracked.double(njetmax),
							 applyNjetCut = cms.untracked.int32(njetCut),
							 applyDphiCut = cms.untracked.int32(dPhiCut),
							 usePrescaleWeight = cms.untracked.int32(usePrescaleWeight),
                      ApplyLumiWeighing = cms.untracked.int32(applyLumiweighing),
                      ApplyEventWeighing = cms.untracked.int32(applyEventweighing)
)


metsigseq = cms.Sequence(
							 metsignomht +
							 metsigmht25 +
							 metsigmht50 +
							 metsigmht75 +
							 metsigmht100 +
							 metsigmht125 +
							 metsigmht150 +
							 metsigmht175 +
							 metsigmht200 +
							 metsigmht350 +
							 metsigmht500 +
							 metsig25nomht +
							 metsig25mht25 +
							 metsig25mht50 +
							 metsig25mht75 +
							 metsig25mht100 +
							 metsig25mht125 +
							 metsig25mht150 +
							 metsig25mht175 +
							 metsig25mht200 +
							 metsig25mht350 +
							 metsig25mht500 +
							 metsig50nomht +
							 metsig50mht25 +
							 metsig50mht50 +
							 metsig50mht75 +
							 metsig50mht100 +
							 metsig50mht125 +
							 metsig50mht150 +
							 metsig50mht175 +
							 metsig50mht200 +
							 metsig50mht350 +
							 metsig50mht500 +
							 metsig75nomht +
							 metsig75mht25 +
							 metsig75mht50 +
							 metsig75mht75 +
							 metsig75mht100 +
							 metsig75mht125 +
							 metsig75mht150 +
							 metsig75mht175 +
							 metsig75mht200 +
							 metsig75mht350 +
							 metsig75mht500 +
							 metsig100nomht +
							 metsig100mht25 +
							 metsig100mht50 +
							 metsig100mht75 +
							 metsig100mht100 +
							 metsig100mht125 +
							 metsig100mht150 +
							 metsig100mht175 +
							 metsig100mht200 +
							 metsig100mht350 +
							 metsig100mht500 +
							 metsig200nomht +
							 metsig200mht25 +
							 metsig200mht50 +
							 metsig200mht75 +
							 metsig200mht100 +
							 metsig200mht125 +
							 metsig200mht150 +
							 metsig200mht175 +
							 metsig200mht200 +
							 metsig200mht350 +
							 metsig200mht500 +
							 metsig300nomht +
							 metsig300mht25 +
							 metsig300mht50 +
							 metsig300mht75 +
							 metsig300mht100 +
							 metsig300mht150 +
							 metsig300mht125 +
							 metsig300mht175 +
							 metsig300mht200 +
							 metsig300mht350 +
							 metsig300mht500
)

