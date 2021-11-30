
# nanoAOD producer for RJPsi analysis

Add your fork of the repository as remote and relevant packages:

```shell
export SCRAM_ARCH=slc7_amd64_gcc700
cmsrel CMSSW_10_6_14
cd CMSSW_10_6_14/src/
cmsenv
git cms-init
git cms-merge-topic -u friti:TransientTracks
git cms-merge-topic -u friti:KinParticleVtxFitter
git cms-merge-topic -u friti:GenParticlesPrecision
git clone git@github.com:friti/RJPsiAnalysis.git  ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
```
