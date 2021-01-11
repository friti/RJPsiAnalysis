#!/bin/bash

set -ex

shopt -s expand_aliases
set +u && source ${CMS_PATH}/cmsset_default.sh; set -u
export SCRAM_ARCH=${SCRAM_ARCHITECTURE}
cmsrel ${CMSSW_RELEASE}
cd ${CMSSW_RELEASE}/src
cmsenv
git cms-init
git cms-merge-topic -u friti:TransientTracks
git cms-merge-topic -u friti:KinParticleVtxFitter
git clone git@github.com:gabriel-rmrz/RJPsiAnalysis.git ./PhysicsTools
git cms-addpkg PhysicsTools/NanoAOD
scram b
