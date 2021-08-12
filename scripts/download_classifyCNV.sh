#!/bin/bash

# install ClassifyCNV
# https://github.com/Genotek/ClassifyCNV

# desired version
version=$1

# download release
curl -LJO https://github.com/Genotek/ClassifyCNV/archive/v"${version}".zip &&

# unzip
unzip ClassifyCNV-"${version}".zip &&

# update ClinGen
cd ClassifyCNV-"${version}" || exit &&
bash update_clingen.sh &&
cd ..



