#!/bin/bash
export MINC_PATH=${MINC_TOOLKIT}

$R CMD INSTALL --build . --configure-args="--with-build-path=${MINC_TOOLKIT}"
