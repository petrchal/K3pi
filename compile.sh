#!/bin/bash

echo compiling in $STAR_VERSION
rsync -va $STAR/StRoot/StEvent ./StRoot/
rsync -va $STAR/StRoot/StMuDSTMaker ./StRoot/

cons

