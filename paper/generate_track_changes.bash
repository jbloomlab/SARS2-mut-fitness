#!/bin/bash

# generate track changes LaTex

latexdiff \
    --exclude-safecmd=citep \
    --exclude-safecmd=cite \
    paper_bioRxiv_v2.tex \
    paper.tex \
    > paper_track_changes.tex
