#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run R script to subtype ATRT 

Rscript 00-ATRT_subtyping.R

