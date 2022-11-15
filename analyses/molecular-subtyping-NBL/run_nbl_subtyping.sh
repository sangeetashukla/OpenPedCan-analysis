#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run notebook to subtype MYCN NBL 

Rscript -e "rmarkdown::render('00-Preprocessing.Rmd')"

Rscript -e "rmarkdown::render('01-Subtyping.Rmd')"
