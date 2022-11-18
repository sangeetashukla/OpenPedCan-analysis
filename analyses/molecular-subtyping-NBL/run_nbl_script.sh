#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run notebook to subtype MYCN NBL 

Rscript -e "rmarkdown::render('00-subset-for-NBL.Rmd')"

Rscript -e "rmarkdown::render('01-find-matched-biospecimen.Rmd')"

Rscript -e "rmarkdown::render('02-find-non-matching-biospecimen.Rmd')"

Rscript -e "rmarkdown::render('03-subtyping.Rmd')"
