#!/bin/bash

set -e
set -o pipefail

# set up running directory
cd "$(dirname "${BASH_SOURCE[0]}")" 

# Run notebook to subtype MYCN NBL 
Rscript 00-select-diagnoses.R

Rscript -e "rmarkdown::render('01-subset-for-NBL.Rmd')"

Rscript -e "rmarkdown::render('02-find-matched-biospecimen.Rmd')"

Rscript -e "rmarkdown::render('03-find-non-matching-biospecimen.Rmd')"

Rscript -e "rmarkdown::render('04-subtyping.Rmd')"

Rscript -e "rmarkdown::render('05-qc-checks.Rmd')"

echo 'NBL Subtyping Completed'

