README.md: README.Rmd
	Rscript --vanilla --verbose -e 'rmarkdown::render("README.Rmd")'

docs:
	Rscript -e 'devtools::document()'

check:
	Rscript -e 'devtools::check()'

install:
	Rscript -e 'devtools::install()'
