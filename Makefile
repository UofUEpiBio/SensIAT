README.md: README.Rmd
	Rscript --vanilla --verbose -e 'rmarkdown::render("README.Rmd")'
