all: README.md README_solutions.md

README.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", output_file = "README.md", params = list(eval = FALSE))'

README_solutions.md: README.Rmd
	R -e 'rmarkdown::render("README.Rmd", output_file = "README_solutions.md", params = list(eval = TRUE))'
