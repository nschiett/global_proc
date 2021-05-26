library(rmarkdown)
library(knitr)

render <- function(target){
  
  knitr::knit(input = paste0("text/", target, ".Rmd"), output =  paste0("text/", target, ".md"))

  pandoc_convert(paste0("text/", target, ".md"), output = "output/text/tmp.tex", wd = I('.'), options = I(c('--bibliography', 'text/model_proc.bib', '--standalone', '--template', 'text/template.latex',  '--pdf-engine', 'xelatex')))

  rearrangePaperAndCleanRefs('output/text/tmp.tex', paste0("output/text/", target, ".tex"))
  # build pdf
  latexBuildClean(paste0("output/text/", target, ".tex"), engine = 'xelatex')
  
  pandoc_convert(paste0("output/text/", target, ".tex"), output = paste0("output/text/", target, "Formatted.md"), wd = I('.'))
  
  #fixFigPathForDocxConversion( paste0("output/text/", target, "Formatted.md"), output =  paste0("output/text/", target, "FormattedWithFigPath.md"))
  
  pandoc_convert(paste0("output/text/", target, "Formatted.md"), output = paste0("output/text/", target, ".docx"), wd = I('.'), options = I('--resource-path=/'))

}

render("suppl_methods")

