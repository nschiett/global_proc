#################
# DOCX CONVERSION
#################
fixFigPathForDocxConversion  <-  function (inputFilePath, outputFilePath) {
  # get's rid of weird Rmd conversions that screw up pdf rendering
  all       <-  readLines(inputFilePath)
  figLines  <-  grep('![](../figures/fig', all, fixed = TRUE)
  writeLines(all[-figLines], outputFilePath)
}

##################
# LATEX CONVERSION
##################
rearrangePaperAndCleanRefs  <-  function (input, output) {
  # cleans some in-text citations, rearrange order of sections 
  all   <-  readLines(input)
  # number of refs and fix reference indentation
  refstag  <-  grep('\\hypertarget{ref-', all, fixed = TRUE)
  all[refstag + 1]   <-  paste0('\\noindent ', all[refstag + 1])
  writeLines(all, output)
}

###########################
# LATEX BUILD
# WRITTEN BY Rich FitzhJohn
###########################
latexBuildClean  <-  function (...) {
  latexBuild(..., clean = TRUE)
}

latexBuild  <-  function (filename, bibliography = NULL, chdir = TRUE, interaction = 'nonstopmode', maxAttempts = 5L, clean = FALSE, engine = 'pdflatex') {
  if (chdir && dirname(filename) != '') {
    owd  <-  setwd(dirname(filename))
    on.exit(setwd(owd))
    filename  <-  basename(filename)
  }
  
  res  <-  runLatex(filename, interaction, engine)
  if (engine == 'xelatex') {
    res <- runLatex(filename, interaction, engine)
  }
  if (!is.null(bibliography)) {
    runBibtex(filename)
    res <- runLatex(filename, interaction, engine)
  }
  
  pat <- c('Rerun to get cross-references right', 'Rerun to get citations correct',
           'Rerun to get outlines right')
  
  isin  <-  function (p, x) {
    any(grepl(p, x))
  }
  
  for (i in seq_len(maxAttempts)) {
    if (any(vapply(pat, isin, logical(1), res))) {
      res  <-  runLatex(filename, interaction, engine)
    } else {
      break
    }
  }
  
  if (clean) {
    latexClean(filename)
  }
  
  invisible(NULL)
}

latexClean <- function (filename) {
  filebase <- sub('.tex$', '', filename)
  exts  <-  c('.log', '.aux', '.bbl', '.blg', '.fls', '.out', '.snm', '.nav', '.tdo',
              '.toc')
  aux  <-  paste0(filebase, exts)
  file.remove(aux[file.exists(aux)])
}

runLatex  <-  function (filename, interaction = 'nonstopmode', engine = 'pdflatex') {
  args  <-  c(paste0('-interaction=', interaction), '-halt-on-error', filename)
  callSystem(SysWhich(engine), args)
}

runBibtex <- function (filename) {
  callSystem(SysWhich('bibtex'), sub('.tex$', '', filename))
}

SysWhich <- function (x) {
  ret  <-  Sys.which(x)
  if (ret == '') {
    stop(sprintf('%s not found in $PATH', x))
  }
  ret
}
##' Function imported from callr package;  makes it easy to call a
##' system command from R and have it behave.
##'
##' This function uses \code{system2} to call a system command fairly
##' portably.  What it adds is a particular way of dealing with
##' errors.  \code{callSystem} runs the command \code{command} with
##' arguments \code{args} (and with optionally set environment
##' variables \code{env}) and hides \emph{all} produced output to
##' stdout and stderr.  If the command fails (currently any nonzero
##' exit code is counted as a failure) then \code{callSystem} will
##' throw an R error giving
##' \itemize{
##' \item the full string of the command run
##' \item the exit code of the command
##' \item any \code{errmsg} attribute that might have been returned
##' \item all output that the program produced to either stdout and
##' stderr
##' }
##'
##' This means that a successful invocation of a program produces no
##' output while the unsuccessful invocation throws an error and
##' prints all information to the screen (though this is delayed until
##' failure happens).
##'
##'
##' \code{callSystem} also returns the contents of both stderr and
##' stdout \emph{invisibly} so that it can be inspected if needed.
##'
##' @title Run a system command, stopping on error
##' @param command The system command to be invoked, as a character
##' string.  \code{\link{Sys.which}} is useful here.
##' @param args A character vector of arguments to \code{command}
##' @param env A character vector of name=value pairs to be set as
##' environment variables (see \code{\link{system2}}).
##' @param maxLines Maximum number of lines of program output to
##' print with the error message.  We may prune further to get the
##' error message under \code{getOption("warn.length")}, however.
##' @param p Fraction of the error message to show from the tail of
##' the output if truncating on error (default is 20\% lines are head,
##' 80\% is tail).
##' @param stdout,stderr Passed to \code{system2}.  Set one of these
##' to \code{FALSE} to avoid capturing output from that stream.  Setting
##' both to \code{FALSE} is not recommended.
##' @export
##' @author Rich FitzJohn
callSystem  <-  function (command, args, env = character(), maxLines = 20,
                          p = 0.8, stdout = TRUE, stderr = TRUE) {
  res  <-  suppressWarnings(system2(command, args, env = env, stdout = stdout, stderr = stderr))
  ok   <-  attr(res, 'status')
  if (!is.null(ok) && ok != 0) {
    maxNc  <-  getOption('warning.length')
    
    cmd     <-  paste(c(env, shQuote(command), args), collapse = ' ')
    msg     <-  sprintf('Running command:\n  %s\nhad status %d', cmd, ok)
    errmsg  <-  attr(cmd, 'errmsg')
    if (!is.null(errmsg)) {
      msg  <-  c(msg, sprintf('%s\nerrmsg: %s', errmsg))
    }
    sep  <-  paste(rep('-', getOption('width')), collapse = '')
    
    ## Truncate message:
    if (length(res) > maxLines) {
      n <- ceiling(maxLines * p)
      res <- c(head(res, ceiling(maxLines - n)),
               sprintf('[[... %d lines dropped ...]]', length(res) - maxLines),
               tail(res, ceiling(n)))
    }
    
    ## compute the number of characters so far, including three new lines:
    nc   <-  (nchar(msg) + nchar(sep) * 2) + 3
    i    <-  max(1, which(cumsum(rev(nchar(res) + 1L)) < (maxNc - nc)))
    res  <-  res[(length(res) - i + 1L):length(res)]
    msg  <-  c(msg, 'Program output:', sep, res, sep)
    stop(paste(msg, collapse = '\n'))
  }
  invisible(res)
}
