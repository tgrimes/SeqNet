## Resubmission

SeqNet was archived on CRAN on 2025-04-23 ("issues were not corrected
despite reminders"). The last version on CRAN, 1.1.3, had two NOTEs in
its check results (no ERRORs or WARNINGs):

* checking Rd cross-references ... NOTE
  Found the following Rd file(s) with Rd \link{} targets missing package
  anchors: plot_modules.Rd, plot_network.Rd, plot_network_diff.Rd
  (all three came from `\pkg{\link{igraph}}` in the shared roxygen docs
  for the `generate_layout` argument; fixed to `\pkg{igraph}`, since
  `\link{}` inside `\pkg{}` doesn't need a target)

* checking dependencies in R code ... NOTE
  Namespace in Imports field not imported from: 'Rdpack'
  (Rdpack is required in Imports because of `RdMacros: Rdpack`, used
  for `@references` formatting, but wasn't otherwise referenced in R
  code; fixed by importing `Rdpack::reprompt`, which is Rdpack's own
  documented convention for this situation)

Both are fixed in this version (1.1.4). No other changes were made to
the package's functionality.

## Test environments

* local macOS 14.6.1 (aarch64-apple-darwin20), R 4.5.2, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  New submission, package was archived on CRAN

  This note is expected for any resubmission of an archived package
  and is not an issue with the package itself.

## Downstream dependencies

There are no downstream dependencies for this package.
