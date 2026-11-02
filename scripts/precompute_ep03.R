# Precompute the output of the "heavy" chunks in episodes/03-community-analysis.Rmd
#
# Why: chunks that use mia and scater need ~170 Bioconductor packages, which
# GitHub Actions cannot install in time. Those chunks are marked
#   ```{r name, eval = FALSE, heavy = TRUE}
# so the website shows their code without running it. This script runs them
# on YOUR computer and saves what they print and plot, and the episode shows
# those saved results in the matching `name-out` chunks.
#
# How: from the top of the lesson folder, with mia, scater, vegan, readr,
# dplyr, tidyr, ggplot2 and knitr installed locally, run
#   Rscript scripts/precompute_ep03.R
# Re-run it whenever you change the episode's code or its data, then commit
#   episodes/data/ep03_precomputed.RData  and  episodes/fig/03-*.png
#
# It knits the episode itself with the heavy chunks switched on, so the code
# that produced the saved output is exactly the code shown in the lesson.

library(knitr)

episode <- "episodes/03-community-analysis.Rmd"
out_rdata <- "episodes/data/ep03_precomputed.RData"
fig_prefix <- "fig/03-"   # relative to episodes/
stopifnot(file.exists(episode))
dir.create("episodes/fig", showWarnings = FALSE)

# Remove old figures so stale images are never left behind.
unlink(Sys.glob("episodes/fig/03-*.png"))

saved <- list(text = list(), figs = list())

# Set up markdown output first, so knit() keeps the hooks we add below.
render_markdown()
opts_chunk$set(
  comment = "",          # save output without "##"; the website adds its own
  error = FALSE,         # stop on the first error instead of saving it
  fig.path = file.path(tempdir(), "live-"),  # figures of normal chunks: throw away
  dpi = 150
)

# heavy = TRUE: run the chunk and save its figures under episodes/fig/
opts_hooks$set(heavy = function(options) {
  if (isTRUE(options$heavy)) {
    options$eval <- TRUE
    options$fig.path <- fig_prefix
  }
  options
})
# saved_output = TRUE: these chunks show old saved results, so skip them here
opts_hooks$set(saved_output = function(options) {
  if (isTRUE(options$saved_output)) options$eval <- FALSE
  options
})

# Record what heavy chunks print and plot.
default_output <- knit_hooks$get("output")
knit_hooks$set(output = function(x, options) {
  if (isTRUE(options$heavy)) {
    saved$text[[options$label]] <<- c(saved$text[[options$label]], x)
  }
  default_output(x, options)
})
default_plot <- knit_hooks$get("plot")
knit_hooks$set(plot = function(x, options) {
  if (isTRUE(options$heavy)) {
    saved$figs[[options$label]] <<- c(saved$figs[[options$label]], x)
  }
  default_plot(x, options)
})

# Knit the whole episode from inside episodes/, as the website build does,
# so "data/..." and "fig/..." paths resolve the same way. The rendered
# markdown is only a by-product.
owd <- setwd("episodes")
invisible(knit(basename(episode), output = file.path(tempdir(), "03-precompute.md"),
               quiet = TRUE))
setwd(owd)

# Check that every display chunk in the episode has something to show.
rmd <- readLines(episode)
wanted <- sub(".*\\{r ([^,}]+)-out,.*", "\\1",
              grep("^```\\{r [^,}]+-out,.*saved_output = TRUE", rmd, value = TRUE))
missing <- setdiff(wanted, c(names(saved$text), names(saved$figs)))
if (length(missing)) stop("No output was saved for: ", paste(missing, collapse = ", "))

precomputed <- saved
save(precomputed, file = out_rdata)
message("Saved ", out_rdata, " (", length(saved$text), " text outputs) and ",
        length(unlist(saved$figs)), " figures in episodes/fig/")
