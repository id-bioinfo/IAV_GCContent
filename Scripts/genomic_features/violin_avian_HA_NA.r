library(ggplot2)
library(ggstatsplot)

avian_ha <- read.table('GC_subtype_avianFamily_data/gc_avianHA_year.txt', header = T, sep = "\t")
plt_ha <- ggbetweenstats(
  data = avian_ha,
  plot.type = "boxviolin",
  x = GeneId,
  y = GCContent,
  results.subtitle = FALSE,
  k = 1,
  ggsignif.args = list(textsize = 1, tip_length = 0.01),
  centrality.plotting = FALSE,
  palette = "Set3"
)

plt_ha


avian_na <- read.table('GC_subtype_avianFamily_data/gc_avianNA_year.txt', header = T, sep = "\t")
plt_na <- ggbetweenstats(
  data = avian_na,
  plot.type = "boxviolin",
  x = GeneId,
  y = GCContent,
  results.subtitle = FALSE,
  k = 1,
  ggsignif.args = list(textsize = 1, tip_length = 0.01),
  centrality.plotting = FALSE,
  palette = "Set3"
)

plt_na
