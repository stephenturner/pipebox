#!/usr/bin/env Rscript

# For interactive development
# args <- "/Users/turner/repos/pipebox/testdata/results/SRR507778"

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  theme_set(theme_bw())
})

# Get the outbase as the first command line argument
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("Usage: pipebox-post <outbase>")
} else {
  outbase <- args[1]
}
rm(args)

# Read in the seqtk results
fqchk <- read_table(paste0(outbase, ".seqtkfqchk.tsv"),
                    skip=3,
                    col_names=c("POS", "bases", "pA", "pC", "pG", "pT", "pN", "avgQ", "errQ", "plow", "phigh"),
                    show_col_types = FALSE)

# Per base quality score plot
p_readqual <-
  fqchk |>
  mutate(Position=(POS)) |>
  mutate(Quality=cut(avgQ, breaks=c(0, 20, 30, Inf), labels = c("<Q20", "Q20-30", ">Q30"))) |>
  ggplot(aes(Position, avgQ)) +
  geom_col(aes(fill=Quality)) +
  theme_classic() +
  scale_fill_manual(values=c("red3", "goldenrod2", "green4"), drop=FALSE) +
  labs(y="Per-base quality score")

# Per base nucleotide composition
p_nuc <-
  fqchk |>
  select(Position=POS, pA:pT) |>
  pivot_longer(!Position, names_to="Nucleotide", values_to="Percent") |>
  mutate(Nucleotide=gsub("^p", "", Nucleotide)) |>
  mutate(Percent=Percent/100) |>
  ggplot(aes(Position, Percent, fill=Nucleotide)) + geom_col() +
  scale_y_continuous(name="Per-base nucleotide composition", labels=scales::percent)

# Read in the samtools stats
samstats <- read_lines(paste0(outbase, ".samtoolsstats.tsv"))

# Get out just the coverage information
samstats_cov <-
  samstats |>
  grep("^COV", x=_, value=TRUE) |>
  I() |>
  read_tsv(col_names=c("COV", "bin", "cov", "n"), show_col_types = FALSE)

# Get a coverage plot
p_cov <-
  samstats_cov |>
  ggplot(aes(cov, n)) +
  geom_point() +
  geom_line() +
  scale_y_log10(labels=scales::number) +
  labs(x="Coverage (X)", y="Number of bases covered at X coverage") +
  coord_cartesian(xlim=c(0,40))

# Read in the samtools stats
bcfstats <- read_lines(paste0(outbase, ".bcfstats.tsv"))

# Get out just the coverage information
bcfstats_snps <-
  bcfstats |>
  grep("^ST", x=_, value=TRUE) |>
  I() |>
  read_tsv(col_names=c("ST", "id", "SNP", "N"), show_col_types = FALSE) |>
  mutate(Type=ifelse(SNP %in% c("A>C", "C>A", "C>T", "T>C"), "Transition", "Transversion"))

# Variant type plot
p_var <-
  bcfstats_snps |>
  ggplot(aes(SNP, N)) + geom_col(aes(fill=Type))

# Combine all the plots with patchwork
p <- (p_readqual+p_nuc)/(p_cov+p_var)
# Add a title
p <- p + plot_annotation(title=basename(outbase))
# Add space between plots
p <- p & theme(plot.margin = unit(c(.5,.5,.5,.5), "cm"), plot.title = element_text(face="bold", hjust = 0.5))
# Save to file
ggsave(paste0(outbase, ".metrics.pdf"), plot = p, width=12, height=8)
ggsave(paste0(outbase, ".metrics.png"), plot = p, width=12, height=8)
