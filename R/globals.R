# =============================================================================
# Global variable declarations
# =============================================================================
# Several functions use non-standard evaluation (dplyr/tidyr verbs) with bare
# column names. Declaring them here prevents spurious R CMD check NOTEs of the
# form "no visible binding for global variable". This has no runtime effect.
# =============================================================================

utils::globalVariables(c(
  # bundled dataset object
  "example_edna",
  # bare column names referenced in dplyr/tidyr NSE
  "x_label",
  "cape_x",
  "SampleID",
  "Species",
  "Counts",
  "TrueCommunity",
  "Relative_frequency",
  "BioRep"
))
