#' ======================================================================
#' variables.R
#'
#' Shared variables: color palettes, cell-type ordering, and tissue colors.
#' Sourced by pipeline and figure scripts to keep palette definitions DRY.
#'
#' Author: Michael B. Thomsen, Ph.D.
#' Created: 2024-07-01
#' ======================================================================


# Class (high-level) colors — 7 classes in the atlas.
class_colors <- c(
  "NEU"   = "#762a83",
  "GLIA"  = "#0072B5",
  "IMM"   = "#CC50A7",
  "FIBRO" = "#C1DB1B",
  "EC"    = "#D7191C",
  "MURAL" = "#EA6337",
  "ERY"   = "#FA242F"
)


# Canonical class-grouped ordering of Major cell types. Matches the labels
# produced by 07_main_namexfer.R. *_CYC labels represent cycling (mitotic)
# subpopulations within each class.
class_ordered_majors <- c(
  "NEU_A", "NEU_C", "NEU_ATF3",                                      # NEU
  "SGC", "MYSC", "NMSC", "G_CYC",                                    # GLIA
  "EF", "EM", "PM", "F_CYC",                                         # FIBRO
  "PERI", "VSMC", "M_CYC",                                           # MURAL
  "VEC_A", "VEC_V", "LEC", "E_CYC",                                  # EC
  "IMG", "MAC", "MONO", "DC", "NT", "EOS", "MAST", "LEUK", "H_CYC",  # IMM
  "ERY"
)


# Major cell type colors.
# Hex codes carried over from figure_1_panels.R; names remapped to the
# current *_CYC convention (GSC→G_CYC, FSC→F_CYC, MSC→M_CYC, VSC→E_CYC, HSC→H_CYC).
major_cell_type_colors <- c(
  "NEU_A"    = "#E8C34A",
  "NEU_C"    = "#E0DF4E",
  "NEU_ATF3" = "#C4E04E",
  "SGC"      = "#DB46E6",
  "MYSC"     = "#A3E6BA",
  "NMSC"     = "#5FE17D",
  "G_CYC"    = "#D1868F",
  "EF"       = "#8546E8",
  "EM"       = "#DE5258",
  "PM"       = "#B6C0DB",
  "F_CYC"    = "#E1A257",
  "PERI"     = "#6D81DE",
  "VSMC"     = "#8AEB44",
  "M_CYC"    = "#6AE1D9",
  "VEC_A"    = "#DDDB8E",
  "VEC_V"    = "#DBDEBE",
  "LEC"      = "#5FE7B1",
  "E_CYC"    = "#D7AB90",
  "IMG"      = "#7894C3",
  "MAC"      = "#67C4E1",
  "MONO"     = "#A2D673",
  "DC"       = "#77A18F",
  "NT"       = "#BE74DA",
  "EOS"      = "#CEA5DD",
  "MAST"     = "#DB4FB4",
  "LEUK"     = "#CBE8E4",
  "H_CYC"    = "#E9CBDB",
  "ERY"      = "#DA77B3"
)


# Major-to-Class mapping (matches 07_main_namexfer.R remap logic).
# Useful for aggregating/coloring Major-level results by high-level class.
major_to_class <- c(
  "NEU_A" = "NEU", "NEU_C" = "NEU", "NEU_ATF3" = "NEU",
  "SGC" = "GLIA", "MYSC" = "GLIA", "NMSC" = "GLIA", "G_CYC" = "GLIA",
  "EF" = "FIBRO", "EM" = "FIBRO", "PM" = "FIBRO", "F_CYC" = "FIBRO",
  "PERI" = "MURAL", "VSMC" = "MURAL", "M_CYC" = "MURAL",
  "VEC_A" = "EC", "VEC_V" = "EC", "LEC" = "EC", "E_CYC" = "EC",
  "IMG" = "IMM", "MAC" = "IMM", "MONO" = "IMM", "DC" = "IMM",
  "NT" = "IMM", "EOS" = "IMM", "MAST" = "IMM", "LEUK" = "IMM", "H_CYC" = "IMM",
  "ERY" = "ERY"
)


# Glia Minor cell type colors — CANONICAL SOURCE OF TRUTH for every glia plot.
# Scheme chosen to match the immunofluorescence figures and the current Figure 3
# legend:  SGC = purples,  nmSC = yellows,  mSC = blues,  G_CYC / G_PROG = greens.
# Every projection / dynamics / pseudotime / velocity script sources THIS vector.
# Do NOT revert to the retired scheme (SGC = blue / nmSC = purple / mSC = teal /
# prog = yellow-green) — that mismatch is what produced the wrong glia colors in the
# label-transfer figures (S8/S10/S11/S13).
glia_cell_type_colors <- c(
  "SGC"        = "#b74ad3",
  "SGC IEG"    = "#7d4ad3",
  "SGC IFN"    = "#9901a3",
  "nmSC"       = "#ffff70",
  "nmSC IEG"   = "#f7f197",
  "nmSC Ngfr+" = "#ffc20a",
  "mSC"        = "#1a85ff",
  "mSC IEG"    = "#005FDB",
  "mSC REP"    = "#00C1FD",
  "G_CYC"      = "#00F7C5",
  "G_PROG"     = "#00D5AF"
)


# Tissue fill colors used across abundance/distribution plots.
tissue_colors <- c(
  "DRG" = "#7AA6DC",
  "SN"  = "#F2AD00"
)
