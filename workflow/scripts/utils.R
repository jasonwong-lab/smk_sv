load_pkg <- function(pkgs) {
  for (pkg in pkgs) {
    if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))) {
      install.packages(pkg)
      if (!suppressWarnings(suppressMessages(require(pkg, character.only = TRUE)))) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
        suppressMessages(require(pkg, character.only = TRUE))
      }
    }
  }
}

load_pkg(
  c("vroom", "tibble", "glue", "dplyr", "tidyr", "purrr", "GenomicRanges", "stringr", "BiocParallel", "parallel")
)

my_vroom <- function(file, na_append = NULL, delim = "\t", comment = "#", col_names = TRUE, guess_lines = 10, n_col = NULL, ...) {
  na_origin <- c("", "NA", "NaN")
  myna <- c(na_origin, na_append)

  if (is.null(n_col)) {
    n_col <- vroom::vroom(
      file,
      delim = delim, n_max = guess_lines, col_types = vroom::cols(),
      na = myna, comment = comment, col_names = col_names, guess_max = Inf, ...
    ) |>
      ncol()
  }

  dat <- vroom::vroom(
    file,
    delim = delim, col_types = paste0(rep_len("c", n_col), collapse = ""),
    na = myna, comment = comment, col_names = col_names, ...
  )

  dat
}

nothing <- function(x) {
  x
}

create_grs <- function(string) {
  if (is.na(string)[1]) return(as("1:1-1", "GRanges"))
  as(string, "GRanges")
}

calculate_overlap_flag <- function(coord, sv) {
  gr <- str_replace_all(coord, "chr", "") |>
    str_split(";") |>
    unlist() |>
    create_grs()
  # pct <- max(width(pintersect(gr, sv)) / width(gr))
  # pct
  pcts_gr <- width(pintersect(gr, sv)) / width(gr)
  pcts_sv <- width(pintersect(gr, sv)) / width(sv)
  any(pcts_gr >= 0.9 & pcts_sv >= 0.8)
}

add_overlap_flag_annotsv <- function(tsv, workers = detectCores()) {
  columns <- c("B_gain_coord", "B_loss_coord", "po_B_loss_allG_coord", "B_ins_coord", "po_B_gain_allG_coord", "B_inv_coord")
  # pct_columns <- c("pct_gain", "pct_loss", "pct_po_loss", "pct_ins", "pct_po_gain", "pct_inv")
  flag_columns <- c("flag_gain", "flag_loss", "flag_po_loss", "flag_ins", "flag_po_gain", "flag_inv")

  # for (i in seq_along(pct_columns)) {
  #   tsv[[pct_columns[i]]] <- rep(0, nrow(tsv))
  # }
  for (i in seq_along(flag_columns)) {
    tsv[[flag_columns[i]]] <- rep(FALSE, nrow(tsv))
  }

  gr_sv <- GRanges(seqnames = as.character(tsv$SV_chrom), ranges = IRanges(as.numeric(tsv$SV_start), as.numeric(tsv$SV_end)))

  register(MulticoreParam(workers = workers))

  for (j in seq_along(columns)) {
    # tsv[[pct_columns[j]]] <- bplapply(seq_len(nrow(tsv)), function(i) {
    #   calculate_overlap(tsv[[columns[j]]][i], gr_sv[i])
    # }) |>
    #   unlist()
    tsv[[flag_columns[j]]] <- bplapply(seq_len(nrow(tsv)), function(i) {
      calculate_overlap_flag(tsv[[columns[j]]][i], gr_sv[i])
    }) |>
      unlist()
  }

  tsv
}

# *--------------------------------------------------------------------------* #
# * QUAL Filtration - Decided not to do this finally                         * #
# * SVision: SV quality of the SV described in this region                   * #
# * Severus: Not found                                                       * #
# * SVIM: Score of each call on a range between 0 and 100 (not phred-scaled) * #
# *       Because the score depends on the number of supporting reads, its   * #
# *       distribution varies with the sequencing coverage in the input.     * #
# *       Therefore, we cannot make a general statement about suitable score * #
# *       cutoffs. For high-coverage datasets (>40x), we would recommend a   * #
# *       threshold of 10-15. For low-coverage datasets, the threshold       * #
# *       should be lower.                                                   * #
# * cuteSV: QUAL of an SV is less than 5 will be marked as q5 and            * #
# *         recommended to filter out.                                       * #
# *         The quality scores are calculated when cuteSV detects the SV     * #
# *         genotypes.                                                       * #
# * Sniffles: This is currently not indicated                                * #
# *--------------------------------------------------------------------------* #

filter_vep <- function(v, maf_vep, t) {
  id_include <- switch(v,
    BND = {
      maf_vep |> nothing()
      # dplyr::filter(IMPACT == "HIGH", BIOTYPE == "protein_coding")
    },
    DEL = {
      maf_vep |> nothing()
      # dplyr::filter(IMPACT == "HIGH", BIOTYPE == "protein_coding")
    },
    INS = {
      maf_vep |> nothing()
      # dplyr::filter(IMPACT == "HIGH", BIOTYPE == "protein_coding")
    },
    INV = {
      maf_vep |> nothing()
      # dplyr::filter(BIOTYPE == "protein_coding", IMPACT %in% c("HIGH", "MODERATE"))
    },
    DUP = {
      maf_vep |> nothing()
      # dplyr::filter(IMPACT %in% c("HIGH", "MODERATE"), BIOTYPE == "protein_coding")
    }
  ) |>
    dplyr::pull(vcf_id) |>
    na.omit() |>
    unique()

  # min_qual <- switch(t, cutesv = 5, severus = 20, sniffles = 20, svim = 5, svision = 20)
  min_qual <- switch(t, cutesv = 0, severus = 0, sniffles = 0, svim = 0, svision = 0, 0)

  id_exclude <- switch(v,
    BND = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    DEL = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    INS = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
            as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 |
            as.numeric(vcf_qual) < min_qual
        )
    },
    INV = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    },
    DUP = {
      maf_vep |>
        dplyr::filter(
          as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 | as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 | as.numeric(vcf_qual) < min_qual
        )
    }
  ) |>
    dplyr::pull(vcf_id) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

filter_snpeff <- function(v, tsv_snpeff, t) {
  id_include <- switch(v,
    BND = {
      tsv_snpeff |> nothing()
      # dplyr::filter(grepl("HIGH", IMPACT), grepl("protein_coding", BIOTYPE))
    },
    DEL = {
      tsv_snpeff |> nothing()
      # dplyr::filter(grepl("HIGH", IMPACT), grepl("protein_coding", BIOTYPE))
    },
    INS = {
      tsv_snpeff |> nothing()
      # dplyr::filter(grepl("HIGH", IMPACT), grepl("protein_coding", BIOTYPE))
    },
    INV = {
      tsv_snpeff |> nothing()
      # dplyr::filter(grepl("protein_coding", BIOTYPE), grepl("HIGH|MODERATE", `ANN[*].IMPACT`))
    },
    DUP = {
      tsv_snpeff |> nothing()
      # dplyr::filter(grepl("HIGH|MODERATE", IMPACT), grepl("protein_coding", `ANN[*].BIOTYPE`))
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  # min_qual <- switch(t, cutesv = 5, severus = 20, sniffles = 20, svim = 5, svision = 20)
  min_qual <- switch(t, cutesv = 0, severus = 0, sniffles = 0, svim = 0, svision = 0, 0)

  id_exclude <- switch(v,
    BND = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    DEL = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    INS = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    INV = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    },
    DUP = {
      tsv_snpeff |>
        dplyr::filter(as.numeric(QUAL) < min_qual)
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

filter_annotsv <- function(v, tsv_annotsv, terms_relative) {
  id_include <- switch(v,
    BND = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
          # // !is.na(P_gain_coord) | !is.na(P_loss_coord) | !is.na(P_ins_coord) | !is.na(po_P_gain_coord) | !is.na(po_P_loss_coord)
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
          # // !is.na(P_loss_coord) | !is.na(P_ins_coord) | !is.na(po_P_gain_coord) | !is.na(po_P_loss_coord)
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
          # // !is.na(P_loss_coord) | !is.na(P_ins_coord) | !is.na(po_P_gain_coord) | !is.na(po_P_loss_coord)
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
          # // !is.na(P_loss_coord) | !is.na(P_ins_coord) | !is.na(po_P_gain_coord) | !is.na(po_P_loss_coord)
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          grepl(terms_relative, GenCC_disease, ignore.case = TRUE) | grepl(terms_relative, OMIM_phenotype, ignore.case = TRUE) | grepl("4|5", ACMG_class) | as.numeric(AnnotSV_ranking_score) >= 0.9
          # // !is.na(P_loss_coord) | !is.na(P_ins_coord) | !is.na(po_P_gain_coord) | !is.na(po_P_loss_coord)
        )
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  id_exclude <- switch(v,
    BND = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_gain_AFmax) >= 0.05 & as.numeric(B_loss_AFmax) >= 0.05 &
          flag_gain & flag_loss
          # & pct_gain >= 0.9 & pct_loss >= 0.9
          # // !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right) | as.numeric(DDD_HI_percent) >= 95
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_loss_AFmax) >= 0.05 & flag_loss) | (!is.na(po_B_loss_allG_coord) & flag_po_loss)
          # (as.numeric(B_loss_AFmax) >= 0.05 & pct_loss >= 0.9) | (!is.na(po_B_loss_allG_coord) & pct_po_loss >= 0.9)
          # // !is.na(po_B_loss_someG_coord) | !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right) | as.numeric(ExAC_delZ) <= 0 | as.numeric(DDD_HI_percent) >= 90
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_ins_AFmax) >= 0.05 & flag_ins) | (!is.na(po_B_gain_allG_coord) & flag_po_gain)
          # (as.numeric(B_ins_AFmax) >= 0.05 & pct_ins >= 0.9) | (!is.na(po_B_gain_allG_coord) & pct_po_gain >= 0.9)
          # // !is.na(po_B_loss_someG_coord) | !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right) | as.numeric(ExAC_cnvZ) <= 0 | as.numeric(DDD_HI_percent) >= 90
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_inv_AFmax) >= 0.05 & flag_inv
          # as.numeric(B_inv_AFmax) >= 0.05 & pct_inv >= 0.9
          # // !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right) | as.numeric(DDD_HI_percent) >= 90
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          (as.numeric(B_gain_AFmax) >= 0.05 & flag_gain) | (!is.na(po_B_gain_allG_coord) & flag_po_gain)
          # (as.numeric(B_gain_AFmax) >= 0.05 & pct_gain >= 0.9) | (!is.na(po_B_gain_allG_coord) & pct_po_gain >= 0.9)
          # // !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right) | as.numeric(DDD_HI_percent) >= 90
        )
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

merge_svid <- function(list_vep, list_snpeff, list_annotsv, stringent = TRUE) {
  if (stringent) {
    ids_include <- unique(
      c(list_vep$id_include, list_snpeff$id_include, list_annotsv$id_include)
    )
    ids_exclude <- unique(
      c(list_vep$id_exclude, list_snpeff$id_exclude, list_annotsv$id_exclude)
    )
    ids <- setdiff(ids_include, ids_exclude) |>
      union(list_annotsv$id_include) |>
      unique()
  } else {
    ids_vep <- setdiff(list_vep$id_include, list_vep$id_exclude)
    ids_snpeff <- setdiff(list_snpeff$id_include, list_snpeff$id_exclude)
    ids <- setdiff(c(ids_vep, ids_snpeff), list_annotsv$id_exclude) |>
      union(list_annotsv$id_include) |>
      unique()
  }

  ids
}

filter_svid <- function(caller, sv_type, id_tab, annotsv_result, vep_result, snpeff_result, terms_relative, ...) {
  if (nrow(vep_result) == 0) {
    list_vep <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_vep <- filter_vep(v = sv_type, maf_vep = vep_result, t = caller)
  }
  if (nrow(snpeff_result) == 0) {
    list_snpeff <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_snpeff <- filter_snpeff(v = sv_type, tsv_snpeff = snpeff_result, t = caller)
  }
  if (nrow(annotsv_result) == 0) {
    list_annotsv <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_annotsv <- filter_annotsv(v = sv_type, tsv_annotsv = annotsv_result, terms_relative = terms_relative)
  }

  ids <- merge_svid(list_vep = list_vep, list_snpeff = list_snpeff, list_annotsv = list_annotsv, ...)

  ids
}
