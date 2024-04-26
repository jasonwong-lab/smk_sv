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
  c(
    "vroom", "tibble", "glue", "dplyr", "tidyr", "purrr"
  )
)

my_vroom <- function(file, na_append = NULL, delim = "\t", comment = "#", col_names = TRUE) {
  na_origin <- c("", "NA", "NaN")
  myna <- c(na_origin, na_append)

  n_col <- vroom::vroom(
    file,
    delim = delim, n_max = 10, col_types = vroom::cols(),
    na = myna, comment = comment, col_names = col_names
  ) |>
    ncol()

  dat <- vroom::vroom(
    file,
    delim = delim, col_types = paste0(rep_len("c", n_col), collapse = ""),
    na = myna, comment = comment, col_names = col_names
  )

  dat
}

nothing <- function(x) {
  x
}

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

  id_exclude <- switch(v,
    BND = {
      if (t != "svim") {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01 |
              as.numeric(vcf_qual) < 20
          )
      } else {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      }
    },
    DEL = {
      if (t != "svim") {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(vcf_qual) < 20 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      } else {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      }
    },
    INS = {
      if (t != "svim") {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(vcf_qual) < 20 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      } else {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      }
    },
    INV = {
      if (t != "svim") {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(vcf_qual) < 20 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      } else {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      }
    },
    DUP = {
      if (t != "svim") {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(vcf_qual) < 20 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      } else {
        maf_vep |>
          dplyr::filter(
            as.numeric(AF) >= 0.01 | as.numeric(EAS_AF) >= 0.01 | as.numeric(ASN_AF) >= 0.01 |
              as.numeric(gnomAD_AF) >= 0.01 | as.numeric(gnomAD_EAS_AF) >= 0.01
          )
      }
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

  id_exclude <- switch(v,
    BND = {
      if (t != "svim") {
        tsv_snpeff |>
          dplyr::filter(as.numeric(QUAL) < 20)
        # tibble::tibble(ID = NA)
      } else {
        tibble::tibble(ID = NA)
      }
    },
    DEL = {
      if (t != "svim") {
        tsv_snpeff |>
          dplyr::filter(as.numeric(QUAL) < 20)
      } else {
        tibble::tibble(ID = NA)
      }
    },
    INS = {
      if (t != "svim") {
        tsv_snpeff |>
          dplyr::filter(as.numeric(QUAL) < 20)
      } else {
        tibble::tibble(ID = NA)
      }
    },
    INV = {
      if (t != "svim") {
        tsv_snpeff |>
          dplyr::filter(as.numeric(QUAL) < 20)
      } else {
        tibble::tibble(ID = NA)
      }
    },
    DUP = {
      if (t != "svim") {
        tsv_snpeff |>
          dplyr::filter(as.numeric(QUAL) < 20)
      } else {
        tibble::tibble(ID = NA)
      }
    }
  ) |>
    dplyr::pull(ID) |>
    na.omit() |>
    unique()

  list(id_include = id_include, id_exclude = id_exclude)
}

filter_annotsv <- function(v, tsv_annotsv) {
  id_include <- switch(v,
    BND = {
      tsv_annotsv |>
        dplyr::filter(
          grepl("leuka?emia", GenCC_disease, ignore.case = TRUE) |
            grepl("leuka?emia", OMIM_phenotype, ignore.case = TRUE) |
            grepl("^4$|^5$", ACMG_class) |
            as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          grepl("leuka?emia", GenCC_disease, ignore.case = TRUE) |
            grepl("leuka?emia", OMIM_phenotype, ignore.case = TRUE) |
            grepl("^4$|^5$", ACMG_class) |
            as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          grepl("leuka?emia", GenCC_disease, ignore.case = TRUE) |
            grepl("leuka?emia", OMIM_phenotype, ignore.case = TRUE) |
            grepl("^4$|^5$", ACMG_class) |
            as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          grepl("leuka?emia", GenCC_disease, ignore.case = TRUE) |
            grepl("leuka?emia", OMIM_phenotype, ignore.case = TRUE) |
            grepl("^4$|^5$", ACMG_class) |
            as.numeric(AnnotSV_ranking_score) >= 0.9
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          grepl("leuka?emia", GenCC_disease, ignore.case = TRUE) |
            grepl("leuka?emia", OMIM_phenotype, ignore.case = TRUE) |
            grepl("^4$|^5$", ACMG_class) |
            as.numeric(AnnotSV_ranking_score) >= 0.9
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
          as.numeric(B_gain_AFmax) >= 0.05 & as.numeric(B_loss_AFmax) >= 0.05 |
            as.numeric(DDD_HI_percent) >= 95
          # !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right)
        )
    },
    DEL = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_loss_AFmax) >= 0.05 |
            as.numeric(ExAC_delZ) <= 0 | as.numeric(DDD_HI_percent) >= 90 |
            !is.na(po_B_loss_allG_coord)
          # !is.na(po_B_loss_someG_coord) | !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right)
        )
    },
    INS = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_ins_AFmax) >= 0.05 |
            as.numeric(ExAC_cnvZ) <= 0 | as.numeric(DDD_HI_percent) >= 90 |
            !is.na(po_B_gain_allG_coord)
          # !is.na(po_B_loss_someG_coord) | !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right)
        )
    },
    INV = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_inv_AFmax) >= 0.05 |
            as.numeric(DDD_HI_percent) >= 90
          # !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right)
        )
    },
    DUP = {
      tsv_annotsv |>
        dplyr::filter(
          as.numeric(B_gain_AFmax) >= 0.05 |
            as.numeric(DDD_HI_percent) >= 90 |
            !is.na(po_B_gain_allG_coord)
          # !is.na(ENCODE_blacklist_left) | !is.na(ENCODE_blacklist_right)
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

filter_svid <- function(caller, sv_type, id_tab, annotsv_result, vep_result, snpeff_result, ...) {
  if (nrow(vep_result) == 0) {
    list_vep <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_vep <- filter_vep(sv_type, vep_result, caller)
  }
  if (nrow(snpeff_result) == 0) {
    list_snpeff <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_snpeff <- filter_snpeff(sv_type, snpeff_result, caller)
  }
  if (nrow(annotsv_result) == 0) {
    list_annotsv <- list(id_include = NULL, id_exclude = NULL)
  } else {
    list_annotsv <- filter_annotsv(sv_type, annotsv_result)
  }

  ids <- merge_svid(list_vep, list_snpeff, list_annotsv, ...)

  ids
}
