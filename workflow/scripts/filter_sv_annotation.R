.libPaths(snakemake@params[["libs_r"]])
snakemake@source("utils.R")

file_annotsv <- snakemake@input[["annotsv"]]
file_maf <- snakemake@input[["maf"]]
file_tsv <- snakemake@input[["tsv"]]
file_tab <- snakemake@input[["tab"]]
s <- snakemake@wildcards[["sample"]]
t <- snakemake@wildcards[["caller"]]
v <- snakemake@wildcards[["type_sv"]]
vcf <- snakemake@input[["vcf"]]
vcf_final <- snakemake@output[["vcf_final"]]
ids_tmp <- snakemake@output[["ids"]]
table_tsv <- snakemake@output[["table"]]


annotsv <- file_annotsv |>
  my_vroom()
maf_vep <- file_maf |>
  my_vroom() |>
  tidyr::drop_na(Start_Position)
tsv_snpeff <- file_tsv |>
  my_vroom()
id_info <- file_tab |>
  my_vroom(col_names = FALSE) |>
  setNames(c("chr", "pos", "ref", "alt", "sniffles", "cutesv", "svim"))
if (nrow(annotsv) != 0) {
  tsv_annotsv <- annotsv |>
    dplyr::filter(ID %in% na.omit(id_info[[t]]))
} else {
  tsv_annotsv <- annotsv
}
ids <- filter_svid(
  caller = t, sv_type = v,
  annotsv_result = tsv_annotsv, vep_result = maf_vep, snpeff_result = tsv_snpeff,
  stringent = FALSE
)

list_annotation <- list(VEP = maf_vep, SnpEff = tsv_snpeff, AnnotSV = tsv_annotsv)
annotation_modified <- lapply(seq_along(list_annotation), function(i) {
  x <- list_annotation[[i]]
  if ("vcf_id" %in% colnames(x)) {
    x <- dplyr::mutate(x, ID = vcf_id)
  }
  if (length(ids) != 0 && nrow(x) != 0) {
    tmp <- x[x$ID %in% ids, ]
    colnames(tmp) <- glue::glue("{names(list_annotation)[i]}.{colnames(tmp)}")
    tmp
  } else {
    NULL
  }
}) |>
  dplyr::bind_rows()

write.table(annotation_modified, table_tsv, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t", na = "")

write.table(ids, ids_tmp, quote = FALSE, col.names = FALSE, row.names = FALSE)
system(glue::glue('bcftools view -i "ID=@{ids_tmp}" {vcf} > {vcf_final}'))

save.image(snakemake@output[["rdata"]])