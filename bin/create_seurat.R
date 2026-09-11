#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
})

parse_args <- function(args) {
    out <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("Unexpected argument: ", key)
        }
        if (i == length(args)) {
            stop("Missing value for ", key)
        }
        out[[substring(key, 3L)]] <- args[[i + 1L]]
        i <- i + 2L
    }
    out
}

required_arg <- function(args, name) {
    value <- args[[name]]
    if (is.null(value) || !nzchar(value)) {
        stop("Missing required argument --", name)
    }
    value
}

read_features <- function(dir) {
    gz_path <- file.path(dir, "features.tsv.gz")
    plain_path <- file.path(dir, "features.tsv")
    path <- if (file.exists(gz_path)) gz_path else plain_path
    if (!file.exists(path)) {
        stop("Missing 10x feature table: ", gz_path, " or ", plain_path)
    }
    con <- if (endsWith(path, ".gz")) gzfile(path, open = "rt") else file(path, open = "rt")
    on.exit(close(con), add = TRUE)
    x <- read.delim(
        con,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE
    )
    if (ncol(x) < 2L) {
        stop("Expected at least two columns in ", path)
    }
    if (ncol(x) < 3L) {
        x[[3L]] <- NA_character_
    }
    names(x)[1:3] <- c("feature_id", "feature_name", "feature_type")
    x$seurat_feature <- make.unique(as.character(x$feature_name))
    rownames(x) <- x$seurat_feature
    x
}

read_10x_matrix <- function(dir) {
    matrix_gz <- file.path(dir, "matrix.mtx.gz")
    matrix_plain <- file.path(dir, "matrix.mtx")
    barcodes_gz <- file.path(dir, "barcodes.tsv.gz")
    barcodes_plain <- file.path(dir, "barcodes.tsv")
    features_gz <- file.path(dir, "features.tsv.gz")
    features_plain <- file.path(dir, "features.tsv")

    matrix_path <- if (file.exists(matrix_gz)) matrix_gz else matrix_plain
    barcodes_path <- if (file.exists(barcodes_gz)) barcodes_gz else barcodes_plain
    features_path <- if (file.exists(features_gz)) features_gz else features_plain

    missing <- c(
        matrix = matrix_path,
        barcodes = barcodes_path,
        features = features_path
    )[!file.exists(c(matrix_path, barcodes_path, features_path))]
    if (length(missing) > 0L) {
        stop("Incomplete 10x MEX directory ", dir, ": missing ", paste(names(missing), collapse = ", "))
    }

    open_text <- function(path) {
        if (endsWith(path, ".gz")) gzfile(path, open = "rt") else file(path, open = "rt")
    }

    matrix_con <- open_text(matrix_path)
    on.exit(close(matrix_con), add = TRUE)
    x <- Matrix::readMM(matrix_con)

    barcode_con <- open_text(barcodes_path)
    on.exit(close(barcode_con), add = TRUE)
    barcodes <- readLines(barcode_con)

    feature_con <- open_text(features_path)
    on.exit(close(feature_con), add = TRUE)
    features <- read.delim(
        feature_con,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE
    )
    if (ncol(features) < 2L) {
        stop("Expected at least two columns in ", features_path)
    }

    if (nrow(x) != nrow(features)) {
        stop("Feature count does not match matrix rows in ", dir)
    }
    if (ncol(x) != length(barcodes)) {
        stop("Barcode count does not match matrix columns in ", dir)
    }

    rownames(x) <- make.unique(as.character(features[[2L]]))
    colnames(x) <- barcodes
    as(x, "dgCMatrix")
}

align_sparse <- function(x, features, cells) {
    x <- as(x, "dgCMatrix")
    sx <- summary(x)
    fi <- match(rownames(x), features)
    ci <- match(colnames(x), cells)

    keep_rows <- !is.na(fi)
    keep_cols <- !is.na(ci)
    keep <- keep_rows[sx$i] & keep_cols[sx$j]

    Matrix::sparseMatrix(
        i = fi[sx$i[keep]],
        j = ci[sx$j[keep]],
        x = sx$x[keep],
        dims = c(length(features), length(cells)),
        dimnames = list(features, cells),
        giveCsparse = TRUE
    )
}

combine_feature_metadata <- function(primary, secondary, features) {
    all_cols <- union(colnames(primary), colnames(secondary))
    out <- as.data.frame(
        matrix(NA_character_, nrow = length(features), ncol = length(all_cols)),
        stringsAsFactors = FALSE
    )
    colnames(out) <- all_cols
    rownames(out) <- features

    for (src in list(primary, secondary)) {
        common_rows <- intersect(rownames(src), rownames(out))
        common_cols <- intersect(colnames(src), colnames(out))
        for (column in common_cols) {
            current <- out[common_rows, column]
            missing <- is.na(current) | current == ""
            if (any(missing)) {
                values <- as.character(src[common_rows, column])
                out[common_rows[missing], column] <- values[missing]
            }
        }
    }
    out
}

productivity_flag <- function(status) {
    out <- rep(NA, length(status))
    out[status == "productive"] <- TRUE
    out[!is.na(status) & nzchar(status) & status != "productive" & status != "unknown_no_v_cds"] <- FALSE
    out
}

add_vdj_metadata <- function(object, calls_path, receptors_path) {
    calls <- read.delim(
        calls_path,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    receptors <- read.delim(
        receptors_path,
        header = TRUE,
        sep = "\t",
        quote = "",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    needed_calls <- c("cell", "recombination_id", "productivity_status")
    needed_receptors <- c("cell", "heavy_recombination_id", "light_recombination_id")
    if (!all(needed_calls %in% colnames(calls))) {
        stop("VDJ calls table lacks required columns: ", paste(setdiff(needed_calls, colnames(calls)), collapse = ", "))
    }
    if (!all(needed_receptors %in% colnames(receptors))) {
        stop("VDJ receptor table lacks required columns: ", paste(setdiff(needed_receptors, colnames(receptors)), collapse = ", "))
    }

    cells <- colnames(object)
    idx <- match(cells, receptors$cell)
    hc <- receptors$heavy_recombination_id[idx]
    lc <- receptors$light_recombination_id[idx]
    hc[is.na(hc) | hc == ""] <- NA_character_
    lc[is.na(lc) | lc == ""] <- NA_character_

    call_key <- paste(calls$cell, calls$recombination_id, sep = "\r")
    status_lookup <- setNames(calls$productivity_status, call_key)

    hc_key <- paste(cells, hc, sep = "\r")
    lc_key <- paste(cells, lc, sep = "\r")
    hc_status <- unname(status_lookup[hc_key])
    lc_status <- unname(status_lookup[lc_key])

    object$HC_recombination_id <- hc
    object$LC_recombination_id <- lc
    object$HC_productive <- productivity_flag(hc_status)
    object$LC_productive <- productivity_flag(lc_status)
    object
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
exonic_dir <- required_arg(args, "exonic")
intronic_dir <- required_arg(args, "intronic")
output <- required_arg(args, "output")
project <- if (is.null(args$project)) "Norn" else args$project

options(Seurat.object.assay.version = "v5")

message("Reading exonic matrix: ", exonic_dir)
exonic <- read_10x_matrix(exonic_dir)
exonic_features <- read_features(exonic_dir)

message("Reading intronic matrix: ", intronic_dir)
intronic <- read_10x_matrix(intronic_dir)
intronic_features <- read_features(intronic_dir)

# The exonic matrix defines the canonical set of cells. Intronic-only cells are
# not introduced as zero-expression cells into the Seurat object.
cells <- colnames(exonic)
features <- union(rownames(exonic), rownames(intronic))

exonic <- align_sparse(exonic, features, cells)
intronic <- align_sparse(intronic, features, cells)

object <- Seurat::CreateSeuratObject(
    counts = exonic,
    assay = "RNA",
    project = project,
    min.cells = 0,
    min.features = 0
)

# Explicit velocity semantics while keeping normal Seurat counts == exonic.
object[["RNA"]]$spliced <- exonic
object[["RNA"]]$unspliced <- intronic

rna_meta <- combine_feature_metadata(exonic_features, intronic_features, features)
rna_meta$gene_id <- rna_meta$feature_id
rna_meta$gene_name <- rna_meta$feature_name
rna_meta <- rna_meta[features, , drop = FALSE]
object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], metadata = rna_meta)

object$cell_id <- colnames(object)

if (!is.null(args$`vdj-calls`) || !is.null(args$`vdj-receptors`)) {
    if (is.null(args$`vdj-calls`) || is.null(args$`vdj-receptors`)) {
        stop("--vdj-calls and --vdj-receptors must be supplied together")
    }
    message("Adding VDJ cell metadata")
    object <- add_vdj_metadata(object, args$`vdj-calls`, args$`vdj-receptors`)
} else {
    object$HC_recombination_id <- NA_character_
    object$LC_recombination_id <- NA_character_
    object$HC_productive <- NA
    object$LC_productive <- NA
}

object@misc$norn <- list(
    created_by = "Norn",
    project = project,
    seurat_version = as.character(utils::packageVersion("Seurat")),
    seurat_object_version = as.character(utils::packageVersion("SeuratObject")),
    RNA_layers = SeuratObject::Layers(object[["RNA"]])
)

message("Writing Seurat object: ", output)
saveRDS(object, file = output, compress = "gzip")

message("Created ", output)
message("  cells: ", ncol(object))
message("  RNA features: ", nrow(object[["RNA"]]))
message("  RNA layers: ", paste(SeuratObject::Layers(object[["RNA"]]), collapse = ", "))
