# ──────────────────────────────────────────────────────────────────────────────
# Dirichlet fallback usage: generate diagnostics (via your helper) and summarize
# Robust to whitespace/casing in glmnet_grouped/reshape_note columns
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(stringr); library(purrr)
    library(tibble); library(here)
})

# ============================== CONFIG =======================================
diag_dir   <- here("output", "dq_diag")                    # where diagnostics live
dictionary_dir <- here("data_raw", "cause-codes")          # holds the two Foreman CSVs
parquet_dir_candidates <- c(                               
    here("data_private", "mcod"),
    here("data_private", "mcod_sample")
)
county_var <- "county_ihme"                                # county column in your data
sample_max_rows <- 50000                                   # cap sample for quick run
out_prefix <- NULL  # e.g., here("output","fallback_summary"); NULL = don't save

# =============== tiny ICD helper for synthetic fallback ======================
canonical_icd <- function(x) stringr::str_remove_all(stringr::str_to_upper(x), "[^A-Z0-9]")

# ===================== summarize fallback (core) =============================
summarize_fallback <- function(diag_df) {
    names(diag_df) <- tolower(trimws(names(diag_df)))
    need <- c("garbage", "n_te", "glmnet_grouped", "reshape_note")
    miss <- setdiff(need, names(diag_df))
    if (length(miss)) stop("Diagnostics missing required columns: ", paste(miss, collapse = ", "))
    
    # robust coercions (trim, lower, map)
    to_logical_robust <- function(x) {
        if (is.logical(x)) return(x)
        sx <- trimws(as.character(x))
        lx <- tolower(sx)
        out <- ifelse(lx %in% c("true","t","1","yes","y"), TRUE,
                      ifelse(lx %in% c("false","f","0","no","n"), FALSE, NA))
        as.logical(out)
    }
    
    diag_df <- diag_df %>%
        mutate(
            garbage = as.character(garbage),
            n_te = suppressWarnings(as.numeric(n_te)),
            glmnet_grouped = to_logical_robust(glmnet_grouped),
            reshape_note = trimws(as.character(reshape_note)),
            reshape_note_l = tolower(reshape_note),
            # Define fallback strictly: only if glmnet_grouped==FALSE OR note=="dirichlet_fallback"
            is_fallback = (glmnet_grouped == FALSE) | (reshape_note_l == "dirichlet_fallback")
        )
    
    diag_df %>%
        mutate(reason = dplyr::case_when(
            reshape_note == "dirichlet_fallback" ~ "glmnet not fit (counts/sparsity)",
            glmnet_grouped == FALSE              ~ "glmnet fit invalid/failed",
            TRUE                                 ~ "ok"
        )) %>%
        select(garbage, n_tr, n_te, k_cand, nz_feat,
               glmnet_grouped, reshape_note, reason)
    
    
    # sanity peek: what got parsed?
    cat("\n[parse] glmnet_grouped value counts (after coercion):\n")
    print(diag_df %>% count(glmnet_grouped, name = "n"))
    cat("\n[parse] reshape_note samples:\n")
    print(diag_df %>% count(reshape_note, name = "n") %>% arrange(desc(n)) %>% head(10))
    
    overall <- diag_df %>%
        summarise(
            chunks = n(),
            chunks_fallback = sum(is_fallback, na.rm = TRUE),
            frac_chunks_fallback = ifelse(chunks > 0, chunks_fallback / chunks, NA_real_),
            rows = sum(n_te, na.rm = TRUE),
            rows_fallback = sum(ifelse(is_fallback, n_te, 0), na.rm = TRUE),
            frac_rows_fallback = ifelse(rows > 0, rows_fallback / rows, NA_real_)
        )
    
    by_bin <- diag_df %>%
        group_by(garbage) %>%
        summarise(
            chunks = n(),
            chunks_fallback = sum(is_fallback, na.rm = TRUE),
            frac_chunks_fallback = ifelse(chunks > 0, chunks_fallback / chunks, NA_real_),
            rows = sum(n_te, na.rm = TRUE),
            rows_fallback = sum(ifelse(is_fallback, n_te, 0), na.rm = TRUE),
            frac_rows_fallback = ifelse(rows > 0, rows_fallback / rows, NA_real_),
            .groups = "drop"
        ) %>%
        arrange(garbage)
    
    list(overall = overall, by_bin = by_bin, raw = diag_df)
}

read_diag_dir <- function(diag_dir) {
    if (!dir.exists(diag_dir)) stop("Directory does not exist: ", diag_dir)
    fs <- list.files(diag_dir, pattern = "^glmnet_prob_diag_.*\\.csv$", full.names = TRUE)
    if (!length(fs)) {
        stop("No diag CSVs found in ", diag_dir, " (expected files named like glmnet_prob_diag_*.csv).")
    }
    purrr::map_dfr(fs, readr::read_csv, show_col_types = FALSE)
}

print_fallback_summary <- function(res) {
    cat("\n================= Dirichlet Fallback Usage =================\n")
    cat("Overall (weighted by n_te rows):\n")
    print(res$overall %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
    cat("\nBy garbage group:\n")
    print(res$by_bin %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
    cat("============================================================\n\n")
}

write_fallback_summary <- function(res, out_prefix = NULL) {
    if (is.null(out_prefix) || !nzchar(out_prefix)) return(invisible(NULL))
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
    readr::write_csv(res$overall, paste0(out_prefix, "_overall.csv"))
    readr::write_csv(res$by_bin,  paste0(out_prefix, "_by_bin.csv"))
    invisible(NULL)
}

# ================ generate diagnostics if missing ============================
ensure_diagnostics_exist <- function() {
    if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
    has_any <- length(list.files(diag_dir, pattern = "^glmnet_prob_diag_.*\\.csv$")) > 0
    if (has_any) {
        message("[ok] Found existing diagnostics in: ", diag_dir)
        return(invisible(TRUE))
    }
    
    message("[info] No diagnostics found in ", diag_dir, " — generating a small run to create them...")
    # 1) source your helper from here("code/helpers")
    helper_path <- here("code", "helpers", "dq_entropy_helper.R")
    if (!file.exists(helper_path)) {
        stop("Helper not found at: ", helper_path, "\nPlease ensure dq_entropy_helper.R exists there.")
    }
    source(helper_path)
    if (!exists("compute_entropy_county_foreman")) {
        stop("compute_entropy_county_foreman() not defined after sourcing helper. Check helper file.")
    }
    
    # 2) find parquet
    pq_dir <- NULL
    for (cand in parquet_dir_candidates) {
        if (dir.exists(cand)) { pq_dir <- cand; break }
    }
    
    # 3) build ds (from parquet if possible, else synthetic tiny ds)
    if (!is.null(pq_dir)) {
        if (!requireNamespace("arrow", quietly = TRUE)) {
            stop("arrow package is required to read parquet. Please install.packages('arrow').")
        }
        files <- list.files(pq_dir, "\\.parquet$", full.names = TRUE)
        if (!length(files)) stop("No parquet files found in: ", pq_dir)
        
        schema_names <- names(arrow::read_parquet(files[1], as_data_frame = FALSE)$schema)
        cols_needed <- intersect(
            c("ucod", county_var, "year", paste0("record", 1:20), paste0("record_", 1:20), "age_years", "ager27"),
            schema_names
        )
        ds <- arrow::read_parquet(files[1], as_data_frame = TRUE, col_select = cols_needed)
        
        # normalize record column names to record_1..record_20
        if (!any(grepl("^record_", names(ds))) && any(grepl("^record[0-9]+$", names(ds)))) {
            ds <- dplyr::rename_with(ds, ~ sub("^record([0-9]+)$","record_\\1", .x), .cols = dplyr::matches("^record[0-9]+$"))
        }
        # add county if absent (dummy single county—ok for diag creation)
        if (!county_var %in% names(ds)) {
            ds[[county_var]] <- "dummy_county"
        }
        # sampling + cleaning
        if (nrow(ds) > sample_max_rows) {
            set.seed(42)
            ds <- ds %>% dplyr::slice_sample(n = sample_max_rows)
        }
        ds <- ds %>% dplyr::mutate(
            ucod = canonical_icd(ucod),
            dplyr::across(dplyr::starts_with("record_"), ~ canonical_icd(.x))
        )
    } else {
        # synthetic tiny ds so helper will still write a diag CSV
        message("[info] No parquet directories found; using a tiny synthetic dataset for diagnostics.")
        set.seed(123)
        n <- 2000
        ds <- tibble(
            ucod = sample(c("I219","J189","A419","X599","Y349","V899","G309"), n, TRUE) |> canonical_icd(),
            !!county_var := sample(sprintf("c%03d", 1:5), n, TRUE),
            year = sample(2015:2020, n, TRUE),
            record_1 = sample(c("I219","J189","A419","T404","Y349","S099"), n, TRUE) |> canonical_icd(),
            record_2 = sample(c("I219","J189","B349","W199","T509","G309"), n, TRUE) |> canonical_icd()
        )
        
        # minimal dictionaries if yours aren't present (so helper can run)
        dir.create(dictionary_dir, recursive = TRUE, showWarnings = FALSE)
        icd_map_path  <- file.path(dictionary_dir, "foreman-icd10-mapping.csv")
        table2_path   <- file.path(dictionary_dir, "foreman-table2-map.csv")
        if (!file.exists(icd_map_path)) {
            tibble::tibble(
                ICD10 = c("I219","J189","A419","X599","Y349","V899","G309","T404","T509","S099","W199","B349"),
                USCOD = c("C_heart","C_resp","C_sepsis","G_1","G_2","G_1","C_neuro","C_drug","C_drug","C_injury","C_injury","C_resp")
            ) %>% readr::write_csv(icd_map_path)
        }
        if (!file.exists(table2_path)) {
            tibble::tibble(
                target_cause = c("C_heart","C_resp","C_neuro","C_drug","C_injury"),
                G_1 = c("TRUE","TRUE","FALSE","TRUE","TRUE"),
                G_2 = c("TRUE","FALSE","TRUE","TRUE","FALSE"),
                G_3 = FALSE, G_4 = FALSE, G_5 = FALSE, G_6 = FALSE, G_7 = FALSE, G_8 = FALSE, G_9 = FALSE
            ) %>% readr::write_csv(table2_path)
        }
    }
    
    # 4) run helper to write diagnostics (exact path your helper uses)
    icd_map_path  <- file.path(dictionary_dir, "foreman-icd10-mapping.csv")
    table2_path   <- file.path(dictionary_dir, "foreman-table2-map.csv")
    if (!file.exists(icd_map_path) || !file.exists(table2_path)) {
        stop("Dictionary files not found. Expected:\n  ",
             icd_map_path, "\n  ", table2_path, "\nAdjust 'dictionary_dir' or create the files.")
    }
    
    options(DQ_DIAG_DIR = diag_dir, DQ_VERBOSE = TRUE, DQ_SEED = 123)
    invisible(compute_entropy_county_foreman(
        ds = ds,
        county_var = county_var,
        dict_dir = dictionary_dir,
        icd_map_path = icd_map_path,
        code_map_path = table2_path
    ))
    message("[ok] Diagnostics generated in: ", diag_dir)
    invisible(TRUE)
}

# ================================ RUN ========================================
ensure_diagnostics_exist()
diag_df <- read_diag_dir(diag_dir)
res <- summarize_fallback(diag_df)

cat("\n================= Dirichlet Fallback Usage =================\n")
cat("Overall (weighted by n_te rows):\n")
print(res$overall %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
cat("\nBy garbage group:\n")
print(res$by_bin %>% mutate(across(where(is.numeric), ~ round(.x, 6))))
cat("============================================================\n\n")

# Optional: write CSVs
if (!is.null(out_prefix) && nzchar(out_prefix)) {
    write_fallback_summary(res, out_prefix)
}
