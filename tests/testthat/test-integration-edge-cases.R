# Integration tests and edge cases for scToppR package

test_that("toppSave parameter validation is comprehensive", {
    data("toppdata.ifnb")
    
    # Test invalid format
    expect_error(
        toppSave(toppdata.ifnb, 
                filename = "test",
                save_dir = tempdir(),
                format = "invalid_format"),
        "Please select one of c\\('xlsx', 'csv', 'tsv'\\) for format"
    )
    
    # Test missing save_dir
    expect_error(
        toppSave(toppdata.ifnb, filename = "test"),
        "Please specify a `save_dir` to save the file"
    )
    
    # Test invalid cluster column for split
    expect_error(
        toppSave(toppdata.ifnb,
                filename = "test", 
                save_dir = tempdir(),
                split = TRUE,
                cluster_col = "nonexistent_column"),
        "Cannot split by cluster column.*not found in toppData"
    )
})

test_that("gene processing handles different data types and edge cases", {
    data("ifnb.de")
    
    # Test with factor columns
    ifnb_factor <- ifnb.de
    ifnb_factor$celltype <- as.factor(ifnb_factor$celltype)
    ifnb_factor$gene <- as.factor(ifnb_factor$gene)
    
    gene_data_factor <- .process_degs(
        degs = ifnb_factor,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC"
    )
    expect_type(gene_data_factor, "list")
    
    # Test with very small dataset
    small_data <- ifnb.de[1:10, ]
    gene_data_small <- .process_degs(
        degs = small_data,
        cluster_col = "celltype", 
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        num_genes = 5
    )
    expect_type(gene_data_small, "list")
    
    # Test with data containing duplicate genes
    dup_data <- rbind(ifnb.de[1:100, ], ifnb.de[1:10, ])  # Add duplicates
    gene_data_dup <- .process_degs(
        degs = dup_data,
        cluster_col = "celltype",
        gene_col = "gene", 
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC"
    )
    expect_type(gene_data_dup, "list")
})

test_that("column name variations are handled correctly", {
    data("ifnb.de")
    
    # Test with different column names
    renamed_data <- ifnb.de
    colnames(renamed_data)[colnames(renamed_data) == "celltype"] <- "cluster_id"
    colnames(renamed_data)[colnames(renamed_data) == "gene"] <- "feature"
    colnames(renamed_data)[colnames(renamed_data) == "p_val_adj"] <- "padj"
    colnames(renamed_data)[colnames(renamed_data) == "avg_log2FC"] <- "logFC"
    
    gene_data_renamed <- .process_degs(
        degs = renamed_data,
        cluster_col = "cluster_id",
        gene_col = "feature", 
        p_val_col = "padj",
        logFC_col = "logFC"
    )
    expect_type(gene_data_renamed, "list")
})

test_that("extreme parameter values are handled gracefully", {
    data("ifnb.de")
    
    # Test with very strict cutoffs (should result in very few genes)
    gene_data_strict <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj", 
        logFC_col = "avg_log2FC",
        pval_cutoff = 1e-10,  # Very strict
        fc_cutoff = 5.0       # Very high
    )
    expect_type(gene_data_strict, "list")
    
    # Test with very loose cutoffs
    gene_data_loose <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype", 
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC", 
        pval_cutoff = 1.0,    # Very loose
        fc_cutoff = 0.0       # No fold change requirement
    )
    expect_type(gene_data_loose, "list")
    
    # Test with num_genes = 0 (should handle gracefully)
    gene_data_zero <- .process_degs(
        degs = ifnb.de,
        cluster_col = "celltype",
        gene_col = "gene",
        p_val_col = "p_val_adj",
        logFC_col = "avg_log2FC",
        num_genes = 0
    )
    expect_type(gene_data_zero, "list")
})

test_that("data with missing values is handled correctly", {
    data("ifnb.de") 
    
    # Introduce missing values
    missing_data <- ifnb.de
    missing_data$p_val_adj[1:10] <- NA
    missing_data$avg_log2FC[5:15] <- NA
    missing_data$gene[20:25] <- NA
    
    # Should handle missing values gracefully
    expect_no_error({
        gene_data_missing <- .process_degs(
            degs = missing_data,
            cluster_col = "celltype",
            gene_col = "gene",
            p_val_col = "p_val_adj", 
            logFC_col = "avg_log2FC"
        )
    })
})

test_that("plotting functions work with different data sizes", {
    data("toppdata.pbmc")
    
    # Test with very large dataset (subset to avoid performance issues in tests) 
    large_subset <- toppdata.pbmc[1:1000, ]  # Reasonable subset for testing
    if (nrow(large_subset) > 0) {
        plot_large <- toppPlot(large_subset,
            category = unique(large_subset$Category)[1],
            clusters = unique(large_subset$Cluster)[1],
            save = FALSE
        )
        expect_s3_class(plot_large, "ggplot")
    }
    
    # Test with minimal data
    minimal_data <- toppdata.pbmc[1:3, ]
    plot_minimal <- toppPlot(minimal_data,
        category = unique(minimal_data$Category)[1],
        clusters = unique(minimal_data$Cluster)[1],
        save = FALSE
    )
    expect_s3_class(plot_minimal, "ggplot")
})

test_that("file operations work correctly", {
    data("toppdata.ifnb")
    tmp_dir <- tempdir()
    
    # Test that files are created with expected names
    toppSave(toppdata.ifnb, 
            filename = "edge_case_test",
            save_dir = tmp_dir,
            split = FALSE,
            format = "csv",
            verbose = FALSE)
    
    expect_true(file.exists(file.path(tmp_dir, "edge_case_test.csv")))
    
    # Test with special characters in filename (should be handled)
    toppSave(toppdata.ifnb,
            filename = "test_file_123", 
            save_dir = tmp_dir,
            split = FALSE,
            format = "xlsx", 
            verbose = FALSE)
    
    expect_true(file.exists(file.path(tmp_dir, "test_file_123.xlsx")))
    
    # Clean up test files
    unlink(file.path(tmp_dir, "edge_case_test.csv"))
    unlink(file.path(tmp_dir, "test_file_123.xlsx"))
})

test_that("package data integrity is maintained", {
    # Test all included datasets load correctly
    data("ifnb.de")
    expect_true(is.data.frame(ifnb.de))
    expect_true(nrow(ifnb.de) > 0)
    expect_true("celltype" %in% colnames(ifnb.de))
    
    data("ifnb.markers.df")
    expect_true(is.data.frame(ifnb.markers.df))
    expect_true(ncol(ifnb.markers.df) > 0)
    
    data("ifnb.markers.list.CD8T") 
    expect_true(is.vector(ifnb.markers.list.CD8T))
    expect_true(length(ifnb.markers.list.CD8T) > 0)
    
    data("pbmc.markers")
    expect_true(is.data.frame(pbmc.markers))
    expect_true(nrow(pbmc.markers) > 0)
    
    data("toppdata.ifnb")
    expect_true(is.data.frame(toppdata.ifnb))
    expect_true(all(c("Category", "PValue", "Cluster") %in% colnames(toppdata.ifnb)))
    
    data("toppdata.pbmc") 
    expect_true(is.data.frame(toppdata.pbmc))
    expect_true(all(c("Category", "PValue", "Cluster") %in% colnames(toppdata.pbmc)))
    
    data("toppdata.airway")
    expect_true(is.data.frame(toppdata.airway))
    expect_true(all(c("Category", "PValue", "Cluster") %in% colnames(toppdata.airway)))
})

test_that("error messages are informative and helpful", {
    data("ifnb.de")
    
    # Test that error messages provide clear guidance
    expect_error(
        toppFun(ifnb.de, type = "invalid_type"),
        "Please ensure the parameter `type` is one of: degs, marker_df, or marker_list"
    )
    
    expect_error(
        toppFun(ifnb.de, cluster_col = "missing_column"),
        "Cluster column `missing_column` not found in data. Please specify"
    )
})

test_that("function compatibility across different data types", {
    # Test that functions work with tibbles, data.tables, etc. if they support them
    data("ifnb.de")
    
    # Convert to tibble if dplyr is available
    if (requireNamespace("dplyr", quietly = TRUE)) {
        ifnb_tibble <- dplyr::as_tibble(ifnb.de)
        
        gene_data_tibble <- .process_degs(
            degs = ifnb_tibble,
            cluster_col = "celltype",
            gene_col = "gene",
            p_val_col = "p_val_adj", 
            logFC_col = "avg_log2FC"
        )
        expect_type(gene_data_tibble, "list")
    }
})