#' Get results from ToppFun
#'
#' @description
#' The toppFun() function takes a data.frame or other tabular data structure and selects genes to use in querying
#' ToppGene.
#'
#' @details The use of data from ToppGene is governed by their Terms of Use: https://toppgene.cchmc.org/navigation/termsofuse.jsp
#'
#' @param input_data A vector of markers or dataframe with columns as cluster labels
#' @param type One of c("degs", "marker_list", or "marker_df). If "degs" is selected,
#' the input_data is assumed to be a data.frame with logfoldchange, pvalue, and gene name columns.
#' If "marker_list" is selected, input_data is assumed to be a list of genes with no other stats,
#' and any thresholds pertaining to "degs" will be ignored.
#' If "marker_df" is selected, the input_data is assumed to be a data.frame with columns as clusters/celltypes,
#' and entries are lists of markers.
#' @param topp_categories A string or vector with specific toppfun categories for the query
#' @param cluster_col Column name for the groups of cells (e.g. cluster or celltype)
#' @param gene_col Column name for genes (e.g. gene or feature)
#' @param logFC_col Column name for the avg log FC column
#' @param p_val_col Column name for the p-value or adjusted p-value (preferred)
#' @param num_genes Number of genes per group to use for toppGene query
#' @param key_type Gene name format
#' @param pval_cutoff (adjusted) P-value cutoff for filtering differentially expressed genes
#' @param fc_filter Include "ALL" genes, or only "UPREG" or "DOWNREG" for each cluster
#' @param fc_cutoff Avg log fold change cutoff for filtering differentially expressed genes
#' @param clusters Which clusters to include in toppGene query
#' @param min_genes Minimum number of genes to match in a query
#' @param max_genes Maximum number of genes to match in a query
#' @param max_results Maximum number of results per cluster
#' @param correction P-value correction method ("FDR" is "BH")
#' @param verbose Verbosity setting, TRUE or FALSE
#' @importFrom stringr str_glue str_c str_subset
#' @importFrom dplyr filter arrange select
#' @importFrom httr2 request req_body_json req_perform resp_body_json
#' @return data.frame
#' @examples
#' data("ifnb.de")
#' toppData <- toppFun(ifnb.de,
#'     topp_categories = NULL,
#'     cluster_col = "celltype",
#'     gene_col = "gene",
#'     p_val_col = "p_val_adj",
#'     logFC_col = "avg_log2FC")
#' @export
toppFun <- function(input_data,
                    type = "degs",
                    topp_categories = NULL,
                    cluster_col = "cluster",
                    gene_col = "gene",
                    p_val_col = "adj_p_val_col",
                    logFC_col = "avg_logFC",
                    num_genes = 1000,
                    pval_cutoff = 0.5,
                    fc_cutoff = 0,
                    fc_filter = "ALL",
                    clusters = NULL,
                    correction ="FDR",
                    key_type = "SYMBOL",
                    min_genes = 2,
                    max_genes = 1500,
                    max_results = 50,
                    verbose = TRUE
                    ) {

  #Print message about the use of ToppGene's data and adhering to their terms of use
  msg <- "This function returns data generated from ToppGene (https://toppgene.cchmc.org/)\n
Any use of this data must be done so under the Terms of Use and citation guide established by ToppGene.\n
Terms of Use: https://toppgene.cchmc.org/navigation/termsofuse.jsp
Citations: https://toppgene.cchmc.org/help/publications.jsp"
  message(msg)

  #parse for degs vs markers
  if (type == "degs") {
    input_data <- as.data.frame(input_data)

    #subset clusters if needed
    if (!(is.null(clusters))) {
      input_data <- input_data |>
        dplyr::filter(!!as.name(cluster_col) %in% clusters)
    }
    if (!(cluster_col %in% colnames(input_data))) {
      stop(paste0("Cluster column `", cluster_col, "` not found in data. Please specify."))
    }
    #parse fc_filter
    if (!(fc_filter %in% c("ALL", "UPREG", "DOWNREG"))){
      stop("please select one of c('ALL', 'UPREG', 'DOWNREG') for fc_filter")
    }
    gene_data <- process_markers(markers=input_data,
                                   cluster_col=cluster_col,
                                   gene_col=gene_col,
                                   p_val_col = p_val_col,
                                   logFC_col = logFC_col,
                                   num_genes=num_genes,
                                   pval_cutoff=pval_cutoff,
                                   fc_cutoff = fc_cutoff,
                                   genes_submit_cutoff=genes_submit_cutoff,
                                   fc_filter=fc_filter)
  } else if (type == "marker_list") {
    gene_data = data.frame(genes=input_data)
  } else if (type == "marker_df") {
    gene_data = input_data
  }

  #parse correction method
  if (!(correction %in% c("none", "FDR", "Bonferroni"))) {
    stop("invalid P-value correction method, please select either 'none', 'FDR', or 'Bonferroni'")
  }

  #parse categories
  if (is.null(topp_categories)) {
    topp_categories = get_ToppCats()
  }

  big_df <- data.frame()
  missing_clusters = c()

  for (col in names(gene_data)) {

    if (!(col %in% c('rank', 'X'))) {
      gene_list = gene_data[[col]]

      if (sum(!(is.na(gene_list))) >= min_genes) {
        if (verbose) {
          cat('Working on cluster:', col, '\n')
        }

        d <- get_topp(gene_list = gene_list,
                      topp_categories = topp_categories,
                      key_type = "SYMBOL",
                      pval_cutoff=pval_cutoff,
                      min_genes=min_genes,
                      max_genes=max_genes,
                      max_results=max_results,
                      correction=correction)

        if (nrow(d) == 0){
          missing_clusters = append(missing_clusters, col)
        } else if (length(names(gene_data)) == 1) {
          big_df <- rbind(big_df, d)
        } else {
          d[['Cluster']] = col
          big_df <- rbind(big_df, d)
        }
      } else {
        missing_clusters = append(missing_clusters, col)
      }
    }
  }
  #report any missing clusters
  if (length(missing_clusters) > 0) {
    write(stringr::str_glue("WARNING: no results found for clusters {stringr::str_c(missing_clusters, collapse=', ')}"),
          stderr())
  }
  return (big_df)
}

#' Convert genes into Entrez format
#'
#' @param genes A list of genes
#' @return a vector of genes in Entrez format
#' @importFrom httr2 request req_body_json req_perform
#' @examples
#' get_Entrez(genes=c("IFNG", "FOXP3"))
#' @export
get_Entrez<- function(genes){
  lookup_url = 'https://toppgene.cchmc.org/API/lookup'
  req <- httr2::request(lookup_url)
  if (length(genes) == 1) {
    resp <- req |>
      httr2::req_body_json(list(Symbols=list(genes))) |>
      httr2::req_perform()
  } else {
    resp <- req |>
      httr2::req_body_json(list(Symbols=genes)) |>
      httr2::req_perform()
  }
  new_gene_list <- c()
  for (g in httr2::resp_body_json(resp)$Genes) {
    new_gene_list <- base::append(new_gene_list, g[['Entrez']])
  }
  return (new_gene_list)
}

##### PROCESS MARKER INPUTS
process_markers <- function (markers, cluster_col, gene_col, p_val_col, logFC_col,
                             num_genes=1000,
                             pval_cutoff=0.5,
                             fc_cutoff=0,
                             fc_filter="ALL",
                             genes_submit_cutoff=1000) {


  #parse columns - tries to find columns if differ from the default TODO

  #Make a list of lists - each sub list has all of the genes (up to num_genes if specified) of the filtered data
  marker_list = list()

  for (cl in unique(markers[[cluster_col]])) {
    all_cl_markers <- markers |>
      dplyr::filter(!!as.name(cluster_col) == cl) |>
      dplyr::filter(!!as.name(p_val_col) < pval_cutoff)
    if (fc_filter == "ALL"){
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(abs(!!as.name(logFC_col)) > fc_cutoff) |>
        dplyr::arrange(-abs(!!as.name(logFC_col))) |>
        dplyr::select(!!as.name(gene_col))
    } else if (fc_filter == "UPREG") {
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(!!as.name(logFC_col) > fc_cutoff) |>
        dplyr::arrange(-!!as.name(logFC_col)) |>
        dplyr::select(!!as.name(gene_col))
    } else if (fc_filter == "DOWNREG") {
      all_cl_markers <- all_cl_markers |>
        dplyr::filter(!!as.name(logFC_col) < fc_cutoff) |>
        dplyr::arrange(!!as.name(logFC_col)) |>
        dplyr::select(!!as.name(gene_col))
    }
    if (!(is.null(num_genes))) {
      if (length(all_cl_markers[[gene_col]]) > num_genes) {
        marker_list[[cl]] <- all_cl_markers[1:num_genes, gene_col] |> unlist() |> as.character()
      } else {
        marker_list[[cl]] <- all_cl_markers[[gene_col]] |> unlist() |> as.character()
      }
    } else {
      marker_list[[cl]] <- all_cl_markers[[gene_col]] |> unlist() |> as.character()
    }
    }
  return (marker_list)
}

get_topp <- function(gene_list,
                      key_type,
                      topp_categories,
                      max_results=10,
                      min_genes=5,
                      max_genes=1500,
                      pval_cutoff=0.05,
                      correction="FDR") {

  #assertions - to add
  #print(gene_list)
  #convert gene names if necessary
  if (key_type != 'ENTREZ') {
    new_gene_list = get_Entrez(gene_list)
  } else{
    new_gene_list = gene_list
  }

  #create payload for POST request
  payload = list()
  payload[['Genes']] = new_gene_list
  category_list = list()
  for (i in 1:length(topp_categories)) {
    cat_dict = list(Type=topp_categories[i],
                    PValue=pval_cutoff,
                    MinGenes=min_genes,
                    MaxGenes=max_genes,
                    MaxResults=max_results,
                    Correction=correction)
    category_list[[i]] <- cat_dict
  }

  #send POST request
  url = 'https://toppgene.cchmc.org/API/enrich'
  req <- request(url)
  resp <- req |>
    httr2::req_body_json(list(Genes=new_gene_list, Categories=category_list)) |>
    httr2::req_perform()

  response_data <- httr2::resp_body_json(resp)[['Annotations']]
  keepers <- c("Category","ID","Name","PValue","QValueFDRBH","QValueFDRBY","QValueBonferroni",
               "TotalGenes","GenesInTerm","GenesInQuery","GenesInTermInQuery","Source","URL")
  results_df <- data.frame()
  for (i in 1:length(response_data)) {
    results_df <- rbind(results_df,response_data[[i]][keepers])
  }
  return (results_df)
}


#' Get a vector of ToppFun categories
#'
#' @return a vector
#' @examples
#' get_ToppCats()
#' @export
get_ToppCats <- function() {
  toppCats <- c("GeneOntologyMolecularFunction",
                "GeneOntologyBiologicalProcess",
                "GeneOntologyCellularComponent",
                "HumanPheno",
                "MousePheno",
                "Domain",
                "Pathway",
                "Pubmed",
                "Interaction",
                "Cytoband",
                "TFBS",
                "GeneFamily",
                "Coexpression",
                "CoexpressionAtlas",
                "ToppCell",
                "Computational",
                "MicroRNA",
                "Drug",
                "Disease")
  return (toppCats)
}

#' Save toppData results (optionally) split by celltype/cluster
#'
#' @param toppData Results from toppFun as a dataframe
#' @param filename filename prefix for each split file
#' @param save_dir the directory to save files
#' @param split Boolean, whether to split the dataframe by celltype/cluster
#' @param format Saved file format, one of c("xlsx", "csv", "tsv")
#' @returns A saved file
#' @importFrom openxlsx write.xlsx
#' @importFrom stringr str_glue
#' @importFrom dplyr filter
#' @importFrom utils write.table
#' @examples
#' data("toppData")
#' toppSave(toppData, filename="toppFun_results", split = TRUE, format = "xlsx")
#' @export
toppSave <- function (toppData,
                      filename = NULL,
                      save_dir = NULL,
                      split = TRUE,
                      format = "xlsx") {
  if (is.null(save_dir)) {
    save_dir = getwd()
  }
  if (isTRUE(split)) {
    for (gr in unique(toppData$Cluster)) {
      tmp_toppData <- toppData |>
        dplyr::filter(Cluster == gr)

      if (format == 'xlsx') {

        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData_{gr}.xlsx")
        } else {
          this_file = stringr::str_glue("{filename}_{gr}.xlsx")
        }

        openxlsx::write.xlsx(tmp_toppData,
                             file = file.path(save_dir, this_file),
                             colNames = TRUE,
                             rowNames = FALSE,
                             borders = "columns",
                             sheetName="toppData",
                             overwrite = TRUE)
        cat("Saving file:", this_file, "\n")

      } else if (format == "csv") {
        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData_{gr}.csv")
        } else {
          this_file = stringr::str_glue("{filename}_{gr}.csv")
        }

        utils::write.table(tmp_toppData,
                  file = file.path(save_dir, this_file),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      } else if (format == "tsv") {
      if (is.null(filename)) {
        this_file = stringr::str_glue("toppData_{gr}.tsv")
      } else {
        this_file = stringr::str_glue("{filename}_{gr}.tsv")
      }

        utils::write.table(tmp_toppData,
                  file = file.path(save_dir, this_file),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      }
    }
  } else {
      if (format == 'xlsx') {

        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData.xlsx")
        } else {
          this_file = stringr::str_glue("{filename}.xlsx")
        }

        openxlsx::write.xlsx(toppData,
                             file = file.path(save_dir, this_file),
                             colNames = TRUE,
                             rowNames = FALSE,
                             borders = "columns",
                             sheetName="toppData",
                             overwrite = TRUE)
        cat("Saving file:", this_file, "\n")

      } else if (format == "csv") {
        if (is.null(filename)) {
          this_file = stringr::str_glue("toppData.csv")
        } else {
          this_file = stringr::str_glue("{filename}.csv")
        }

        utils::write.table(toppData,
                  file = file.path(save_dir, this_file),
                  sep = ",",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
        cat("Saving file:", this_file, "\n")
      } else if (format == "tsv") {
      if (is.null(filename)) {
        this_file = stringr::str_glue("toppData.tsv")
      } else {
        this_file = stringr::str_glue("{filename}.tsv")
      }

      utils::write.table(toppData,
                  file = file.path(save_dir, this_file),
                  sep = "\t",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
      cat("Saving file:", this_file, "\n")
    }

  }

}
