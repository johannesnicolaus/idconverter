
#' Convert ID using biomaRt (ensembl database)
#'
#' @param genelist List of gene IDs in the form of character vector.
#' @param organism Organism of the input gene ID (can be checked using check_organisms_biomart()).
#' @param supplied_id The key of the gene ID supplied (e.g. "ensembl_gene_id" or
#' "external_gene_name"). Can be checked using check_attributes_biomart().
#' @param attribs Gene ID to convert the input ID to. Can be checked using 
#' check_attributes_biomart().
#' @param host Host server for Ensembl database, closest server will automatically
#' selected. However if there is an error, it can be set to the following:
#' "uswest.ensembl.org", "asia.ensembl.org", "useast.ensembl.org".
#' @param unique_rows Removes duplicated rows. 
#' @param replace_extgenename_with_ensembl When converting ensemblid and there 
#' is no corresponding gene name, replace with ensembl gene id.
#' @param remove_geneid_version Remove gene id version from ensembl ID (e.g. ENSGALG00000025996.4 to ENSGALG00000025996).
#'
#' @return data frame consisting of gene ids.
#' @export
#'
#' @examples convert_biomart(c("ENSGALG00000025996", "ENSGALG00000017312.5"), organism = "ggallus")
convert_biomart <- function(genelist,
                            organism = "hsapiens",
                            supplied_id = "ensembl_gene_id",
                            attribs = c('ensembl_gene_id', 'external_gene_name', 'description', 'mmusculus_homolog_ensembl_gene'),
                            host = "www.ensembl.org", 
                            unique_rows = FALSE,
                            replace_extgenename_with_ensembl = T, 
                            remove_geneid_version = T){
  
  if(remove_geneid_version == T){
    genelist <- gsub(x = genelist, pattern = "\\..*", replacement = "")
  }
  
  speciesmart <- biomaRt::useMart("ensembl", dataset = paste0(organism ,"_gene_ensembl"), host = host)
  
  converted_genes <- biomaRt::getBM(attributes=attribs, 
                           filters = supplied_id, 
                           values = genelist, 
                           mart = speciesmart,
                           uniqueRows = unique_rows
  )
  
  # fill blanks of external gene name with ensembl
  if (replace_extgenename_with_ensembl == T) {
    converted_genes$external_gene_name[converted_genes$external_gene_name == ""] <- converted_genes$ensembl_gene_id[converted_genes$external_gene_name == ""]
  }
  
  return(converted_genes)
  
}



#' Check organism biomart (Ensembl database)
#'
#' @return List of organism that are available on the Ensembl database.
#' @export
#'
#' @examples check_organisms_biomart()
check_organisms_biomart <- function(){
  mart <- biomaRt::useMart('ensembl')
  biomaRt::listDatasets(mart)
}


#' Check available attributes (key types) for selected organism
#'
#' @param organism Organism to check attributes from. Can be checked using check_organisms_biomart()
#' @param host Host server for Ensembl database, closest server will automatically
#' selected. However if there is an error, it can be set to the following:
#'
#' @return List of attributes()
#' @export
#'
#' @examples check_attributes_biomart("ggallus", host = "asia.ensembl.org")
check_attributes_biomart <- function(organism, host = "www.ensembl.org"){
  speciesmart <- biomaRt::useMart("ensembl", dataset = paste0(organism ,"_gene_ensembl"), host = host)
  return(biomaRt::listAttributes(speciesmart))
  
}
