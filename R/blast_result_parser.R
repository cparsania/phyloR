
#' Wrapper function around taxize::genbank2uid.
#'
#' Given a genBank accession alphanumeric string or a gi numeric string \code{(x)}, it returns tibble of taxid, name and other related columns.
#' @param x vector of genBank accession alphanumeric string, or a gi numeric string \code{(x)}.
#' @param ... other parameters to be passed to \code{taxize::genbank2uid}
#'
#' @return a tbl with colnames x, taxid, class, match, multiple_matches, pattern_match, uri, name
#' @export
#' @importFrom taxize genbank2uid
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @importFrom purrr map_df
#' @examples
#' \dontrun{
#' x <- c("XP_022900619.1", "XP_022900618.1", "XP_018333511.1", "XP_018573075.1")
#' genbank2uid_tbl(x = x)
#' }
genbank2uid_tbl <- function(x , ...){

        start_time <- lubridate::now()
        uid_list <- taxize::genbank2uid(x ,  ...)
        uid_tbl <- tibble::tibble(x = x, taxid = unlist(uid_list)) %>%
                dplyr::bind_cols( purrr::map_df(uid_list , attributes))
        time_taken <- start_time - lubridate::now()
        cat_green_tick("done. ", " Time taken " , time_taken)
        return(uid_tbl)

}




#' Get phylogenetic rank for a given ncbi taxonomy id.
#'
#'
#' @param x a vector of valid ncbi taxonomy id.
#' @param rank a character string indicating required phylogenetic rank. Default "kingdom". Can be one of the following
#' \enumerate{
#' \item superkingdom
#' \item kingdom
#' \item phylum
#' \item subphylum
#' \item class
#' \item subclass
#' \item infraclass
#' \item cohort
#' \item order
#' \item suborder
#' \item infraorder
#' \item superfamily
#' \item family
#' \item subfamily
#' \item genus
#' \item species
#' \item tribe
#' \item no rank
#' }
#' @return a tbl with the below column.
#' \enumerate{
#' \item query_taxon
#' \item kingdom
#' \item kingdom_id
#' \item rank
#' }
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom rlang arg_match
#' @importFrom taxizedb classification
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @importFrom lubridate now
#' @examples
#' \dontrun{
#' x <- c(166361, 166361, 224129, 217634)
#' get_taxon_rank(x)
#' }
#'
get_taxon_rank <-  function(x , rank = "kingdom"){
        #x <- tt
        x <- tibble::tibble(query_taxon = unique(x) %>% as.character())
        rlang::arg_match(rank  , c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"))

        if(length(rank) !=1 ) stop("argument `rank` must be of length 1")

        start_time <- lubridate::now()
        cli::cat_bullet("Starting rank search...",col = "red")
        result_ranks <- taxizedb::classification(x$query_taxon)  %>%
                tibble::tibble(taxid = names(.) , all_ranks = .) %>%
                tidyr::unnest(cols = all_ranks) %>%
                dplyr::filter(.$rank == !!rank)
        time_taken <- start_time - lubridate::now()
        cat_green_tick("done. ", " Time taken " , time_taken)
        rank_id_col_name <- paste(rank , "id" ,sep = "_")
        x %>% dplyr::left_join(result_ranks, by = c("query_taxon" = "taxid")) %>%
                ## rename cols
                dplyr::rename(!!rank_id_col_name := id , !!rank := name) %>%
                ## arrange cols
                dplyr::select(query_taxon,!!rank ,!!rank_id_col_name,rank)

}


#' Filter blast hits
#'
#' Given a tbl of blast output format 7 hits can be filtered by \code{evalue}, \code{bit_score}, \code{query_cov}, \code{identity} and  \code{query_length}
#'
#' @param blast_tbl an object of class tbl from blast output format 7. It must contains below columns.
#' \enumerate{
#' \item query_acc_ver
#' \item subject_acc_ver
#' \item identity
#' \item alignment_length
#' \item mismatches
#' \item gap_opens
#' \item q_start
#' \item q_end
#' \item s_start
#' \item s_end
#' \item evalue
#' \item bit_score
#' \item positives
#' }
#'
#' Column names can be assigned using the function \code{get_blast_outformat_7_colnames()}
#'
#' @param evalue a numeric value; default NULL. Only hits having evalue <= than given will be retained.
#' @param bit_score a numeric value; default NULL. Only hits having bit_score >= than given will be retained.
#' @param query_cov a numeric value between 0 and 100 ; default NULL. Only hits having query_cov >= than given will be retained. query_cov will be calculated from the column \code{q_start} \code{q_end} and param \code{query_length}
#' @param identity a numeric value between 0 and 100;default NULL. Only hits having identity >= than given will be retained.
#' @param query_length a numeric value indicating length of query sequence.
#'
#' @return a tbl of filtered blast hits.
#' @importFrom tibble is_tibble
#' @importFrom dplyr filter
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' \dontrun{
#'  f <- system.file("extdata","example_blast_output_01.txt", package = "phyloR")
#'  d <- readr::read_delim(f, delim ="\t" , col_names = F)
#'  colnames(d) <- phyloR::get_blast_outformat_7_colnames()
#'  filtered <- filter_blast_hits(d , query_cov = 90 , identity = 40, query_length = 249, evalue = 1e-6)
#' }
#' @seealso \link{get_blast_outformat_7_colnames}
#'
filter_blast_hits <- function(blast_tbl ,
                              evalue = NULL,
                              bit_score = NULL,
                              query_cov = NULL,
                              identity = NULL,
                              query_length = NULL ){

        #validate args
        stopifnot(tibble::is_tibble(blast_tbl))
        null_or_numeric(evalue)
        null_or_numeric(bit_score)
        null_or_numeric(query_cov)
        null_or_numeric(identity)
        null_or_numeric(query_length)


        blast_tbl_cols <- colnames(blast_tbl)

        ## filter by evalue
        if(!is.null(evalue)){
                if( ! ("evalue" %in% blast_tbl_cols) ) {
                        stop("blast_tbl must contain column 'evalue'")
                }

                blast_tbl <- blast_tbl %>%
                        dplyr::filter(evalue <= !!evalue)
        }

        ## filter by bitscore
        if(!is.null(bit_score)){
                if( ! ("bit_score" %in% blast_tbl_cols) ) {
                        stop("blast_tbl must contain column 'bit_score'")
                }

                blast_tbl <- blast_tbl %>%
                        dplyr::filter(bit_score >= !!bit_score)
        }

        ## filter by query_cov
        if(!is.null(query_cov)){
                if( ! ("q_start" %in% blast_tbl_cols) ) {
                        stop("blast_tbl must contain column 'q_start'")
                }

                if( ! ("q_end" %in% blast_tbl_cols)) {
                        stop("blast_tbl must contain column 'q_end'")
                }

                if(is.null(query_length)) {
                        stop("query_length can not be NULL when query_cov cutoff is applied.")
                }

                blast_tbl <- blast_tbl %>%
                        dplyr::filter(get_query_cov(qstart = q_start , qend = q_end ,qlen = query_length) >= !!query_cov)
        }

        ## filter by identity
        if(!is.null(identity)){
                if( ! ("identity" %in% blast_tbl_cols) ) {
                        stop("blast_tbl must contain column 'identity'")
                }

                blast_tbl <- blast_tbl %>%
                        dplyr::filter(identity >= !!identity)
        }

        return(blast_tbl)

}


#' Get blast output format 7 column names
#'
#' @return A character vector with values below
#' \itemize{
#' \item query_acc_ver
#' \item subject_acc_ver
#' \item identity
#' \item alignment_length
#' \item mismatches
#' \item gap_opens
#' \item q_start
#' \item q_end
#' \item s_start
#' \item s_end
#' \item evalue
#' \item bit_score
#' \item positives
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' get_blast_outformat_7_colnames()
#' }
#'
#'
get_blast_outformat_7_colnames <- function(){

        c("query_acc_ver",
          "subject_acc_ver",
          "identity",
          "alignment_length",
          "mismatches",
          "gap_opens",
          "q_start",
          "q_end",
          "s_start",
          "s_end",
          "evalue",
          "bit_score",
          "positives")

}



#' Calculate sequence coverage
#'
#'
#' @param start A numeric denoting start
#' @param end A numeric denoting end
#' @param len A numeric denoting length
#' @seealso get_query_cov
#' @seealso get_subj_cov
#' @keywords internal
#' @return A numeric
#'
get_cov <- function(start , end , len){

        ( (end - start +1 )/ len) * 100

}


#' Calculate \% query coverage.
#'
#' Given a sequence start, end and length coverage can be calculated.
#'
#' @param qstart A numeric denoting query start
#' @param qend A numeric denoting query end
#' @param qlen A numeric denoting query length
#' @export
#' @seealso \link{get_subj_cov}
#' @return A numeric
#'
get_query_cov <- function(qstart , qend , qlen){
        get_cov(start = qstart , end = qend , len = qlen)
}


#' Calculate \% subject coverage.
#'
#' Given a sequence start, end and length coverage can be calculated.
#'
#' @param sstart A numeric denoting subject start
#' @param send A numeric denoting subject end
#' @param slen A numeric denoting subject length
#' @export
#' @seealso \link{get_query_cov}
#' @return A numeric
#'
get_subj_cov <- function(sstart , send , slen){
        get_cov(start = sstart , end = send , len = slen)

}



#' Remove redundant hits from blast tabular output
#'
#' @description Given the output of blast tabular format in a tbl, it removes redundancy subject hit accession. For each redundant hit longest hit will be returned.
#'
#' @param blast_output_tbl an object of class tbl containing blast tabular output
#' @param subject_acc_colname a string denoting a column of subject hits in a \code{blast_output_tbl}. Default "subject_acc_ver"
#' @param alignment_length_colname a string denoting a column of alignment length in a \code{blast_output_tbl}. Default "alignment_length"
#'
#' @return a tbl
#' @export
#' @importFrom tibble is_tibble
#' @importFrom rlang sym
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @examples
remove_redundant_hits <-  function(blast_output_tbl ,
                               subject_acc_colname = "subject_acc_ver" ,
                               alignment_length_colname = "alignment_length"){

        if(!tibble::is_tibble(blast_output_tbl)){
                stop("blast_output_tbl must be a class of tibble")
        }

        subject_acc_colname = rlang::sym(subject_acc_colname)
        alignment_length_colname = rlang::sym(alignment_length_colname)

        blast_output_tbl %>%
                dplyr::group_by(!!subject_acc_colname) %>%
                dplyr::arrange(!!alignment_length_colname) %>%
                dplyr::slice(1)

}




## given the blast hits in tabular format add taxonomy annotations

#' Add taxonomy annotations to blast tabular output
#'
#' @description Given the blast output (\code{-outfmt 7}) in a tabular format add taxonomic annotations to each subject hit.
#' @param blast_output_tbl a tbl of blast output format 7 (\code{-outfmt 7})
#' @param subject_acc_colname a string (default : "subject_acc_ver")denoting column name of subject_acc_colname.
#' @param ncbi_acc_key user specific ENTREZ api key. Get one via \code{taxize::use_entrez()}
#' @param taxonomy_level a string indicating level of taxonomy to be mapped. Can be one of the following
#' \enumerate{
#' \item superkingdom
#' \item kingdom
#' \item phylum
#' \item subphylum
#' \item class
#' \item subclass
#' \item infraclass
#' \item cohort
#' \item order
#' \item suborder
#' \item infraorder
#' \item superfamily
#' \item family
#' \item subfamily
#' \item genus
#' \item species
#' \item tribe
#' \item no rank
#' }
#' @param map_superkindom logical (default TRUE). map superkingdom if kingdom not found. Valid only when taxonomy_level == "kingdom".
#'
#' @return
#' @export
#' @importFrom tibble is_tibble
#' @importFrom rlang arg_match
#' @importFrom cli cat_rule cli_alert_info cli_alert
#' @importFrom dplyr pull select bind_rows select rename left_join right_join
#' @importFrom purrr map
#' @importFrom TidyWrappers tbl_keep_rows_NA_any
#' @examples
add_taxonomy_columns <- function(blast_output_tbl,
                                 subject_acc_colname = "subject_acc_ver",
                                 ncbi_acc_key =NULL,
                                 taxonomy_level = "kingdom",
                                 map_superkindom = TRUE){

        ## validate user inputs
        ## blast_output_tbl must be tbl
        if(!tibble::is_tibble(blast_output_tbl)){
                stop("blast_output_tbl must be a class of tibble")
        }

        ## taxonomy_level must be one one of these
        rlang::arg_match(taxonomy_level  ,
                         c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"))


        ## blast_output_tbl mus not contain column names 'taxid' and taxonomy_level
        if(any(colnames(blast_output_tbl) %in% c(taxonomy_level))) {

                stop("Column '", taxonomy_level ,"'" , " must not present in 'blast_output_tbl'.",
                     " Either remove it or rename it")

        }

        is_taxid_present <- FALSE
        if(any(colnames(blast_output_tbl) %in% "taxid")) {
                is_taxid_present <- TRUE
                cli::cat_rule("WARNING")
                cli::cli_alert_info("As column 'taxid' present in blast output tbl, same will be used to map taxonomy level.")
                cli::cli_alert_info("To perform new 'taxid' search either remove  or rename columnn 'taxid'." )
                cli::cat_rule("WARNING ENDS")
        }

        ## get subject accession
        subject_acc_id <- blast_output_tbl %>%
                dplyr::pull(subject_acc_ver) %>%
                unique()


        if(is_taxid_present){
                taxon_data_sub <- blast_output_tbl %>%
                        dplyr::select(subject_acc_ver , taxid)
        } else{
                ## split subject accession vector in chunks
                ## splitting long vector of subject id in chunks is necessary as ping many ids at a time to NCBI frequently throws error
                chunk_max_length = 20 ## each chunk can be maximum of length "chunk_max_length"
                chunk_length <-   length(subject_acc_id) / chunk_max_length
                subject_acc_id_splts <- split(subject_acc_id , f = rep(1:chunk_length, length.out = length(subject_acc_id)))

                cli::cli_rule("'taxid' mapping starts")
                cli::cli_alert(paste("Number of iterations to take place : "  , chunk_length ,sep = " "))

                ## get taxid for each subject hit
                taxon_data <- purrr::map(subject_acc_id_splts ,
                                         ~ phyloR::genbank2uid_tbl(..1,key = ncbi_acc_key) )

                cat_green_tick("Mapping done.")

                ## list of tbl to tbl
                taxon_data_sub <- taxon_data %>%
                        dplyr::bind_rows() %>%
                        dplyr::select(1,2) %>%
                        dplyr::rename("subject_acc_ver" = "x")
        }


        ## get kingdom
        taxid_to_taxon <- phyloR::get_taxon_rank(taxon_data_sub$taxid , rank = taxonomy_level) %>%
                dplyr::select(1,2)

        ##if kingdom is NA map superkingdom
        if(map_superkindom && taxonomy_level=="kingdom"){
                if(nrow(taxid_to_taxon %>% TidyWrappers::tbl_keep_rows_NA_any()) >=1 ){
                        skingdom_query <- taxid_to_taxon %>% TidyWrappers::tbl_keep_rows_NA_any() %>% dplyr::pull(query_taxon)
                        taxid_to_skingdom <- phyloR::get_taxon_rank(skingdom_query , rank = "superkingdom") %>% dplyr::select(1,2)
                        taxid_to_taxon <- taxid_to_taxon %>% dplyr::left_join(taxid_to_skingdom , by = "query_taxon") %>%
                                dplyr::mutate(kingdom = dplyr::coalesce(kingdom, superkingdom)) %>% dplyr::select(1,2)
                }
        }

        ## add column taxid if it is not present in original data
        if(!is_taxid_present){
                blast_output_tbl <- blast_output_tbl %>%
                        dplyr::left_join(taxon_data_sub , by = c("subject_acc_ver"))
        }

        ## add user asked taxonomy_level
        blast_output_tbl <- blast_output_tbl %>%
                dplyr::left_join(taxid_to_taxon , by = c("taxid" = "query_taxon"))

        return(blast_output_tbl)
}


subset_fasta <- function(fasta_file){

}


#' format fasta headers
#' @description format the sequence headers for fasta file obtained from NCBI blast output
#' @param fasta_file string denoting full path of a fasta file.
#' @param keep_alignemnt_coord logical (default : TRUE) decides whether alignment coordinates to keep in headers or not. When TRUE header must contains alignment coordinates in this format \code{:([:digit:]+-[:digit:]+)}. (e.g. CEJ90625.1:1-252)
#' @importFrom Biostrings readBStringSet
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @importFrom stringr str_match_all
#' @importFrom tidyr unnest_wider
#' @importFrom dplyr select rename mutate_all mutate
#' @importFrom stringr str_match str_trim str_remove str_replace_all
#' @importFrom glue glue_data
#' @return sequences as an object of class \code{Biostrings}
#' @export
#'
#' @examples
format_fasta_headers <- function(fasta_file = NULL, keep_alignemnt_coord= TRUE){

        ## validate inputs
        if(!file.exists(fasta_file)){
            stop("fasta_file does not exist.")
        }

        fa_seq <- Biostrings::readBStringSet(fasta_file)
        fa_headers <- names(fa_seq)

        fa_headers_dislocate <- fa_headers %>% tibble::tibble(fa_headers = .) %>%
                dplyr::mutate(header_elems = purrr::map(fa_headers ,
                                                        ~( stringr::str_match_all(string = ..1 , pattern = "([^\\s]+)\\s(.*)\\[(.*)\\]") %>%
                                                                                unlist() ))) %>%
                tidyr::unnest_wider(col = header_elems) %>%
                dplyr::select(-...1) %>%
                dplyr::rename(subject_id = "...2" ,
                              desc = "...3",
                              species = "...4") %>%
                dplyr::mutate_all(stringr::str_trim) %>%
                dplyr::mutate(align_coord = stringr::str_match(string = subject_id,".*:([:digit:]+-[:digit:]+)") %>% .[,2] ) %>% ## add alignment coord column
                dplyr::mutate(subject_id = stringr::str_remove(string = subject_id,":.*"))  ## remove alignment coord from subject id

        if(keep_alignemnt_coord){
                fa_new_headers <-fa_headers_dislocate %>% glue::glue_data("{subject_id}__{align_coord}__{desc}__{species}")
        } else {
                fa_new_headers <-fa_headers_dislocate %>% glue::glue_data("{subject_id}__{desc}__{species}")
        }

        ## Replace any special char with '_'
        fa_new_headers <- fa_new_headers %>% stringr::str_replace_all(pattern = "\\W+" , "_") ## replace special char _

        names(fa_seq) <- fa_new_headers
        return(fa_seq)

}


## subset fasta

#' get records from BStringset
#'
#' @description an object of class BStringset stores sequences and headers from FASTA or FASTQ files. See \link{Biostrings::readBStringSet()}. This function helps users to filter sequences by headers with full or partial match.
#' @param x character vector of query ids
#' @param y an object of class \code{BStringset} from which sequences to be filtered
#' @param partial_match logical (default TRUE).
#' \enumerate {
#' \item TRUE : return sequences if element of x matches anywhere in the headers of y.
#' \item FALSE : return sequences only of exact match occur between element of x and headers of y.
#' }
#'
#' @return an object of class BStringset
#' @export
#' @importFrom purrr map
#' @importFrom stringr str_which fixed
#' @examples
subset_bstringset <- function(x , y,  partial_match  = TRUE){

        ## validate inputs
        if(! is(x , "character")){
                stop("'x' must be an object of class 'character'" )
        }

        if(! is(y , "BStringSet")){
             stop("'y' must be an object of class 'BStringset'" )
        }

        fa_headers <- names(y)

        if(partial_match){
                match_index <- purrr::map(x , ~ stringr::str_which( fa_headers , stringr::fixed(..1) ) )
        } else {
                match_index <- purrr::map(x , ~ match(..1 , fa_headers))
        }
        match_index <-  unlist(match_index) %>% unique()
        match_index <- match_index[! match_index %>% is.na() ]

        if(length(match_index) == 0 ){
                stop("No match found")
        } else{
                match_seq <- y[match_index]
        }

        return(match_seq)
}




