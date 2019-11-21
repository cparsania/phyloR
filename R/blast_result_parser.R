
#' Wrapper function around taxize::genbank2uid.
#'
#' Given a genBank accession alphanumeric string, or a gi numeric string \code{(x)}, it returns tibble of taxid, name and other columns.
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

        start_time <- Sys.time()
        uid_list <- taxize::genbank2uid(x ,  ...)
        uid_tbl <- tibble::tibble(x = x, taxid = unlist(uid_list)) %>%
                dplyr::bind_cols( purrr::map_df(uid_list , attributes))
        time_taken <- start_time - Sys.time()
        message(time_taken)
        return(uid_tbl)

}




#' Title
#'
#' @param x a vector of valid ncbi taxid
#' @param rank default "kingdom". Can be one of the "no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"
#'
#' @return a tbl
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
#' @param blast_tbl an object of class tbl from blast output format 7. Column names must be identical to the given below
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
#' It can be assigned using the function \code{get_blast_outformat_7_colnames()}
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
#'
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
#' @return A character vector
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



#' Calculate coverage
#'
#' @param start A numeric denoting start
#' @param end A numeric denoting end
#' @param len A numeric denoting length
#'
#' @return A numeric
#'
get_cov <- function(start , end , len){

        ( (end - start +1 )/ len) * 100

}


#' Calculate query coverage.
#'
#' @param qstart A numeric denoting query start
#' @param qend A numeric denoting query end
#' @param qlen A numeric denoting query length
#' @export
#' @return A numeric
#'
get_query_cov <- function(qstart , qend , qlen){
        get_cov(start = qstart , end = qend , len = qlen)
}


#' Calculate subject coverage.
#'
#' @param sstart A numeric denoting subject start
#' @param send A numeric denoting subject end
#' @param slen A numeric denoting subject length
#' @export
#' @return A numeric
#'
get_subj_cov <- function(sstart , send , slen){
        get_cov(start = sstart , end = send , len = slen)

}


