
#' A Wrapper function around \link{taxize::genbank2uid}.
#'
#' Given a ncbi gene bank accession alphanumeric string or a gi numeric string \code{(x)}, it returns a tibble contaning ncbi taxonomy id, name and other related columns.
#' @param x vector of ncbi gene bank accession alphanumeric string, or a gi numeric string \code{(x)}.
#' @param ... other parameters to be passed to \code{taxize::genbank2uid}.
#'
#' @return a tbl with colnames x, taxid, class, match, multiple_matches, pattern_match, uri, name
#' @export
#' @importFrom lubridate now
#' @importFrom taxize genbank2uid
#' @importFrom tibble tibble is_tibble
#' @importFrom dplyr bind_cols filter left_join rename select group_by mutate arrange slice ungroup pull bind_rows coalesce mutate_all
#' @importFrom purrr map_df map
#' @importFrom rlang arg_match sym as_name as_string
#' @importFrom cli cat_bullet cat_rule cli_alert_info
#' @importFrom taxizedb classification
#' @importFrom tidyr unnest unnest_wider
#' @importFrom glue glue glue_data
#' @importFrom TidyWrappers tbl_keep_rows_NA_any
#' @importFrom Biostrings readBStringSet
#' @importFrom stringr str_match_all str_trim str_match str_remove str_replace_all str_which fixed
#' @importFrom methods is
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
        time_taken <-  (lubridate::now() - start_time) %>% round(2)
        cat_green_tick("Done. ", " Time taken " , time_taken)
        cat_rule()
        return(uid_tbl)

}




#' Get ncbi phylogeny level for a given ncbi taxonomy id.
#'
#'
#' @param x a vector of valid ncbi taxonomy id.
#' @param rank a string denoting required ncbi phylogenetic level. Default "kingdom". An input can be one of the followings
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

        all_ranks <- id <- name <- query_taxon <- NULL

        x <- tibble::tibble(query_taxon = unique(x) %>% as.character())
        rlang::arg_match(rank  , c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"))

        if(length(rank) !=1 ) stop("argument `rank` must be of length 1")

        start_time <- lubridate::now()
        cli::cat_bullet("Rank search begins...",col = "red")
        cli::cat_rule()
        result_ranks <- taxizedb::classification(x$query_taxon)  %>%
                tibble::tibble(taxid = names(.) , all_ranks = .) %>%
                purrr::map(~.) %>%
                tidyr::unnest(cols = all_ranks) %>%
                dplyr::filter(.$rank == !!rank)
        time_taken <- (lubridate::now() - start_time) %>% round(2)
        cat_green_tick("Done. ", " Time taken " , time_taken)
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
#'  f <- system.file("extdata","blast_output_01.txt" ,package = "phyloR")
#'  d <- readr::read_delim(f, delim ="\t" , col_names = F , comment = "#")
#'  colnames(d) <- phyloR::get_blast_outformat_7_colnames()
#'  filtered <- filter_blast_hits(d ,
#'  query_cov = 90 ,
#'  identity = 40,
#'  query_length = 249,
#'  evalue = 1e-6)
#' }
#' @seealso \link{get_blast_outformat_7_colnames}
#'
filter_blast_hits <- function(blast_tbl ,
                              evalue = NULL,
                              bit_score = NULL,
                              query_cov = NULL,
                              identity = NULL,
                              query_length = NULL ){

        # global variable declaration
        q_start <- q_end <- NULL

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
#' @description Given the output of blast tabular format in a tbl, it removes redundant hits by subject accession.
#' For each redundant hit longest hit will be kept.
#'
#' @param blast_output_tbl an object of class tbl containing blast tabular output.
#' @param subject_start_colname a string denoting a column of subject start. Default "s_start".
#' @param subject_end_colname a string denoting a column of subject end. Default "s_end".
#' @param keep_length logical, default FALSE, indicates whether to keep column of subject length.
#' @param subject_acc_colname a string denoting a column of subject hits in a \code{blast_output_tbl}. Default "subject_acc_ver".
#'
#' @return a tbl
#' @export
#' @importFrom tibble is_tibble
#' @importFrom rlang sym
#' @importFrom dplyr group_by
#' @importFrom dplyr arrange
#' @importFrom dplyr slice
#' @examples
#' \dontrun{
#'  f <- system.file("extdata","blast_output_01.txt" ,package = "phyloR")
#'  d <- readr::read_delim(f, delim ="\t" , col_names = F , comment = "#")
#'  colnames(d) <- phyloR::get_blast_outformat_7_colnames()
#'  remove_redundant_hits(d)
#' }
remove_redundant_hits <-  function(blast_output_tbl,
                                   subject_acc_colname = "subject_acc_ver" ,
                                   subject_start_colname = "s_start",
                                   subject_end_colname = "s_end" ,
                                   keep_length = FALSE){

        if(!tibble::is_tibble(blast_output_tbl)){
                stop("blast_output_tbl must be a class of tibble")
        }

        subject_acc_colname = rlang::sym(subject_acc_colname)
        subject_start_colname = rlang::sym(subject_start_colname)
        subject_end_colname = rlang::sym(subject_end_colname)

        blast_output_tbl22 <- blast_output_tbl %>%
                dplyr::group_by(!!subject_acc_colname) %>%
                dplyr::mutate(subject_aligned_length = !!subject_end_colname - !!subject_start_colname + 1) %>%
                dplyr::arrange(desc(subject_aligned_length)) %>%
                dplyr::slice(1) %>%
                dplyr::ungroup()

        if(keep_length){
                blast_output_tbl22
        } else {
                blast_output_tbl22 %>% dplyr::select(-subject_aligned_length)
        }



}




## Given the blast hits in tabular format add taxonomy annotations

#' Add ncbi taxonomy levels to ncbi protein accession
#'
#' @description Given a tbl with a column of valid ncbi protein accession,
#' the function assigns the ncbi taxonomy levels to each ncbi protein accession.
#'
#' @param tbl an object of class tbl
#' @param ncbi_accession_colname a string (default : "ncbi_accession") denoting column name of ncbi accession.
#' @param ncbi_acc_key user specific ENTREZ api key. Get one via \code{taxize::use_entrez()}
#' @param taxonomy_level a string indicating level of ncbi taxonomy to be assigned to each ncbi protein accession. An input can be one of the followings
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
#' @param map_superkindom logical (default FALSE). Assign superkingdom if kingdom is not found. Valid only when taxonomy_level == "kingdom".
#' @param batch_size The number of queries to submit at a time.
#' @details The aim of this function is to assign the specific level of ncbi taxonomy to the ncbi accession (protein).
#' To do so, it requires a tibble with at least one column of ncbi (protein) accession.
#' Returned taxonomy columns will be added on input tibble object keeping original columns as they were.
#' Internally, first, it finds the ncbi taxonomy id for each ncbi accession and then it maps required taxonomy level.
#' Assigning taxonomy id to each ncbi accession may take time depending upon number of input ncbi accessions.
#' On subsequent runs or in a first run you may supply taxonomy column ('taxid') in input tibble,
#' which will reduce the time to find taxonomy ids and directly
#' assign the taxonomy level to given taxonomy id.
#' To map taxonomy levels for large number of ncbi accession one may choose parallel processing approach as shown in the example.
#' @return a tbl.
#' @export
#' @importFrom tibble is_tibble
#' @importFrom rlang arg_match
#' @importFrom cli cat_rule cli_alert_info cli_alert
#' @importFrom dplyr pull select bind_rows select rename left_join right_join
#' @importFrom purrr map
#' @importFrom TidyWrappers tbl_keep_rows_NA_any
#' @examples
#' \dontrun{
#' f <- system.file("extdata","blast_output_01.txt" ,package = "phyloR")
#' d <- readr::read_delim(f, delim ="\t" , col_names = F , comment = "#")
#' colnames(d) <- phyloR::get_blast_outformat_7_colnames()
#'
#' ## add kingdom
#' with_kingdom <- d %>%
#'         dplyr::slice(1:50) %>%
#'         add_taxonomy_columns(ncbi_accession_colname ="subject_acc_ver" )
#'
#' ## add species
#' with_kingdom_and_species <- with_kingdom %>%
#'         add_taxonomy_columns(ncbi_accession_colname ="subject_acc_ver",taxonomy_level = "species")
#' dplyr::glimpse(with_kingdom_and_species)
#'
#' #------------------------------------
#' ## using parallel processing approach
#'
#' library(furrr)
#' num_of_splits <- 10
#' d <- d %>% dplyr::slice(1:100)
#' split_vec <- rep(1:num_of_splits , length.out = nrow(d))
#' qq_split <- d %>% dplyr::mutate(split_vec = split_vec)  %>%
#' dplyr::group_by(split_vec) %>% dplyr::group_split()
#' future::plan("multiprocess")
#' out <- qq_split[1:num_of_splits] %>%
#'         future_map( ~ phyloR::add_taxonomy_columns(tbl = ..1 ,
#'  taxonomy_level = "species" ,map_superkindom = F,
#'  ncbi_accession_colname = "subject_acc_ver" , batch_size = 20,
#'  ncbi_acc_key = "64c65ab9c52e0312bbcf4c32d3056cbcaa09"),
#'                    .progress = TRUE) %>%
#'         dplyr::bind_rows()
#' }
add_taxonomy_columns <- function(tbl,
                                 ncbi_accession_colname = "ncbi_accession",
                                 ncbi_acc_key =NULL,
                                 taxonomy_level = "kingdom",
                                 map_superkindom = FALSE,
                                 batch_size = 20 ){

        ## globle variable declaration

        taxid <- query_taxon <- kingdom <- superkingdom <- NULL

        ## validate user inputs
        ## tbl must be tbl
        if(!tibble::is_tibble(tbl)){
                stop("blast_output_tbl must be a class of tibble")
        }

        ## taxonomy_level must be one one of these
        rlang::arg_match(taxonomy_level  ,
                         c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"))


        ## blast_output_tbl mus not contain column names 'taxid' and taxonomy_level
        if(any(colnames(tbl) %in% c(taxonomy_level))) {

                stop("Column '", taxonomy_level ,"'" , " must not present in 'blast_output_tbl'.",
                     " Either remove it or rename it")

        }
        ncbi_accession_colname <- rlang::sym(ncbi_accession_colname)

        ## ncbi_accession_colname must present in the tbl.
        if(!rlang::as_name(ncbi_accession_colname) %in% colnames(tbl)){
                stop(glue::glue("'tbl' must contains column '{ncbi_accession_colname}'. Use argument 'ncbi_accession_colname' to specify a column containing valid ncbi accession in a 'tbl'"))
        }


        is_taxid_present <- FALSE
        if(any(colnames(tbl) %in% "taxid")) {
                is_taxid_present <- TRUE
                cli::cat_rule("WARNING")
                cli::cli_alert_info("As column 'taxid' present in blast output tbl, same will be used to map taxonomy level.")
                cli::cli_alert_info("To perform new 'taxid' search either remove  or rename columnn 'taxid'." )
                cli::cat_rule("WARNING ENDS")
        }

        ## get subject accession
        ncbi_acc <- tbl %>%
                dplyr::pull(!!ncbi_accession_colname) %>%
                unique()

        if(is_taxid_present){
                taxon_data_sub <- tbl %>%
                        dplyr::select(!!ncbi_accession_colname , taxid)
        } else{
                ## get taxid for each subject hit
                taxon_data <- phyloR::genbank2uid_tbl(ncbi_acc,key = ncbi_acc_key, batch_size = batch_size)

                ## list of tbl to tbl
                taxon_data_sub <- taxon_data %>%
                        dplyr::bind_rows() %>%
                        dplyr::select(1,2) %>%
                        dplyr::rename(!!ncbi_accession_colname := "x")
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
                tbl <- tbl %>%
                        dplyr::left_join(taxon_data_sub , by = rlang::as_name(ncbi_accession_colname))
        }

        ## add user asked taxonomy_level
        tbl <- tbl %>%
                dplyr::left_join(taxid_to_taxon , by = c("taxid" = "query_taxon"))

        return(tbl)
}




#' Format fasta headers
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
#' @importFrom rlang .data
#' @return sequences as an object of class \code{Biostrings}
#' @export
#'
#' @examples
#' \dontrun{
#'  f <- system.file("extdata","blast_output_01.fasta" ,package = "phyloR")
#'  ## existing headers
#'  Biostrings::readAAStringSet(f) %>% names() %>% head()
#'  ## new  headers
#'  f %>% format_fasta_headers() %>% names() %>% head()
#' }
format_fasta_headers <- function(fasta_file = NULL, keep_alignemnt_coord= TRUE){

        # global variable

        header_elems <- subject_id <- NULL

        ## validate inputs
        if(!file.exists(fasta_file)){
                stop("fasta_file does not exist.")
        }

        fa_seq <- Biostrings::readBStringSet(fasta_file)
        fa_headers <- names(fa_seq)

        fa_headers_dislocate <- fa_headers %>%
                tibble::tibble(fa_headers = .) %>%
                dplyr::mutate(header_elems = purrr::map(fa_headers ,
                                                        ~( stringr::str_match_all(string = ..1 , pattern = "([^\\s]+)\\s(.*)\\[(.*)\\]") %>%
                                                                   unlist() ))) %>%
                tidyr::unnest_wider(col = header_elems) %>%
                dplyr::select(-.data$...1) %>%
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

#' Get records from BStringset
#'
#' @description an object of class BStringset stores sequences and headers from FASTA or FASTQ files. See \code{Biostrings::readBStringSet()}. This function helps users to filter sequences by headers with full or partial match.
#' @param x character vector of query ids
#' @param y an object of class \code{BStringset} from which sequences to be filtered
#' @param partial_match logical (default TRUE).
#' \enumerate{
#' \item TRUE : return sequences if element of x matches anywhere in the headers of y.
#' \item FALSE : return sequences only of exact match occur between element of x and headers of y.
#' }
#'
#' @return an object of class BStringset
#' @export
#' @importFrom purrr map
#' @importFrom stringr str_which fixed
#' @examples
#' \dontrun{
#' x <- c("EIT75269.1",
#'   "TGO19408.1",
#'   "KAF2153260.1",
#'   "OAA41719.1",
#'   "OSS52177.1",
#'   "XP_018252424.1",
#'   "XP_008598593.1",
#'   "KXN65110.1",
#'   "XP_018147989.1",
#'   "XP_022493698.1",
#'   "RII05464.1",
#'   "XP_018703519.1",
#'   "RZR67285.1",
#'   "OLY78428.1",
#'   "XP_007819064.1",
#'   "PQK17331.1",
#'   "KXN66278.1",
#'   "CRK21695.1",
#'   "CVK85925.1",
#'   "KID81639.1")
#' y_file <- system.file("extdata" ,"blast_output_01.fasta" , package = "phyloR")
#' y <- Biostrings::readBStringSet(y_file)
#' ## subset fasta by partial match
#' sub <- phyloR::subset_bstringset(x = x , y = y ,partial_match = T)
#' subset_bstringset(x , y,  partial_match  = TRUE)
#' }
subset_bstringset <- function(x , y,  partial_match  = TRUE){

        ## validate inputs
        if(! methods::is(x , "character")){
                stop("'x' must be an object of class 'character'" )
        }

        if(! methods::is(y , "BStringSet")){
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


#' Assign (NCBI) taxonomy levels to phylo convertible tibble.
#'
#' @param tree_data a tbl containing minimum columns required to create an object of class phylo.
#' @param ncbi_accession_colname a string (default : "ncbi_accession") denoting column name of ncbi accession.
#' @param taxonomy_levels a character vector containing levels of ncbi taxonomy.
#' For each of these levels corresponding values will be mapped for ncbi accession.
#'
#' @description An object of class \code{phylo} can be converted to an object of class \code{tibble} and vice versa.
#' The given function adds columns of ncbi taxonomy levels to the tibble created from class \code{phylo}.
#' Resultant tibble can be passed to \code{ggtree::ggtree()} to visualize phylogenetic tree.
#' Added taxonomy columns can be used for aesthetics such as tip color, node color etc.
#' @return a tbl containing all columns from \code{tree_data} + columns of assigned taxonomy levels.
#' @export
#' @examples
#' \dontrun{
#' ## read tree from newick string
#' tree_string <- "((XP_005187699_1__Musca_domestica:0.070627277,(XP_019893806_1__Musca_domestica:0.071069674,((XP_013113221_1__Stomoxys_calcitrans:0.1494662042,ACB98719_1__Glossina_morsitans_morsitans:0.3489851076)67.4/100:0.0470213767,XP_013102958_1__Stomoxys_calcitrans:0.1794878827)98.1/100:0.0959227604)88.2/99:0.0323598861)93/99:0.0435291148,((XP_017472861_1__Rhagoletis_zephyria:0.0049337059,XP_017472862_1__Rhagoletis_zephyria:0.0112391294)97.3/100:0.0860969479,(XP_020713236_1__Ceratitis_capitata:0.2642805176,(XP_014102010_1__Bactrocera_oleae:0.1183517872,XP_018784523_1__Bactrocera_latifrons:0.1137567198)29.6/88:0.0758551876)99.9/100:0.247740081)92/100:0.0716529011)34.3/66:2.487103817;"
#' tree_objct <- read.tree(text = tree_string)

#' tree_tbl <- tree_objct %>% ggtree::fortify()

#' tree_tbl <- tree_tbl  %>%
#'         dplyr::mutate( seqid =  dplyr::case_when(isTip ~ stringr::str_replace(label , pattern = "__.*","" ) %>%  ## split by '__'
#'                                                         stringr::str_replace(pattern = "_\\d$" , ""), ## remove trailing digits from seqid
#'                                                  TRUE ~ label
#'         )
#'        )
#' ## add taxonomy
#' tree_tbl_with_taxonomy <- tree_tbl %>%
#'        phyloR::tidy_taxonomy_tree(ncbi_accession_colname = "seqid",taxonomy_levels = c("species" ,"kingdom","family"))
#'
#' ## visualize  tree
#'
#' # tips colored by species
#' tree_tbl_with_taxonomy %>% ggtree::ggtree() + ggtree::geom_tiplab(aes(color = species))
#' # tips colored by family
#' tree_tbl_with_taxonomy %>% ggtree::ggtree() + ggtree::geom_tiplab(aes(color = family))
#' # tips colored by kingdom
#' tree_tbl_with_taxonomy %>% ggtree::ggtree() + ggtree::geom_tiplab(aes(color = kingdom))
#' }
tidy_taxonomy_tree <- function(tree_data,
                               ncbi_accession_colname = "ncbi_accession",
                               taxonomy_levels = c("species" ,"kingdom")){


        stopifnot(all(taxonomy_levels %in% c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe")))

        ## get accession column name
        accession_col_name = rlang::sym(ncbi_accession_colname)

        ## prepare accession column tibble
        accession_tbl <- tree_data %>%
                dplyr::filter(isTip) %>%
                dplyr::select(!!accession_col_name)

        ## map taxonomy levels
        for (i in seq_len(length(taxonomy_levels))){
                accession_tbl <- accession_tbl %>%
                        phyloR::add_taxonomy_columns(ncbi_accession_colname = rlang::as_string(accession_col_name),
                                                     taxonomy_level = taxonomy_levels[i])
        }
        ## add taxonomy to tree data
        tree_data_with_taxonomy <- tree_data %>%
                dplyr::left_join(accession_tbl , by = rlang::as_string(accession_col_name))


}

