
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
        #x <- c(166361, 166361, 224129, 217634)

        x <- tibble::tibble(query_taxon = unique(x) %>% as.character())
        rlang::arg_match(rank  , c("no rank", "superkingdom", "kingdom", "phylum", "subphylum", "class", "subclass", "infraclass", "cohort", "order", "suborder", "infraorder", "superfamily", "family", "subfamily", "genus", "species", "tribe"))

        if(length(rank) !=1 ) stop("argument `rank` must be of length 1")

        result_ranks <- taxizedb::classification(x$query_taxon)  %>%
                tibble::tibble(taxid = names(.) , all_ranks = .) %>%
                tidyr::unnest(cols = all_ranks) %>%
                dplyr::filter(.$rank == !!rank)

        rank_id_col_name <- paste(rank , "id" ,sep = "_")

        x %>% dplyr::left_join(result_ranks, by = c("query_taxon" = "taxid")) %>%
                ## rename cols
                dplyr::rename(!!rank_id_col_name := id , !!rank := name) %>%
                ## arrange cols
                dplyr::select(query_taxon,!!rank ,!!rank_id_col_name,rank)

}


