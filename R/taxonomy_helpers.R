
#' Wrapper function around taxize::genbank2uid.
#'
#' Given a genBank accession alphanumeric string, or a gi numeric string \code{(x)}, it returns tibble of taxid, name and other columns.
#' @param x vector of genBank accession alphanumeric string, or a gi numeric string \code{(x)}.
#' @param ... other parameters to be passed to taxize::genbank2uid
#'
#' @return a tbl with colnames x, taxid, class, match, multiple_matches, pattern_match, uri, name
#' @export
#' @importFrom taxize genbank2uid
#' @importFrom tibble tibble
#' @importFrom dplyr bind_cols
#' @examples
#' \dontrun{
#' x <- c("XP_022900619.1", "XP_022900618.1", "XP_018333511.1", "XP_018573075.1")
#' genbank2uid_tbl(x)
#' }
genbank2uid_tbl <- function(x , ...){

        start_time <- Sys.time()
        uid_list <- taxize::genbank2uid(x , ...)
        uid_tbl <- tibble::tibble(x = x, taxid = unlist(cc)) %>%
                dplyr::bind_cols( map_df(cc , attributes))
        time_taken <- tt - Sys.time()
        message(time_taken)
        return(uid_tbl)

}
