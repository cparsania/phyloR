
#' Print green tick on terminal. Wrapper function on  cli::cat_bullet
#'
#' @param ... Other arguments to be passed to \code{cli::cat_bullet}
#' @keywords internal
#' @return
#'
cat_green_tick <- function (...){

        cli::cat_bullet(..., bullet = "tick", bullet_col = "green")
}


#' Check value NULL or Numeric
#'
#' @param x a value to be checked for numeric or NULL
#' @keywords internal
#' @return
#'
null_or_numeric <- function(x){

        arg  <- rlang::enexpr(x)
        if(!is.null(x)){
                if(!is.numeric(x)){
                        stop(arg ," must be either NULL or numeric ")
                }
        }
}

#' Check value NULL or character
#'
#' @param x a value to be checked for character or NULL
#' @keywords internal
#' @return
#'
null_or_character <- function(x){

        arg  <- rlang::enexpr(x)
        if(!is.null(x)){
                if(!is.character(x)){
                        stop(arg ," must be either NULL or character ")
                }
        }
}

#' Check value NULL or LOGICAL
#'
#' @param x a value to be checked for logical or NULL
#' @keywords internal
#' @return
#'
null_or_logical <- function(x){

        arg  <- rlang::enexpr(x)
        if(!is.null(x)){
                if(!is.logical(x)){
                        stop(arg ," must be either NULL or logical ")
                }
        }
}


