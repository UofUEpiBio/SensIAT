match.names <- function(x, names, allow.multiple = NA){
    i <- pmatch(names(x), names, duplicates.ok = TRUE)
    if(!isTRUE(allow.multiple) && anyDuplicated(i, NA_integer_)){
        if(is.na(allow.multiple))
            warning("Multiple matches found")
        else
            stop("Multiple matches found")
    }
    names(x)[!is.na(i)] <- names[i[!is.na(i)]]
    return(x)
}
