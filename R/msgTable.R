msgTableCreate <- function() {
  t <- data.frame(id = character(0), msg = character(0), stringsAsFactors = FALSE)
  class(t) <- c("msgTable", "data.frame")
  t
}

msgTableAppend <- function(...) {
  u <- rbind(...)
  class(u) <- c("msgTable", "data.frame")
  u
}

msgTableAdd <- function(t, id, msg) {
  newlines <- data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  msgTableAppend(t, newlines)
}

msgTableSingleton <- function(id, msg) {
  msgTableAdd(msgTableCreate(), id, msg)
}

msgTableIsEmpty <- function(x)
  dim(x)[1] == 0
  
#' Print \code{msgTables} objects
#' 
#' Print in the REPL the \code{msgTables}
#' 
#' @param x an object of class \code{msgTables}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return Print in the REPL the \code{msgTables}
#' 
#' @export
print.msgTable <- function(x, ...) {
  if (msgTableIsEmpty(x)) {
    cat("Correct format\n")
  }
  else {
    cat("Message(s):\n")
    for (m in x$msg) {
      cat(paste("\t",m,"\n",sep=""))
    }
  }
}
