ntable <- function(x,y=NULL,percent=1,digits=2,row=FALSE,col=FALSE){
    if (is.null(y)){
        cat(rownames(table(x)),"\n",sep="\t")
        cat(round(table(x)/sum(table(x)),digits),"\n",sep="\t")
        cat(table(x),"\n",sep="\t")
    }
    else{
        if(length(x)!=length(y))
            stop("x and y do not have same length\n")
        tab <- table(x,y)
        pt <- round(prop.table(tab,percent),digits)
        k <- !is.na(x) & !is.na(y)
        if (row) {
            pt <- cbind(pt,round(table(x)/sum(k),digits))
            tab <- cbind(tab,table(x[k]))
        }
        cat("",colnames(tab),"\n",sep="\t")
        j <- 0
        for (i in rep(c(0,1),nrow(tab))) {
            j <- j + 1 - i
            if (i==0) {
                cat(row.names(pt)[j],pt[j,],"\n",sep="\t")
            }
            else {
                cat("",tab[j,],"\n",sep="\t")
            }
        }
        if (col) {
            cat("\n")
            cat("",round(table(y[!is.na(x)])/sum(k),digits),"\n",sep="\t")
            cat("",table(y[!is.na(x)]),sum(k),"\n",sep="\t")}
    }

    invisible(NULL)
}
