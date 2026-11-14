ref <- commandArgs(TRUE)[1]
stopifnot(nzchar(ref))
oracle <- new.env(parent=globalenv())
oracle$signature <- methods::signature
oracle$setMethod <- function(f, signature, definition, ...) {
    if (unname(signature['x']) == 'matrix') oracle[[f]] <- definition
}
for (f in c('helper.R', 'RUVg-methods.R', 'RUVr-methods.R', 'RUVs-methods.R'))
    sys.source(file.path(ref, 'R', f), oracle)
set.seed(42)
x <- matrix(rpois(30*16, 20), 30, 16)
x[,9:16] <- 2*x[,9:16]
e <- matrix(rnorm(30*16), 30, 16)
g <- oracle$makeGroups(rep(letters[1:4], each=4))
write.table(x,'data/input.csv',sep=',',row.names=FALSE,col.names=FALSE)
write.table(e,'data/residuals.csv',sep=',',row.names=FALSE,col.names=FALSE)
write.table(g,'data/groups_r.csv',sep=',',row.names=FALSE,col.names=FALSE)
cases <- list(r=list(fun='RUVr', args=list(residuals=e)),
              r_log=list(fun='RUVr',args=list(residuals=e,isLog=TRUE,center=FALSE)),
              s=list(fun='RUVs',args=list(scIdx=g)),
              s_log=list(fun='RUVs',args=list(scIdx=g,isLog=TRUE)))
for (name in names(cases)) {
    item <- cases[[name]]
    input <- if (isTRUE(item$args$isLog)) log(x+1) else x
    result <- do.call(oracle[[item$fun]],c(list(x=input,cIdx=1:12,k=2,round=FALSE),item$args))
    write.table(result$normalizedCounts,paste0('data/',name,'_corrected.csv'),sep=',',row.names=FALSE,col.names=FALSE)
    w <- result$W
    projection <- w %*% solve(t(w)%*%w,t(w))
    write.table(projection,paste0('data/',name,'_projection.csv'),sep=',',row.names=FALSE,col.names=FALSE)
}
