library(gplots)

outm <- c()

for ( i in 30:30 ) {
    iname <- sprintf('sign_vec_len%d.csv', i)
    oname <- sprintf('heatmap_len%d.pdf', i)
    m <- read.csv( iname, sep=' ', header=TRUE)
    m1 <- as.matrix(m[2:65])
    pdf(file=oname,width=20,height=10)
    for(line in m ) {line;break}
    rownames(m1) <- line
    # cexRow = 1/log10(1000)
    # cexCol = 1/log10(10)
    # heatmap(m1,Rowv=cexRow, Colv=cexCol,col = heat.colors(256), scale="column", margins=c(5,10))
    my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 200)
    # heatmap.2(m1,col = heat.colors(256), scale="column", margins=c(10,25), trace='none')
    cluster <- heatmap.2(m1,col = my_palette, scale="column", margins=c(10,25), trace='none')
    # heatmap(m1, col = heat.colors(256))
    dev.off()
    clu1 <- gsub("list","",toString(cluster$rowDendrogram[1]))
    clu2 <- gsub("list","",toString(cluster$rowDendrogram[2]))
    clu1 <- gsub("[[:punct:]]", "", clu1)
    clu2 <- gsub("[[:punct:]]", "", clu2)
    clu1 <- lapply(strsplit(clu1, " "), as.numeric)[[1]]
    clu2 <- lapply(strsplit(clu2, " "), as.numeric)[[1]]
    len <- length(clu1) + length(clu2)
    if ( len %in% clu1 ) {
        outm <- c(outm, toString(c(sprintf('len_%d', i),clu2)))
    } else {
        outm <- c(outm, toString(c(sprintf('len_%d', i),clu1)))
    }
}

write(outm,file='cluster_info.txt')
