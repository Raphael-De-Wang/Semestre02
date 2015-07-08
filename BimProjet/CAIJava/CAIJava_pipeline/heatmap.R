library(gplots)

m <- read.csv('sign_vec.csv',sep=' ',header=TRUE)
m1 <- as.matrix(m[2:65])
pdf(file='heatmap.pdf',width=20,height=10)
for(line in m ) {line;break}
rownames(m1) <- line

# cexRow = 1/log10(1000)
# cexCol = 1/log10(10)
# heatmap(m1,Rowv=cexRow, Colv=cexCol,col = heat.colors(256), scale="column", margins=c(5,10))
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 200)
# heatmap.2(m1,col = heat.colors(256), scale="column", margins=c(10,25), trace='none')
heatmap.2(m1,col = my_palette, scale="column", margins=c(10,25), trace='none')

# heatmap(m1, col = heat.colors(256))
dev.off()

