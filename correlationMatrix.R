# requires corrplot library
library(corrplot)
# filepath for input tabulated text file
filename <- "//mrc-smb.lifesci.dundee.ac.uk/mrc-group-folder/ALESSI/Toan/For Golgitag_Paper/For_Pearson_Corr_02.txt"
df <- read.table(filename, header = TRUE, sep="\t")
df <- df[colnames(df)[1:which(colnames(df) == "HA.WCL_06")]]
cor_mat <- cor(as.matrix(df), use="everything")
# filepath for output pdf file
pdf("//mrc-smb.lifesci.dundee.ac.uk/mrc-group-folder/ALESSI/Toan/For Golgitag_Paper/For_Pearson_Corr_02.txt.pdf")
corrplot(cor_mat, order="hclust", type="lower", method="ellipse")
dev.off()