library(ggplot2)

# Read the CSV file
CSV<-read.csv("C:\\Users\\skyeb\\OneDrive\\Desktop\\newclass\\course-code-repo-main\\bioinformatics_(BIOL130-230)\\genomics_gene-parsing-processing\\RGCAULevels.csv")

# Extract chromosome numbers
chromosome<-substr(CSV$Gene,2,2)

# Group Polar and NonPolar by chromosome
final_datafile<-data.frame(Gene=chromosome, NonPolarPercent=CSV$Non.Polar., PolarPercent=CSV$Polar., Length=CSV$Length)

oneway.test(GCContent~Gene,final_datafile)
oneway.test(AUContent~Gene,final_datafile)

# Calculate mean Polar% and Non-Polar% and make plot
polar_plot <- ggplot(final_datafile, aes(x = Gene, y = PolarPercent)) +
  geom_bar(stat = "summary", fun = "mean", fill = "blue") +  
  geom_errorbar(stat = "summary", fun.data = "mean_cl_normal", width = 0.2) +  
  labs(title = "% Polar AA by Chromosome", x = "Chromosome", y = "% Polar AA")

nonpolar_plot <- ggplot(final_datafile, aes(x = Gene, y = NonPolarPercent)) +
  geom_bar(stat = "summary", fun = "mean", fill = "red") +  
  geom_errorbar(stat = "summary", fun.data = "mean_cl_normal", width = 0.2) +  
  labs(title = "% NonPolar AA by Chromosome", x = "Chromosome", y = "% NonPolar AA")

# save to PDF
ggsave("C:\\Users\\skyeb\\OneDrive\\Desktop\\newclass\\course-code-repo-main\\bioinformatics_(BIOL130-230)\\genomics_gene-parsing-processing\\PlotOutput.pdf", polar_plot)
ggsave("C:\\Users\\skyeb\\OneDrive\\Desktop\\newclass\\course-code-repo-main\\bioinformatics_(BIOL130-230)\\genomics_gene-parsing-processing\\PlotOutput.pdf", nonpolar_plot)
