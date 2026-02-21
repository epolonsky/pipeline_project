library(sleuth)
library(dplyr)

stab <- data.frame(sample = c("SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"), condition = c("2dpi", "6dpi", "2dpi", "6dpi"), donor = c("1", "1", "3", "3"), path = c( "results/kallisto/SRR5660030", "results/kallisto/SRR5660033", "results/kallisto/SRR5660044", "results/kallisto/SRR5660045"), stringsAsFactors = FALSE) # make a table, "stab", describing samples and kallisto output

so <- sleuth_prep(stab, ~ donor + condition) # initialize sleuth object using sleuth_prep function from sleuth library
so <- sleuth_fit(so, ~ donor, 'reduced') # donor only
so <- sleuth_fit(so, ~ donor + condition, 'full')  # donor + condition 
so <- sleuth_lrt(so, 'reduced', 'full')   # test condition effect 

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) # extract the test results from the sleuth object
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) # filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant_table <- data.frame(target_id = sleuth_significant$target_id, test_stat = sleuth_significant$test_stat, pval = sleuth_significant$pval, qval = sleuth_significant$qval) # extract only the target_id, test_stat, pval, and qval columns from the results

write.table(sleuth_significant_table, file="results/sleuth_significant.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, append = TRUE) # write FDR < 0.05 transcripts to file
