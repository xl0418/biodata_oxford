if __name__ == '__main__':
	import ChromAnalysis3 
	import numpy as np

	input_file = "Merged_BiallelicOnly_NoMissing_annot.txt"   # origional data
	filter_outfile = "filter.txt" # store the data that satisfies the length criteria
	total_reads_file = "total.txt" # further substract the total reads and save for futher calculation of means and stds
	result_file = "result.txt" # output file names for the support and quality

	test = ChromAnalysis3.ChromAnalysis(input_file=input_file, 
		filter_outfile = filter_outfile, 
		total_reads_file = total_reads_file,
		filter_length=10000,
		result_file = result_file, 
		cores=10)
	test.historgram()
