if __name__ == '__main__':
	import ChromAnalysis3 
	import numpy as np

	input_file = "Merged_BiallelicOnly_NoMissing_annot.txt"   # origional data
	filter_outfile = "filter2.txt" # store the data that satisfies the length criteria
	total_reads_file = "total2.txt" # further substract the total reads and save for futher calculation of means and stds
	result_file = "result2.txt" # output file names for the support and quality

	test = ChromAnalysis3.ChromAnalysis(input_file=input_file, 
		filter_outfile = filter_outfile, 
		total_reads_file = total_reads_file,
		filter_length=10000,
		result_file = result_file, 
		cores=20)
	test.reduce_data()
	# test.filter_rawdata(chunksize = 1024 * 1024 * 10)
	# test.write_summary_stats(chunksize = 1024 * 1024 * 10)
	# test.compute_summary(chunksize = 1024 * 1024 * 10)
	# test.historgram()

	# lower = np.zeros(len(test.individualnames))
	# lower.fill(5)
	# upper = test.total_means + test.total_stds

	# test.data_analysis(lower = lower, upper=upper, mode="support", chunksize = 1024 * 1024 * 10)
	# test.data_analysis(lower = lower, upper=upper, mode="quality", chunksize = 1024 * 1024 * 10)

	
	# Merged_BiallelicOnly_NoMissing_annot