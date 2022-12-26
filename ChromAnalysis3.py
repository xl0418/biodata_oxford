import numpy as np
import os
import multiprocessing as mp
from itertools import repeat
import logging  
import sys
import matplotlib.pyplot as plt
import time  

logging.basicConfig(filename="ChromAnalysis.log", filemode="a",
	format='%(asctime)s - %(message)s', level=logging.INFO)
logger = logging.getLogger("CALog")
#logger.addHandler(logging.StreamHandler(sys.stdout))
class ChromAnalysis:
	# initialize input and output file names
	def __init__(self, input_file, filter_length,filter_outfile, result_file, total_reads_file, cores = 4):

		self.input_file = input_file
		self.result_file = result_file
		self.contigs_id = []
		self.filter_outfile = filter_outfile
		self.total_reads_file = total_reads_file
		self.cores = cores
		self.contigs_reduce = []
		self.reduced_file = "reduced_" + self.input_file

		# substract the contigs with the length larger than the filter criteria and the column names of all individual samples. 
		with open(self.input_file) as f:
			for line in f:
				if line.startswith("##contig"):
					temp_contigs = line.split(",")
					temp_contig_length = int(temp_contigs[1].split("=")[1])
					if temp_contig_length > filter_length:
						self.contigs_id.append(temp_contigs[0].split("=")[2])
					else:
						self.contigs_reduce.append(temp_contigs[0].split("=")[2])

				if not line.startswith('##') and line.startswith('#'):
					colnames = line.split("\t") # grab the column names
					self.ncol = len(colnames)
					col_remove_ids = [0, 1, 3, 4, 5, 6, 7, 8]
					for index in sorted(col_remove_ids, reverse=True):
						del colnames[index]
					break

		self.individualnames = colnames[1:]
		self.individualnames[-1] = self.individualnames[-1].split("\n")[0]
		first_row = '#ID'
		# create new headings for the filtered output.
		for indiv_name in self.individualnames:
			first_row = first_row + '\t' + indiv_name

		self.filtered_concise_data = [first_row]

		logger.info(f"Reading the source data {self.input_file}...")
		logger.info(f"The filtered data is stored in {self.filter_outfile}")
		logger.info(f"The total reads of the filtered data is stored in {self.total_reads_file}" )
		logger.info(f"The final output is stored in Support_/Quality_{self.result_file}")
		logger.info(f"In total, {len(self.contigs_id)} contigs are selected with the criteria of the length larger than {filter_length}.")
		logger.info(f"{len(self.individualnames)} individuals/samples are assessed.")
		logger.info("---------------------------------------------------------")

	def process_filter(self, chunkStart, chunkSize):
		# Read lines from the data file.
		with open(self.input_file) as f:
			f.seek(chunkStart)
			lines = f.read(chunkSize).splitlines()
			Reads_line = []
			for line in lines:
				if not line.startswith("#"):
					split_string = line.split("\t")
					chrom_pos_id = split_string[0]
					if chrom_pos_id in self.contigs_id:
						pos_info = split_string[2] # ID
						for individual_id in range(9, self.ncol):
							ATinfo, treads = split_string[individual_id].split(":")[1:3]
							pos_info = pos_info + '\t' + ATinfo + ':' + treads
						Reads_line.append(pos_info)
				else:
					continue
		return Reads_line
	
	def filter_rawdata(self, chunksize = 1024 * 1024):
		logger.info("Start to filter the data...")
		pool = mp.Pool(self.cores)
		jobs = []

		#create jobs
		for chunkStart,chunkSize in self.chunkify(filename = self.input_file, size = chunksize):
		    jobs.append( pool.apply_async(self.process_filter,(chunkStart,chunkSize)) )

		#wait for all jobs to finish
		total_reads_list = []
		for job in jobs:
			total_reads_list.extend(job.get())  

		#clean up
		pool.close()
		with open(self.filter_outfile, 'w') as fsw:
			for item in total_reads_list:
				fsw.write("%s\n" % item)
		logger.info(f"Filtering data is done! -- {self.filter_outfile}")


	# read chunks from the txt file for parallel computation
	def chunkify(self,filename, size): 
		fileEnd = os.path.getsize(filename)
		with open(filename,'rb') as f:
			chunkEnd = f.tell()
			while True:
				chunkStart = chunkEnd
				f.seek(size,1)
				f.readline()
				chunkEnd = f.tell()
				yield chunkStart, chunkEnd - chunkStart
				if chunkEnd > fileEnd:
					break

	def reduce_data(self, chunksize = 1024 * 1024):
		print('copying "{}" --> "{}"'.format(self.input_file, self.reduced_file))
		logging.info('copying "{}" --> "{}"'.format(self.input_file, self.reduced_file))

		if os.path.exists(self.input_file) is False:   
			print('ERROR: file does not exist: "{}"'.format(self.input_file))
			sys.exit(1)
		if os.path.exists(self.reduced_file) is True:
			os.remove(self.reduced_file)
		if os.path.exists(self.reduced_file) is True:   
			print('ERROR: file exists, cannot overwrite it: "{}"'.format(self.reduced_file))
			sys.exit(1)

		with open(self.input_file) as f:
			with open(self.reduced_file, 'w') as ofp:
				for line in f:
					if line.startswith("##"):
						if line.startswith("##contig"):
							temp_contigs = line.split(",")
							temp_contig_id = temp_contigs[0].split("=")[2]
							if temp_contig_id in self.contigs_reduce:
								pass
							else:
								ofp.write(line)
						else:
							ofp.write(line)
					elif not line.startswith('#'):
						temp_contig_id = line.split("\t")[0] 
						if temp_contig_id in self.contigs_reduce:
							continue
						else:
							ofp.write(line)
					else:
						ofp.write(line)
		

		logger.info(f"Reduced the data by removing the small contigs! -- {self.reduced_file}")


	# Extract total reads from chunks
	def process_wrapper(self, chunkStart, chunkSize):
		with open(self.filter_outfile) as f:
			f.seek(chunkStart)
			lines = f.read(chunkSize).splitlines()
			Reads_line = []
			for line in lines:
				if not line.startswith("#"):
					split_string = line.split("\t")
					pos_info = split_string[0]  # ID
					for individual_id in range(1, len(self.individualnames) + 1):
						treads = split_string[individual_id].split(":")[1]
						pos_info = pos_info + '\t' + treads
					Reads_line.append(pos_info)
				else:
					continue
		return Reads_line



    # write the total reads in disk for future use
	def write_summary_stats(self, chunksize = 1024 * 1024):
		logger.info("--------")

		logger.info("Start to write the summary stats to the file in parallel...")

		if not os.path.isfile(self.filter_outfile): 
			return print("Please run filter_rawdata function first or specifiy the correct filtered data file.")
		else:
			pool = mp.Pool(self.cores)
			jobs = []

			#create jobs
			for chunkStart,chunkSize in self.chunkify(filename = self.filter_outfile, size = chunksize):
			    jobs.append( pool.apply_async(self.process_wrapper,(chunkStart,chunkSize)) )

			#wait for all jobs to finish
			total_reads_list = []
			for job in jobs:
				total_reads_list.extend(job.get())  

			#clean up
			pool.close()

			with open(self.total_reads_file, 'w') as fsw:
				for item in total_reads_list:
					fsw.write("%s\n" % item)

			logger.info(f"Writing the summary stats to the file is done! -- {self.total_reads_file}")


    # convert the total reads list to ndarray for summary computing
	def list_to_array(self, chunkStart, chunkSize):
		with open(self.total_reads_file) as f:
			f.seek(chunkStart)
			lines = f.read(chunkSize).splitlines()
			Reads_line = []
			for line in lines:
				split_string = [int(x) for x in line.split("\t")[1:]]
				Reads_line.extend(split_string)
		return Reads_line

	# compute summaries to set upper and lower boundary of the total reads for each individual
	def compute_summary(self, chunksize = 1024*1024):
		logger.info("--------")

		logger.info("Start to compute the summary stats to the file in parallel...")

		pool = mp.Pool(self.cores)
		jobs = []

		#create jobs
		for chunkStart,chunkSize in self.chunkify(filename = self.total_reads_file, size = chunksize):
		    jobs.append( pool.apply_async(self.list_to_array,(chunkStart,chunkSize)) )

		total_reads_list = []
		for job in jobs:
			total_reads_list.extend(job.get())  

		#clean up
		pool.close()

		total_reads_array = np.array(total_reads_list).reshape(int(len(total_reads_list)/len(self.individualnames)), len(self.individualnames))
		self.total_means = np.mean(total_reads_array,axis=0)
		self.total_stds = np.std(total_reads_array,axis=0)

		logger.info(f"Computing the summary stats is done!")

		logger.info(f"In total, {total_reads_array.shape[0]} positions are seleted to compute the means and std of reads depth.")

		logger.info(f"The means of each individual sample is {self.total_means} ...")

		logger.info(f"The standard derivation of each individual sample is {self.total_stds} ...")


	# Extract total reads from chunks
	def process_support(self, chunkStart, chunkSize, lower, upper, mode):
		with open(self.filter_outfile) as f:
			f.seek(chunkStart)
			lines = f.read(chunkSize).splitlines()
			output_support = []
			output_quality = []

			for line in lines:
				if not line.startswith("#"):
					temp_line = [""]
					split_string = line.split("\t")
					individual_count = 0
					for individual_id in range(1, len(self.individualnames) + 1):
						TReind = int(split_string[individual_id].split(":")[1])
						if mode == "quality":
							if TReind >= lower[individual_id - 1] and TReind <= upper[individual_id - 1]:
								# Good info
								Geind = "G"
							elif TReind >  upper[individual_id - 1]:  
								Geind = "H"
							else:
								Geind = "B"
							temp_line[0] = temp_line[0] + '\t' + Geind
						elif mode == "support":
							# Support info
							ATeind = split_string[individual_id].split(":")[0]
							# Reference / Alternative support
							Aeind, Teind = ATeind.split(",")[0:2]
							if TReind < lower[individual_id - 1] or TReind >upper[individual_id - 1]: 
								# Result of the mark: Missing info "-"; Reference support "1"; Alternative support "0"
								Reind = "-"
							else:
								if int(Aeind) == TReind:
									Reind = "1"
								elif int(Teind) == TReind:
									Reind = "0"
								else:
									Reind = "-"
							temp_line[0] = temp_line[0] + '\t' + Reind
						else:
							break
					if mode == "support":
						output_support.append(temp_line)
					else:
						output_quality.append(temp_line)
				else:
					continue
		if mode == "support":
			return output_support
		else:
			return output_quality

	def assign_to_output(self, individual_id, result_row, mode):
		for nrow in range(len(self.analysis_result)):
			#print(analysis_result[l])
			if mode == "support":
				result_row = result_row + self.analysis_result[nrow][0].split("\t")[individual_id+1]
			elif mode == "quality":
				result_row = result_row + self.analysis_result[nrow][0].split("\t")[individual_id+1]
			else:
				break
		return result_row

	def historgram(self, histplot = False, chunksize = 1024 * 1024):
		logger.info("--------")
		logger.info("Start plotting historgrams...")

		pool = mp.Pool(self.cores)
		jobs = []

		#create jobs
		for chunkStart,chunkSize in self.chunkify(filename = self.total_reads_file, size = chunksize):
		    jobs.append( pool.apply_async(self.list_to_array,(chunkStart,chunkSize)) )

		total_reads_list = []
		for job in jobs:
			total_reads_list.extend(job.get())  

		#clean up
		pool.close()

		total_reads_array = np.array(total_reads_list).reshape(int(len(total_reads_list)/len(self.individualnames)), len(self.individualnames))

		if histplot:
			for individual_id in range(len(self.individualnames)):
				# the histogram of the data
				plt.hist(total_reads_array[:,individual_id], bins = 'auto', facecolor='r', alpha=0.75)
				plt.xlabel('Reads')
				plt.ylabel('Counts')
				plt.title('Histogram of Reads')
				plt.grid(True)
				plt.savefig(f"Hist_sample_{self.individualnames[individual_id]}.png")

		with open("Quantiles.txt","w") as f:
		    f.write("\n".join(" ".join(map(str, x)) for x in (
		    	self.individualnames,
				np.quantile(total_reads_array, 0.1, axis = 0),
				np.quantile(total_reads_array, 0.25, axis = 0),
				np.quantile(total_reads_array, 0.5, axis = 0),
				np.quantile(total_reads_array, 0.75, axis = 0),
				np.quantile(total_reads_array, 0.9, axis = 0)
		    	)))

		logger.info("Histograms are done!")



	def data_analysis(self, lower, upper, mode, chunksize = 1024*1024):
		logger.info("The lower boundary of reads of samples are: {}".format(' '.join(map(str, lower.tolist()))))
		
		logger.info("The upper boundary of reads of samples are: {}".format(' '.join(map(str, upper.tolist()))))

		logger.info("--------")

		logger.info("Start to assign Support and Quality index to the file in parallel...")

		pool = mp.Pool(self.cores)
		jobs = []

		#create jobs
		for chunkStart,chunkSize in self.chunkify(filename = self.filter_outfile, size = chunksize):
			jobs.append( pool.apply_async(self.process_support,(chunkStart,chunkSize, lower, upper, mode)) )
		#clean up
		pool.close()

		#wait for all jobs to finish
		self.analysis_result = []
		for job in jobs:
			self.analysis_result.extend(job.get())  

		logger.info(f"{len(self.analysis_result)} positions are processed in the data analysis step.")


		result = []
		# individual names for results files.
		for indiv_name in self.individualnames:  
			result.append(indiv_name + '\t')	


		with mp.Pool(self.cores) as pool:
			result_out = pool.starmap(self.assign_to_output, zip(range(len(self.individualnames)),
				result, repeat(mode)))

		if mode == "support":
			with open('Support_' + self.result_file, 'w') as fs:  
				for item in result_out:
					fs.write("%s\n" % item)
		elif mode == "quality":
			with open('Quality_' + self.result_file, 'w') as fq:  
				for item in result_out:
					fq.write("%s\n" % item)
		else:
			print("Please sepcify the mode, i.e. support or quality.")

		logger.info(f"All programms are done! The results are stored in {'Support_' + self.result_file}.")

