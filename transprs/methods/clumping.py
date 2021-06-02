from transprs.methods.utils import tmp_extract
import subprocess
import time
import datetime
import pandas as pd

def clumping(
    processor,
    clump_p1=1,
    clump_r2=0.5,
    clump_kb=250,
    ):
    
	start_time = time.time()
	print("Extracting data...")
	tmp_extract(processor)
	print("Done extract data!")
	print("Clumping is running...") 
	subprocess.call("""
		plink \
		    --bfile tmp \
		    --clump-p1 %s \
		    --clump-r2 %s \
		    --clump-kb %s \
		    --clump tmp_ss \
		    --clump-snp-field SNP \
		    --clump-field P \
		    --out tmp_out
		""" % (str(clump_p1),str(clump_r2),str(clump_kb)), shell=True)


    
	subprocess.call("""
    awk 'NR!=1{print $3}' tmp_out.clumped >  tmp_out.valid.snp
		""", shell=True)
    
    
	print("Done clumping!")
	valid_snps = list(pd.read_table("tmp_out.valid.snp",header=None)[0])
	processor.adjusted_ss["clumping"] =  processor.sumstats[processor.sumstats.SNP.isin(valid_snps)]

	print("The clumping result stores in .adjusted_ss['clumping']!")
    
	subprocess.call("""
    rm ./tmp*
		""", shell=True)
    
	print("--- Done in %s ---" % (str(datetime.timedelta(seconds=round(time.time() - start_time)))))
    