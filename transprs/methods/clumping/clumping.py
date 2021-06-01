from transprs.methods.utils import tmp_extract
import subprocess


def clumping(processor):
    tmp_extract(processor)
    subprocess.call(
        """
		plink \
		    --bfile tmp \
		    --clump-p1 1 \
		    --clump-r2 0.1 \
		    --clump-kb 250 \
		    --clump tmp_ss \
		    --clump-snp-field SNP \
		    --clump-field P \
		    --out tmp_out
		""",
        shell=True,
    )
