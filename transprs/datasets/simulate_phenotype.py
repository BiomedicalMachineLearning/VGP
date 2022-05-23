import time
import datetime

from rpy2 import robjects


def simulate_phenotype(list_genotype_prefix, gene_correlation_matrix, nsnp, h2):

    start_time = time.time()

    print("Phenotype simulation is running...")

    path = os.path.dirname(transprs.__file__)
    list_snps_path = path + "/datasets/list_snps"

    robjects.globalenv["ethniclist"] = list_genotype_prefix
    robjects.globalenv["gencov"] = gene_correlation_matrix
    robjects.globalenv["num_snp"] = nsnp
    robjects.globalenv["h2"] = list_genotype_prefix
    robjects.globalenv["list_snps_path"] = list_snps_path

    robjects.r(
        """
        library(snpStats)
        library(PhenotypeSimulator)

        num_ethnics = length(ethniclist)


        list_snps = read.table(list_snps_path, header=F, stringsAsFactors=F)
        v1 = seq(1:nrow(list_snps))
        v2 = sample(v1, num_snp)

        write.table(list_snps[v2,], file="tmp_causalsnp.txt", quote=F, row.names=F, col.names=F)

        h2factor = h2/num_snp
        for (i in 1:num_ethnics)
            for (j in 1:num_ethnics) gencov[i, j] = gencov[i, j] * h2factor

        library(MASS)

        snp_weight = mvrnorm(num_snp, rep(0,num_ethnics), gencov) ###

        for (ethnic_i in 1:num_ethnics) {
            bed_file = ethniclist[ethnic_i]
            genotype = read.plink(bed=bed_file)
            print(bed_file)
            geno = as(genotype$genotypes, "numeric")
            print(nrow(geno))
            geno_scale = standardiseGenotypes(geno)
            # geno_scale1 = scale(geno)

            pheno = scale(geno_scale[,v2] %*% snp_weight[,ethnic_i] + rnorm(nrow(geno_scale), 0, sqrt(1-opt$h2)))
            out = data.frame(rownames(pheno), rownames(pheno), pheno)
            out_file = paste0("pheno_", ethniclist[ethnic_i], "_ncausal", num_snp, ".txt")
            write.table(out, file=out_file, quote = F, col.name = F, row.name = F)
            print(out_file)
        }


    """
    )

    print(
        "--- Done in %s ---"
        % (str(datetime.timedelta(seconds=round(time.time() - start_time))))
    )
