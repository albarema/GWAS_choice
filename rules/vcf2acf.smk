configfile: 'config.yaml'
CHR = range(1,23)

rule all:
    input:
        expand("temp_vcf_chr{chrom}.acf", chrom=CHR),
        "paneldir/1000gp_inds.acf.gz"
	
rule filter_snps:
    """
    Filter the vcf file - keep only biallelic sites and MAF>5%.
    Further filtering of samples should be done here (e.g.: --samples-file path/to/sampleList.txt)
    """
    input:
        vcf=config['path_2vcf'],
    output:
        vcf="vcf/chr{chrom}_filtered.vcf.gz"
    shell:
        """
        bcftools -m2 -M2 -v snps -q 0.05:minor {input.vcf} | bgzip -c > {output.vcf} ;
        bcftools index -t {output.vcf}
        """

rule get_acf:
    input:
        vcf="vcf/chr{chrom}_filtered.vcf.gz",
        epofile=config['epo_file'],
        faifile=config['fai_file'],
        lbdfile=config['lbd_file'],
    output:
        tmp1=temp("temp_vcf_chr{chrom}.acf"),
    params:
        glPath=config['galactools_path']
    shell:
         """
         {params.glPath}/glactools vcfm2acf --epo {input.epofile} --fai {input.faifile} <(bcftools view {input.vcf}) > {output.tmp1};
         """
		
rule get_panel_acf:
    input:
        tmp1="temp_vcf_chr{chrom}.acf",
	panelfile=config['pops_panel']
    output:
	tmp2=temp("temp2_vcf_chr{chrom}.acf")
    params:
        glPath=config['galactools_path']
    shell:
         """
         {params.glPath}/glactools meld -f {input.panelfile} {input.tmp1} > {output.tmp2}
         """

rule merge_galact:
    input:
        tmp2=expand("temp2_vcf_chr{chrom}.acf", chrom=CHR,allow_missing=True)
    output:
        acftp=temp("temp2_vcf_allchr.acf"),
        popfile=config['acf_file']
    params:
        glPath=config['galactools_path']
    shell:
         """
        {params.glPath}/glactools cat {input.tmp2} > {output.acftp}
        cat <({params.glPath}/glactools view -h {output.acftp} | head -1) <({params.glPath}/glactools view {output.acftp}| sort -k1,1 -k2,2g) | bgzip -c > {output.popfile}
        tabix -s 1 -b 2 -e 2 {output.popfile}
        """


