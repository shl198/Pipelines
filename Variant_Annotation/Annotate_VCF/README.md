This is the pipeline to annotate the vcf file. Start from the gzip vcf file. Steps are:

1. Make sure nextflow is installed and in the pathway.
2. Pull the image using the command: **singularity pull docker://shl198/vep_loftee:101.0_gnomad_pLI**, this will download image file **vep_loftee_101.0_gnomad_pLI.sif**.
3. Configure the parameters in file **p01_annotate_vcf_parameters.config**.
4. Run the following command:
     
		bsub -o log.txt -q short -n 16 -M 41457280 -R "span[ptile=16]" "nextflow run p01_annotate_vcf.nf -c p01_annotate_vcf_Parameters.config  -resume"
5. Files startswith m are for further downstream analysis after vep annotation.


Attention: In the configure file, you can see all file paths have prefix /media/, that's because for singularity we set a parameter -B /:/media. this means the singularity mount the folder / in your local computer to /media in the container, so all the full path of the files in local computer would have a /media prefix in the container.