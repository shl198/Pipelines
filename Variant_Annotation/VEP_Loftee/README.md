* Dockerfile: commands to build the docker image that have VEP and plugins including LOFTEE.
* gnomad.v2.1.1.oe_lof.by_gene.txt: gnomad recommend to replace ExACpLI loss of
function intolerant score with observed/expected (oe) metric. Details are here:https://macarthurlab.org/2018/10/17/gnomad-v2-1/. So I downloaded the constraint file from gnomad and extract two columns (gene symbol and oe_lof).
* Dockerfile_v1: commands to build the docker image that have VEP and plugins including LOFTEE, replace ExACpLI with gnomad. Better to use this compaired to Dockerfile.

If you don't want to build the container by yourself, you can directly download from docker://shl198/vep_loftee:101.0_gnomad_pLI.

Run the vep annotation pipeline in HPC using singularity
--------------------------------------------------------

1. Pull the image using the command: **singularity pull docker://shl198/vep_loftee:101.0_gnomad_pLI**, this will download image file **vep_loftee_101.0_gnomad_pLI.sif**.
2. Download GRCh38 reference genome from ensembl(http://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz). <br />
3. Check this website 'https://github.com/konradjk/loftee/tree/grch38', download the loftee software. Then follow the instructions in readme to download file **human_ancestor_fa.gz**, and bigwig file **gerp_conservation_scores.homo_sapiens.GRCh38.bw**, and PhyloCSF file **phylocsf_gerp.sql**, finaly put the files into a loftee folder <br />.
4. Download the vep cache folder from 'http://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_vep_101_GRCh38.tar.gz'.
5. Download this folder <br />
   Example command of running vep: bolded parameters are the ones you need to change <br /> 

        bsub -q short -o log.txt "singularity run -B /:/media /path/to/aingularity/image/file \
        vep -i **/media/path/to/sample.vcf** \
        --plugin LoF,loftee_path:**/media/path/to/loftee**,human_ancestor_fa:**/media/path/to/human_ancestor_fa.gz**,gerp_bigwig:**/media/path/to/gerp_conservation_scores.homo_sapiens.GRCh38.bw**,conservation_file:**/media/path/to/loftee_grch38/phylocsf_gerp.sql** \
		--dir_plugins **/media/path/to/loftee_grch38** \
		--custom **/media/hpc/grid/wip_drm_targetsciences/projects/gnomAD/gnomad_v3/gnomad.genomes.r3.0.sites.chr16.vcf.bgz,gnomADg,vcf,exact,0,AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_nfe,AF_oth** \
		--plugin Carol \
		--plugin Condel,/opt/.vep/Plugins/config/Condel/config,b \
		--plugin ExACpLI,/opt/gnomad.v2.1.1.oe_lof.by_gene.txt \
		--plugin LoFtool,/opt/.vep/Plugins/LoFtool_scores.txt \
		--plugin CADD,**/media/path/to/CADDv1_6/gnomad.genomes.r3.0.indel.tsv.gz \
		-o /media/path/to/output.vcf.gz** \
		--cache --force_overwrite --buffer_size 10000 \
		--species homo_sapiens --assembly GRCh38 \
		--dir /media/path/to/vep_db_101_GRCH38 \
		--offline --variant_class --fork 1 --hgvs -e \
		--fa /media/path/to/GRCh38_vep101.fa \
		--minimal  --compress_output gzip \
		--allele_number --check_existing --vcf" <br />

There's one parameter for singularity that I need to explain here, the -B /:/media, this means the singularity mount the folder / in your local computer to /media in the container, so all the full path of the files in local computer would have a /media prefix in the container.