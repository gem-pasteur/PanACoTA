oriname=$1
gemname=$2
cpus=$3

source /local/gensoft2/adm/etc/profile.d/modules.sh

module load hmmer/3.1b1 aragorn/1.2.36 barrnap minced/0.1.6 blast+ prodigal infernal/1.1 ncbi_toolbox/20151127 signalp
module load prokka/1.11

prokka --outdir $oriname-prokka11Resbis --cpus $cpus $oriname"_genomic.fna"
# prokka --outdir $oriname-prokka11Res --cpus $cpus --prefix $gemname $oriname"_genomic.fna"
