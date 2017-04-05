
if [ "$#" != 4 ]; then
    echo "Usage : $0 <oriname> <db_path> <gemname> <cpus>

    <oriname>: name of original file (in dbpath)
    <dbpath>: absolute path to the folder containing oriname
    <gemname>: name to give to the genome after annotation
    <cpus>: number of cpus to run"
    exit 1
fi

oriname=$1
dbpath=$2
gemname=$3
cpus=$4

qsub -q gem -wd $dbpath -pe thread $cpus run_prokka-SGE.sh $oriname $gemname $cpus
