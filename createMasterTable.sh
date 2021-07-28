#!/bin/bash
#$ -N qsubcreateMasterTable
#$ -cwd
usage() { echo "Generates MasterTable.csv in bed format. Usage: cmd [-w] wildtype bam(s) [-t] treatment bam(s) [-c] parclip bam [-g] gtf [-d] TRUE for paired-end [-o] existing output dir" 1>&2; exit 1; }
while getopts ":hw:t:c:g:d:o:" opt; do
  case ${opt} in
    h ) usage
      ;;
    w ) multiw+=(\""$OPTARG"\"\,)
      ;;
    t ) multit+=(\""$OPTARG"\"\,)
      ;;
    c ) parclip=(\""$OPTARG"\")
      ;;
    g ) gtf=(\""$OPTARG"\")
      ;;
    d ) paired=("$OPTARG")
      ;;
    o ) STR="$OPTARG"
        length=${#STR}
        last_char=${STR:length-1:1}
        [[ $last_char != "/" ]] && STR="$STR/"; :
        output=("$STR")
      ;;    
  esac
done
shift $((OPTIND -1))

multiwall=$(echo ${multiw[@]})
wlength=${#multiwall}
multiwall=${multiwall:0:wlength-1}

multitall=$(echo ${multit[@]})
tlength=${#multitall}
multitall=${multitall:0:tlength-1}

WT=$(echo "myWTpaths = c(${multiwall})")
Treatment=$(echo "myTreatmentpaths = c(${multitall})")
PARCLIP=$(echo "myPARCLIPpath = ${parclip}")
GTF=$(echo "mygtfpath = ${gtf}")
PAIRED=$(echo myispairedendread = ${paired})
OUTPUTDIR=$(echo "outputwd = \"${output}\"")


cp /home/manzouro/bin/ggCDPbamv1run.R ${output}ggCDPbamv1run.R
sed -i "s|^myWTpaths.*|${WT}|" ${output}ggCDPbamv1run.R
sed -i "s|^myTreatmentpaths.*|${Treatment}|" ${output}ggCDPbamv1run.R
sed -i "s|^myPARCLIPpath.*|${PARCLIP}|" ${output}ggCDPbamv1run.R
sed -i "s|^mygtfpath.*|${GTF}|" ${output}ggCDPbamv1run.R
sed -i "s|^myispairedendread.*|${PAIRED}|" ${output}ggCDPbamv1run.R
sed -i "s|^outputwd.*|${OUTPUTDIR}|" ${output}ggCDPbamv1run.R

module load R/3.5.3
Rscript ${output}ggCDPbamv1run.R
