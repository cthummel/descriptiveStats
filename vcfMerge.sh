#!/bin/bash

###################################################################################
# Use as follows:                                                                 #
# bash ~/bin/vcfMerge -p PATH/TO/VCF/FILES                                        # 
###################################################################################

## BEGIN SCRIPT
usage()
{
    cat << EOF

usage: $0 OPTIONS

OPTIONS can be:
    -h      Show this message
    -p      Path to files
    -s      Target File Suffix
EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

PATH=
# Check options passed in.
while getopts "h p:" OPTION
do
    case $OPTION in
        h)
            usage
            exit 1
            ;;
        p)
            PATH=$OPTARG
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

FILES=()
for filename in $PATH*.FINAL.vcf.gz; do
    FILES+=("$filename")
done

module load bcftools
bcftools merge "${FILES[@]}" 