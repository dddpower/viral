val1=$(basename "$1")
seqkit rmdup -s < $1 > unique.$val1
