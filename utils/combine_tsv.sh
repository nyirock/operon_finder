awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' *.tsv >all.tsv
