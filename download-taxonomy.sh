mkdir -p dashboard
cd dashboard


if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2024-06-01.zip
    unzip taxdmp_2024-06-01.zip
fi
