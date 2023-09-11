mkdir -p dashboard
cd dashboard


if [ ! -e names.dmp ] ; then
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-12-01.zip
    unzip taxdmp_2022-12-01.zip
fi
