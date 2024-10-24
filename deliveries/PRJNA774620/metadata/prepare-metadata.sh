cat raw_metadata.tsv | grep "Untreated faeces" | awk '{print $1}' > metadata.tsv
