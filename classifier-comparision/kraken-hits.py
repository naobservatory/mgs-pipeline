import os
import json

def classify(sample, read_id, cdata):
    hvr = cdata["hvreads"].get(read_id, None)
    if not hvr:
        return False
    
    for kraken_token in hvr[1].split():
        taxid, kmers = kraken_token.split(":")
        if taxid.isdigit() and int(taxid) in cdata[
                "human_viruses"]:
            return True
    return False
