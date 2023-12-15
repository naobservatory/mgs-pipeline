import os
import json

def classify(sample, read_id, cdata):
    hvr = cdata["hvreads"].get(read_id, None)
    if not hvr:
        return False
    return hvr[0] in cdata["human_viruses"]
