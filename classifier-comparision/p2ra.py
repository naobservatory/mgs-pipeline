import re
import os
import json
import math

def classify(sample, read_id, cdata):
    hvr = cdata["hvreads"].get(read_id, None)
    if not hvr:
        return False
    if hvr[0] not in cdata["human_viruses"]:
        return False

    hv_al = cdata["alignments"]["hv"].get(read_id, None)
    if not hv_al:
        return False

    score = 0
    length = 0
    for read_part in hv_al:
        _, _, _, _, _, part_score, part_length = read_part
        score += int(part_score)
        length += int(part_length)

    return score / math.log(length) >= 22
