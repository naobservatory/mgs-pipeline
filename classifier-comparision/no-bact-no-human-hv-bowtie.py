import re
import os
import json

def classify(sample, read_id, cdata):
    human_al = cdata["alignments"]["human"].get(read_id, None)
    if human_al:
        return False

    kraken_result = cdata["processed"]
    taxid, = re.findall(".*taxid ([0-9]+).*", kraken_result[2])
    taxid = int(taxid)

    if taxid == 9606:  # Kraken assigned human
        return False
    if cdata["is_bacterial"](taxid):
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

    bend1_length = 40
    bend1_score = 60

    bend2_length = 80
    bend2_score = 95

    max_length = 200
    max_score = 100

    min_length_cutoff = 28
    uncollapsed_score = 269

    if length < min_length_cutoff:
        return False
    elif not read_id.startswith("M_"):
        return score >= uncollapsed_score
    elif length < bend1_length:
        return score >= bend1_score
    elif length < bend2_length:
        return score >= (
            (length - bend1_length) / (bend2_length - bend1_length)
            * (bend2_score - bend1_score)) + bend1_score
    else:
        return score >= (
            (length - bend2_length) / (max_length - bend2_length)
            * (max_score - bend2_score)) + bend2_score
