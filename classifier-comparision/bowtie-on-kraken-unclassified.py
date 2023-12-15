import re
import os
import json

def classify(sample, read_id, cdata):
    notes = {}
    notes["collapsed"] = read_id.startswith("M_")

    hv_al = cdata["alignments"]["hv"].get(read_id)
    if not hv_al:
        notes["non-aligned"] = True
    else:
        score = 0
        length = 0
        for read_part in hv_al:
            _, _, _, _, _, part_score, part_length = read_part
            score += int(part_score)
            length += int(part_length)

        notes["hv_score"] = score
        notes["hv_length"] = length
    
    human_al = cdata["alignments"]["human"].get(read_id)
    if human_al:
        notes["human"] = True

    kraken_result = cdata["processed"]
    taxid, = re.findall(".*taxid ([0-9]+).*", kraken_result[2])
    taxid = int(taxid)

    if taxid in cdata["human_viruses"]:
        notes["kraken_assigned_hv"] = True
    elif taxid != 0:
        notes["kraken_assigned_nonhv"] = True

    return eval_notes(notes), notes

def eval_notes(notes):
    if (notes.get("human") or
        notes.get("non-aligned") or
        notes.get("kraken_assigned_nonhv")):
        return False

    score = notes["hv_score"]
    length = notes["hv_length"]
    
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
    elif not notes["collapsed"]:
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
