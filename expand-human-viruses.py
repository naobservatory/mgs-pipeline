#!/usr/bin/env python3
import sys

# For human-viruses.tsv we want a list of all viruses that infect humans, but
# human-viruses-raw.tsv from the Virus Host DB doesn't always include all
# taxonomic children.  Make an expanded list by adding the missing ones.

in_fname, out_fname = sys.argv[1:]

raw_hv = set()
with open(in_fname) as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        raw_hv.add(int(taxid))

children = {}
with open("dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)

        if parent_taxid not in children:
            children[parent_taxid] = []
        children[parent_taxid].append(child_taxid)

hv = set()


def add_children(taxid):
    hv.add(taxid)
    for child in children.get(taxid, []):
        add_children(child)


for taxid in raw_hv:
    add_children(taxid)


# taxid -> [name]
# first name is scientific name
taxonomic_names = {}
with open("dashboard/names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)

        if taxid in hv:
            if taxid not in taxonomic_names or name_class == "scientific name":
                taxonomic_names[taxid] = name

with open(out_fname, "w") as outf:
    for taxid in sorted(hv):
        outf.write("%s\t%s\n" % (taxid, taxonomic_names[taxid]))
