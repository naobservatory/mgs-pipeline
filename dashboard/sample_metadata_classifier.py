def interpret(project, papers, bits):
    if project in papers["Rothman 2021"]["projects"]:
        sample, date, wtp, is_enriched = bits
        if wtp == "JW":
            # Rothman confirmed over email that JW = JWPCP.
            wtp = "JWPCP"

        return sample, dict(
            date=date,
            country="United States",
            state="California",
            location="Los Angeles",
            county={
                # Hyperion
                "HTP": "Los Angeles County",
                # San Jose Creek
                "SJ": "Los Angeles County",
                # Joint Water Pollution Control Plant
                "JWPCP": "Los Angeles County",
                # Orange County
                "OC": "Orange County",
                # Point Loma
                "PL": "San Diego County",
                # South Bay
                "SB": "San Diego County",
                # North City
                "NC": "San Diego County",
                # Escondido Hale Avenue Resource Recovery Facility
                "ESC": "San Diego County",
            }[wtp],
            fine_location=wtp,
            enrichment="panel" if is_enriched == "1" else "viral")
    elif project in papers["Crits-Christoph 2021"]["projects"]:
        sample, municipality, date, method, sequencing = bits
        return sample, dict(
            date=date,
            country="United States",
            state="California",
            location="San Francisco",
            county={
                "Berkeley": "Alameda County",
                "Marin": "Marin County",
                "Oakland": "Alameda County",
                "SF": "San Francisco County",
            }[municipality],
            fine_location=municipality,
            method=method,
            enrichment="panel" if sequencing == "enriched" else "viral")
    elif project in papers["Brumfield 2022"]["projects"]:
        sample, na_type, date = bits
        return sample, dict(
            date=date,
            country="United States",
            location="Maryland",
            fine_location="Manhole",
            na_type=na_type)
    elif project in papers["Bengtsson-Palme 2016"]["projects"]:
        sample, location, site = bits
        return sample, dict(
            date="2012-09",
            country="Sweden",
            location=location,
            fine_location=site)
    elif project in papers["Brinch 2020"]["projects"]:
        sample, loc, date = bits
        return sample, dict(
            date=date,
            country="Denmark",
            location="Copenhagen",
            fine_location=loc)
    elif project in papers["Spurbeck 2023"]["projects"]:
        sample, loc, date = bits
        return sample, dict(
            date=date,
            country="United States",
            state="Ohio",
            location="Ohio",
            # https://github.com/naobservatory/mgs-pipeline/issues/9
            county={
                "A": "Summit County",
                "B": "Trumbull County",
                "C": "Lucas County",
                "D": "Lawrence County",
                "E": "Sandusky County",
                "F": "Franklin County",
                "G": "Licking County",
                "H": "Franklin County",
                "I": "Greene County",
                "J": "Montgomery County",
            }[loc],
            fine_location=loc,
            enrichment="viral",
            method={
                "A": "AB",
                "B": "AB",
                "C": "C",
                "D": "D",
                "E": "EFGH",
                "F": "EFGH",
                "G": "EFGH",
                "H": "EFGH",
                "I": "IJ",
                "J": "IJ",
            }[loc])

    elif project in papers["Munk 2022"]["projects"]:
        sample, country, location, date = bits
        return sample, dict(
            date=date,
            country=country,
            location=location)
    elif project in papers["Petersen 2015"]["projects"]:
        sample, country, city = bits
        return sample, dict(
            country=country,
            location=city,
            # Per Supplementary Table 7 they're all one of 23-08-2013,
            # 27-06-2013, 29-08-2013, 24-08-2013.  But the mapping
            # between samples and dates doesn't seem to be in the
            # paper.
            date="Summer 2013")
    elif project in papers["Maritz 2019"]["projects"]:
        sample = bits[0]
        return sample, dict(
            country="United States",
            location="New York City",
            # Paper has "17 raw sewage samples collected from 14 DEP
            # wastewater treatment plants from the five NYC boroughs in
            # November 2014".
            date="2014-11")
    elif project in papers["Fierer 2022"]["projects"]:
        sample = bits[0]
        return sample, dict(
            country="United States",
            location="Boulder, CO",
            # I can't find metadata on which samples are from which
            # days or which spots on campus.
            date="2020-09")
    elif project in papers["Ng 2019"]["projects"]:
        sample, stage, date = bits
        return sample, dict(
            country="Singapore",
            location="Singapore",
            date=date,
            fine_location={
                "Effluent from Membrane Biorector (MBR)": "MBR",
                "Effluent from Primary Settling Tank (PST)": "PST",
                "Effluent from Secondary Settling Tank (SST)": "SST",
                "Effluent from Wet Well (WW)": "WW",
                "Influent": "Influent",
                "Sludge (SLUDGE)": "Sludge",
            }[stage])
    elif project in papers["Hendriksen 2019"]["projects"]:
        sample, date, cluster = bits
        return sample, dict(
            country="Kenya",
            location="Kibera",
            fine_location=cluster,
            date=date)
    elif project in papers["Yang 2020"]["projects"]:
        sample, city = bits
        return sample, dict(
            country="China",
            location=city,
            date="2018",
            enrichment="viral")
    elif project in papers["Wang 2022"]["projects"]:
        sample, date, hospital = bits
        return sample, dict(
            country="Saudi Arabia",
            location="Jeddah",
            date=date,
            fine_location=hospital)
    elif project in papers["Johnson 2023"]["projects"]:
        sample, _, _ = bits
        return sample, dict(
            country="United States",
            location="Municipal",
            date="2022",
            enrichment="viral")
    elif project in papers["Cui 2023"]["projects"]:
        sample, = bits
        return sample, dict(
            country="China",
            # Also possible this was Changchun
            location="Harbin",
            # Also possible this was 2020-10-15
            date="2022-10-19")
    elif project in papers["McCall 2023"]["projects"]:
        sample, reads, method, designation = bits
        designation = designation.replace("Rice_ViroCap_", "")
        fine_location, sample_letter = designation.split(" ")

        targeted_DE = sample_letter.endswith("_targetedDE")
        sample_letter = sample_letter.split("_")[0]

        fine_location = {
            "ECC": "Child Care",
            "HS": "High School",
            "Jail": "Jail",
            "Shel": "Shelter",
            "NH": "Nursing Home",
            "WWTP": "WTP",
        }[fine_location] + "-" + sample_letter

        if targeted_DE:
            fine_location += "-DE"

        enrichment = {
            "Targeted-Capture": "panel",
            "WGS": "viral",
        }[method]

        return sample, dict(
            country="United States",
            state="Texas",
            county="Harris County",
            location="Houston",
            fine_location=fine_location,
            # Paper doesn't seem to list a date.  Guess 2022 while we find out.
            # https://twist.com/a/197793/ch/619193/t/4469435/
            date="2022",
            enrichment=enrichment)
    else:
        raise Exception("Metadata format for %s unknown" % project)
