import pysam
import re
import math

bam_file = "/path/to/.sorted.aligned.bam"

with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam.fetch("chrX", 19990922, 19990973):

        # # Skip unmapped or missing qualities
        if read.query_qualities is None:
            continue
        
        read_seq = read.get_forward_sequence()
        # read_seq = read.query_sequence 
        if read_seq is None:
            continue
        
        # Haplotype tag (HP) if available
        try:
            haplotype = read.get_tag("HP")
        except KeyError:
            haplotype = "NA"

        base_qualities = read.query_qualities
        mapq = read.mapping_quality
        mean_bq = sum(base_qualities) / len(base_qualities)

        
        # CCA repeats
        
        if read.is_reverse:
            strand = '-'
            cca_matches = list(re.finditer(r"(TGG){15,}", read_seq))
        else:
            strand = '+'
            cca_matches = list(re.finditer(r"(CCA){15,}", read_seq))

        if cca_matches:
            longest_cca = max(cca_matches, key=lambda m: m.end() - m.start())
            span_cca = longest_cca.span()

            cca_length = (span_cca[1] - span_cca[0]) // 3

            region_bq = base_qualities[span_cca[0]:span_cca[1]]
            mean_bq_cca = sum(region_bq) / len(region_bq) if len(region_bq) > 0 else math.nan
        else:
            cca_length = 0
            mean_bq_cca = math.nan

        
        # CCG repeats
        
        if read.is_reverse:
            ccg_matches = list(re.finditer(r"(CGG){15,}", read_seq))
        else:
            ccg_matches = list(re.finditer(r"(CCG){15,}", read_seq))

        if ccg_matches:
            longest_ccg = max(ccg_matches, key=lambda m: m.end() - m.start())
            span_ccg = longest_ccg.span()

            ccg_length = (span_ccg[1] - span_ccg[0]) // 3

            region_bq = base_qualities[span_ccg[0]:span_ccg[1]]
            mean_bq_ccg = sum(region_bq) / len(region_bq) if len(region_bq) > 0 else math.nan
        else:
            ccg_length = 0
            mean_bq_ccg = math.nan

        
        # Output
        
        if cca_matches or ccg_matches:
            print(
                f"{read.query_name}\t"
                f"strand={strand}\tHP={haplotype}\t"
                f"MAPQ={mapq}\t"
                f"MeanBQ={mean_bq:.2f}\t"
                f"MeanBQ_CCA={mean_bq_cca:.2f}\tLongestCCA={cca_length}\t"
                f"MeanBQ_CCG={mean_bq_ccg:.2f}\tLongestCCG={ccg_length}"
            )
