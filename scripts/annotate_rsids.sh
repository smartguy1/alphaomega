#!/bin/bash
set -e

CLINVAR_VCF="db/clinvar_grch37.vcf.gz"
RAW_VCF="profiles/Martin/Martin.raw.vcf.gz"
TAB_RAW="/tmp/rsid_map.tab"
TAB_SORTED="/tmp/rsid_map_sorted.tab"

echo "=== Step 1: Extracting rsIDs from ClinVar GRCh37 VCF ==="
python3 -c "
import gzip
count = 0
with gzip.open('${CLINVAR_VCF}', 'rt', errors='replace') as fin:
    with open('${TAB_RAW}', 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            chrom = fields[0]
            pos = fields[1]
            info = fields[7]
            rs_num = None
            for item in info.split(';'):
                if item.startswith('RS='):
                    rs_num = item[3:]
                    break
            if rs_num and rs_num != '.':
                fout.write(chrom + '\t' + pos + '\t' + 'rs' + rs_num + '\n')
                count += 1
print(f'Extracted {count} rsID mappings')
"

echo "=== Step 2: Sorting ==="
sort -k1,1V -k2,2n "$TAB_RAW" > "$TAB_SORTED"

echo "=== Step 3: Compressing and indexing ==="
bgzip -f "$TAB_SORTED"
tabix -s1 -b2 -e2 "${TAB_SORTED}.gz"

echo "=== Step 4: Annotating VCF with rsIDs ==="
bcftools annotate \
    -a "${TAB_SORTED}.gz" \
    -c CHROM,POS,ID \
    "$RAW_VCF" \
    -Oz -o profiles/Martin/Martin.tmp.vcf.gz

mv profiles/Martin/Martin.tmp.vcf.gz "$RAW_VCF"
bcftools index -f "$RAW_VCF"

echo "=== Step 5: Counting rsIDs ==="
RSID_COUNT=$(bcftools query -f '%ID\n' "$RAW_VCF" | grep -c '^rs' || true)
echo "SUCCESS: $RSID_COUNT variants now have rsIDs!"
