import sqlite3
import sys

conn = sqlite3.connect('db/snpedia.sqlite')
conn.row_factory = sqlite3.Row

# 1. Check rs149106583
row = conn.execute("SELECT rsid, content FROM snps WHERE rsid='Rs149106583'").fetchone()
if row:
    print(f"--- Rs149106583 FOUND ---")
    print(row['content'].split('\n')[:10])
else:
    print("Rs149106583 not in SNPedia.")

# 2. Find any Y-markers that are also in SNPedia
print("\nSearching for Y-markers in SNPedia...")
rows = conn.execute("SELECT rsid, content FROM snps WHERE content LIKE '%Chromosome=Y%' OR content LIKE '%[[Haplogroup R]]%' LIMIT 50").fetchall()
for r in rows:
    # Look for the 'tree=' or 'clade=' tags
    content = r['content'].lower()
    if 'haplogroup' in content:
        print(f"Found {r['rsid']}")
        # Extract the clade
        for line in r['content'].split('\n'):
            if 'clade_haplogroup=' in line.lower() or 'derived_haplogroup=' in line.lower():
                print(f"  {line.strip()}")
