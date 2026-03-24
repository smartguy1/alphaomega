"""
AlphaOmega — ClinVar Local Mirror

Downloads the ClinVar VCF from NCBI FTP and provides a query interface.
ClinVar contains clinical significance annotations for human variants,
including pathogenicity, disease associations, and review status.

The VCF is parsed into a local SQLite database for efficient querying.
"""

import gzip
import logging
import os
import re
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# ClinVar FTP URLs
CLINVAR_BASE = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/"
CLINVAR_FILENAME = "clinvar.vcf.gz"
CLINVAR_URL = f"{CLINVAR_BASE}{CLINVAR_FILENAME}"

DEFAULT_DB_DIR = "db"
DEFAULT_DB_NAME = "clinvar.sqlite"


@dataclass
class ClinVarEntry:
    """A single ClinVar annotation."""
    rsid: str
    chromosome: str
    position: int
    ref: str
    alt: str
    clinical_significance: str    # e.g., "Pathogenic", "Benign"
    review_status: str            # e.g., "criteria_provided,_single_submitter"
    review_stars: int             # 0-4
    conditions: list[str]         # Disease names
    gene: str
    molecular_consequence: str
    origin: str                   # e.g., "germline"
    clinvar_id: str

    @property
    def is_pathogenic(self) -> bool:
        return "pathogenic" in self.clinical_significance.lower()

    @property
    def is_benign(self) -> bool:
        return "benign" in self.clinical_significance.lower()

    @property
    def is_uncertain(self) -> bool:
        return "uncertain" in self.clinical_significance.lower()


# Review status → star rating mapping
REVIEW_STARS = {
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_interpretation_for_the_single_variant": 0,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}


class ClinVarMirror:
    """
    Interface to the local ClinVar database.

    Usage:
        mirror = ClinVarMirror(db_dir="db")
        mirror.ensure_downloaded()
        entries = mirror.lookup("rs12345")
    """

    def __init__(self, db_dir: str = DEFAULT_DB_DIR, build: str = "GRCh38"):
        self.db_dir = Path(db_dir)
        self.db_path = self.db_dir / DEFAULT_DB_NAME
        self.vcf_path = self.db_dir / CLINVAR_FILENAME
        self.build = build
        self._conn: Optional[sqlite3.Connection] = None

    @property
    def is_available(self) -> bool:
        """Check if the database has been built."""
        return self.db_path.exists()

    def ensure_downloaded(self, force: bool = False) -> bool:
        """
        Download ClinVar VCF and build local SQLite database.
        Returns True if download/build was performed.
        """
        if self.is_available and not force:
            logger.info(f"ClinVar database already exists: {self.db_path}")
            return False

        self.db_dir.mkdir(parents=True, exist_ok=True)

        # Download VCF
        url = CLINVAR_URL
        if self.build == "GRCh37":
            url = url.replace("GRCh38", "GRCh37")

        logger.info(f"Downloading ClinVar VCF ({self.build})...")
        logger.info(f"URL: {url}")

        try:
            response = requests.get(url, stream=True, timeout=300)
            response.raise_for_status()

            total = int(response.headers.get("content-length", 0))
            logger.info(f"Download size: {total / 1024 / 1024:.1f} MB")

            with open(self.vcf_path, "wb") as f:
                downloaded = 0
                for chunk in response.iter_content(chunk_size=65536):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total > 0 and downloaded % (5 * 1024 * 1024) < 65536:
                        pct = (downloaded / total) * 100
                        logger.info(f"  Downloaded {pct:.0f}%")

            logger.info("Download complete. Building SQLite database...")
            self._build_database()
            return True

        except requests.RequestException as e:
            logger.error(f"Download failed: {e}")
            raise RuntimeError(f"Failed to download ClinVar: {e}") from e

    def _build_database(self):
        """Parse the ClinVar VCF and build the SQLite database."""
        conn = sqlite3.connect(str(self.db_path))

        conn.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                rsid TEXT,
                chromosome TEXT,
                position INTEGER,
                ref TEXT,
                alt TEXT,
                clinical_significance TEXT,
                review_status TEXT,
                review_stars INTEGER,
                conditions TEXT,
                gene TEXT,
                molecular_consequence TEXT,
                origin TEXT,
                clinvar_id TEXT
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_rsid ON variants(rsid)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_gene ON variants(gene)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_sig ON variants(clinical_significance)")

        count = 0
        batch = []
        batch_size = 5000

        logger.info(f"Parsing ClinVar VCF: {self.vcf_path}")

        open_func = gzip.open if str(self.vcf_path).endswith(".gz") else open

        with open_func(self.vcf_path, "rt", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue

                try:
                    entry = self._parse_vcf_line(line)
                    if entry is None:
                        continue

                    batch.append((
                        entry.rsid, entry.chromosome, entry.position,
                        entry.ref, entry.alt, entry.clinical_significance,
                        entry.review_status, entry.review_stars,
                        "|".join(entry.conditions), entry.gene,
                        entry.molecular_consequence, entry.origin,
                        entry.clinvar_id,
                    ))

                    if len(batch) >= batch_size:
                        conn.executemany(
                            "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                            batch,
                        )
                        conn.commit()
                        count += len(batch)
                        batch.clear()
                        if count % 50000 == 0:
                            logger.info(f"  Parsed {count} variants...")

                except Exception as e:
                    logger.debug(f"Skipping malformed line: {e}")
                    continue

        # Flush remaining
        if batch:
            conn.executemany(
                "INSERT INTO variants VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                batch,
            )
            conn.commit()
            count += len(batch)

        conn.close()
        logger.info(f"ClinVar database built: {count} variants")

    def _parse_vcf_line(self, line: str) -> Optional[ClinVarEntry]:
        """Parse a single VCF line into a ClinVarEntry."""
        fields = line.strip().split("\t")
        if len(fields) < 8:
            return None

        chrom, pos, var_id, ref, alt, qual, filt, info_str = fields[:8]

        # Parse INFO field
        info = self._parse_info(info_str)

        # Extract rsID
        rsid = ""
        rs_val = info.get("RS", "")
        if rs_val:
            rsid = f"rs{rs_val}"
        elif var_id.startswith("rs"):
            rsid = var_id

        if not rsid:
            return None  # Skip variants without rsID

        # Clinical significance
        clnsig = info.get("CLNSIG", "not_provided")
        clnsig = clnsig.replace("_", " ").replace("/", ", ")

        # Review status
        review = info.get("CLNREVSTAT", "no_assertion_provided")
        stars = REVIEW_STARS.get(review, 0)

        # Conditions
        conditions_raw = info.get("CLNDN", "not_provided")
        conditions = [
            c.replace("_", " ")
            for c in conditions_raw.split("|")
            if c != "not_provided"
        ]

        # Gene
        gene_info = info.get("GENEINFO", "")
        gene = gene_info.split(":")[0] if gene_info else ""

        # Molecular consequence
        mc = info.get("MC", "")
        consequence = ""
        if mc:
            parts = mc.split("|")
            if len(parts) >= 2:
                consequence = parts[1] if len(parts) > 1 else parts[0]

        # Origin
        origin_code = info.get("ORIGIN", "")
        origin_map = {"1": "germline", "2": "somatic", "4": "inherited"}
        origin = origin_map.get(origin_code, origin_code)

        # ClinVar ID
        clinvar_id = info.get("CLNID", info.get("CLNACC", var_id))

        return ClinVarEntry(
            rsid=rsid,
            chromosome=chrom,
            position=int(pos),
            ref=ref,
            alt=alt,
            clinical_significance=clnsig,
            review_status=review,
            review_stars=stars,
            conditions=conditions,
            gene=gene,
            molecular_consequence=consequence,
            origin=origin,
            clinvar_id=str(clinvar_id),
        )

    @staticmethod
    def _parse_info(info_str: str) -> dict[str, str]:
        """Parse a VCF INFO field into a dict."""
        info = {}
        for item in info_str.split(";"):
            if "=" in item:
                key, val = item.split("=", 1)
                info[key] = val
            else:
                info[item] = ""
        return info

    def _get_connection(self) -> sqlite3.Connection:
        """Get or create database connection."""
        if self._conn is None:
            if not self.is_available:
                raise RuntimeError(
                    "ClinVar database not found. Run ensure_downloaded() first."
                )
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def lookup(self, rsid: str) -> list[ClinVarEntry]:
        """Look up all ClinVar entries for a given rsID."""
        conn = self._get_connection()

        rsid = rsid.lower().strip()
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"

        rows = conn.execute(
            "SELECT * FROM variants WHERE rsid = ?", (rsid,)
        ).fetchall()

        results = []
        for row in rows:
            row_dict = dict(row)
            conditions = row_dict["conditions"].split("|") if row_dict["conditions"] else []
            results.append(ClinVarEntry(
                rsid=row_dict["rsid"],
                chromosome=row_dict["chromosome"],
                position=row_dict["position"],
                ref=row_dict["ref"],
                alt=row_dict["alt"],
                clinical_significance=row_dict["clinical_significance"],
                review_status=row_dict["review_status"],
                review_stars=row_dict["review_stars"],
                conditions=conditions,
                gene=row_dict["gene"],
                molecular_consequence=row_dict["molecular_consequence"],
                origin=row_dict["origin"],
                clinvar_id=row_dict["clinvar_id"],
            ))

        return results

    def lookup_gene(self, gene: str) -> list[ClinVarEntry]:
        """Look up all ClinVar entries for a given gene."""
        conn = self._get_connection()
        rows = conn.execute(
            "SELECT * FROM variants WHERE gene = ?", (gene.upper(),)
        ).fetchall()
        return [self._row_to_entry(dict(r)) for r in rows]

    def lookup_pathogenic(
        self, min_stars: int = 1
    ) -> list[ClinVarEntry]:
        """Find all pathogenic variants above a review threshold."""
        conn = self._get_connection()
        rows = conn.execute(
            "SELECT * FROM variants "
            "WHERE clinical_significance LIKE '%pathogenic%' "
            "AND review_stars >= ?",
            (min_stars,)
        ).fetchall()
        return [self._row_to_entry(dict(r)) for r in rows]

    def get_stats(self) -> dict:
        """Return database statistics."""
        conn = self._get_connection()
        total = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        pathogenic = conn.execute(
            "SELECT COUNT(*) FROM variants WHERE clinical_significance LIKE '%pathogenic%'"
        ).fetchone()[0]
        benign = conn.execute(
            "SELECT COUNT(*) FROM variants WHERE clinical_significance LIKE '%benign%'"
        ).fetchone()[0]

        return {
            "db_path": str(self.db_path),
            "db_size_mb": self.db_path.stat().st_size / 1024 / 1024 if self.db_path.exists() else 0,
            "total_variants": total,
            "pathogenic": pathogenic,
            "benign": benign,
            "uncertain": total - pathogenic - benign,
        }

    def _row_to_entry(self, row_dict: dict) -> ClinVarEntry:
        """Convert a database row to ClinVarEntry."""
        conditions = row_dict["conditions"].split("|") if row_dict["conditions"] else []
        return ClinVarEntry(
            rsid=row_dict["rsid"],
            chromosome=row_dict["chromosome"],
            position=row_dict["position"],
            ref=row_dict["ref"],
            alt=row_dict["alt"],
            clinical_significance=row_dict["clinical_significance"],
            review_status=row_dict["review_status"],
            review_stars=row_dict["review_stars"],
            conditions=conditions,
            gene=row_dict["gene"],
            molecular_consequence=row_dict["molecular_consequence"],
            origin=row_dict["origin"],
            clinvar_id=row_dict["clinvar_id"],
        )

    def close(self):
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
