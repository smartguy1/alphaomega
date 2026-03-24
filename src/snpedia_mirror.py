"""
AlphaOmega — SNPedia Local Mirror

Downloads the SNPedia Zenodo SQLite dump and provides a query interface.
The Zenodo dump (record 16053572) is a pre-built SQLite database (~24MB zipped)
containing ~111,728 SNPs with genotype annotations.

Since the dump schema may vary, this module:
  1. Downloads and extracts the ZIP
  2. Discovers the schema automatically
  3. Provides a normalized query interface regardless of underlying schema
  4. Can export rsID lists for VCF filtering
"""

import io
import logging
import os
import sqlite3
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

# Zenodo download URL
ZENODO_RECORD_ID = "16053572"
ZENODO_DOWNLOAD_URL = (
    f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files/SNPedia2025.zip?download=1"
)

DEFAULT_DB_DIR = "db"
DEFAULT_DB_NAME = "snpedia.sqlite"


@dataclass
class SNPediaGenotype:
    """A single genotype annotation from SNPedia."""
    rsid: str
    genotype: str
    magnitude: float
    repute: str          # "good", "bad", "neutral", or empty
    summary: str
    details: str = ""
    frequency: Optional[float] = None


@dataclass
class SNPediaEntry:
    """Full SNPedia entry for a given rsID."""
    rsid: str
    chromosome: str = ""
    position: int = 0
    gene: str = ""
    orientation: str = ""
    summary: str = ""
    genotypes: list[SNPediaGenotype] = field(default_factory=list)


class SNPediaMirror:
    """
    Interface to the local SNPedia SQLite database.

    Usage:
        mirror = SNPediaMirror(db_dir="db")
        mirror.ensure_downloaded()
        entry = mirror.lookup("rs1234")
        genotype = mirror.lookup_genotype("rs1234", "AG")
    """

    def __init__(self, db_dir: str = DEFAULT_DB_DIR):
        self.db_dir = Path(db_dir)
        self.db_path = self.db_dir / DEFAULT_DB_NAME
        self._conn: Optional[sqlite3.Connection] = None
        self._schema_info: Optional[dict] = None

    @property
    def is_available(self) -> bool:
        """Check if the database has been downloaded."""
        return self.db_path.exists()

    def ensure_downloaded(self, force: bool = False) -> bool:
        """
        Download the SNPedia database if not present.
        Returns True if download was performed, False if already existed.
        """
        if self.is_available and not force:
            logger.info(f"SNPedia database already exists: {self.db_path}")
            return False

        self.db_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Downloading SNPedia database from Zenodo...")
        logger.info(f"URL: {ZENODO_DOWNLOAD_URL}")

        try:
            response = requests.get(ZENODO_DOWNLOAD_URL, stream=True, timeout=120)
            response.raise_for_status()

            # Get total size for progress reporting
            total = int(response.headers.get("content-length", 0))
            logger.info(f"Download size: {total / 1024 / 1024:.1f} MB")

            # Download to memory and extract
            content = io.BytesIO()
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                content.write(chunk)
                downloaded += len(chunk)
                if total > 0 and downloaded % (1024 * 1024) == 0:
                    pct = (downloaded / total) * 100
                    logger.info(f"  Downloaded {pct:.0f}%")

            content.seek(0)

            # Extract ZIP
            logger.info("Extracting database...")
            with zipfile.ZipFile(content) as zf:
                # Find the SQLite file inside the ZIP
                sqlite_files = [
                    n for n in zf.namelist()
                    if n.endswith((".sqlite", ".db", ".sqlite3"))
                ]

                if not sqlite_files:
                    # Try extracting all and looking for the largest file
                    zf.extractall(str(self.db_dir))
                    logger.info(f"Extracted all files to {self.db_dir}")
                    # Find the extracted DB
                    self._find_and_rename_db()
                else:
                    # Extract the SQLite file directly
                    src_name = sqlite_files[0]
                    logger.info(f"Extracting: {src_name}")

                    with zf.open(src_name) as src, open(self.db_path, "wb") as dst:
                        dst.write(src.read())

            logger.info(f"SNPedia database ready: {self.db_path}")
            logger.info(f"Database size: {self.db_path.stat().st_size / 1024 / 1024:.1f} MB")

            # Discover schema
            self._discover_schema()
            return True

        except requests.RequestException as e:
            logger.error(f"Download failed: {e}")
            raise RuntimeError(f"Failed to download SNPedia database: {e}") from e

    def _find_and_rename_db(self):
        """Find the extracted database file and rename to standard name."""
        for f in self.db_dir.iterdir():
            if f.suffix in (".sqlite", ".db", ".sqlite3") and f.name != DEFAULT_DB_NAME:
                f.rename(self.db_path)
                logger.info(f"Renamed {f.name} → {DEFAULT_DB_NAME}")
                return

        # If no DB file found, check if there's a large file
        for f in sorted(self.db_dir.iterdir(), key=lambda x: x.stat().st_size, reverse=True):
            if f.stat().st_size > 1_000_000 and f.suffix not in (".zip",):
                try:
                    # Test if it's a valid SQLite file
                    conn = sqlite3.connect(str(f))
                    conn.execute("SELECT 1")
                    conn.close()
                    f.rename(self.db_path)
                    logger.info(f"Found and renamed database: {f.name}")
                    return
                except sqlite3.Error:
                    continue

        raise RuntimeError("Could not find SQLite database in extracted files")

    def _get_connection(self) -> sqlite3.Connection:
        """Get or create database connection."""
        if self._conn is None:
            if not self.is_available:
                raise RuntimeError(
                    "SNPedia database not found. Run ensure_downloaded() first."
                )
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def _discover_schema(self) -> dict:
        """
        Discover the database schema.
        Returns a dict mapping table names to their columns.
        """
        if self._schema_info is not None:
            return self._schema_info

        conn = self._get_connection()
        cursor = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        )
        tables = [row["name"] for row in cursor.fetchall()]

        schema = {}
        for table in tables:
            cursor = conn.execute(f"PRAGMA table_info({table})")
            columns = [row["name"] for row in cursor.fetchall()]
            schema[table] = columns

            # Log row count
            count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            logger.info(f"  Table '{table}': {count} rows, columns: {columns}")

        self._schema_info = schema
        return schema

    def get_schema(self) -> dict:
        """Return the discovered schema as {table_name: [column_names]}."""
        return self._discover_schema()

    def get_stats(self) -> dict:
        """Return database statistics."""
        conn = self._get_connection()
        schema = self._discover_schema()

        stats = {
            "db_path": str(self.db_path),
            "db_size_mb": self.db_path.stat().st_size / 1024 / 1024,
            "tables": {},
        }

        for table in schema:
            count = conn.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            stats["tables"][table] = count

        return stats

    def lookup(self, rsid: str) -> Optional[SNPediaEntry]:
        """
        Look up a SNP by rsID.
        Handles schema variations gracefully.
        """
        conn = self._get_connection()
        schema = self._discover_schema()

        # Normalize rsID
        rsid = rsid.lower().strip()
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"

        entry = SNPediaEntry(rsid=rsid)

        # Try to find the SNP in any table that has rsid-like columns
        for table, columns in schema.items():
            # Find rsID column (could be 'rsid', 'name', 'snp', 'id', etc.)
            rsid_col = self._find_column(columns, ["rsid", "name", "snp_name", "id", "snp"])
            if rsid_col is None:
                continue

            try:
                # Also try uppercase
                row = conn.execute(
                    f"SELECT * FROM {table} WHERE {rsid_col} = ? OR {rsid_col} = ?",
                    (rsid, rsid.upper())
                ).fetchone()

                if row is None:
                    # Try without 'rs' prefix
                    row_dict_keys = dict(row) if row else {}
                    continue

                row_dict = dict(row)

                # Extract fields by guessing column names
                entry.chromosome = str(
                    row_dict.get("chromosome", row_dict.get("chr", row_dict.get("chrom", "")))
                )
                entry.position = int(
                    row_dict.get("position", row_dict.get("pos", row_dict.get("bp", 0)))
                )
                entry.gene = str(
                    row_dict.get("gene", row_dict.get("gene_name", ""))
                )
                entry.summary = str(
                    row_dict.get("summary", row_dict.get("description", row_dict.get("text", "")))
                )

            except sqlite3.Error as e:
                logger.debug(f"Error querying table {table}: {e}")
                continue

        # Try to find genotype data
        entry.genotypes = self._lookup_genotypes(rsid)

        return entry if (entry.summary or entry.genotypes) else None

    def _lookup_genotypes(self, rsid: str) -> list[SNPediaGenotype]:
        """Find all genotype annotations for a given rsID."""
        conn = self._get_connection()
        schema = self._discover_schema()
        genotypes = []

        for table, columns in schema.items():
            # Look for tables with genotype-related columns
            has_genotype = self._find_column(columns, ["genotype", "geno", "alleles"])
            has_rsid = self._find_column(columns, ["rsid", "name", "snp_name", "snp", "id"])

            if has_genotype is None or has_rsid is None:
                continue

            try:
                rows = conn.execute(
                    f"SELECT * FROM {table} WHERE {has_rsid} = ? OR {has_rsid} = ?",
                    (rsid, rsid.upper())
                ).fetchall()

                for row in rows:
                    row_dict = dict(row)

                    mag_col = self._find_column(
                        list(row_dict.keys()), ["magnitude", "mag", "importance"]
                    )
                    rep_col = self._find_column(
                        list(row_dict.keys()), ["repute", "reputation", "effect"]
                    )
                    sum_col = self._find_column(
                        list(row_dict.keys()), ["summary", "description", "text", "detail"]
                    )

                    genotypes.append(SNPediaGenotype(
                        rsid=rsid,
                        genotype=str(row_dict.get(has_genotype, "")),
                        magnitude=float(row_dict.get(mag_col, 0) or 0) if mag_col else 0.0,
                        repute=str(row_dict.get(rep_col, "") or "") if rep_col else "",
                        summary=str(row_dict.get(sum_col, "") or "") if sum_col else "",
                    ))

            except sqlite3.Error as e:
                logger.debug(f"Error querying genotypes in {table}: {e}")
                continue

        return genotypes

    def lookup_genotype(
        self, rsid: str, genotype: str
    ) -> Optional[SNPediaGenotype]:
        """Look up a specific genotype for a given rsID."""
        entry = self.lookup(rsid)
        if entry is None:
            return None

        genotype_normalized = "".join(sorted(genotype.upper()))
        for gt in entry.genotypes:
            gt_normalized = "".join(sorted(gt.genotype.upper().replace("(", "").replace(")", "").replace(";", "")))
            if gt_normalized == genotype_normalized:
                return gt

        return None

    def export_rsid_list(self, output_path: str) -> int:
        """
        Export all rsIDs to a text file (one per line).
        Used for VCF filtering in variant_caller.
        Returns the number of rsIDs exported.
        """
        conn = self._get_connection()
        schema = self._discover_schema()
        rsids = set()

        for table, columns in schema.items():
            rsid_col = self._find_column(columns, ["rsid", "name", "snp_name", "snp", "id"])
            if rsid_col is None:
                continue

            try:
                rows = conn.execute(f"SELECT DISTINCT {rsid_col} FROM {table}").fetchall()
                for row in rows:
                    val = str(row[0]).strip().lower()
                    if val.startswith("rs") and val[2:].isdigit():
                        rsids.add(val)
            except sqlite3.Error:
                continue

        # Write to file
        os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
        with open(output_path, "w") as f:
            for rsid in sorted(rsids):
                f.write(f"{rsid}\n")

        logger.info(f"Exported {len(rsids)} rsIDs to {output_path}")
        return len(rsids)

    def build_position_index(self) -> dict[tuple[str, int], str]:
        """
        Build a (chromosome, position) → rsid lookup dict.
        Used by the annotator to resolve rsIDs for variants that
        bcftools outputs without IDs (just '.').

        Returns dict like {("1", 752566): "rs3094315", ...}
        """
        conn = self._get_connection()
        schema = self._discover_schema()
        index: dict[tuple[str, int], str] = {}

        for table, columns in schema.items():
            rsid_col = self._find_column(columns, ["rsid", "name", "snp_name", "snp", "id"])
            chr_col = self._find_column(columns, ["chromosome", "chr", "chrom"])
            pos_col = self._find_column(columns, ["position", "pos", "bp"])

            if not all([rsid_col, chr_col, pos_col]):
                continue

            try:
                rows = conn.execute(
                    f"SELECT {rsid_col}, {chr_col}, {pos_col} FROM {table} "
                    f"WHERE {pos_col} IS NOT NULL AND {pos_col} > 0"
                ).fetchall()

                for row in rows:
                    rsid = str(row[0]).strip().lower()
                    chrom = str(row[1]).strip().replace("chr", "")
                    try:
                        pos = int(row[2])
                    except (ValueError, TypeError):
                        continue

                    if rsid.startswith("rs") and pos > 0:
                        index[(chrom, pos)] = rsid

            except sqlite3.Error as e:
                logger.debug(f"Error building position index from {table}: {e}")
                continue

        logger.info(f"Built position index: {len(index)} SNPs with coordinates")
        return index

    def search_by_category(self, category: str) -> list[str]:
        """Search for rsIDs by category/trait keyword."""
        conn = self._get_connection()
        schema = self._discover_schema()
        results = []

        for table, columns in schema.items():
            cat_col = self._find_column(
                columns,
                ["category", "categories", "trait", "condition", "tag", "topic"]
            )
            rsid_col = self._find_column(columns, ["rsid", "name", "snp_name", "snp", "id"])

            if cat_col is None or rsid_col is None:
                continue

            try:
                rows = conn.execute(
                    f"SELECT DISTINCT {rsid_col} FROM {table} "
                    f"WHERE {cat_col} LIKE ?",
                    (f"%{category}%",)
                ).fetchall()

                for row in rows:
                    results.append(str(row[0]))
            except sqlite3.Error:
                continue

        return results

    @staticmethod
    def _find_column(
        columns: list[str], candidates: list[str]
    ) -> Optional[str]:
        """
        Find a column name that matches one of the candidates.
        Case-insensitive matching.
        """
        col_lower = {c.lower(): c for c in columns}
        for candidate in candidates:
            if candidate.lower() in col_lower:
                return col_lower[candidate.lower()]
        return None

    def close(self):
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()
