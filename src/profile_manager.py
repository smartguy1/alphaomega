"""
AlphaOmega — Profile Manager

Hybrid SQLite + file storage for managing family member profiles.
  - Large artifacts (VCF, reports) stored as files in profiles/{name}/
  - Annotations and metadata indexed in db/profiles.sqlite for cross-family queries
"""

import json
import logging
import os
import sqlite3
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Optional

from src.annotator import AnnotatedVariant
from src.trait_clusterer import Domain

logger = logging.getLogger(__name__)

DEFAULT_DB_DIR = "db"
DEFAULT_PROFILES_DIR = "profiles"
PROFILES_DB_NAME = "profiles.sqlite"


@dataclass
class Profile:
    """A family member profile."""
    id: int
    name: str
    bam_path: str
    detected_build: str
    bam_hash: str
    created_at: str
    last_analyzed: str


class ProfileManager:
    """
    Manages family member profiles with hybrid SQLite + file storage.

    Usage:
        pm = ProfileManager(db_dir="db", profiles_dir="profiles")
        pm.initialize()
        pm.create_profile("martin", bam_path="/path/to/bam")
        pm.store_annotations("martin", annotated_variants)
        shared = pm.find_shared_variants(["martin", "spouse"])
    """

    def __init__(
        self,
        db_dir: str = DEFAULT_DB_DIR,
        profiles_dir: str = DEFAULT_PROFILES_DIR,
    ):
        self.db_dir = Path(db_dir)
        self.profiles_dir = Path(profiles_dir)
        self.db_path = self.db_dir / PROFILES_DB_NAME
        self._conn: Optional[sqlite3.Connection] = None

    def initialize(self):
        """Create database tables and directories if they don't exist."""
        self.db_dir.mkdir(parents=True, exist_ok=True)
        self.profiles_dir.mkdir(parents=True, exist_ok=True)

        conn = self._get_connection()

        conn.execute("""
            CREATE TABLE IF NOT EXISTS profiles (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE NOT NULL,
                bam_path TEXT,
                detected_build TEXT,
                bam_hash TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                last_analyzed TIMESTAMP
            )
        """)

        conn.execute("""
            CREATE TABLE IF NOT EXISTS profile_annotations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                profile_id INTEGER NOT NULL,
                rsid TEXT NOT NULL,
                genotype TEXT,
                chromosome TEXT,
                position INTEGER,
                snpedia_magnitude REAL,
                snpedia_repute TEXT,
                snpedia_summary TEXT,
                clinvar_significance TEXT,
                clinvar_condition TEXT,
                clinvar_review_stars INTEGER,
                pharmcat_gene TEXT,
                pharmcat_phenotype TEXT,
                domain TEXT,
                FOREIGN KEY (profile_id) REFERENCES profiles(id) ON DELETE CASCADE
            )
        """)

        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_annotations_rsid "
            "ON profile_annotations(rsid)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_annotations_domain "
            "ON profile_annotations(domain)"
        )
        conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_annotations_profile "
            "ON profile_annotations(profile_id)"
        )

        conn.commit()
        logger.info(f"Profile database initialized: {self.db_path}")

    def _get_connection(self) -> sqlite3.Connection:
        """Get or create database connection."""
        if self._conn is None:
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
            self._conn.execute("PRAGMA foreign_keys = ON")
        return self._conn

    # ── Profile CRUD ─────────────────────────────────────────────────────

    def create_profile(
        self,
        name: str,
        bam_path: str = "",
        detected_build: str = "",
        bam_hash: str = "",
    ) -> int:
        """
        Create a new profile. Returns the profile ID.
        Also creates the profile directory on disk.
        """
        conn = self._get_connection()

        # Create profile directory
        profile_dir = self.profiles_dir / name
        profile_dir.mkdir(parents=True, exist_ok=True)

        try:
            cursor = conn.execute(
                "INSERT INTO profiles (name, bam_path, detected_build, bam_hash) "
                "VALUES (?, ?, ?, ?)",
                (name, bam_path, detected_build, bam_hash),
            )
            conn.commit()
            profile_id = cursor.lastrowid
            logger.info(f"Created profile '{name}' (id={profile_id})")
            return profile_id

        except sqlite3.IntegrityError:
            # Profile already exists, update it
            conn.execute(
                "UPDATE profiles SET bam_path=?, detected_build=?, bam_hash=? "
                "WHERE name=?",
                (bam_path, detected_build, bam_hash, name),
            )
            conn.commit()
            row = conn.execute(
                "SELECT id FROM profiles WHERE name=?", (name,)
            ).fetchone()
            logger.info(f"Updated existing profile '{name}'")
            return row["id"]

    def get_profile(self, name: str) -> Optional[Profile]:
        """Get a profile by name."""
        conn = self._get_connection()
        row = conn.execute(
            "SELECT * FROM profiles WHERE name=?", (name,)
        ).fetchone()

        if row is None:
            return None

        return Profile(
            id=row["id"],
            name=row["name"],
            bam_path=row["bam_path"] or "",
            detected_build=row["detected_build"] or "",
            bam_hash=row["bam_hash"] or "",
            created_at=str(row["created_at"] or ""),
            last_analyzed=str(row["last_analyzed"] or ""),
        )

    def list_profiles(self) -> list[Profile]:
        """List all profiles."""
        conn = self._get_connection()
        rows = conn.execute("SELECT * FROM profiles ORDER BY name").fetchall()
        return [
            Profile(
                id=r["id"], name=r["name"],
                bam_path=r["bam_path"] or "",
                detected_build=r["detected_build"] or "",
                bam_hash=r["bam_hash"] or "",
                created_at=str(r["created_at"] or ""),
                last_analyzed=str(r["last_analyzed"] or ""),
            )
            for r in rows
        ]

    def delete_profile(self, name: str) -> bool:
        """Delete a profile and its annotations (files are NOT deleted)."""
        conn = self._get_connection()
        profile = self.get_profile(name)
        if profile is None:
            return False

        conn.execute("DELETE FROM profile_annotations WHERE profile_id=?", (profile.id,))
        conn.execute("DELETE FROM profiles WHERE id=?", (profile.id,))
        conn.commit()
        logger.info(f"Deleted profile '{name}' (files in profiles/{name}/ preserved)")
        return True

    # ── Annotations ──────────────────────────────────────────────────────

    def store_annotations(
        self,
        profile_name: str,
        variants: list[AnnotatedVariant],
        domain_map: Optional[dict[str, str]] = None,
    ):
        """
        Store annotated variants for a profile.
        Clears existing annotations first (full replace).

        Args:
            profile_name: Profile name
            variants: List of annotated variants
            domain_map: Optional {rsid: domain_value} mapping from trait clusterer
        """
        conn = self._get_connection()
        profile = self.get_profile(profile_name)

        if profile is None:
            raise ValueError(f"Profile '{profile_name}' not found")

        # Clear existing annotations
        conn.execute(
            "DELETE FROM profile_annotations WHERE profile_id=?",
            (profile.id,),
        )

        # Insert new annotations
        batch = []
        for v in variants:
            domain = ""
            if domain_map and v.rsid in domain_map:
                domain = domain_map[v.rsid]

            batch.append((
                profile.id,
                v.rsid,
                v.genotype,
                v.chromosome,
                v.position,
                v.snpedia.magnitude if v.snpedia else None,
                v.snpedia.repute if v.snpedia else None,
                v.snpedia.summary if v.snpedia else None,
                v.clinvar.significance if v.clinvar else None,
                "|".join(v.clinvar.conditions) if v.clinvar else None,
                v.clinvar.review_stars if v.clinvar else None,
                v.pharmcat.gene if v.pharmcat else None,
                v.pharmcat.phenotype if v.pharmcat else None,
                domain,
            ))

        conn.executemany(
            "INSERT INTO profile_annotations "
            "(profile_id, rsid, genotype, chromosome, position, "
            "snpedia_magnitude, snpedia_repute, snpedia_summary, "
            "clinvar_significance, clinvar_condition, clinvar_review_stars, "
            "pharmcat_gene, pharmcat_phenotype, domain) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            batch,
        )

        # Update last_analyzed timestamp
        conn.execute(
            "UPDATE profiles SET last_analyzed=? WHERE id=?",
            (datetime.now().isoformat(), profile.id),
        )

        conn.commit()
        logger.info(
            f"Stored {len(variants)} annotations for profile '{profile_name}'"
        )

    def get_annotations(
        self,
        profile_name: str,
        domain: Optional[str] = None,
        min_magnitude: float = 0.0,
    ) -> list[dict]:
        """
        Retrieve annotations for a profile.
        Optionally filter by domain or minimum magnitude.
        """
        conn = self._get_connection()
        profile = self.get_profile(profile_name)
        if profile is None:
            return []

        query = "SELECT * FROM profile_annotations WHERE profile_id=?"
        params: list = [profile.id]

        if domain:
            query += " AND domain=?"
            params.append(domain)

        if min_magnitude > 0:
            query += " AND (snpedia_magnitude >= ? OR clinvar_review_stars >= ?)"
            params.extend([min_magnitude, int(min_magnitude)])

        query += " ORDER BY COALESCE(snpedia_magnitude, 0) DESC"

        rows = conn.execute(query, params).fetchall()
        return [dict(r) for r in rows]

    def get_annotation_count(self, profile_name: str) -> int:
        """Get the number of annotations for a profile."""
        conn = self._get_connection()
        profile = self.get_profile(profile_name)
        if profile is None:
            return 0

        row = conn.execute(
            "SELECT COUNT(*) as cnt FROM profile_annotations WHERE profile_id=?",
            (profile.id,),
        ).fetchone()
        return row["cnt"]

    # ── Cross-family queries ─────────────────────────────────────────────

    def find_shared_variants(
        self, profile_names: list[str]
    ) -> list[dict]:
        """
        Find variants shared across multiple profiles.
        Returns variants that appear in ALL specified profiles.
        """
        conn = self._get_connection()

        profiles = [self.get_profile(n) for n in profile_names]
        profiles = [p for p in profiles if p is not None]

        if len(profiles) < 2:
            return []

        profile_ids = [p.id for p in profiles]
        placeholders = ",".join("?" * len(profile_ids))

        rows = conn.execute(
            f"SELECT rsid, genotype, snpedia_summary, clinvar_significance, "
            f"domain, COUNT(DISTINCT profile_id) as profile_count "
            f"FROM profile_annotations "
            f"WHERE profile_id IN ({placeholders}) "
            f"GROUP BY rsid "
            f"HAVING profile_count = ?",
            (*profile_ids, len(profile_ids)),
        ).fetchall()

        return [dict(r) for r in rows]

    def find_unique_variants(
        self, profile_name: str, compare_to: list[str]
    ) -> list[dict]:
        """
        Find variants unique to one profile compared to others.
        """
        conn = self._get_connection()

        target = self.get_profile(profile_name)
        if target is None:
            return []

        others = [self.get_profile(n) for n in compare_to]
        others = [p for p in others if p is not None]

        if not others:
            return self.get_annotations(profile_name)

        other_ids = [p.id for p in others]
        placeholders = ",".join("?" * len(other_ids))

        rows = conn.execute(
            f"SELECT * FROM profile_annotations "
            f"WHERE profile_id = ? "
            f"AND rsid NOT IN ("
            f"  SELECT DISTINCT rsid FROM profile_annotations "
            f"  WHERE profile_id IN ({placeholders})"
            f") ORDER BY COALESCE(snpedia_magnitude, 0) DESC",
            (target.id, *other_ids),
        ).fetchall()

        return [dict(r) for r in rows]

    def get_variant_across_profiles(self, rsid: str) -> list[dict]:
        """
        Look up a specific variant across all profiles.
        Useful for checking if family members share a trait.
        """
        conn = self._get_connection()
        rows = conn.execute(
            "SELECT pa.*, p.name as profile_name "
            "FROM profile_annotations pa "
            "JOIN profiles p ON pa.profile_id = p.id "
            "WHERE pa.rsid = ?",
            (rsid.lower(),),
        ).fetchall()

        return [dict(r) for r in rows]

    # ── File management ──────────────────────────────────────────────────

    def get_profile_dir(self, profile_name: str) -> Path:
        """Get the file directory for a profile."""
        return self.profiles_dir / profile_name

    def get_profile_file(self, profile_name: str, filename: str) -> Path:
        """Get a specific file path within a profile directory."""
        return self.get_profile_dir(profile_name) / filename

    # ── Cleanup ──────────────────────────────────────────────────────────

    def close(self):
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self):
        self.initialize()
        return self

    def __exit__(self, *args):
        self.close()
