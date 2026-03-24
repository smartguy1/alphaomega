"""
AlphaOmega — Command-Line Interface

Click-based CLI for the genomic analysis pipeline.

Commands:
  setup    — Download reference genomes and databases
  analyze  — Run full pipeline: BAM → VCF → annotate → report
  report   — Generate report from existing profile
  compare  — Compare two profiles (shared/unique/relatedness)
  status   — Show database versions, profiles, disk usage
  export   — Export profile for FragmentOfUs biographer
"""

import json
import logging
import os
import sys

import click
import yaml

logger = logging.getLogger("alphaomega")


def _setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )


def _load_config(config_path: str) -> dict:
    """Load YAML configuration."""
    if os.path.exists(config_path):
        with open(config_path) as f:
            return yaml.safe_load(f) or {}
    return {}


@click.group()
@click.option(
    "--config", "-c",
    default="config/default_config.yaml",
    help="Path to configuration file",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.pass_context
def cli(ctx: click.Context, config: str, verbose: bool):
    """🧬 AlphaOmega — Local Genomic Analysis Pipeline"""
    _setup_logging(verbose)
    ctx.ensure_object(dict)
    ctx.obj["config"] = _load_config(config)
    ctx.obj["verbose"] = verbose


# ─── Setup ───────────────────────────────────────────────────────────────────

@cli.command()
@click.option(
    "--grch37", is_flag=True, help="Also download GRCh37 reference",
)
@click.option(
    "--skip-reference", is_flag=True, help="Skip reference genome download",
)
@click.pass_context
def setup(ctx: click.Context, grch37: bool, skip_reference: bool):
    """Download reference genomes and knowledge databases."""
    config = ctx.obj["config"]

    click.echo("🧬 AlphaOmega Setup")
    click.echo("=" * 40)

    # Download knowledge databases
    click.echo("\n📦 Setting up knowledge databases...")

    from src.snpedia_mirror import SNPediaMirror
    from src.clinvar_mirror import ClinVarMirror
    from src.pharmcat_runner import PharmCATRunner

    db_dir = config.get("paths", {}).get("database_dir", "db")

    # SNPedia
    click.echo("\n[1/3] SNPedia database...")
    with SNPediaMirror(db_dir=db_dir) as mirror:
        if mirror.ensure_downloaded():
            stats = mirror.get_stats()
            click.echo(f"  ✅ Downloaded ({stats['db_size_mb']:.1f} MB)")
        else:
            click.echo("  ✅ Already up to date")

    # ClinVar
    click.echo("\n[2/3] ClinVar database...")
    with ClinVarMirror(db_dir=db_dir) as mirror:
        if mirror.ensure_downloaded():
            stats = mirror.get_stats()
            click.echo(f"  ✅ Downloaded ({stats['total_variants']} variants)")
        else:
            click.echo("  ✅ Already up to date")

    # PharmCAT
    click.echo("\n[3/3] PharmCAT...")
    runner = PharmCATRunner(db_dir=db_dir)
    if runner.ensure_installed():
        click.echo("  ✅ Downloaded")
    else:
        click.echo("  ✅ Already installed")

    if not runner.java_available:
        click.echo("  ⚠️  Java not found — PharmCAT requires Java 17+")

    # Reference genomes
    if not skip_reference:
        click.echo("\n🧬 Reference genomes:")
        click.echo("  Run: bash scripts/setup_reference.sh" +
                    (" --both" if grch37 else ""))
        click.echo("  (This downloads ~3GB per genome build)")

    # Initialize profile database
    from src.profile_manager import ProfileManager
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")
    with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir):
        click.echo("\n📁 Profile database initialized")

    click.echo("\n✅ Setup complete!")


# ─── Analyze ─────────────────────────────────────────────────────────────────

@cli.command()
@click.option("--bam", required=True, type=click.Path(exists=True), help="Path to BAM file")
@click.option("--profile", "-p", required=True, help="Profile name")
@click.option("--skip-calling", is_flag=True, help="Skip variant calling (use cached VCF)")
@click.option("--force", is_flag=True, help="Force re-analysis even if cached")
@click.pass_context
def analyze(ctx: click.Context, bam: str, profile: str, skip_calling: bool, force: bool):
    """Run full analysis pipeline: BAM → VCF → annotate → report."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")
    ref_dir = config.get("paths", {}).get("reference_dir", "reference")
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")
    output_dir = config.get("paths", {}).get("output_dir", "output")
    vc_config = config.get("variant_calling", {})
    ann_config = config.get("annotation", {})

    click.echo(f"🧬 Analyzing: {bam}")
    click.echo(f"   Profile: {profile}")

    # Step 1: Variant Calling
    if not skip_calling:
        click.echo("\n[1/4] Calling variants...")
        from src.variant_caller import call_variants

        result = call_variants(
            bam_path=bam,
            output_dir=os.path.join(profiles_dir, profile),
            reference_dir=ref_dir,
            profile_name=profile,
            threads=vc_config.get("threads", 4),
            min_quality=vc_config.get("min_quality", 20),
            min_depth=vc_config.get("min_depth", 5),
            db_dir=db_dir,
            force=force,
        )

        if not result.success:
            click.echo(f"  ❌ Variant calling failed: {result.error}")
            sys.exit(1)

        click.echo(f"  ✅ {result.total_variants} variants called ({result.build.value})")
        vcf_path = result.vcf_path
        build = result.build.value
        bam_hash = result.bam_hash
    else:
        # Use cached VCF
        vcf_path = os.path.join(profiles_dir, profile, f"{profile}.filtered.vcf.gz")
        if not os.path.exists(vcf_path):
            vcf_path = os.path.join(profiles_dir, profile, f"{profile}.raw.vcf.gz")
        if not os.path.exists(vcf_path):
            click.echo("  ❌ No cached VCF found. Run without --skip-calling first.")
            sys.exit(1)
        build = ""
        bam_hash = ""
        click.echo(f"  ✅ Using cached VCF: {vcf_path}")

    # Step 2: Create/update profile
    click.echo("\n[2/4] Setting up profile...")
    from src.profile_manager import ProfileManager

    with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir) as pm:
        pm.create_profile(profile, bam_path=bam, detected_build=build, bam_hash=bam_hash)
        click.echo(f"  ✅ Profile '{profile}' ready")

        # Step 3: Annotate
        click.echo("\n[3/4] Annotating variants...")
        from src.snpedia_mirror import SNPediaMirror
        from src.clinvar_mirror import ClinVarMirror
        from src.annotator import Annotator
        from src.avatar import AVATAR_MARKERS

        with SNPediaMirror(db_dir=db_dir) as snpedia:
            with ClinVarMirror(db_dir=db_dir) as clinvar:
                annotator = Annotator(
                    snpedia=snpedia if snpedia.is_available else None,
                    clinvar=clinvar if clinvar.is_available else None,
                    min_magnitude=ann_config.get("snpedia_min_magnitude", 0),
                    min_review_stars=ann_config.get("clinvar_min_review_stars", 0),
                    force_include_rsids=AVATAR_MARKERS,
                )

                variants = annotator.annotate_vcf(vcf_path)
                click.echo(f"  ✅ {len(variants)} variants annotated")

                # Cluster
                from src.trait_clusterer import TraitClusterer
                clusterer = TraitClusterer()
                clusters = clusterer.cluster(variants)

                domain_map = {}
                for domain, cluster in clusters.items():
                    for v in cluster.variants:
                        domain_map[v.rsid] = domain.value

                # Store annotations
                pm.store_annotations(profile, variants, domain_map=domain_map)

        # Step 4: Generate report
        click.echo("\n[4/4] Generating report...")
        from src.report_generator import ReportGenerator

        gen = ReportGenerator()
        report_path = gen.generate_report(
            profile_name=profile,
            variants=variants,
            output_path=os.path.join(profiles_dir, profile, "report.html"),
            build=build,
        )
        click.echo(f"  ✅ Report saved: {report_path}")

        # Also save JSON summary
        summary = gen.generate_summary_json(profile, clusters)
        json_path = os.path.join(profiles_dir, profile, "summary.json")
        gen.save_summary_json(summary, json_path)

    click.echo(f"\n🎉 Analysis complete for '{profile}'!")
    click.echo(f"   Open: {report_path}")


# ─── Report ──────────────────────────────────────────────────────────────────

@cli.command()
@click.option("--profile", "-p", required=True, help="Profile name")
@click.option("--output", "-o", default=None, help="Output path (default: profiles/{name}/report.html)")
@click.pass_context
def report(ctx: click.Context, profile: str, output: str):
    """Generate report from existing profile annotations."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")

    from src.profile_manager import ProfileManager
    from src.report_generator import ReportGenerator
    from src.annotator import AnnotatedVariant, SNPediaAnnotation, ClinVarAnnotation
    from src.trait_clusterer import TraitClusterer

    with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir) as pm:
        p = pm.get_profile(profile)
        if p is None:
            click.echo(f"❌ Profile '{profile}' not found")
            sys.exit(1)

        annots = pm.get_annotations(profile)
        if not annots:
            click.echo(f"❌ No annotations found for '{profile}'. Run 'analyze' first.")
            sys.exit(1)

        # Reconstruct AnnotatedVariant objects from DB rows
        variants = []
        for a in annots:
            v = AnnotatedVariant(
                rsid=a["rsid"],
                genotype=a.get("genotype", ""),
                chromosome=a.get("chromosome", ""),
                position=a.get("position", 0),
            )
            if a.get("snpedia_summary"):
                v.snpedia = SNPediaAnnotation(
                    magnitude=a.get("snpedia_magnitude", 0) or 0,
                    repute=a.get("snpedia_repute", "") or "",
                    summary=a.get("snpedia_summary", "") or "",
                )
            if a.get("clinvar_significance"):
                conditions = (a.get("clinvar_condition", "") or "").split("|")
                v.clinvar = ClinVarAnnotation(
                    significance=a.get("clinvar_significance", ""),
                    conditions=[c for c in conditions if c],
                    review_stars=a.get("clinvar_review_stars", 0) or 0,
                )
            variants.append(v)

        if output is None:
            output = os.path.join(profiles_dir, profile, "report.html")

        gen = ReportGenerator()
        path = gen.generate_report(
            profile_name=profile,
            variants=variants,
            output_path=output,
            build=p.detected_build,
        )
        click.echo(f"✅ Report generated: {path}")


# ─── Compare ─────────────────────────────────────────────────────────────────

@cli.command()
@click.option("--profiles", "-p", required=True, help="Comma-separated profile names")
@click.option("--relatedness", "-r", is_flag=True, help="Run filiation/relatedness check")
@click.pass_context
def compare(ctx: click.Context, profiles: str, relatedness: bool):
    """Compare variants across family profiles."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")
    names = [n.strip() for n in profiles.split(",")]

    if len(names) < 2:
        click.echo("❌ Need at least 2 profiles to compare")
        sys.exit(1)

    from src.profile_manager import ProfileManager

    with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir) as pm:
        # Shared variants
        click.echo(f"🔬 Comparing profiles: {', '.join(names)}")
        shared = pm.find_shared_variants(names)
        click.echo(f"\n📊 Shared variants: {len(shared)}")

        if shared:
            click.echo("\nTop shared findings:")
            for v in shared[:10]:
                summary = v.get("snpedia_summary", "") or v.get("clinvar_significance", "")
                click.echo(f"  {v['rsid']:12s} {summary[:60]}")

        # Unique per profile
        for name in names:
            others = [n for n in names if n != name]
            unique = pm.find_unique_variants(name, others)
            click.echo(f"\n🔹 Unique to {name}: {len(unique)} variants")

        # Relatedness analysis
        if relatedness and len(names) == 2:
            click.echo(f"\n🧬 Relatedness analysis: {names[0]} ↔ {names[1]}")
            from src.relatedness import analyze_relatedness

            db_path = os.path.join(db_dir, "profiles.sqlite")
            result = analyze_relatedness(names[0], names[1], db_path)

            click.echo(f"   Relationship: {result.relationship.value}")
            click.echo(f"   Confidence:   {result.confidence:.1%}")
            click.echo(f"   Shared SNPs:  {result.shared_snps}")
            click.echo(f"   IBS0/1/2:     {result.ibs0:.1%} / {result.ibs1:.1%} / {result.ibs2:.1%}")
            click.echo(f"   Mendelian errors: {result.mendelian_errors} ({result.mendelian_error_rate:.2%})")
            click.echo(f"   Details: {result.details}")


# ─── Export ──────────────────────────────────────────────────────────────────

@cli.command(name="export")
@click.option("--profile", "-p", required=True, help="Profile name")
@click.option("--output", "-o", default=None, help="Output JSON path")
@click.pass_context
def export_cmd(ctx: click.Context, profile: str, output: str):
    """Export profile for FragmentOfUs AI biographer."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")

    from src.profile_manager import ProfileManager
    from src.annotator import AnnotatedVariant, SNPediaAnnotation
    from src.export import export_for_biographer, save_export

    with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir) as pm:
        annots = pm.get_annotations(profile)
        if not annots:
            click.echo(f"❌ No annotations for '{profile}'. Run 'analyze' first.")
            sys.exit(1)

        # Reconstruct variants
        variants = []
        for a in annots:
            v = AnnotatedVariant(rsid=a["rsid"], genotype=a.get("genotype", ""))
            if a.get("snpedia_summary"):
                v.snpedia = SNPediaAnnotation(
                    magnitude=a.get("snpedia_magnitude", 0) or 0,
                    repute=a.get("snpedia_repute", "") or "",
                    summary=a.get("snpedia_summary", "") or "",
                )
            variants.append(v)

        data = export_for_biographer(profile, variants, profile_manager=pm)

        if output is None:
            output = os.path.join(profiles_dir, profile, "biographer_export.json")

        path = save_export(data, output)
        click.echo(f"✅ Exported for FragmentOfUs: {path}")
        click.echo(f"   {len(data.get('genomic_highlights', {}))} domains, "
                    f"{sum(d['finding_count'] for d in data.get('genomic_highlights', {}).values())} findings")


# ─── Status ──────────────────────────────────────────────────────────────────

@cli.command()
@click.pass_context
def status(ctx: click.Context):
    """Show database versions, profiles, and disk usage."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")
    profiles_dir = config.get("paths", {}).get("profiles_dir", "profiles")

    click.echo("🧬 AlphaOmega Status")
    click.echo("=" * 40)

    # Databases
    click.echo("\n📦 Knowledge Databases:")
    from src.snpedia_mirror import SNPediaMirror
    from src.clinvar_mirror import ClinVarMirror
    from src.pharmcat_runner import PharmCATRunner

    mirror = SNPediaMirror(db_dir=db_dir)
    if mirror.is_available:
        stats = mirror.get_stats()
        click.echo(f"  SNPedia:  ✅ {stats['db_size_mb']:.1f} MB")
    else:
        click.echo("  SNPedia:  ❌ Not downloaded")

    clinvar = ClinVarMirror(db_dir=db_dir)
    if clinvar.is_available:
        stats = clinvar.get_stats()
        click.echo(f"  ClinVar:  ✅ {stats['total_variants']} variants ({stats['db_size_mb']:.1f} MB)")
    else:
        click.echo("  ClinVar:  ❌ Not downloaded")

    runner = PharmCATRunner(db_dir=db_dir)
    click.echo(f"  PharmCAT: {'✅' if runner.is_available else '❌'} "
               f"{'Java ' + ('✅' if runner.java_available else '❌') if runner.is_available else 'Not downloaded'}")

    # Profiles
    click.echo("\n👤 Profiles:")
    from src.profile_manager import ProfileManager
    try:
        with ProfileManager(db_dir=db_dir, profiles_dir=profiles_dir) as pm:
            profiles = pm.list_profiles()
            if profiles:
                for p in profiles:
                    count = pm.get_annotation_count(p.name)
                    click.echo(f"  {p.name:15s} {count:>5d} annotations  "
                               f"(build={p.detected_build or '?'}, "
                               f"analyzed={p.last_analyzed or 'never'})")
            else:
                click.echo("  No profiles found")
    except Exception:
        click.echo("  Profile database not initialized (run 'setup' first)")

    # Disk usage
    click.echo("\n💾 Disk Usage:")
    for dirname in [db_dir, profiles_dir, "reference"]:
        if os.path.exists(dirname):
            total = sum(
                f.stat().st_size
                for f in __import__("pathlib").Path(dirname).rglob("*")
                if f.is_file()
            )
            click.echo(f"  {dirname + '/':15s} {total / 1024 / 1024:.1f} MB")
        else:
            click.echo(f"  {dirname + '/':15s} (not created)")


# ─── Update DB ───────────────────────────────────────────────────────────────

@cli.command(name="update-db")
@click.pass_context
def update_db(ctx: click.Context):
    """Refresh SNPedia and ClinVar databases."""
    config = ctx.obj["config"]
    db_dir = config.get("paths", {}).get("database_dir", "db")

    from src.snpedia_mirror import SNPediaMirror
    from src.clinvar_mirror import ClinVarMirror

    click.echo("🔄 Updating databases...")

    with SNPediaMirror(db_dir=db_dir) as mirror:
        mirror.ensure_downloaded(force=True)
        click.echo("  ✅ SNPedia updated")

    with ClinVarMirror(db_dir=db_dir) as mirror:
        mirror.ensure_downloaded(force=True)
        click.echo("  ✅ ClinVar updated")

    click.echo("✅ Databases updated!")


def main():
    """Entry point."""
    cli()


if __name__ == "__main__":
    main()
