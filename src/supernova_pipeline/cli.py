from __future__ import annotations

import argparse
from pathlib import Path

from rich.console import Console
from rich.table import Table

from .archive_matrix import build_archive_products
from .candidate_followup import run_candidate_followup
from .candidate_ledger import build_detection_queue
from .difference_upgrade import run_difference_followup, run_difference_upgrade, run_supernova_benchmark
from .extensions import run_extensions
from .galaxy_catalog import GalaxyBuildConfig, build_galaxy_master
from .observational_closure import run_observational_closure
from .pixel_search import run_pixel_search
from .quantum_precheck import run_quantum_precheck
from .quantum_submit import submit_quantum_pilots
from .utils import ensure_dir


console = Console()


def bootstrap_dirs(root_dir: Path) -> None:
    for relative in [
        "data/cache",
        "catalogs",
        "archive",
        "candidates",
        "candidate_packets",
        "closure",
        "difference_upgrade",
        "extensions",
        "figures",
        "notes",
        "outputs",
        "sed",
        "logs",
        "src",
    ]:
        ensure_dir(root_dir / relative)


def command_build_galaxies(args: argparse.Namespace) -> dict[str, Path]:
    config = GalaxyBuildConfig(
        max_distance_mpc=args.max_distance_mpc,
        max_av_mag=args.max_av_mag,
        top_n=args.top_n,
        gladeplus_radius_deg=args.gladeplus_radius_deg,
        gladeplus_z_window=args.gladeplus_z_window,
    )
    return build_galaxy_master(root_dir=args.root_dir, config=config)


def command_build_archive(args: argparse.Namespace) -> dict[str, Path]:
    return build_archive_products(
        root_dir=args.root_dir,
        galaxy_master_path=args.galaxy_master_path,
        archive_top_n=args.archive_top_n,
        radius_deg=args.radius_deg,
        min_baseline_days=args.min_baseline_days,
        max_workers=args.max_workers,
    )


def command_init_candidates(args: argparse.Namespace) -> dict[str, Path]:
    return build_detection_queue(args.epoch_pairs_path, root_dir=args.root_dir, top_per_galaxy=args.top_per_galaxy)


def command_run_pixel_search(args: argparse.Namespace) -> dict[str, Path]:
    compatibilities = tuple(part.strip() for part in str(args.compatibilities).split(",") if part.strip())
    return run_pixel_search(
        root_dir=args.root_dir,
        epoch_pairs_path=args.epoch_pairs_path,
        observation_matrix_path=args.observation_matrix_path,
        max_pairs=args.max_pairs,
        per_galaxy=args.per_galaxy,
        compatibilities=compatibilities,
        include_cross_collection=not args.same_collection_only,
        max_workers=args.max_workers,
    )


def command_run_followup(args: argparse.Namespace) -> dict[str, Path]:
    statuses = tuple(part.strip() for part in str(args.statuses).split(",") if part.strip())
    return run_candidate_followup(
        root_dir=args.root_dir,
        candidate_ledger_path=args.candidate_ledger_path,
        observation_matrix_path=args.observation_matrix_path,
        galaxy_master_path=args.galaxy_master_path,
        statuses=statuses,
        max_candidates=args.max_candidates,
        max_workers=args.max_workers,
    )


def command_run_extensions(args: argparse.Namespace) -> dict[str, Path]:
    return run_extensions(
        root_dir=args.root_dir,
        candidate_ledger_path=args.candidate_ledger_path,
        observation_matrix_path=args.observation_matrix_path,
        galaxy_master_path=args.galaxy_master_path,
    )


def command_run_observational_closure(args: argparse.Namespace) -> dict[str, Path]:
    return run_observational_closure(
        root_dir=args.root_dir,
        candidate_ledger_path=args.candidate_ledger_path,
    )


def command_run_difference_upgrade(args: argparse.Namespace) -> dict[str, Path]:
    compatibilities = tuple(part.strip() for part in str(args.compatibilities).split(",") if part.strip())
    return run_difference_upgrade(
        root_dir=args.root_dir,
        epoch_pairs_path=args.epoch_pairs_path,
        observation_matrix_path=args.observation_matrix_path,
        output_dir=args.output_dir,
        pair_subset_path=args.pair_subset_path,
        max_pairs=args.max_pairs,
        per_galaxy=args.per_galaxy,
        compatibilities=compatibilities,
        include_cross_collection=not args.same_collection_only,
        sign_mode=args.sign_mode,
        resume=bool(args.resume),
    )


def command_run_difference_followup(args: argparse.Namespace) -> dict[str, Path]:
    statuses = tuple(part.strip() for part in str(args.statuses).split(",") if part.strip())
    return run_difference_followup(
        root_dir=args.root_dir,
        detections_path=args.detections_path,
        observation_matrix_path=args.observation_matrix_path,
        output_dir=args.output_dir,
        statuses=statuses,
        max_candidates=args.max_candidates,
    )


def command_run_supernova_benchmark(args: argparse.Namespace) -> dict[str, Path]:
    return run_supernova_benchmark(
        root_dir=args.root_dir,
        epoch_pairs_path=args.epoch_pairs_path,
        observation_matrix_path=args.observation_matrix_path,
        output_dir=args.output_dir,
        resume=bool(args.resume),
    )


def command_run_quantum_precheck(args: argparse.Namespace) -> dict[str, Path]:
    return run_quantum_precheck(root_dir=args.root_dir, output_dir=args.output_dir)


def command_submit_quantum_pilots(args: argparse.Namespace) -> dict[str, Path]:
    return submit_quantum_pilots(root_dir=args.root_dir, output_dir=args.output_dir)


def command_run_pilot(args: argparse.Namespace) -> None:
    bootstrap_dirs(args.root_dir)
    galaxy_paths = command_build_galaxies(args)
    archive_paths = build_archive_products(
        root_dir=args.root_dir,
        galaxy_master_path=galaxy_paths["galaxy_master"],
        archive_top_n=args.archive_top_n,
        radius_deg=args.radius_deg,
        min_baseline_days=args.min_baseline_days,
        max_workers=args.max_workers,
    )
    candidate_paths = build_detection_queue(archive_paths["epoch_pairs"], root_dir=args.root_dir, top_per_galaxy=args.top_per_galaxy)
    pixel_paths = run_pixel_search(
        root_dir=args.root_dir,
        epoch_pairs_path=archive_paths["epoch_pairs"],
        observation_matrix_path=archive_paths["observation_matrix"],
        max_pairs=args.pixel_max_pairs,
        per_galaxy=args.pixel_per_galaxy,
        compatibilities=tuple(part.strip() for part in str(args.pixel_compatibilities).split(",") if part.strip()),
        include_cross_collection=not args.pixel_same_collection_only,
        max_workers=args.pixel_max_workers,
    )

    table = Table(title="SUPERNOVA Pilot Outputs")
    table.add_column("Artifact")
    table.add_column("Path")
    for label, path in [
        ("galaxy_master", galaxy_paths["galaxy_master"]),
        ("galaxy_scores", galaxy_paths["galaxy_scores"]),
        ("observation_matrix", archive_paths["observation_matrix"]),
        ("epoch_pairs", archive_paths["epoch_pairs"]),
        ("detection_queue", candidate_paths["detection_queue"]),
        ("candidate_ledger", candidate_paths["candidate_ledger"]),
        ("pixel_pair_summary", pixel_paths["pair_summary"]),
        ("scanned_candidate_ledger", pixel_paths["candidate_ledger"]),
    ]:
        table.add_row(label, str(path))
    console.print(table)


def build_parser() -> argparse.ArgumentParser:
    root_dir = Path("/home/primary/SUPERNOVA")
    parser = argparse.ArgumentParser(description="SUPERNOVA pipeline CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    bootstrap = sub.add_parser("bootstrap", help="Create the standard workspace directories.")
    bootstrap.set_defaults(func=lambda args: bootstrap_dirs(args.root_dir))

    build_galaxies = sub.add_parser("build-galaxies", help="Build the nearby-galaxy master catalog.")
    build_galaxies.add_argument("--top-n", type=int, default=130)
    build_galaxies.add_argument("--max-distance-mpc", type=float, default=40.0)
    build_galaxies.add_argument("--max-av-mag", type=float, default=0.5)
    build_galaxies.add_argument("--gladeplus-radius-deg", type=float, default=0.25)
    build_galaxies.add_argument("--gladeplus-z-window", type=float, default=0.003)
    build_galaxies.set_defaults(func=command_build_galaxies)

    build_archive = sub.add_parser("build-archive", help="Query MAST and build the observation matrix.")
    build_archive.add_argument("--galaxy-master-path", type=Path, default=root_dir / "catalogs" / "galaxy_master.parquet")
    build_archive.add_argument("--archive-top-n", type=int, default=25)
    build_archive.add_argument("--radius-deg", type=float, default=0.15)
    build_archive.add_argument("--min-baseline-days", type=float, default=30.0)
    build_archive.add_argument("--max-workers", type=int, default=4)
    build_archive.set_defaults(func=command_build_archive)

    init_candidates = sub.add_parser("init-candidates", help="Initialize the detection queue and candidate ledger.")
    init_candidates.add_argument("--epoch-pairs-path", type=Path, default=root_dir / "archive" / "epoch_pairs.parquet")
    init_candidates.add_argument("--top-per-galaxy", type=int, default=5)
    init_candidates.set_defaults(func=command_init_candidates)

    pixel_search = sub.add_parser("run-pixel-search", help="Download calibrated images and run the fading-source scan.")
    pixel_search.add_argument("--epoch-pairs-path", type=Path, default=root_dir / "archive" / "epoch_pairs.parquet")
    pixel_search.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    pixel_search.add_argument("--max-pairs", type=int, default=12)
    pixel_search.add_argument("--per-galaxy", type=int, default=2)
    pixel_search.add_argument("--compatibilities", type=str, default="exact,very_similar")
    pixel_search.add_argument("--max-workers", type=int, default=1)
    pixel_search.add_argument("--same-collection-only", action="store_true")
    pixel_search.set_defaults(func=command_run_pixel_search)

    followup = sub.add_parser("run-followup", help="Run all-epoch forced photometry and SED follow-up for surviving candidates.")
    followup.add_argument("--candidate-ledger-path", type=Path, default=root_dir / "candidates" / "candidate_ledger.parquet")
    followup.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    followup.add_argument("--galaxy-master-path", type=Path, default=root_dir / "catalogs" / "galaxy_master.parquet")
    followup.add_argument("--statuses", type=str, default="PASS,REVIEW")
    followup.add_argument("--max-candidates", type=int, default=0)
    followup.add_argument("--max-workers", type=int, default=1)
    followup.set_defaults(func=command_run_followup)

    extensions = sub.add_parser("run-extensions", help="Run the extended novel search methods on the current dataset.")
    extensions.add_argument("--candidate-ledger-path", type=Path, default=root_dir / "candidates" / "candidate_ledger.parquet")
    extensions.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    extensions.add_argument("--galaxy-master-path", type=Path, default=root_dir / "catalogs" / "galaxy_master.parquet")
    extensions.set_defaults(func=command_run_extensions)

    closure = sub.add_parser("run-observational-closure", help="Adjudicate disappearing-star candidates against export, dust, and systematic explanations.")
    closure.add_argument("--candidate-ledger-path", type=Path, default=root_dir / "candidates" / "candidate_ledger.parquet")
    closure.set_defaults(func=command_run_observational_closure)

    difference_upgrade = sub.add_parser("run-difference-upgrade", help="Run the empirically registered, difference-first search upgrade.")
    difference_upgrade.add_argument("--epoch-pairs-path", type=Path, default=root_dir / "archive" / "epoch_pairs.parquet")
    difference_upgrade.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    difference_upgrade.add_argument("--output-dir", type=Path, default=None)
    difference_upgrade.add_argument("--pair-subset-path", type=Path, default=None)
    difference_upgrade.add_argument("--max-pairs", type=int, default=60)
    difference_upgrade.add_argument("--per-galaxy", type=int, default=2)
    difference_upgrade.add_argument("--compatibilities", type=str, default="exact,very_similar")
    difference_upgrade.add_argument("--sign-mode", type=str, choices=["fade", "brighten", "both"], default="fade")
    difference_upgrade.add_argument("--same-collection-only", action="store_true")
    difference_upgrade.add_argument("--resume", action="store_true")
    difference_upgrade.set_defaults(func=command_run_difference_upgrade)

    difference_followup = sub.add_parser("run-difference-followup", help="Run registered difference-image forced photometry follow-up for upgraded candidates.")
    difference_followup.add_argument("--detections-path", type=Path, default=root_dir / "difference_upgrade" / "fade_candidates.parquet")
    difference_followup.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    difference_followup.add_argument("--output-dir", type=Path, default=None)
    difference_followup.add_argument("--statuses", type=str, default="PASS,REVIEW")
    difference_followup.add_argument("--max-candidates", type=int, default=0)
    difference_followup.set_defaults(func=command_run_difference_followup)

    benchmark = sub.add_parser("run-supernova-benchmark", help="Run the blinded known-supernova recovery benchmark with the upgraded detector.")
    benchmark.add_argument("--epoch-pairs-path", type=Path, default=root_dir / "archive" / "epoch_pairs.parquet")
    benchmark.add_argument("--observation-matrix-path", type=Path, default=root_dir / "archive" / "observation_matrix.parquet")
    benchmark.add_argument("--output-dir", type=Path, default=None)
    benchmark.add_argument("--resume", action="store_true")
    benchmark.set_defaults(func=command_run_supernova_benchmark)

    quantum_precheck = sub.add_parser("run-quantum-precheck", help="Design and locally preflight the pilot Braket quantum tests without submitting tasks.")
    quantum_precheck.add_argument("--output-dir", type=Path, default=None)
    quantum_precheck.set_defaults(func=command_run_quantum_precheck)

    quantum_submit = sub.add_parser("submit-quantum-pilots", help="Submit the prechecked pilot quantum tests to live Braket hardware.")
    quantum_submit.add_argument("--output-dir", type=Path, default=None)
    quantum_submit.set_defaults(func=command_submit_quantum_pilots)

    run_pilot = sub.add_parser("run-pilot", help="Run the enacted pilot workflow end to end.")
    run_pilot.add_argument("--top-n", type=int, default=130)
    run_pilot.add_argument("--archive-top-n", type=int, default=25)
    run_pilot.add_argument("--top-per-galaxy", type=int, default=5)
    run_pilot.add_argument("--max-distance-mpc", type=float, default=40.0)
    run_pilot.add_argument("--max-av-mag", type=float, default=0.5)
    run_pilot.add_argument("--gladeplus-radius-deg", type=float, default=0.25)
    run_pilot.add_argument("--gladeplus-z-window", type=float, default=0.003)
    run_pilot.add_argument("--radius-deg", type=float, default=0.15)
    run_pilot.add_argument("--min-baseline-days", type=float, default=30.0)
    run_pilot.add_argument("--max-workers", type=int, default=4)
    run_pilot.add_argument("--pixel-max-pairs", type=int, default=12)
    run_pilot.add_argument("--pixel-per-galaxy", type=int, default=2)
    run_pilot.add_argument("--pixel-compatibilities", type=str, default="exact,very_similar")
    run_pilot.add_argument("--pixel-max-workers", type=int, default=1)
    run_pilot.add_argument("--pixel-same-collection-only", action="store_true")
    run_pilot.set_defaults(func=command_run_pilot)

    for subparser in [bootstrap, build_galaxies, build_archive, init_candidates, pixel_search, followup, extensions, closure, difference_upgrade, difference_followup, benchmark, quantum_precheck, quantum_submit, run_pilot]:
        subparser.add_argument("--root-dir", type=Path, default=root_dir)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
