"""PWA submission command-line interface.

This script provides a modern CLI for submitting PWA fits with YAML configuration files.
"""

import argparse
import sys
from pathlib import Path

from .submission_manager import SubmissionManager


def create_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser.

    Returns:
        argparse.ArgumentParser: Configured argument parser.
    """
    parser = argparse.ArgumentParser(
        description="Submit PWA fits with YAML configuration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
            Examples:
            # Submit using YAML configuration
            submit my_config.yaml
            
            # Create a template configuration file
            submit --template my_template.yaml
            
            # Test mode (don't actually submit jobs)
            submit my_config.yaml --test
        """,
    )

    # Main operation mode
    parser.add_argument("config", nargs="?", help="YAML configuration file path")

    parser.add_argument(
        "--template",
        metavar="OUTPUT_PATH",
        help="Create a template YAML configuration file at the specified path",
    )

    parser.add_argument(
        "--test",
        action="store_true",
        help=(
            "Test mode: show what would be submitted without actually submitting."
            " Saves AmpTools cfg file and slurm submission file in current directory."
        ),
    )

    return parser


def main() -> int:
    """Main entry point for the CLI.

    Returns:
        int: Exit code (0 for success, non-zero for error).
    """
    parser = create_parser()
    args = parser.parse_args()

    manager = SubmissionManager()

    try:
        # Handle template creation
        if args.template:
            manager.config_manager.create_template_config(args.template)
            print(f"Template configuration created at: {args.template}")
            return 0

        # Handle YAML configuration mode
        if args.config:
            if not Path(args.config).exists():
                print(f"Error: Configuration file not found: {args.config}")
                return 1

            # Load PWA config from YAML file
            config = manager.config_manager.load_yaml_config(args.config)

            # Override test mode from command line
            if args.test:
                config.compute.test = True

            # Submit jobs
            job_ids = manager.submit_fits(config)

            if job_ids:
                print(f"\nSuccessfully submitted {len(job_ids)} jobs:")
                for job_id in job_ids:
                    print(f"  - {job_id}")
            elif config.compute.test:
                print("\nTest mode: No jobs were actually submitted.")
            else:
                print("\nNo jobs were submitted (possibly waiting for data files).")

            return 0

        else:
            parser.print_help()
            return 1

    except Exception as e:
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
