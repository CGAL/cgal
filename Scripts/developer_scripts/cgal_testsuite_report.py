#!/usr/bin/env python3
import os
import json
from typing import Dict, List
from dataclasses import dataclass
from datetime import datetime
import subprocess
import re
import requests

CGAL_SERVER_URL = "https://cgal.geometryfactory.com/CGAL"
LATEST_VERSION_URL = f"{CGAL_SERVER_URL}/Releases/LATEST"
JSON_DATA_URL_TEMPLATE = f"{
    CGAL_SERVER_URL}/testsuite/CGAL-{{version}}/search_index.json"
TESTSUITE_URL_TEMPLATE = f"{
    CGAL_SERVER_URL}/testsuite/results-{{version}}.shtml"
TIMEOUT_DURATION = 10


@dataclass
class TPLInfo:
    name: str
    version: str
    status: str


@dataclass
class PlatformInfo:
    name: str
    debug: str
    os: str
    tester: str
    compiler: str
    tpl_info: List[TPLInfo]


def fetch_data_from_url(url: str) -> str:
    """Fetch data from a given URL."""
    response = requests.get(url, timeout=TIMEOUT_DURATION)
    response.raise_for_status()
    return response.text.strip()


def get_latest_version() -> str:
    """Return latest CGAL version from LATEST (CGAL-<version>.tar.gz)"""
    tarball_name = fetch_data_from_url(LATEST_VERSION_URL)
    match = re.match(r'CGAL-([^.]+\.[^-]+-[^-]+-\d+)', tarball_name)
    if not match:
        raise ValueError(f"Unexpected tarball name format: {tarball_name}")
    return match.group(1)


def fetch_json_data(version: str) -> Dict:
    """Fetch JSON data for the given CGAL testsuite."""
    url = JSON_DATA_URL_TEMPLATE.format(version=version)
    json_data = fetch_data_from_url(url)
    return json.loads(json_data)


def analyze_tpl_data(json_data: Dict) -> List[PlatformInfo]:
    """Analyze TPL data from JSON and return a list of PlatformInfo."""
    platforms_info = []
    for platform in json_data.get('platforms', []):
        tpl_list = [
            TPLInfo(
                name=item.get('name', 'Unknown'),
                version=item.get('version', 'N/A'),
                status=item.get('status', 'unknown')
            )
            for item in platform.get('third_party_libs', [])
        ]
        platform_info = PlatformInfo(
            name=platform.get('platform_name', 'Unknown Platform'),
            debug=platform.get('debug', '-'),
            os=platform.get('operating_system', '-'),
            tester=platform.get('tester_name', '-'),
            compiler=platform.get('compiler', '-'),
            tpl_info=tpl_list
        )
        platforms_info.append(platform_info)
    return platforms_info


def get_docker_images() -> Dict[str, List[str]]:
    """
    Get Docker image information by calling `list_test_runner_machines`.
    Returns a dictionary with machine names as keys and lists of images as values.
    """
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        result = subprocess.run(
            [os.path.join(script_dir, 'list_test_runner_machines'), '--plain'],
            capture_output=True,
            text=True,
            check=True
        )
        output = result.stdout.strip()

        machines_info = {}
        current_machine = None
        parsing_images = False

        for line in output.splitlines():
            if line.startswith("## "):
                current_machine = line.strip("# ").strip()
                machines_info[current_machine] = []
                parsing_images = False

            elif line.startswith("Tested images:"):
                parsing_images = True

            elif parsing_images and (line.startswith("cgal/testsuite-docker:") or line.startswith("docker.io/cgal/testsuite-docker:")):
                machines_info[current_machine].append(line.strip())

        return machines_info

    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"Error running `list_test_runner_machines`: {e}") from e
    except Exception as e:
        raise RuntimeError(f"Error parsing Docker information: {e}") from e


def add_docker_summary(report: List[str], machines_info: Dict[str, List[str]]):
    """Add a summary of Docker images used on each machine to the report."""
    report.append("\n## Docker Test Summary\n")
    for machine, images in machines_info.items():
        report.append(f"\n### Machine: {machine} ({len(images)} images)")
        report.append("\n#### Tested Images:")
        for image in images:
            report.append(f"- {image}")
        report.append("")


def generate_markdown_report(platforms_info: List[PlatformInfo], version: str) -> str:
    """Generate a markdown report from the platforms information."""
    machines_info = get_docker_images()
    report = []
    report.append("# TestSuite Report")
    report.append(f"\nGenerated on: {
                  datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    url = TESTSUITE_URL_TEMPLATE.format(version=version)
    report.append(f"\nCGAL Version: [{version}]({url})\n")
    add_docker_summary(report, machines_info)
    report.append("\n## Platforms Summary\n")
    report.append("| Platform | Debug | OS | Tester | Compiler |")
    report.append("|----------|-------|----|--------|----------|")
    for platform in platforms_info:
        report.append(
            f"| {platform.name} | {platform.debug} | {platform.os} | "
            f"{platform.tester} | {platform.compiler} |"
        )
    report.append("\n## Detailed Third-party Libraries")
    for platform in platforms_info:
        report.append(f"\n### Platform: {platform.name}\n")
        tpl_list = sorted(platform.tpl_info, key=lambda x: x.name)
        report.append("| Library Name | Version | Status |")
        report.append("|--------------|---------|--------|")
        for tpl in tpl_list:
            version_str = str(tpl.version) if tpl.version else "N/A"
            status_str = "❌" if tpl.version == "not found" else "✅"
            report.append(f"| {tpl.name} | {version_str} | {status_str} |")
        found_tpls = sum(1 for tpl in tpl_list if tpl.version != "not found")
        total_tpls = len(tpl_list)
        report.append(
            f"\n**Summary**: found {found_tpls} third-party libraries out of {total_tpls}")
    return "\n".join(report)


def main():
    """Main function to generate the testsuite report."""
    try:
        version = get_latest_version()
        json_data = fetch_json_data(version)
        platforms_info = analyze_tpl_data(json_data)
        markdown_report = generate_markdown_report(platforms_info, version)
        print(markdown_report)
    except requests.RequestException as e:
        print(f"**Error fetching data:**\n\n```\n{str(e)}\n```\n")
        raise
    except json.JSONDecodeError as e:
        print(f"**Error: Invalid JSON data**\n\n```\n{str(e)}\n```")
        print(f"\nFile:\n\n```json\n{e.doc}\n```")
        raise
    except Exception as e:
        print(f"**Error processing data:**\n\n```\n{str(e)}\n```\n")
        raise


if __name__ == "__main__":
    main()
