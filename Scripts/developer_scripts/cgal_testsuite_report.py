import json
from typing import Dict, List
from dataclasses import dataclass
from datetime import datetime
import re
import requests

LATEST_VERSION_URL = "https://cgal.geometryfactory.com/CGAL/Releases/LATEST"
JSON_DATA_URL_TEMPLATE = "https://cgal.geometryfactory.com/CGAL/testsuite/CGAL-{version}/search_index.json"
TESTSUITE_URL_TEMPLATE = "https://cgal.geometryfactory.com/CGAL/testsuite/results-{version}.shtml"
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
            for item in platform.get('tpl', [])
        ]
        platform_info = PlatformInfo(
            name=platform.get('name', 'Unknown Platform'),
            debug=platform.get('debug', '-'),
            os=platform.get('os', '-'),
            tester=platform.get('tester', '-'),
            compiler=platform.get('compiler', '-'),
            tpl_info=tpl_list
        )
        platforms_info.append(platform_info)
    return platforms_info

def fragment_name(platform: PlatformInfo) -> str:
    """Return a fragment name from a given platform."""
    return f"platform-{platform.name.lower().replace(' ', '-')}"

def generate_markdown_report(platforms_info: List[PlatformInfo], version: str) -> str:
    """Generate a markdown report from the platforms information."""
    report = []
    report.append("# TestSuite Report")
    report.append(f"\nGenerated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    url = TESTSUITE_URL_TEMPLATE.format(version=version)
    report.append(f"\nCGAL Version: [{version}]({url})\n")
    report.append("## Platforms Summary\n")
    report.append("| Platform | Debug | OS | Tester | Compiler |")
    report.append("|----------|-------|----|--------|----------|")
    for platform in platforms_info:
        report.append(
            f"| [{platform.name}](#{fragment_name(platform)}) | {platform.debug} | {platform.os} | "
            f"{platform.tester} | {platform.compiler} |"
        )
    report.append("\n## Detailed Third-party Libraries\n")
    for platform in platforms_info:
        report.append(f"### Platform: {platform.name}")
        tpl_list = sorted(platform.tpl_info, key=lambda x: x.name)
        report.append("\n| Library Name | Version | Status |")
        report.append("|----------|----------|---------|")
        for tpl in tpl_list:
            version_str = str(tpl.version) if tpl.version else "N/A"
            status_str = "✅" if tpl.status == "found" else "❌"
            report.append(f"| {tpl.name} | {version_str} | {status_str} |")
        found_tpls = sum(1 for tpl in tpl_list if tpl.status == "found")
        total_tpls = len(tpl_list)
        report.append(f"**Summary**: found {found_tpls} third-party libraries out of {total_tpls}")
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
