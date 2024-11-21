import json
import sys
from typing import Dict, List
from dataclasses import dataclass
from datetime import datetime
import re
import requests

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

def get_latest_version() -> str:
    url = "https://cgal.geometryfactory.com/CGAL/Releases/LATEST"
    response = requests.get(url)
    response.raise_for_status()
    version_text = response.text.strip()
    match = re.match(r'CGAL-([^.]+\.[^-]+-[^-]+-\d+)', version_text)
    if not match:
        raise ValueError(f"Unexpected version format: {version_text}")
    return match.group(1)

def fetch_json_data(version: str) -> Dict:
    url = f"https://cgal.geometryfactory.com/CGAL/testsuite/CGAL-{version}/search_index.json"
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

def analyze_tpl_data(json_data: Dict) -> List[PlatformInfo]:
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

def generate_markdown_report(platforms_info: List[PlatformInfo], version: str) -> str:
    report = []
    report.append("# TestSuite Report")
    report.append(f"Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    report.append(f"CGAL Version: {version}\n")
    report.append("## Platforms Summary")
    report.append("\n| Platform | Debug | OS | Tester | Compiler |")
    report.append("|----------|--------|----|---------| ---------|")
    for platform in platforms_info:
        report.append(
            f"| [{platform.name}](#platform-{platform.name.lower().replace(' ', '-')}) | {platform.debug} | {platform.os} | "
            f"{platform.tester} | {platform.compiler} |"
        )
    report.append("\n")
    report.append("## Detailed TPL\n")
    for platform in platforms_info:
        report.append(f"### Platform: {platform.name} <a name='{platform.name.lower().replace(' ', '-')}'></a>")
        tpl_list = sorted(platform.tpl_info, key=lambda x: x.name)
        report.append("\n| TPL Name | Version | Status |")
        report.append("|----------|----------|---------|")
        for tpl in tpl_list:
            version_str = str(tpl.version) if tpl.version else "N/A"
            status_str = "✅" if tpl.status == "found" else "❌"
            report.append(f"| {tpl.name} | {version_str} | {status_str} |")
        report.append("\n")
        found_tpls = sum(1 for tpl in tpl_list if tpl.status == "found")
        total_tpls = len(tpl_list)
        report.append(f"**Summary**: {found_tpls}/{total_tpls} TPLs found")
    return "\n".join(report)

def main():
    try:
        version = get_latest_version()
        json_data = fetch_json_data(version)
        platforms_info = analyze_tpl_data(json_data)
        markdown_report = generate_markdown_report(platforms_info, version)
        print(markdown_report)
    except requests.RequestException as e:
        print(f"Error fetching data: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError:
        print("Error: Invalid JSON data", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing data: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
