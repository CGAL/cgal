import sys
from urllib.parse import urljoin, urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from bs4 import BeautifulSoup
import requests
import os

# Set to keep track of already visited URLs
visited = set()
# Lock to synchronize access to 'visited' set between threads
visited_lock = threading.Lock()
# List to store broken links
broken_links = []

def is_internal(url, domain):
    """
    Check if a URL is internal (belongs to the same domain).
    """
    parsed = urlparse(url)
    return parsed.netloc == "" or parsed.netloc == domain

def check_and_crawl(url, referrer, domain):
    """
    Check if the URL has already been visited, then try to fetch it.
    If the link is broken, add it to 'broken_links'.
    Return new internal links found on the page.
    """
    with visited_lock:
        if url in visited:
            return []
        visited.add(url)

    try:
        response = requests.get(url, timeout=10)
        if response.status_code >= 400:
            print(f"[BROKEN] {url} (found on: {referrer}) - Status: {response.status_code}")
            broken_links.append((url, response.status_code, referrer))
            return []

        print(f"[OK] {url}")
        soup = BeautifulSoup(response.text, "html.parser")
        new_links = []
        # Find all <a href="..."> tags
        for tag in soup.find_all("a", href=True):
            href = tag["href"]
            full_url = urljoin(url, href)
            # Ignore anchors and external links
            if "#" in href or not is_internal(full_url, domain):
                continue
            new_links.append((full_url, url, domain))
        return new_links
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] {url} (found on: {referrer}) - {e}")
        broken_links.append((url, str(e), referrer))
        return []

def threaded_crawl(start_url, max_workers=20):
    """
    Start crawling from the start URL using multiple threads.
    """
    domain = urlparse(start_url).netloc
    to_visit = [(start_url, "ROOT", domain)]
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        while to_visit:
            # Launch checks in parallel
            futures = [executor.submit(check_and_crawl, url, ref, dom) for url, ref, dom in to_visit]
            to_visit = []
            for future in as_completed(futures):
                new_links = future.result()
                to_visit.extend(new_links)

def save_broken_links_report(filename="broken_links_report.txt", links=None):
    """
    Save the broken links report to a text file.
    """
    with open(filename, "w", encoding="utf-8") as f:
        f.write("--- Final Broken Links Report ---\n")
        for url, error, referrer in links:
            f.write(f"{url} (found on: {referrer}) => {error}\n")

def filter_broken_links_by_packages(allowed_packages_file):
    """
    Filter broken_links to keep only those in 'Manual' or in allowed packages.
    """
    with open(allowed_packages_file, "r", encoding="utf-8") as f:
        allowed = set(line.strip() for line in f if line.strip())

    filtered = []
    for url, error, referrer in broken_links:
        path_parts = urlparse(url).path.strip("/").split("/")
        if len(path_parts) < 3:
            continue
        # path_parts = ['8186', 'v12', 'Manual', ...] or ['8186', 'v12', 'Alpha_shapes_3', ...]
        section = path_parts[2]
        if section == "Manual" or section in allowed:
            filtered.append((url, error, referrer))
    return filtered

if __name__ == "__main__":
    args = sys.argv[1:]
    versions_to_check = []
    allowed_packages_file = None
    if args:
        if os.path.isfile(args[0]) or args[0].endswith(".txt"):
            print("ERROR: First argument should be a URL, not a package file.")
            sys.exit(1)
        versions_to_check = [args[0]]
        print(f"Checking custom documentation URL: {args[0]}")
        if len(args) > 1:
            allowed_packages_file = args[1]
            print(f"Using allowed package list from: {allowed_packages_file}")
    else:
        # Default list of documentation versions to check
        versions_to_check = [
            "https://doc.cgal.org/latest/Manual/index.html", # LATEST
            "https://cgal.geometryfactory.com/CGAL/doc/master/Manual/index.html", # master
        ]
        print(f"Checking {len(versions_to_check)} CGAL documentation versions:")
        for version in versions_to_check:
            print(f" - {version}")

    for version_url in versions_to_check:
        print(f"\nStarting crawl for version: {version_url}")
        threaded_crawl(version_url, max_workers=40)

    print("\n--- Final Broken Links Report ---")

    if allowed_packages_file:
        final_broken_links = filter_broken_links_by_packages(allowed_packages_file)
    else:
        final_broken_links = broken_links
    for url, error, referrer in broken_links:
        print(f"{url} (found on: {referrer}) => {error}")

    # Save the report to a file
    save_broken_links_report(links=final_broken_links)
