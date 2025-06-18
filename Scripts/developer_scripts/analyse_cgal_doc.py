import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

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
    except Exception as e:
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

def save_broken_links_report(filename="broken_links_report.txt"):
    """
    Save the broken links report to a text file.
    """
    with open(filename, "w") as f:
        f.write("--- Final Broken Links Report ---\n")
        for url, error, referrer in broken_links:
            f.write(f"{url} (found on: {referrer}) => {error}\n")

if __name__ == "__main__":
    # List of documentation versions to check
    all_versions = [
        "https://doc.cgal.org/latest/Manual/index.html", # LATEST
        "https://cgal.geometryfactory.com/CGAL/doc/master/Manual/index.html", # master
    ]

    print(f"Checking {len(all_versions)} CGAL documentation versions:")
    for version in all_versions:
        print(f" - {version}")

    for version_url in all_versions:
        print(f"\n Starting crawl for version: {version_url}")
        threaded_crawl(version_url, max_workers=40)

    print("\n--- Final Broken Links Report ---")
    for url, error, referrer in broken_links:
        print(f"{url} (found on: {referrer}) => {error}")

    # Save the report to a file
    save_broken_links_report()
