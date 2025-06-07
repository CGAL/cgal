import hashlib, glob, json, sys, os

def hash_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()

def hash_sources(file_patterns):
    h = hashlib.sha256()
    files = []
    for pattern in file_patterns:
        for filename in sorted(glob.glob(pattern, recursive=True)):
            files.append(filename)
            with open(filename, 'rb') as f:
                h.update(f.read())
    return h.hexdigest(), files

if __name__ == "__main__":
    exe_path = sys.argv[1]
    source_patterns = sys.argv[2:]  # e.g., "*.cpp" "*.h"
    exe_hash = hash_file(exe_path)
    source_hash, source_files = hash_sources(source_patterns)
    metadata = {
        "exe_hash": exe_hash,
        "source_hash": source_hash,
        "source_files": source_files,
    }
    with open(exe_path + ".metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    print(f"Wrote metadata to {exe_path}.metadata.json")
