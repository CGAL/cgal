import hashlib, glob, json, sys, os, traceback

def hash_file(path):
    h = hashlib.sha256()
    try:
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                h.update(chunk)
        return h.hexdigest()
    except FileNotFoundError:
        print(f"[metadata] Warning: file not found: {path}", file=sys.stderr)
        return ""
    except PermissionError as e:
        print(f"[metadata] Warning: permission denied reading {path}: {e}", file=sys.stderr)
        return ""
    except Exception as e:
        print(f"[metadata] Warning: failed to hash file {path}: {e}", file=sys.stderr)
        return ""

def hash_sources(file_patterns):
    h = hashlib.sha256()
    files = []
    for pattern in file_patterns:
        matched = sorted(glob.glob(pattern, recursive=True))
        if not matched:
            # Allow plain file path without glob metacharacters as well
            if os.path.isfile(pattern):
                matched = [pattern]
            else:
                print(f"[metadata] Info: no files matched pattern '{pattern}'", file=sys.stderr)
        for filename in matched:
            files.append(filename)
            try:
                with open(filename, 'rb') as f:
                    # Read in chunks to avoid loading very large files in memory
                    for chunk in iter(lambda: f.read(4096), b""):
                        h.update(chunk)
            except Exception as e:
                print(f"[metadata] Warning: failed to read source file {filename}: {e}", file=sys.stderr)
                continue
    # If no files were hashed, return empty hash string for clarity
    return (h.hexdigest() if files else ""), files

if __name__ == "__main__":
    try:
        if len(sys.argv) < 2:
            print("[metadata] Error: missing arguments. Usage: python generate_source_exec_metadata.py <exe_path> [<source_glob> ...]", file=sys.stderr)
            sys.exit(0)

        exe_path = sys.argv[1]
        source_patterns = sys.argv[2:]  # e.g., "*.cpp" "*.h"

        if not os.path.isfile(exe_path):
            print(f"[metadata] Warning: target binary does not exist: {exe_path}", file=sys.stderr)

        exe_hash = hash_file(exe_path)
        source_hash, source_files = hash_sources(source_patterns)

        metadata = {
            "exe_hash": exe_hash,
            "source_hash": source_hash,
            "source_files": source_files,
        }

        out_path = exe_path + ".metadata.json"
        try:
            with open(out_path, "w") as f:
                json.dump(metadata, f, indent=2)
            print(f"[metadata] Wrote metadata to {out_path}")
        except Exception as e:
            print(f"[metadata] Warning: failed to write metadata file {out_path}: {e}", file=sys.stderr)
            # Do not fail the build if we cannot write the metadata
        # Always succeed to avoid breaking the build on non-critical metadata issues
        sys.exit(0)
    except SystemExit:
        # Allow explicit sys.exit(0) to pass through
        raise
    except Exception:
        # Last-resort catch: never fail the build because of metadata generation
        print("[metadata] Unexpected error during metadata generation:", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(0)
