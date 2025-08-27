import sys
import re
import os
import collections
import logging
from itertools import islice

CONFIG_REGEX = re.compile(
    r"(.*Configuring (examples|demo|test) *in *(test/|examples/|demo/))((?!done)\w+)"
)
DEMO_REGEX = re.compile(r".*in demo/")
EXAMPLES_REGEX = re.compile(r'.*in examples/')
SEPARATOR = "------------------------------------------------------------------"

def find_third_separator(contents):
    """Find the position of the third separator line in the contents.
    If there are less than 3 separators, then return the position where the third separator
    should be inserted.
    """
    separator_positions = (
        i for i, line in enumerate(contents) if line.strip() == SEPARATOR
    )
    return next(islice(separator_positions, 2, None), len(contents) + 2)

def last(iterator):
    """Return the last item of an iterator or None if empty."""
    return collections.deque(iterator, maxlen=1).pop()

def find_last(contents, query_string):
    """Find the number of the last line matching the query string."""
    position, _ = last(filter(lambda x: x[1].strip() == query_string, enumerate(contents)))
    return position

def read_file_lines(file_path):
    """Read the lines of a file and return them as a list."""
    try:
        with open(file_path, "r", encoding="utf-8") as file:
            return file.readlines()
    except IOError as e:
        print(f"Error opening file {file_path}: {e}")
        return []

def write_file_lines(file_path, contents):
    """Write the contents to a file. The contents should be a list of strings."""
    try:
        with open(file_path, "w", encoding="utf-8") as file:
            file.write("".join(contents))
    except IOError as e:
        print(f"Error writing to file {file_path}: {e}")

def mark_package_as_missing_requirements(global_report_file_name, name):
    """Mark a package as missing requirements in the global report file."""
    try:
        with open(global_report_file_name, "a+", encoding="utf-8") as global_report:
            print(f"{name} r", file=global_report)
    except IOError as e:
        print(f"Error opening global report file {global_report_file_name}: {e}")

def handle_end_of_package(package_name, report_file_name, lines_to_write):
    """Handle the end of a package by inserting the lines to write into the report file."""
    if not lines_to_write:
        return

    file_path = f"{package_name}/{report_file_name}"
    contents = read_file_lines(file_path)
    position = find_third_separator(contents)

    if not any(re.search("- CMake Results .*", content) for content in contents):
        lines_to_write.insert(0, f"""
{SEPARATOR}
- CMake Results for {package_name}
{SEPARATOR}

""")
    contents[position:position] = lines_to_write

    write_file_lines(file_path, contents)


SCM_BRANCH_FILE_CONTENT = read_file_lines(f"{os.getcwd()}/../../../../../.scm-branch")

def handle_new_package__is_ignored(name, report_file_name, cmake_logs):
    """Handle new package creation or update logs if package already exists."""
    if not os.path.isdir(name):
        os.mkdir(name)
        write_file_lines(f"{name}/{report_file_name}", SCM_BRANCH_FILE_CONTENT)
        return True
    else:
        file_path = f"{name}/{report_file_name}"
        contents = read_file_lines(file_path)
        position = find_third_separator(contents)

        if not any(re.search("- CMake Logs .*", content) for content in contents):
            contents.insert(
                position - 1,
                SEPARATOR + "\n- CMake Logs from Installation \n" + SEPARATOR + "\n\n",
            )
        for log in cmake_logs:
            contents.insert(position, log)
            position += 1

        write_file_lines(file_path, contents)
        return False

def retrieve_cmake_logs(file_path):
    """Retrieve the CMake logs from a file and return them as a list."""
    logging.debug("Opening file %s", file_path)
    contents = read_file_lines(file_path)

    position_begin = 1 + find_last(contents, SEPARATOR)
    position_end = 1 + find_last(contents, "== Generating build files for tests ==")

    cmake_logs = contents[position_begin:position_end]

    logging.debug("CMake log beginning is at line %d", position_begin)
    logging.debug("CMake log       end is at line %d", position_end)
    logging.debug("Length of contents is %d", len(contents))
    logging.debug("Length of installation CMake logs is %d", len(cmake_logs))
    logging.debug("Installation CMake logs are %s", "".join(cmake_logs))
    return cmake_logs

def main():
    """Main function that processes the input report file and performs necessary operations."""
    input_report_file_name = sys.argv[1]
    report_file_name = sys.argv[2]
    global_report_file_name = sys.argv[3]

    cmake_logs = retrieve_cmake_logs(f"Installation/{report_file_name}")

    package_name = ""
    lines_to_write = []

    for line in read_file_lines(input_report_file_name):

        line_matches_new_package = CONFIG_REGEX.match(line)
        if package_name and line_matches_new_package:
            handle_end_of_package(package_name, report_file_name, lines_to_write)
            lines_to_write = []
            package_name = ""

        if line_matches_new_package:
            logging.debug("Found new package %s", line_matches_new_package.group(0))
            logging.debug("          group 1 %s", line_matches_new_package.group(1))
            new_package_name = line_matches_new_package.group(0).replace(
                line_matches_new_package.group(1), ""
            )
            logging.debug("Setting package name to %s", new_package_name)
            package_name = new_package_name
            if DEMO_REGEX.match(line):
                package_name = f"{package_name}_Demo"
            elif EXAMPLES_REGEX.match(line):
                package_name = f"{package_name}_Examples"

            if package_name == "incomplete":
                package_name = ""
                continue
            else:
                is_ignored = handle_new_package__is_ignored(
                    package_name, report_file_name, cmake_logs
                )
                logging.debug("Is package %s ignored? %s", package_name, is_ignored)
                if is_ignored:
                    mark_package_as_missing_requirements(global_report_file_name, package_name)

        if package_name and not line_matches_new_package and line.strip() != "":
            lines_to_write.append(line)

if __name__ == "__main__":
    main()
