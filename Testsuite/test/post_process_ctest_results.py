import sys
import re
import os
import collections

CONFIG_REGEX = re.compile(
    r"(.*Configuring (examples|demo|test)*( in )*(test/|examples/|demo/)*)((?!done)\w+)"
)
DEMO_REGEX = re.compile(r".*in demo/")
EXAMPLES_REGEX = re.compile(r'.*in examples/')
SEPARATOR = "------------------------------------------------------------------"


# open the Installation report
# For each NAME, check if NAME is a directory. If not, create one, create a
# text report, and write everything that is in the report until the next NAME
# in it. Then, add 'NAME r' in the global report. This should allow to get all
# the NOTICE and other info explaining why the configuration is skipped.

def find_third_separator(contents):
    separator_count = 0
    for j, line in enumerate(contents):
        if line.strip() == SEPARATOR:
            separator_count += 1
            if separator_count == 3:
                return j
    return len(contents) + 2

def last(iterator):
    return collections.deque(iterator, maxlen=1).pop()

def find_last_separator(contents):
    position, _ = last(filter(lambda x: x[1].strip() == SEPARATOR, enumerate(contents)))
    return position

def process_report(input_report_file_name, report_file_name, global_report_file_name):
    name = ""
    is_writing = False
    is_ignored = False
    position = 0
    lines_to_write = []
    installation_cmake_logs = []

    file_path = f"Installation/{report_file_name}"
    print(f"Debug: Opening file {file_path}")
    with open(file_path, "r", encoding="utf-8") as file:
        contents = file.readlines()
    position = find_last_separator(contents)
    print(f"Debug: Position is {position}")
    print(f"Debug: Length of contents is {len(contents)}")
    for i, line in enumerate(contents):
        if i > position:
            if line.strip() == "== Generating build files for tests ==":
                print(f"Debug: Found the line {line} at position {i}")
                break
            installation_cmake_logs.append(line)
    contents = []
    print(f"Debug: Length of installation CMake logs is {len(installation_cmake_logs)}")
    print(f"Debug: Installation CMake logs are {"".join(installation_cmake_logs)}")
    global_report = open(global_report_file_name, "a+", encoding="utf-8")
    with open(input_report_file_name, "rt", encoding="utf-8") as input_report_file:
        for line in input_report_file:
            match = CONFIG_REGEX.match(line)
            if is_writing:
                if match:
                    is_writing = False
                    if is_ignored:
                        print(f"{name} r", file=global_report)
                        is_ignored = False
                    if lines_to_write:
                        file_path = f"{name}/{report_file_name}"
                        if os.path.exists(file_path):
                            with open(file_path, "r", encoding="utf-8") as file:
                                contents = file.readlines()
                        else:
                            contents = []

                        position = find_third_separator(contents)

                        if not any(
                            re.search("- CMake Results .*", content)
                            for content in contents
                        ):
                            lines_to_write.insert(
                                0,
                                f"{SEPARATOR}\n- CMake Results for {name}\n{SEPARATOR}\n\n",
                            )
                        lines_to_write.insert(0, "\n")
                        contents[position:position] = lines_to_write

                        with open(file_path, "w", encoding="utf-8") as file:
                            file.write("".join(contents))

                        lines_to_write = []
                    if is_ignored:
                        is_ignored = False
                else:
                    if line.strip() != "":
                        lines_to_write.append(line)
            if not is_writing:
                if match:
                    name = match.group(0).replace(match.group(1), "")
                    print(f"Debug: Found name {name}")
                    if DEMO_REGEX.match(line):
                        name = f"{name}_Demo"
                    elif EXAMPLES_REGEX.match(line):
                        name = f"{name}_Examples"
                    elif name == "libCGAL":
                        name = "libCGAL_shared"
                    elif name == "libCGAL_Core":
                        name = "libCGALCore_shared"
                    elif name == "libCGAL_ImageIO":
                        name = "libCGALimageIO_shared"
                    elif name == "libCGAL_Qt6":
                        name = "libCGALQt6_shared"
                    if name == "incomplete":
                        is_writing = False
                        is_ignored = False
                        continue
                    else:
                        if not os.path.isdir(name):
                            is_ignored = True
                            os.mkdir(name)
                            with open(
                                f"{os.getcwd()}/../../../../../.scm-branch",
                                "r",
                                encoding="utf-8",
                            ) as scm_branch_file:
                                scm_branch_content = scm_branch_file.read()

                            with open(
                                f"{name}/{report_file_name}", "w+", encoding="utf-8"
                            ) as report_file_handle:
                                report_file_handle.write(scm_branch_content)
                        else:
                            is_ignored = False
                            file_path = f"{name}/{report_file_name}"
                            if os.path.exists(file_path):
                                with open(file_path, "r", encoding="utf-8") as file:
                                    contents = file.readlines()
                            else:
                                contents = []

                            position = find_third_separator(contents)

                            if not any(
                                re.search("- CMake Logs .*", content)
                                for content in contents
                            ):
                                contents.insert(
                                    position - 1,
                                    SEPARATOR
                                    + "\n- CMake Logs from Installation \n"
                                    + SEPARATOR
                                    + "\n\n",
                                )
                            for log in installation_cmake_logs:
                                contents.insert(position, log)
                                position += 1

                            with open(file_path, "w", encoding="utf-8") as file:
                                file.write("".join(contents))

                        is_writing = True

    if is_writing:
        is_writing = False
        if is_ignored:
            print(f"{name} r", file=global_report)
            is_ignored = False
    global_report.close()

def main():
    input_report_file_name = sys.argv[1]
    report_file_name = sys.argv[2]
    global_report_file_name = sys.argv[3]
    process_report(input_report_file_name, report_file_name, global_report_file_name)

if __name__ == "__main__":
    main()
