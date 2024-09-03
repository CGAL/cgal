import sys
import re
import os

input_report_file_name=sys.argv[1]
report_file_name=sys.argv[2]
global_report_file_name=sys.argv[3]
config_regex=re.compile(r'(.*Configuring (examples|demo|test)*( in )*(test/|examples/|demo/)*)((?!done)\w+)')
demo_regex=re.compile(r'.*in demo/')
examples_regex=re.compile(r'.*in examples/')
Separator = "------------------------------------------------------------------"


#open the Installation report
#For each NAME, check if NAME is a directory. If not, create one, create a
#text report, and write everything that is in the report until the next NAME
#in it. Then, add 'NAME r' in the global report. This should allow to get all
#the NOTICE and other info explaining why the configuration is skipped.

name=""
is_writing=False
is_ignored=False
position = 0
lines_to_write = []
installation_cmake_logs = []



def find_third_separator(inner_contents):
    separator_count = 0
    for j, inner_line in enumerate(inner_contents):
        if inner_line.strip() == Separator:
            separator_count += 1
            if separator_count == 3:
                return j
    return len(inner_contents) + 2

def find_last_separator(inner_contents):
    for j, inner_line in enumerate(inner_contents):
        if inner_line.strip() == Separator:
            inner_position = j
    return inner_position

with open ("{dir}/{file}".format(dir="Installation",file=report_file_name), "r", encoding="utf-8") as file:
    contents = file.readlines()
position = find_last_separator(contents)
for i, line in enumerate(contents):
    if i > position:
        installation_cmake_logs.append(line)
    if line.strip() == "== Generating build files for tests ==":
        break
contents = []

global_report = open(global_report_file_name, "a+", encoding="utf-8")
with open(input_report_file_name, "rt", encoding="utf-8") as input_report_file:
    for myline in input_report_file:
        match = config_regex.match(myline)
        if is_writing:
            if match:
                is_writing = False
                if is_ignored:
                    print("{label} {result}".format(label=name, result='r'), file=global_report)
                    is_ignored = False
                if lines_to_write:
                    file_path = f"{name}/{report_file_name}"
                    if os.path.exists(file_path):
                        with open(file_path, "r", encoding="utf-8") as file:
                            contents = file.readlines()
                    else:
                        contents = []

                    position = find_third_separator(contents)

                    if not any(re.search("- CMake Results .*", content) for content in contents):
                        lines_to_write.insert(0, f"{Separator}\n- CMake Results for {name}\n{Separator}\n\n")
                    lines_to_write.insert(0, "\n")
                    contents[position:position] = lines_to_write

                    with open(file_path, "w", encoding="utf-8") as file:
                        file.write("".join(contents))

                    lines_to_write = []
                if is_ignored:
                    is_ignored = False
            else:
                if myline.strip() != "":
                    lines_to_write.append(myline)
        if not is_writing:
            if match:
                name=match.group(0).replace(match.group(1), "")
                if demo_regex.match(myline):
                    name="{str}_Demo".format(str=name)
                elif examples_regex.match(myline):
                    name="{str}_Examples".format(str=name)
                elif name == "libCGAL":
                    name="libCGAL_shared"
                elif name == "libCGAL_Core":
                    name="libCGALCore_shared"
                elif name == "libCGAL_ImageIO":
                    name="libCGALimageIO_shared"
                elif name == "libCGAL_Qt6":
                    name="libCGALQt6_shared"
                if name=="incomplete":
                    is_writing=False
                    is_ignored=False
                    continue
                else:
                    if not os.path.isdir(name):
                        is_ignored = True
                        os.mkdir(name)
                        with open("{}/../../../../../.scm-branch".format(os.getcwd()), 'r', encoding="utf-8") as scm_branch_file:
                            scm_branch_content = scm_branch_file.read()

                        with open("{dir}/{file}".format(dir=name, file=report_file_name), "w+", encoding="utf-8") as report_file_handle:
                            report_file_handle.write(scm_branch_content)
                    else:
                        is_ignored = False
                        file_path = "{dir}/{file}".format(dir=name, file=report_file_name)
                        if os.path.exists(file_path):
                            with open(file_path, "r", encoding="utf-8") as file:
                                contents = file.readlines()
                        else:
                            contents = []

                        position = find_third_separator(contents)

                        if not any(re.search("- CMake Logs .*", content) for content in contents):
                            contents.insert(position - 1, Separator + "\n- CMake Logs from Installation \n" + Separator + "\n\n")
                        for log in installation_cmake_logs:
                            contents.insert(position, log)
                            position += 1

                        with open(file_path, "w", encoding="utf-8") as file:
                            file.write("".join(contents))

                    is_writing = True

if is_writing:
    is_writing=False
    if is_ignored:
        print("{label} {result}".format(label=name, result='r'), file=global_report)
        is_ignored=False
global_report.close()
