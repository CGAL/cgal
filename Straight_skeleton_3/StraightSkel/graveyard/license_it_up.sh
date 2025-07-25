#!/bin/bash

# Set the path to the license header file
LICENSE_FILE="license.txt"

# Function to process a single file
process_file()
{
  local file="$1"
  local temp_file="$file.tempitempo"

  echo $file

  # Add the license header to the temporary file
  cat "$LICENSE_FILE" "$file" > "$temp_file"

  # Replace the original file with the temporary file
  mv "$temp_file" "$file"
}

# Function to process files in a directory recursively
process_directory()
{
  local directory="$1"

  # Process .h files in the current directory
  for file in "$directory"/*.h; do
    if [ -f "$file" ]; then
      process_file "$file"
    fi
  done

  # Process .cpp files in the current directory
  for file in "$directory"/*.cpp; do
    if [ -f "$file" ]; then
      process_file "$file"
    fi
  done

  # Recursively process subdirectories
  for subdir in "$directory"/*/; do
    if [ -d "$subdir" ]; then
      process_directory "$subdir"
    fi
  done
}

# Start processing from the current directory
process_directory "."

echo "License headers added successfully."
