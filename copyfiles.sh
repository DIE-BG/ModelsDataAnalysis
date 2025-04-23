#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <input_file> <destination_folder>"
  exit 1
fi

# Assign the input file and destination folder to variables
input_file="$1"
destination_folder="$2"

# Check if the input file exists and is readable
if [ ! -r "$input_file" ]; then
  echo "Error: Input file '$input_file' does not exist or is not readable."
  exit 1
fi

# Check if the destination folder exists
if [ ! -d "$destination_folder" ]; then
  echo "Error: Destination folder '$destination_folder' does not exist."
  exit 1
fi

# Read each line from the input file and remove carriage returns
while IFS= read -r file_path; do
  file_path="${file_path//$'\r'/}" # Remove carriage return

  # Check if the file exists before attempting to copy
  if [ -f "$file_path" ]; then
    echo "Copying '$file_path' to '$destination_folder'..."
    cp "$file_path" "$destination_folder"
    if [ "$?" -ne 0 ]; then
      echo "Error copying '$file_path'."
    fi
  else
    echo "Warning: File '$file_path' not found."
  fi
done < "$input_file"

echo "File copying process completed."