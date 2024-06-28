import argparse
import os
import re

def sort_fasta_by_score(fasta_file, output_file, num_records=10):
  """Sorts a FASTA file by score in ascending order and outputs the lowest num_records scores."""
  # Open the FASTA file for reading
  with open(fasta_file, "r") as handle:
    # Initialize variables to store entries
    entries = []
    current_entry = None
    for line in handle:
      line = line.rstrip()
      # Check for header line
      if line.startswith(">"):
        # Create a new entry if needed
        if current_entry:
          entries.append(current_entry)
        current_entry = {"header": line[1:], "sequence": ""}
      else:
        # Append sequence lines to the current entry
        current_entry["sequence"] += line

  # Define a function to extract score from the sequence ID
  def get_score(header):
    score_str = re.search(r" score=(\d+\.\d+)", header).group(1)
    return float(score_str)

  # Check if any entries were parsed
  if not entries:
    print(f"No FASTA entries found in {fasta_file}!")
    return

  # Sort entries by score in ascending order
  entries.sort(key=lambda entry: get_score(entry["header"]))

  # Open the output file for writing
  with open(output_file, "w") as handle:
    # Write the first num_records entries (lowest scores) to the output file
    #print(len(entries))
    top_seq = int(num_records)*len(entries)/100
    for entry in entries[:int(top_seq)]:
      handle.write(f">{entry['header']}\n{entry['sequence']}\n")



def process_fasta_file(filepath):
  """Reads a FASTA file and calculates the consensus sequence."""
  sequences = []
  with open(filepath, "r") as handle:
    current_seq = ""
    for line in handle:
      if line.startswith(">"):  # Header line
        if current_seq:
          sequences.append(current_seq)  # Save previous sequence
        current_seq = ""  # Start a new sequence
      else:
        current_seq += line.strip()  # Append sequence line
    if current_seq:  # Add the last sequence
      sequences.append(current_seq)

  # Calculate consensus sequence
  consensus = ""
  for i in range(min(len(seq) for seq in sequences)):  # Iterate to the shortest length
    column = [seq[i] for seq in sequences]
    most_common = max(set(column), key=column.count)  # Find most common base
    consensus += most_common

  return consensus


def main():
  """Processes all FASTA files and combines consensus sequences into a single file."""
  # Parse arguments
  parser = argparse.ArgumentParser(description="Calculate consensus sequences from FASTA files")
  parser.add_argument("input_folder", help="Path to the folder containing FASTA files")
  parser.add_argument("output_file", help="Path to the output file for combined consensus sequences")
  parser.add_argument("--binder_sequence", help="Sequence of binder", required=False)
  parser.add_argument("--top", help="Best x %")
  args = parser.parse_args()

  consensus_sequences = []
  for filename in os.listdir(args.input_folder):
    if filename.endswith(".fa"):
      
      fasta_file = os.path.join(args.input_folder, filename)
      output_file = os.path.join(args.input_folder, f"sorted_{filename}")
      sort_fasta_by_score(fasta_file, output_file, args.top)
      
      consensus_sequence = process_fasta_file(output_file)
      if args.binder_sequence:
        binder_string=args.binder_sequence
        binder_string=binder_string.replace("Z","\n")
        #consensus_sequences.append(f">A\n{binder_string}:{consensus_sequence}\n")  # Format with filename
        consensus_sequences.append(f">A\n{consensus_sequence}:{binder_string}\n")
      else:
        binder_string=""
        consensus_sequences.append(f">A\n{consensus_sequence}{binder_string}\n")


  # Write combined consensus sequences to output file
  with open(args.output_file, "w") as handle:
    handle.writelines(consensus_sequences)
  print(f"Combined consensus sequences written to {args.output_file}.")




if __name__ == "__main__":
  main()

