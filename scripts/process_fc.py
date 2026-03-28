import os, sys
import pandas as pd

def extract_sample_name(bam_path):
    base = os.path.basename(bam_path)
    if "Aligned.sortedByCoord.out.bam" in base:
        return base.replace("Aligned.sortedByCoord.out.bam", "")
    base = base.replace(".bam", "")
    if "_S" in base:
        return base.split("_S")[0]
    return base

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_fc.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)
    input_file, output_file = sys.argv[1], sys.argv[2]
    fc = pd.read_csv(input_file, sep="\t", header=1)
    fc = fc.drop(columns=["Chr","Start","End","Strand"], errors="ignore")
    fc.columns = [c if c in ("Geneid","Length") else extract_sample_name(c) for c in fc.columns]
    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
    fc.to_csv(output_file, sep="\t", index=False)
    print(f"Processed {len(fc)} genes × {len(fc.columns)-2} samples → {output_file}")

if __name__ == "__main__":
    main()
