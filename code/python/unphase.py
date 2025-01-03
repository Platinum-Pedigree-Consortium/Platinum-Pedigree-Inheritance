import argparse

def unphase_vcf(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):  # Write header lines unchanged
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            for i in range(9, len(fields)):  # Genotype fields start at column 10
                genotype = fields[i]
                if "|" in genotype:  # Phased genotype
                    alleles = genotype.split("|")
                    if "." in alleles:  # Handle no-calls explicitly
                        fields[i] = "./." if alleles[0] == "." and alleles[1] == "." else f"./{alleles[1]}" if alleles[0] == "." else f"./{alleles[0]}"
                    else:  # Sort the alleles numerically
                        sorted_alleles = sorted(alleles, key=lambda x: int(x))
                        fields[i] = f"{sorted_alleles[0]}/{sorted_alleles[1]}"
                elif "/" in genotype:  # Already unphased, ensure sorting
                    alleles = genotype.split("/")
                    if "." in alleles:  # Handle no-calls explicitly
                        fields[i] = "./." if alleles[0] == "." and alleles[1] == "." else f"./{alleles[1]}" if alleles[0] == "." else f"./{alleles[0]}"
                    else:  # Sort numerically
                        sorted_alleles = sorted(alleles, key=lambda x: int(x))
                        fields[i] = f"{sorted_alleles[0]}/{sorted_alleles[1]}"
            outfile.write("\t".join(fields) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Unphase a VCF file and sort genotype alleles.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Path to the output VCF file")
    
    args = parser.parse_args()
    
    unphase_vcf(args.input, args.output)
    print(f"Unphased VCF written to {args.output}")

if __name__ == "__main__":
    main()
