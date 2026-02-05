import csv
import os

class FusionReporter:
    def __init__(self, output_path):
        self.output_path = output_path

    def save_tsv(self, results):
        if not results:
            # Create an empty file with just the header if no fusions found
            print(f"[!] No fusions found for this sample.")
        
        # Expanded fieldnames to match your new Assembler output
        fieldnames = [
            'sample', 'gene_a', 'gene_b', 'chrom_a','pos_a', 
            'chrom_b', 'pos_b', 'support'
        ]
        
        with open(self.output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in results:
                # Ensure we only write keys that exist in our fieldnames
                filtered_row = {k: v for k, v in row.items() if k in fieldnames}
                writer.writerow(filtered_row)
                
        print(f"[*] Report saved: {self.output_path}")