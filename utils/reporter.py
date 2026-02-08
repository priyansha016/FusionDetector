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
                # Report positions in 1-based coordinates (genomics standard; assembler uses 0-based)
                out = {k: v for k, v in row.items() if k in fieldnames}
                if 'pos_a' in out and out['pos_a'] is not None:
                    out['pos_a'] = int(out['pos_a']) + 1
                if 'pos_b' in out and out['pos_b'] is not None:
                    out['pos_b'] = int(out['pos_b']) + 1
                writer.writerow(out)
                
        print(f"[*] Report saved: {self.output_path}")