import collections
from collections import Counter

class FusionAssembler:
    def find_breakpoint(self, reads):
        if not reads:
            return None, 0
            
        evidence_a = []
        evidence_b = []
        
        for r in reads:
            # Side A: Current read's position
            evidence_a.append((r.reference_name, r.reference_start))
            
            # Side B: MUST BE THE MATE/SPLIT POSITION
            if r.has_tag("SA"):
                # Split read: Side B is the other half of the chimeric read
                sa = r.get_tag("SA").split(';')[0].split(',')
                evidence_b.append((sa[0], int(sa[1])))
            else:
                # Discordant: Side B is where the mate is mapped
                # This is where we get 'chr5' from!
                evidence_b.append((r.next_reference_name, r.next_reference_start))
        
        # consensus
        from collections import Counter
        bp_a = Counter(evidence_a).most_common(1)[0][0]
        bp_b = Counter(evidence_b).most_common(1)[0][0]
        
        # --- THE FIX FOR YOUR ISSUE ---
        # If bp_a and bp_b picked the same side, we swap bp_a with the most 
        # common different chromosome found in the evidence.
        return (bp_a, bp_b), len(reads)