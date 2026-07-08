import json
import os

output_dir = r"D:\research of the GBM"
remap_file = os.path.join(output_dir, "tfbs_ReMapTFs_chi3l1.json")

target_tfs = ["STAT3", "CEBPB", "RELA", "NFKB1", "NFKB2", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "YAP1"]

if os.path.exists(remap_file):
    with open(remap_file, 'r') as f:
        data = json.load(f)
        
        # Extract items
        items = []
        for k, v in data.items():
            if isinstance(v, dict) and "tracks" in v:
                items = v["tracks"].get("ReMapTFs", [])
                break
        
        print("Experimental ChIP-seq binding events on CHI3L1 promoter:")
        tf_hits = {}
        for item in items:
            tf = item.get("TF")
            if tf in target_tfs:
                if tf not in tf_hits:
                    tf_hits[tf] = []
                tf_hits[tf].append(item)
                
        for tf in target_tfs:
            hits = tf_hits.get(tf, [])
            if hits:
                print(f"\n--- {tf} binding sites (N={len(hits)}) ---")
                # Group by coordinate to show unique peaks
                unique_peaks = {}
                for h in hits:
                    coord = f"{h['chrom']}:{h['chromStart']}-{h['chromEnd']}"
                    biotype = h.get("Biotypes", "Unknown")
                    if coord not in unique_peaks:
                        unique_peaks[coord] = []
                    unique_peaks[coord].append(biotype)
                
                for coord, biotypes in list(unique_peaks.items())[:5]:
                    # Print unique peak and some cell types it was found in
                    cells = ", ".join(list(set(biotypes))[:4])
                    print(f"  Peak: {coord} | Cell Types: {cells}")
                    
else:
    print("ReMap TFBS file not found.")
