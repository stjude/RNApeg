package org.stjude.compbio.rnapeg;
// TO DO:
// - track exact genes hit

import java.io.*;
import java.util.*;
import java.util.zip.*;

public class IntronCache {
  //
  // load refGene-style records and cache intron positions by reference sequence.
  // should be more efficient since we are loading a very large number of these.
  //
  //  HashMap<String,HashSet<String>> cache;
  // reference sequence -> junctions
  HashMap<String,HashMap<String,HashSet<UCSCRefGene>>> cache;
  // reference sequence -> junction -> matching UCSCRefGene objects
  private String strand_filter = null;
  
  public IntronCache (UCSCRefGeneReader r) throws IOException {
    setup(r);
  }

  public IntronCache (UCSCRefGeneReader r, String strand_filter) throws IOException {
    this.strand_filter = strand_filter;
    setup(r);
  }

  private void setup (UCSCRefGeneReader r) throws IOException {
    //    System.err.println("setup() start");  // debug
    cache = new HashMap<String,HashMap<String,HashSet<UCSCRefGene>>>();

    long start = System.currentTimeMillis();
    int skipped = 0;

    for (UCSCRefGene rg : r) {

      boolean usable = true;

      if (strand_filter != null) {
	//	System.err.println("transcript: " + rg.name + " strand:" + rg.strand);  // debug
	if (!rg.strand.equals(strand_filter)) {
	  //	  System.err.println("skipping " + rg.name + " on " + rg.chrom + " " + rg.strand);  // debug
	  usable = false;
	}
      }

      if (usable) {
	//	System.err.println("using " + rg.name + " on " + rg.chrom + " " + rg.strand);  // debug
	Transcript t = new Transcript(rg, true);
	for (TranscriptIntron ti : t.introns) {
	  String ref_name = Chromosome.standardize_name(ti.get_reference());
	  if (ref_name.indexOf("chr") != 0)
	    System.err.println("WARNING: non-chr* reference sequence " + ref_name);

	  HashMap<String,HashSet<UCSCRefGene>> junction_bucket = cache.get(ref_name);
	  if (junction_bucket == null) {
	    //	    System.err.println("create bucket for " + ref_name);  // debug
	    cache.put(ref_name, junction_bucket = new HashMap<String,HashSet<UCSCRefGene>>());
	  }
	  String range = ti.get_range_digest().toString();

	  HashSet<UCSCRefGene> rg_bucket = junction_bucket.get(range);
	  if (rg_bucket == null) junction_bucket.put(range, rg_bucket = new HashSet<UCSCRefGene>());

	  rg_bucket.add(rg);
	}
      } else {
	skipped++;
      }
    }

    if (skipped > 0) System.err.println("transcripts skipped on undesired strand: " + skipped);  // debug


    System.err.println("intron cache:");  // debug
    int total = 0;
    for (String key : cache.keySet()) {
      int count = cache.get(key).size();
      System.err.println("  " + key + ": " + count);
      total += count;
    }
    System.err.println("  total="+total);  // debug
    System.err.println("flatfile db load: " + (System.currentTimeMillis() - start) + " ms");  // debug
  }

  public HashSet<UCSCRefGene> find_exon_junction (SplicedReadInfo sri) {
    HashSet<UCSCRefGene> result = null;
    String ref_name = Chromosome.standardize_name(sri.reference_name);
    HashMap<String,HashSet<UCSCRefGene>> junction_bucket = cache.get(ref_name);
    if (junction_bucket == null) {
      if (Chromosome.valueOfString(sri.reference_name) != null) 
	System.err.println("no gene annotations for " + sri.reference_name);
      // happens in mm9 for chrM (no gene annotations).
      // "mostly harmless"
    } else {
      String key = (sri.segment_1_end + 1) +  "-" +
      (sri.segment_2_start - 1);
      // convert from exon bases to intron bases
      //      System.err.println("key="+key + " bucket="+ bucket.size());  // debug
      result = junction_bucket.get(key);
    }
    return result;
  }

}
