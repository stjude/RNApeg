package org.stjude.compbio.rnapeg;

import java.util.*;
import java.awt.Rectangle;
import java.security.MessageDigest;

public class Transcript {
  // raw/user-specified fields:
  public String sample_name;
  // e.g. from source file
  public String transcript_id;
  public String strand;
  public ArrayList<TranscriptExon> exons;

  // cooked/generated fields:
  public ArrayList<TranscriptIntron> introns;
  private boolean is_translocation = false;
  private HashSet<String> all_genes;
  // FIX ME: case sensitivity?
  public RangeDigest screen_position;
  // X coordinates postprocessed for screen display purposes
  public boolean is_reference = false;
  public int product_size = 0;
  public int span_size = 0;
  public int total_intron_size = 0;
  MD5Util md5 = null;
  public Transcript duplicate = null;

  public Transcript() {
    exons = new ArrayList<TranscriptExon>();
  }

  public Transcript (UCSCRefGene rg, boolean fix_exon_starts) {
    exons = new ArrayList<TranscriptExon>();
    sample_name = rg.name;
    transcript_id = rg.name;
    strand = rg.strand;
    for (int i = 0; i < rg.exonCount; i++) {
      TranscriptExon te = new TranscriptExon();
      te.transcript = this;
      te.feature_number = i + 1;
      te.reference = rg.chrom;
      te.gene_name = rg.name2;
      te.start = rg.exonStarts[i] + (fix_exon_starts ? 1 : 0);
      // UCSC's exonStarts seems to be in 0-based space,
      // and exonEnds in 1-based space.
      // tweak by converting exonStart to 1-based space.
      te.end = rg.exonEnds[i];
      exons.add(te);
    }
    setup();
  }

  public void set_gene (String gene) {
    // hack to manually set gene symbol for transcript after the fact
    // (e.g. if set from isoform db rather than refGene table)
    all_genes = new HashSet<String>();
    all_genes.add(gene);
    for (TranscriptExon te : exons) {
      te.gene_name = gene;
    }
  }

  public void setup() {
    //
    // call after transcript has been fully loaded/populated
    //
    all_genes = new HashSet<String>();
    for (TranscriptExon te : exons) {
      te.raw_start = te.start;
      te.raw_end = te.end;
      // back up raw start and end, as we may tweak these
      // first 1st/last exon to match reference
      all_genes.add(te.gene_name);
    }
    is_translocation = all_genes.size() > 1;

    product_size = 0;
    for (TranscriptExon e : exons) {
      product_size += e.get_size();
    }

    span_size = (exons.get(exons.size() - 1).end - exons.get(0).start) + 1;

    generate_introns();

    total_intron_size = 0;
    for (TranscriptIntron i : introns) {
      total_intron_size += i.get_size();
    }

  }

  public boolean is_translocation() {
    return is_translocation;
  }

  private void generate_introns() {
    //
    // generate intron objects based on exon sites
    //
    introns = new ArrayList<TranscriptIntron>();
    int size = exons.size();
    for (int i=0; i < size - 1; i++) {
      introns.add(new TranscriptIntron(
				       this,
				       exons.get(i),
				       exons.get(i + 1)));
    }
  }

  public void reset_exon_positions() {
    for (TranscriptExon te : exons) {
      te.start = te.raw_start;
      te.end = te.raw_end;
    }
  }

  public boolean contains_gene (String gene) {
    return all_genes.contains(gene);
  }

  public boolean has_null_genes() {
    return all_genes.size() == 1 && contains_gene(TranscriptLoader.GENE_UNKNOWN);
  }
  
  public Collection<String> get_genes() {
    return all_genes;
  }

  public void generate_screen_position() {
    screen_position = new RangeDigest();
    for (TranscriptExon te : exons) {
      screen_position.add(te.screen_position);
    }
  }

  public RangeDigest get_range_digest() {
    RangeDigest rd = new RangeDigest();
    for (TranscriptExon te : exons) {
      rd.add_start(te.start);
      rd.add_end(te.end);
    }
    return rd;
  }

  public boolean is_plus_strand() {
    //    System.err.println("strand="+strand);  // debug
    return strand.equals("+");
  }

  public void remap_ends_to_reference (ArrayList<Transcript> t_others) {
    // many transcripts appear to be missing precise (a) start positions
    // for the first exon and (b) end positions for the last exon.
    // This may be due to datasets built from intra-exonic 
    // splice junctions only.
    //
    // first exon: modify start only if:
    //  - exon does not already perfectly match a known exon, and
    //  - exon end exactly matches a reference exon end
    TranscriptExon e_first = exons.get(0);
    boolean already_perfect = false;
    for (Transcript t_other : t_others) {
      for (TranscriptExon e_other : t_other.exons) {
	if (e_first.start == e_other.start &&
	    e_first.end == e_other.end) {
	  already_perfect = true;
	  break;
	}
      }
    }
    if (!already_perfect) {
      for (Transcript t_other : t_others) {
	for (TranscriptExon e_other : t_other.exons) {
	  if (e_first.end == e_other.end) {
	    e_first.start = e_other.start;
	    break;
	  }
	}
      }
    }

    // last exon: modify start only if:
    //  - exon does not already perfectly match a known exon, and
    //  - exon start exactly matches a reference exon start
    TranscriptExon e_last = exons.get(exons.size() - 1);
    already_perfect = false;
    for (Transcript t_other : t_others) {
      for (TranscriptExon e_other : t_other.exons) {
	if (e_last.start == e_other.start &&
	    e_last.end == e_other.end) {
	  already_perfect = true;
	  break;
	}
      }
    }
    if (!already_perfect) {
      for (Transcript t_other : t_others) {
	for (TranscriptExon e_other : t_other.exons) {
	  if (e_last.start == e_other.start) {
	    e_last.end = e_other.end;
	    break;
	  }
	}
      }
    }
  }

  public void compare_with (ArrayList<Transcript> t_others, HashSet<String> all_reference_introns) {
    if (t_others.size() == 0) return;

    //
    //  intron comparison:
    //
    for (TranscriptIntron ti : introns) {
      ti.is_perfect_reference_match = all_reference_introns.contains(ti.get_range_digest().toString());
    }

    //
    //  exon comparison:
    //
    for (TranscriptExon te : exons) {
      HashSet<ExonMatchType> matches = new HashSet<ExonMatchType>();
      te.is_putative_span = false;
      te.is_perfect_reference_match = false;

      for (Transcript t_other : t_others) {

	TranscriptExon hit_start = null;
	TranscriptExon hit_end = null;

	for (TranscriptExon te_other : t_other.exons) {
	  if (!(te.end < te_other.start || te.start > te_other.end)) {
	    // some overlap
	    ExonMatchType match_type;
	    if (te.start == te_other.start && te.end == te_other.end) {
	      match_type = ExonMatchType.PERFECT;
	      te.is_perfect_reference_match = true;
	    } else if (te.start == te_other.start) {
	      match_type = ExonMatchType.SAME_START_OR_END;
	      hit_start = te_other;
	    } else if (te.end == te_other.end) {
	      match_type = ExonMatchType.SAME_START_OR_END;
	      hit_end = te_other;
	    } else {
	      match_type = ExonMatchType.PARTIAL;
	    }
	    matches.add(match_type);
	  }
	}

	if (hit_start != null && hit_end != null) {
	  te.is_putative_span = true;
	  te.putative_span_end = hit_end;
	}

      }

      if (matches.size() == 0) {
	te.reference_match_type = ExonMatchType.NO_MATCH;
      } else {
	if (matches.size() > 1) matches.remove(ExonMatchType.NO_MATCH);
	if (matches.size() > 1) matches.remove(ExonMatchType.PARTIAL);
	if (matches.size() > 1) matches.remove(ExonMatchType.SAME_START_OR_END);
	ArrayList<ExonMatchType> results = new ArrayList<ExonMatchType>(matches);
	te.reference_match_type = results.get(0);
	//	System.err.println("final="+te.reference_match_type);  // debug

      }


    }
  }

  public int get_first_exon_start() {
    return exons.get(0).start;
  }

  public int get_last_exon_end() {
    return exons.get(exons.size() - 1).end;
  }

  public int get_exon_count() {
    return exons.size();
  }

  private void build_md5() {
    if (md5 == null) {
      try {
	md5 = new MD5Util();
	md5.reset();
	ArrayList<TranscriptFeature> features = new ArrayList<TranscriptFeature>();
	features.addAll(exons);
	features.addAll(introns);
	for (TranscriptFeature tf : features) {
	  //	System.err.println("prep: " + tf.start + "-" + tf.end);  // debug
	  md5.update(tf.start);
	  md5.update(tf.end);
	}
      } catch (Exception e) {
	System.err.println("ERROR: " + e);  // debug
      }
    }
  }

  public byte[] get_md5 () {
    // of limited utility since entries frequently consist of only
    // partial matches, skip exons, etc.
    // Still, maybe some evidence of duplicates for longer exon
    // counts w/same starts, e.g. 316528 vs 367812 in MYB
    build_md5();
    return md5.get_digest();
  }

  public String get_md5_hex () {
    // printable
    build_md5();
    return md5.get_hex_digest();
  }


}
