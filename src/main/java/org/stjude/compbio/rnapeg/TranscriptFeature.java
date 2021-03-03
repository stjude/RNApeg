package org.stjude.compbio.rnapeg;

import java.awt.Rectangle;

public class TranscriptFeature {
  Transcript transcript;
  // parent
  TranscriptFeatureType type;
  int feature_number;
  // e.g. exon or intron number

  String reference, gene_name;
  // reference sequence name (i.e. chrom), gene symbol
  // track these at this level in case of translocations
  int start, end;

  boolean is_perfect_reference_match;

  //
  // cooked data:
  //
  Rectangle screen_position;
  // X coordinates postprocessed for screen display purposes.
  //
  // FIX ME: kind of clunky to store here, but most convenient.
  // alternative: copy objects and modify coordinates there?

  public RangeDigest get_range_digest() {
    RangeDigest rd = new RangeDigest();
    rd.add_start(start);
    rd.add_end(end);
    return rd;
  }

  public int get_size() {
    return (end - start) + 1;
  }

  public String get_basic_hash_string() {
    return type.toString() + ":" + reference + ":" + start + "-" + end;
  }

  public String get_reference() {
    return reference;
  }

  public boolean overlaps (int start, int end) {
    System.err.println("overlaps " + this.start + "->" + this.end + " " + start + "->" + end);  // debug

    //    return !(end <= this.start || start >= this.end);
    boolean result;
    if (end < this.start) {
      System.err.println("1");  // debug
      result = false;
    } else if (start > this.end) {
      System.err.println("2");  // debug
      result = false;
    } else {
      result = true;
    }
    System.err.println("result="+result);  // debug
    System.exit(1);
    return result;

  }


}