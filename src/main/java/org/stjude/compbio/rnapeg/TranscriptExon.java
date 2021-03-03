package org.stjude.compbio.rnapeg;

import java.awt.Rectangle;

public class TranscriptExon extends TranscriptFeature {

  int raw_start, raw_end;

  ExonMatchType reference_match_type;
  boolean is_putative_span;
  // exon start matches the start of one exon and the end of another.
  TranscriptExon putative_span_end;
  int drawn_right_y;
  int next_junction_coverage;
  // temporary

  public TranscriptExon () {
    type = TranscriptFeatureType.EXON;
    reference_match_type = null;
    is_putative_span = false;
    putative_span_end = null;
    next_junction_coverage = -1;
  }

}
