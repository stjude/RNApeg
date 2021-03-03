package org.stjude.compbio.rnapeg;

import java.awt.Rectangle;

public class TranscriptIntron extends TranscriptFeature {
  int junction_coverage;
  static final int COVERAGE_NULL = -1;

  public TranscriptIntron () {
    type = TranscriptFeatureType.INTRON;
    junction_coverage = COVERAGE_NULL;
  }

  public TranscriptIntron (Transcript transcript, TranscriptExon te1, TranscriptExon te2) {
    type = TranscriptFeatureType.INTRON;
    this.transcript = transcript;
    this.reference = te1.reference;
    feature_number = te1.feature_number;
    start = te1.end + 1;
    end = te2.start - 1;
    junction_coverage = COVERAGE_NULL;
    //    System.err.println("new intron: " + intron_number + " " + start + " " + end);  // debug
    
  }

}
