package org.stjude.compbio.rnapeg;

public enum ExonMatchType {
  PERFECT,
    // perfectly matches reference exon
  SAME_START_OR_END,
    // overlaps and shares a start or end site with a reference
  PARTIAL,
    // overlaps, but without same start/end
  NO_MATCH;
    // novel
}