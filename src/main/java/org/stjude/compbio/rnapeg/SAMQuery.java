package org.stjude.compbio.rnapeg;

import htsjdk.samtools.*;

// simple wrapper for querying:
//   - an entire BAM
//   - reads mapped to a reference
//   - reads mapped to a reference, start and end
// automatically:
//   - disambiguates reference names
//   - looks up reference lengths
//
// MNE 4/2013

public class SAMQuery {
  SamReader sfr;
  ChromosomeDisambiguator cd;

  public SAMQuery (SamReader sfr) {
    this.sfr = sfr;
    cd = new ChromosomeDisambiguator(sfr);
  }

  public SAMRecordIterator query (SAMRegion region) {
    SAMRecordIterator result = null;
    if (region == null || region.tname == null) {
      //
      // all BAM reads:
      //
      result = sfr.iterator();
    } else if (region.isValid()) {
      //
      // reads mapped to a reference, overlapping a start and end:
      //
      result = sfr.queryOverlapping(cd.find(region.tname),
				    region.range.start,
				    region.range.end);
    } else {
      //
      // reads mapped to an entire reference:
      //
      String ref_name = cd.find(region.tname);
      result = sfr.queryOverlapping(ref_name,
				    1,
				    cd.get_length(ref_name));
    }

    return result;
  }

}
