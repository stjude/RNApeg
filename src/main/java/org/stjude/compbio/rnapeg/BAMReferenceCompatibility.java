package org.stjude.compbio.rnapeg;
// check compatibility between a BAM header and a reference sequence.
// MNE 5/2013
//
// TO DO:
// - return array of available reference sequence indexes
// - differentiate between problems it may be possible to work around
//   (e.g. missing secondary sequences) vs. severe errors (e.g. sequence
//   length mismatches)
// 

import htsjdk.samtools.*;
import java.util.*;
import java.io.IOException;

public class BAMReferenceCompatibility {
  private SAMFileHeader sfh;
  private ReferenceSequence rs;
  boolean HAS_MISSING_SEQUENCES, HAS_LENGTH_MISMATCHES;
  private ArrayList<String> errors;
  private ArrayList<String> compatible_sequences;

  public BAMReferenceCompatibility (SAMFileHeader sfh, ReferenceSequence rs) throws IOException {
    this.sfh = sfh;
    this.rs = rs;
    HAS_MISSING_SEQUENCES = HAS_LENGTH_MISMATCHES = false;
    setup();
  }

  private void setup() throws IOException {
    errors = new ArrayList<String>();
    compatible_sequences = new ArrayList<String>();
    if (rs.supports_sequence_list()) {
      //
      // get list of sequences supported by ReferenceSequence implementation:
      //
      ArrayList<String> ref_names = rs.get_sequence_names();

      System.err.println("reference sequence names: " + ref_names);  // debug

      ChromosomeDisambiguator cd = new ChromosomeDisambiguator(ref_names);

      //
      // for each reference sequence in the BAM,
      //  - does it exist in the reference set? (correct for ambiguity)
      //  - does the length match?
      //
      SAMSequenceDictionary ssd = sfh.getSequenceDictionary();
      for (SAMSequenceRecord ssr : ssd.getSequences()) {
	String ref_name_bam = ssr.getSequenceName();
	String ref_name_ref = cd.find(ref_name_bam);
	System.err.println("ref_name_bam:" + ref_name_bam + " ref_name_ref:" + ref_name_ref);  // debug


	if (ref_name_ref == null) {
	  // sequence in BAM header is not present in ReferenceSequence
	  HAS_MISSING_SEQUENCES = true;
	  errors.add("no local reference sequence found for BAM reference name " + ref_name_bam);
	} else {
	  // compare sequence lengths
	  int ref_length_ref = rs.get_length(ref_name_ref);
	  if (ssr.getSequenceLength() == ref_length_ref) {
	    compatible_sequences.add(ref_name_bam);
	  } else {
	    HAS_LENGTH_MISMATCHES = true;
	    errors.add("reference length mismatch for " + ref_name_bam + " bam=" + ssr.getSequenceLength() + " reference:" + ref_length_ref);
	  }
	}

      }
    } else {
      throw new IOException("reference sequence doesn't support sequence list");
    }
  }

  public boolean has_any_incompatibility() {
    // broadest measure of incompatibility
    return HAS_MISSING_SEQUENCES || HAS_LENGTH_MISMATCHES;
  }

  public void report_errors () {
    for (String msg : errors) {
      System.err.println("ERROR: " + msg);  // debug
    }
  }

  public ArrayList<String> get_compatible_reference_sequences() {
    return compatible_sequences;
  }

}
