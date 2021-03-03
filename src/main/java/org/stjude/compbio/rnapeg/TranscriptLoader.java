package org.stjude.compbio.rnapeg;

import java.util.*;

public abstract class TranscriptLoader {
  // wrapper for a collection of Transcript objects for a single sample
  // from various possible sources:
  //
  // - to start: Michael Rusch's UCSC GTF tracks
  // - maybe later:
  //   - simplified flatfile, etc.
  //   - database

  //  FIX ME:
  //   - option to parse in new thread??
  //

  public static final String GENE_UNKNOWN = "unknown_or_span";

  protected HashMap<String,Transcript> id2transcript;
  protected ArrayList<Transcript> all_transcripts;

  public TranscriptLoader () {
    {
      // initializer block: should apply to all constructors
      // (i.e. even for subclasses)
      id2transcript = new HashMap<String,Transcript>();
      all_transcripts = new ArrayList<Transcript>();
    }
  }

  protected void transcript_setup() {
    // shared initialization after raw data loaded
    for (Transcript t : all_transcripts) {
      t.setup();
    }

  }

  public ArrayList<Transcript> get_transcripts() {
    return all_transcripts;
  }

}

