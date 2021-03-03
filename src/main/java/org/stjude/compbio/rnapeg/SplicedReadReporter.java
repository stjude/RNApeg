package org.stjude.compbio.rnapeg;

import htsjdk.samtools.*;
import java.io.*;
import java.util.*;

public class SplicedReadReporter {
  File bam;
  HashMap<String,SplicedReadInfo> events_by_name;
  IntronCache ic;
  //  int minimum_observations_to_report = 5;
  int minimum_observations_to_report = 1;
  // previously we required a minimum of 5 supporting reads to
  // support a junction.  However:
  //
  // 1. raw junctions are uncorrected: after ambiguity correction,
  //    the final count may (greatly) increase.  Thus apparently
  //    marginal junctions may be salvagable and shouldn't be discarded
  //    at this early step.
  // 2. for reference junctions we want to be very conservative in any case
  //    (and the same situation applies: a low reference junction count
  //    may be supplemented by correctable novel junctions)
  //

  int flush_check_interval = 100000;
  boolean wrote_header;
  HashCounter junction_counter;

  boolean READ_NAME_RESTRICT_MODE = false;
  //    static String READ_NAME_RESTRICT_STRING = "HWUSI-EAS1714_15:8:119:8487:14583";

  boolean PRIMARY_READ_MODE = false;
  // whether to process only reads with primary mappings 
  // (see SAM spec)

  boolean EXCLUDE_DUPLICATES = false;
  // whether to exclude optical/PCR duplicates

  boolean IGNORE_INCOMPATIBLE_REFERENCE = false;

  String READ_NAME_RESTRICT_STRING = null;

  int RESTRICT_JUNCTION_START = 0;
  int RESTRICT_JUNCTION_END = 0;

  boolean novel_only = true;
  boolean write_header = true;

  WorkingFile wf = null;
  PrintStream ps = null;

  Reporter rpt_reads = null;

  boolean VERBOSE = false;
  boolean SORT = false;
  // still problems w/sorting

  boolean WRITE_BED = false;
  boolean WRITE_DELIMITED = false;
  boolean ANNOTATE_MODE = false;

  //  String RGB_NOVEL = "255,0,0";
  String RGB_NOVEL = "192,0,0";

  //  String RGB_KNOWN = "0,255,0";
  String RGB_KNOWN = "0,128,0";

  // TO DO:
  // maybe full brightness if covered in both directions,
  // dimmer value if only in 1 direction?

  int QC_MIN_FLANKING_NT_READ_EDGE = 18;
  // minimum nt of mapped sequence between a junction and a read start/end.
  // while this might be legitimate we want to avoid calling junctions
  // where the ONLY evidence is very limited as these short sequences 
  // imply ambiguous/unreliable mappings due to insufficient anchoring.
  // In other words if we're too near the read start/end the data are
  // too sparse for a reliable mapping.
  // Use a value that has a decent chance of being genomically unique
  // (e.g. a PCR primer size)

  int QC_MIN_FLANKING_NT_INTERNAL = 10;
  // some reads contain evidence for multiple junctions.
  // junctions internal to a sequence are allowed more leeway
  // because (a) they are anchored by additional sequence beyond
  // that immediately flanking the junction and (b) some exons
  // are legimately small.

  int MIN_QUALITY_FOR_MISMATCH_CHECK = 15;

  double MAX_JUNK_RATIO_FOR_CLEAN_CALL = 0.05;

  ReferenceSequence reference_sequence = null;

  boolean JUNK_REPORT = false;
  double JUNK_REPORT_MIN_RATIO_TO_COMPLAIN = 0.10;
  double JUNK_REPORT_MIN_LENGTH_TO_COMPLAIN = 15;

  //  SAMRegion query_region = null;
  ArrayList<SAMRegion> query_regions;

  int MIN_SPAN = 0;

  public SplicedReadReporter() {
    ic = null;
    //    query_region = new SAMRegion();
    query_regions = new ArrayList<SAMRegion>();
  }

  public SAMRegion get_query_region() {
    //    return query_region;
    return query_regions.get(query_regions.size() - 1);
  }

  public SAMRegion add_query_region() {
    query_regions.add(new SAMRegion());
    return get_query_region();
  }

  public void set_known_rgb (String rgb) {
    RGB_KNOWN = rgb;
  }

  public void set_novel_rgb (String rgb) {
    RGB_NOVEL = rgb;
  }
  
  public void set_reference_sequence (ReferenceSequence rs) {
    reference_sequence = rs;
  }

  public void set_primary_read_mode (boolean v) {
    PRIMARY_READ_MODE = v;
  }

  public void set_exclude_duplicates (boolean v) {
    EXCLUDE_DUPLICATES = v;
  }

  public void set_ignore_incompatible (boolean v) {
    IGNORE_INCOMPATIBLE_REFERENCE = v;
  }

  public void set_restrict_read_name (String name) {
    READ_NAME_RESTRICT_MODE = true;
    READ_NAME_RESTRICT_STRING = name;
  }

  public void set_restrict_junction_start (int start) {
    RESTRICT_JUNCTION_START = start;
  }

  public void set_restrict_junction_end (int end) {
    RESTRICT_JUNCTION_END = end;
  }

  public void set_bam (File bam) {
    this.bam = bam;
  }

  public void set_intron_cache (IntronCache ic) {
    this.ic = ic;
  }

  public void set_write_header (boolean v) {
    write_header = v;
  }

  public void set_novel_only(boolean v) {
    novel_only = v;
  }

  public void set_annotate_mode (boolean v) {
    novel_only = false;
    ANNOTATE_MODE = true;
    WRITE_DELIMITED = true;
  }

  public void set_outfile (String outfile) throws FileNotFoundException,IOException {
    wf = new WorkingFile(outfile);
    ps = wf.getPrintStream();
  }

  public void set_min_span (int min_span) {
    MIN_SPAN = min_span;
  }

  public void set_read_report (String rr_fn) throws FileNotFoundException,IOException {
    rpt_reads = new Reporter();
    rpt_reads.set_output_filename(rr_fn);
    rpt_reads.add_header("read_name");
    rpt_reads.add_header("strand");
    rpt_reads.add_header("reference");
    rpt_reads.add_header("junction_start");
    rpt_reads.add_header("junction_end");
    rpt_reads.add_header("flanking_nt_left");
    rpt_reads.add_header("flanking_nt_right");
    rpt_reads.add_header("mapping_perfect");
    rpt_reads.add_header("mapping_clean");
  }

  public void set_minimum_observations_to_report (int count) {
    minimum_observations_to_report = count;
  }

  public void report() throws IOException{
    if (bam == null) throw new IOException("no bam specified (-bam)");
    if (wf == null) throw new IOException("no outfile specified (-of)");
    if (reference_sequence == null) throw new IOException("no reference sequence specified: use -2bit FILE or -fasta FILE (.fai indexed)");

    wrote_header = write_header ? false : true;
    junction_counter = new HashCounter();

    events_by_name = new HashMap<String,SplicedReadInfo>();
    SamReader sfr = SamReaderFactory.makeDefault().open(bam);

    //
    // check compatibility between reference sequence and BAM:
    //
    System.err.print("checking BAM/reference compatibility...");  // debug
    BAMReferenceCompatibility brc = new BAMReferenceCompatibility(
								  sfr.getFileHeader(),
								  reference_sequence
								  );
    System.err.println("done");  // debug
    if (brc.has_any_incompatibility()) {
      brc.report_errors();
      String msg = "BAM header not fully compatible with specified reference sequence";
      if (IGNORE_INCOMPATIBLE_REFERENCE) {
	if (query_regions.size() > 0) throw new IOException("can't manually specify query regions with partially incompatible BAM header");
	ArrayList<String> compatible = brc.get_compatible_reference_sequences();
	if (compatible.size() > 0) {
	  for (String name : compatible) {
	    System.err.println("compatible sequence: " + name);  // debug
	    query_regions.add(new SAMRegion(name));
	  }
	} else {
	  throw new IOException("no compatible reference sequences in BAM header");
	}
      } else {
	throw new IOException(msg);
      }
    }

    CigarOperator co;
    int len;
    int current_reference_i = -1;
    String current_reference_name = null;
    int checkpoint = 0;
    byte[] refseq = null;

    if (READ_NAME_RESTRICT_MODE) {
      System.err.println("DEBUG: restricting to read " + READ_NAME_RESTRICT_STRING);  // debug
    }

    int broken_edge_reads = 0;
    int broken_read_index = 0;
    int broken_qual_index = 0;

    if (query_regions.size() == 0) add_query_region();
    // if no query ranges set by user, add a blank region,
    // which will query the entire BAM.

    //    SAMRecordIterator query = sfr.iterator();

    // summarize:
    System.err.println("extracting primary reads only?: " + PRIMARY_READ_MODE);  // debug
    System.err.println("excluding optical/PCR duplicates?: " + EXCLUDE_DUPLICATES);  // debug

    final boolean VERBOSE = false;
    // debug

    for (SAMRegion query_region : query_regions) {
      SAMQuery sq = new SAMQuery(sfr);
      System.err.println("query region: " + query_region);  // debug

      SAMRecordIterator query = sq.query(query_region);

      while (query.hasNext()) {
	SAMRecord sr = query.next();

	if (sr.getReadUnmappedFlag()) continue;

	if (PRIMARY_READ_MODE && sr.getNotPrimaryAlignmentFlag()) continue;

	if (EXCLUDE_DUPLICATES && sr.getDuplicateReadFlag()) continue;

	if (READ_NAME_RESTRICT_MODE &&
	    !(sr.getReadName().equals(READ_NAME_RESTRICT_STRING))) continue;

	if (sr.getReferenceIndex() != current_reference_i) {
	  flush_check(true, null);
	  current_reference_i = sr.getReferenceIndex();
	  current_reference_name = Chromosome.standardize_name(sr.getReferenceName());
	  // e.g. .bed format requires "chr1", not "1"

	  //	if (VERBOSE) System.err.println("new reference: " + current_reference_name);  // debug
	  System.err.println("processing reference: " + current_reference_name);  // debug

	  if (reference_sequence != null) {
	    System.err.print("load reference for " + current_reference_name + "...");
	    refseq = reference_sequence.get_all(current_reference_name);
	    // FIX ME: may need disambiguation/lookup
	    for (int i = 0; i < refseq.length; i++) {
	      refseq[i] = (byte) Character.toUpperCase(refseq[i]);
	    }
	    System.err.println("done");  // debug
	  }

	} else if (++checkpoint % flush_check_interval == 0) {
	  flush_check(false, sr);
	}

	int ref_base = sr.getUnclippedStart();

	if (VERBOSE) System.err.println("new read:" + sr.getReadName() + " unclipped_start:" + ref_base);

	Cigar c = sr.getCigar();

	HashMap<SplicedReadInfo,SplicedReadFlankingInfo> sri2left = new HashMap<SplicedReadInfo,SplicedReadFlankingInfo>();
	HashMap<SplicedReadInfo,SplicedReadFlankingInfo> sri2right = new HashMap<SplicedReadInfo,SplicedReadFlankingInfo>();
      
	SplicedReadFlankingInfo fi = new SplicedReadFlankingInfo();

	byte[] read_bases = sr.getReadBases();
	byte[] read_qualities = sr.getBaseQualities();
	int read_i = 0;
	// for mismatch detection of course it's desirable to use
	// SAMRecord.getAlignmentBlocks() but then these have to be mapped
	// to the appropriate flanking regions.

	int i;
	byte base_ref;

	ArrayList<SplicedReadInfo> sri_ordered = new ArrayList<SplicedReadInfo>();
	// we could use LinkedHashMap for sri2left etc., but that
	// doesn't have get() methods

	for (CigarElement ce : c.getCigarElements()) {
	  co = ce.getOperator();
	  len = ce.getLength();
	  //	System.err.println("op " + co + " " + len);  // debug
	  if (co.equals(CigarOperator.MATCH_OR_MISMATCH)) {
	    for (i = 0; i < len; i++, ref_base++, read_i++) {
	      fi.last_read_i = read_i;
	      if (ref_base > refseq.length) {
		System.err.println("WARNING: read mapped beyond end of reference: " + sr.getReadName() + " at " + sr.getReferenceName() + "." + sr.getAlignmentStart());  // debug
		broken_edge_reads++;
		break;
	      }
	      base_ref = refseq[ref_base - 1];

	      if (read_i >= read_bases.length) {
		System.err.println("WARNING: read index beyond end of read bases: " + sr.getReadName() + " at " + sr.getReferenceName() + "." + sr.getAlignmentStart());
		broken_read_index++;
		break;
	      }
	      if (read_i >= read_qualities.length) {
		System.err.println("WARNING: read index beyond end of quality array: " + sr.getReadName() + " at " + sr.getReferenceName() + "." + sr.getAlignmentStart());
		broken_qual_index++;
		break;
	      }

	      if (
		  base_ref != read_bases[read_i] &&
		  read_qualities[read_i] >= MIN_QUALITY_FOR_MISMATCH_CHECK &&
		  (read_bases[read_i] == 'A' ||
		   read_bases[read_i] == 'C' ||
		   read_bases[read_i] == 'G' ||
		   read_bases[read_i] == 'T') &&
		  (base_ref == 'A' ||
		   base_ref == 'C' ||
		   base_ref == 'G' ||
		   base_ref == 'T')
		  ) {
		fi.count_aligned_mismatched_bases++;
		if (VERBOSE) System.err.println("mismatch base");  // debug
	      }
	    }
	    fi.count_aligned_bases += len;
	  } else if (co.equals(CigarOperator.SOFT_CLIP)) {
	    ref_base += len;
	    read_i += len;
	    fi.count_soft_clip_bases += len;
	  } else if (co.equals(CigarOperator.SKIPPED_REGION)) {
	    int seq1_end = ref_base - 1;
	    // base # before skip
	    int seq2_start = ref_base + len;
	    // base # after skip
	    SplicedReadInfo sri = track_site(current_reference_name, seq1_end, seq2_start, sr);
	    //	  System.err.println("money! splice=" + seq1_end + "-" + seq2_start + " insertion:" + has_insertion + " del:"+has_deletion + " read="+ sr.getReadName());

	    sri_ordered.add(sri);

	    sri2left.put(sri, fi);
	    // the current accumulator is shared, it provides right-flanking
	    // information for the previous skip (if any) as well as 
	    // left-flanking information for this skip

	    fi = new SplicedReadFlankingInfo();
	    sri2right.put(sri, fi);
	    // start new accumulator for right-flanking info for this skip

	    ref_base += len;
	    // read index: since a skip, no change
	  } else if (co.equals(CigarOperator.INSERTION)) {
	    // insertion does not affect reference sequence location
	    fi.count_inserted_bases += len;
	    read_i += len;
	    // inserted nucleotides appear in the read but not the reference
	  } else if (co.equals(CigarOperator.DELETION)) {
	    fi.count_deleted_bases += len;
	    ref_base += len;
	    // deleted nucleotides appear in the reference but not the read
	  } else if (co.equals(CigarOperator.HARD_CLIP)) {
	    // hard clips are in the reference (because we are using unclipped start)
	    // but not the sequence.
	    ref_base += len;
	  } else {
	    System.err.println("unhandled CIGAR operator " + co);  // debug
	    System.exit(1);
	  }
	}

	// CIGAR parsed:
	// evaluate left and right flanking nt and set QC flag

	SplicedReadInfo sri_first = null;
	SplicedReadInfo sri_last = null;
	if (sri_ordered.size() > 0) {
	  sri_first = sri_ordered.get(0);
	  sri_last = sri_ordered.get(sri_ordered.size() - 1);
	}

	for (SplicedReadInfo sri : sri2left.keySet()) {
	  //
	  //  foreach junction supported by read...
	  //
	  SplicedReadFlankingInfo fi_left = sri2left.get(sri);
	  SplicedReadFlankingInfo fi_right = sri2right.get(sri);

	  double junk_ratio_l = fi_left.get_junk_ratio();
	  double junk_ratio_r = fi_right.get_junk_ratio();

	  //
	  //  debug output for troubleshooting / pipeline QC:
	  //
	  boolean complain_l = JUNK_REPORT &&
	    junk_ratio_l >= JUNK_REPORT_MIN_RATIO_TO_COMPLAIN &&
	    fi_left.count_aligned_bases >= JUNK_REPORT_MIN_LENGTH_TO_COMPLAIN;

	  boolean complain_r = JUNK_REPORT &&
	    junk_ratio_r >= JUNK_REPORT_MIN_RATIO_TO_COMPLAIN &&
	    fi_right.count_aligned_bases >= JUNK_REPORT_MIN_LENGTH_TO_COMPLAIN;

	  if (complain_l || complain_r) {

	    System.err.print(sri.get_name() + 
			     " read=" + sr.getReadName() + 
			     " CIGAR=" + sr.getCigar());

	    System.err.print(" l_mm_count=" + fi_left.count_aligned_mismatched_bases);  // debug

	    if (complain_l) {
	      System.err.print(" LEFT:" +
			       " aligned=" + fi_left.count_aligned_bases + 
			       " junk=" + fi_left.get_junk_ratio() + 
			       " mm=" + fi_left.get_junk_ratio_mismatches() +
			       " sc=" + fi_left.get_junk_ratio_soft_clip() +
			       " indel=" + fi_left.get_junk_ratio_indel() 
			       );
	    }
	    if (complain_r) {
	      System.err.print(
			       " RIGHT:" +
			       " aligned=" + fi_right.count_aligned_bases + 
			       " junk=" + fi_right.get_junk_ratio() + 
			       " mm=" + fi_right.get_junk_ratio_mismatches() +
			       " sc=" + fi_right.get_junk_ratio_soft_clip() +
			       " indel=" + fi_right.get_junk_ratio_indel() 
			       );
	    }
	    System.err.println("");  // debug
	  }

	  boolean mapping_perfect = false;
	  boolean mapping_clean = false;

	  if (junk_ratio_l == 0 && junk_ratio_r == 0) {
	    mapping_perfect = true;
	    sri.counter_perfect++;
	    //	  sri.counter_clean++;
	    // should we also count perfect reads as "clean" (good-enough) reads?
	  } else if (junk_ratio_l <= MAX_JUNK_RATIO_FOR_CLEAN_CALL &&
		     junk_ratio_r <= MAX_JUNK_RATIO_FOR_CLEAN_CALL) {
	    mapping_clean = true;
	    sri.counter_clean++;
	  }

	  //	boolean flank_qc = (fi_left.count_aligned_bases >= MIN_FLANKING_NT_FOR_QC &&
	  //			    fi_right.count_aligned_bases >= MIN_FLANKING_NT_FOR_QC);
	  boolean flank_qc = get_flank_qc(fi_left, sri.equals(sri_first)) &&
	    get_flank_qc(fi_right, sri.equals(sri_last));

	  //	System.err.println(sri.get_name() + " L=" + flank_l + " R=" + flank_r + " qc=" + flank_qc);  // debug
	  if (flank_qc) sri.increment_flanking_qc();

	  if (rpt_reads != null) {
	    //
	    // write read-level report for each supported junction:
	    //
	    rpt_reads.set_value("read_name", sr.getReadName());
	    rpt_reads.set_value("strand", sr.getReadNegativeStrandFlag() ? "-" : "+");
	    rpt_reads.set_value("reference", Chromosome.standardize_name(sri.reference_name));
	    rpt_reads.set_value("junction_start", Integer.toString(sri.segment_1_end));
	    rpt_reads.set_value("junction_end", Integer.toString(sri.segment_2_start));
	    rpt_reads.set_value("flanking_nt_left", Integer.toString(fi_left.count_aligned_bases));
	    rpt_reads.set_value("flanking_nt_right", Integer.toString(fi_right.count_aligned_bases));
	    rpt_reads.set_value("mapping_perfect", mapping_perfect ? "1" : "0");
	    rpt_reads.set_value("mapping_clean", mapping_clean ? "1" : "0");

	    boolean usable = true;
	    if (RESTRICT_JUNCTION_START > 0 && sri.segment_1_end != RESTRICT_JUNCTION_START) usable = false;
	    if (RESTRICT_JUNCTION_END > 0 && sri.segment_2_start != RESTRICT_JUNCTION_END) usable = false;

	    if (usable) rpt_reads.end_row();
	  }

	}  // SplicedReadInfo

      }  // query.hasNext()

      flush_check(true, null);
      // finish this query

      query.close();

    }

    wf.finish();
    if (rpt_reads != null) rpt_reads.close();

    if (broken_edge_reads > 0) {
      System.err.println("reads mapped beyond end of reference: " + broken_edge_reads);  // debug
    }

    if (broken_read_index > 0) {
      System.err.println("reads where CIGAR goes past read end: " + broken_read_index);  // debug
    }
    if (broken_qual_index > 0) {
      System.err.println("reads where CIGAR goes past quality array end: " + broken_qual_index);  // debug
    }

  }

  private void flush_check (boolean force, SAMRecord sr) {
    if (VERBOSE) System.err.println("flush check, size=" + events_by_name.size());  // debug
    ArrayList<SplicedReadInfo> to_report = new ArrayList<SplicedReadInfo>();
    if (force) {
      to_report.addAll(events_by_name.values());
    } else {
      int current_start = sr.getAlignmentStart();
      if (VERBOSE) System.err.println("current start="+current_start);  // debug

      ArrayList<SplicedReadInfo> all = new ArrayList<SplicedReadInfo>();
      ArrayList<SplicedReadInfo> candidates = new ArrayList<SplicedReadInfo>();
      
      for (SplicedReadInfo sri : events_by_name.values()) {
	all.add(sri);
	if (current_start > sri.segment_2_start) {
	  //	  System.err.println("FLUSH 1");  // debug
	  candidates.add(sri);
	}
      }

      if (SORT) {
	for (SplicedReadInfo candidate : candidates) {
	  int overlaps = 0;
	  //	System.err.println("candidate=" + candidate.get_name() + " cs=" + current_start);  // debug

	  for (SplicedReadInfo sri : all) {
	    if (candidate.overlaps(sri)) {
	      //	    System.err.println(" overlap with " + sri.get_name());  // debug
	      overlaps++;
	    } else {
	      //	    System.err.println(" NO overlap with " + sri.get_name());  // debug
	    }
	  }
	  //	System.err.println("  overlaps="+overlaps);  // debug

	  if (overlaps > 1) {
	    // overlaps a still-pending span
	    // allow one overlap (for this record)
	    //	  System.err.println("KEEP2!");  // debug
	  } else {
	    // flushable
	    to_report.add(candidate);
	  }
	}
      } else {
	to_report.addAll(candidates);
      }
    }

    Collections.sort(to_report);

    if (false) {
      if (to_report.size() > 0) {
	System.err.println("sorted block:");  // debug
	for (SplicedReadInfo sri : to_report) {
	  System.err.println("  " + sri.get_name() + " count=" + sri.get_count());  // debug

	}
      }
    }

    boolean usable;
    boolean known;
    for (SplicedReadInfo sri : to_report) {
      int count = sri.get_count();
      usable = count >= minimum_observations_to_report;

      if (usable && sri.get_span_size() < MIN_SPAN) usable = false;
      // minimum span size

      HashSet<UCSCRefGene> rgs = null;
      if (usable && ic != null) {
	// if a junction file is specified
	rgs = ic.find_exon_junction(sri);
	known = rgs != null && rgs.size() > 0;
	if (ANNOTATE_MODE) {
	  usable = true;
	} else {
	  usable = novel_only ? !known : known;
	  // only junctions which are either (a) novel or (b) match a reference
	}
      }

      if (usable) {
	if (WRITE_BED) {
	  write_bed(sri, rgs);
	} else if (WRITE_DELIMITED) {
	  write_delimited(sri, rgs);
	} else {
	  write_counts(sri);
	}
      }

      events_by_name.remove(sri.get_name());
    }
  }

  public static void main (String[] argv) {
    SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
    // STFU

    SplicedReadReporter srr = new SplicedReadReporter();

    try {
      UCSCRefGeneReader rg = null;
      String strand_filter = null;

      for (int i=0; i < argv.length; i++) {
	if (argv[i].equals("-bam")) {
	  srr.set_bam(new File(argv[++i]));
	} else if (argv[i].equals("-refflat")) {
	  rg = new UCSCRefGeneReader(argv[++i]);
	} else if (argv[i].equals("-strand")) {
	  strand_filter = argv[++i];
	  if (rg == null) {
	    System.err.println("ERROR: specify -refflat before -strand");  // debug
	    System.exit(1);
	  } else if (strand_filter.equals("+") ||
		     strand_filter.equals("-")) {
	    // OK
	  } else {
	    System.err.println("ERROR: strand must be + or -");  // debug
	    System.exit(1);
	  }
	} else if (argv[i].equals("-t2g")) {
	  // transcript->gene mapping file (e.g. for ENSEMBL)
	  if (rg == null) {
	    System.err.println("ERROR: specify -refflat before -t2g");  // debug
	    System.exit(1);
	  } else {
	    rg.set_transcript2gene(argv[++i]);
	  }
	} else if (argv[i].equals("-min-reads")) {
	  srr.set_minimum_observations_to_report(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-v")) {
	  srr.VERBOSE = true;
	} else if (argv[i].equals("-bed")) {
	  srr.WRITE_BED = true;
	} else if (argv[i].equals("-no-header")) {
	  srr.set_write_header(false);
	} else if (argv[i].equals("-annotate")) {
	  srr.set_annotate_mode(true);
	} else if (argv[i].equals("-of")) {
	  srr.set_outfile(argv[++i]);
	} else if (argv[i].equals("-no-sort")) {
	  srr.SORT = false;
	} else if (argv[i].equals("-reference-only")) {
	  // extract only junctions matching a reference
	  srr.set_novel_only(false);
	  srr.WRITE_DELIMITED = true;
	  // change default output format to tab-delimited
	} else if (argv[i].equals("-rgb-known")) {
	  srr.set_known_rgb(argv[++i]);
	} else if (argv[i].equals("-rgb-novel")) {
	  srr.set_novel_rgb(argv[++i]);
	} else if (argv[i].equals("-2bit")) {
	  srr.set_reference_sequence(new TwoBitFile(argv[++i]));
	} else if (argv[i].equals("-fasta")) {
	  srr.set_reference_sequence(new FASTAIndexedFAI(argv[++i]));
	} else if (argv[i].equals("-restrict-read-name")) {
	  // restrict processing to a specific read name
	  srr.set_restrict_read_name(argv[++i]);
	} else if (argv[i].equals("-restrict-start")) {
	  // currently only affects read report
	  srr.set_restrict_junction_start(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-restrict-end")) {
	  // currently only affects read report
	  srr.set_restrict_junction_end(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-name")) {
	  // single region query: specify reference name
	  // This must now be specified BEFORE -start and -end.
	  srr.add_query_region().tname = new String(argv[++i]);
	  // create a new region
	} else if (argv[i].equals("-start")) {
	  // query region start
	  srr.get_query_region().set_start(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-end")) {
	  // query region end
	  srr.get_query_region().set_end(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-region")) {
	  // specify a region to query (may be used multiple times)
	  srr.add_query_region().parse(argv[++i]);
	} else if (argv[i].equals("-read-report")) {
	  // write read-level report
	  srr.set_read_report(argv[++i]);
	} else if (argv[i].equals("-primary-only")) {
	  srr.set_primary_read_mode(true);
	} else if (argv[i].equals("-no-duplicates")) {
	  srr.set_exclude_duplicates(true);
	} else if (argv[i].equals("-min-span")) {
	  srr.set_min_span(Integer.parseInt(argv[++i]));
	} else if (argv[i].equals("-ignore-incompatible")) {
	  srr.set_ignore_incompatible(true);
	} else {
	  System.err.println("ERROR: unknown parameter " + argv[i]);  // debug
	  System.exit(1);
	}
      }

      if (rg != null) {
	IntronCache ic;
	if (strand_filter == null) {
	  ic = new IntronCache(rg);
	} else {
	  ic = new IntronCache(rg, strand_filter);
	}
	srr.set_intron_cache(ic);
      }

      srr.report();
    } catch (Exception e) {
      System.err.println("ERROR: " + e);  // debug
      e.printStackTrace();
    }
  }

  private SplicedReadInfo track_site (String ref_name, int seq1_end, int seq2_start, SAMRecord sr) {
    if (VERBOSE) System.err.println("tracking " + ref_name + " " + seq1_end + " " + seq2_start + " read="+sr.getReadName());  // debug

    SplicedReadInfo sri = new SplicedReadInfo(ref_name, seq1_end, seq2_start);
    String name = sri.get_name();
    // bleh; better way to do this?
    if (events_by_name.containsKey(name)) {
      sri = events_by_name.get(name);
    } else {
      events_by_name.put(name, sri);
    }
    sri.increment_counter();

    if (sr.getReadNegativeStrandFlag()) {
      sri.increment_minus();
    } else {
      sri.increment_plus();
    }

    return sri;
  }

  private void write_delimited (SplicedReadInfo sri, HashSet<UCSCRefGene> rgs) {
    if (!wrote_header) {
      ArrayList<String> fields = new ArrayList<String>();
      fields.add("junction");
      fields.add("count");
      if (ANNOTATE_MODE) fields.add("type");
      if (!novel_only) {
	fields.add("genes");
	fields.add("transcripts");
      }
      if (ANNOTATE_MODE) {
	fields.add("qc_flanking");
	fields.add("qc_plus");
	fields.add("qc_minus");
	fields.add("qc_perfect_reads");
	fields.add("qc_clean_reads");
      }
      ps.println(Str.join("\t", fields));
      wrote_header = true;
    }

    ArrayList<String> fields = new ArrayList<String>();
    fields.add(format_junction(sri));
    fields.add(Integer.toString(sri.counter));
    if (!novel_only) {
      //
      //  reference only genes: report exact gene hit
      //
      if (rgs == null || rgs.size() == 0) {
	if (ANNOTATE_MODE) fields.add("novel");
	fields.add("");
	fields.add("");
      } else {
	//
	//  TO DO: sort UCSCRefGene entries by transcript name!
	//  output is currently stochastic which makes diffs difficult
	//

	ArrayList<String> genes = new ArrayList<String>();
	ArrayList<String> transcripts = new ArrayList<String>();
	for (UCSCRefGene rg : rgs) {
	  if (rg.name2 == null || rg.name2.length() == 0) {
	    genes.add("");
	  } else {
	    genes.add(rg.name2);
	  }
	  transcripts.add(rg.name);
	}
	//	if (ANNOTATE_MODE) fields.add(transcripts.size() > 0 ? "known" : "novel");
	if (ANNOTATE_MODE) fields.add("known");
	fields.add(Str.join(",", genes));
	fields.add(Str.join(",", transcripts));
	// creates duplicate entries for the same gene,
	// but each index is mappable to the transcript column 
      }
      if (ANNOTATE_MODE) {
	//	fields.add(sri.get_flanking_qc_flag() ? "1" : "0");
	fields.add(Integer.toString(sri.counter_flanking_qc));
	//	fields.add(sri.counter_plus > 0 ? "1" : "0");
	fields.add(Integer.toString(sri.counter_plus));
	//	fields.add(sri.counter_minus > 0 ? "1" : "0");
	fields.add(Integer.toString(sri.counter_minus));
	fields.add(Integer.toString(sri.counter_perfect));
	fields.add(Integer.toString(sri.counter_clean));
      }

    }
    ps.println(Str.join("\t", fields));
  }

  private String format_junction (SplicedReadInfo sri) {
    return sri.reference_name + ":" + sri.segment_1_end + ":+," +
      //      sri.reference_name + ":" + sri.segment_2_start + ":+,";
      sri.reference_name + ":" + sri.segment_2_start + ":+";
    // not sure why that trailing comma was in there  :/
  }

  private void write_counts (SplicedReadInfo sri) {
    ps.println(
	       format_junction(sri) + 
	       "\t" +
	       sri.counter
	       );
  }

  private void write_bed (SplicedReadInfo sri, HashSet<UCSCRefGene> rgs) {
    // .bed format
    // http://genome.ucsc.edu/FAQ/FAQformat.html#format1
    if (!wrote_header) {
      ps.println("track name=junctions description=\"BAM junctions\" visibility=3 itemRgb=\"On\"");
      // 3 = pack (so label is displayed immediately to left of each junction)
      wrote_header = true;
    }
    ArrayList<String> fields = new ArrayList<String>();
    fields.add(sri.reference_name);
    // 1. reference/chrom name
    fields.add(Integer.toString(sri.segment_1_end - 1));
    // 2. start base number (0-based)
    fields.add(Integer.toString(sri.segment_2_start));
    // 3. end base number: base is not counted, so don't need to adjust  :/

    String gene_name = null;
    if (rgs != null) {
      for (UCSCRefGene rg : rgs) {
	gene_name = rg.name2;
	if (gene_name != null) break;
      }
    }
    if (gene_name == null) gene_name = "junction";

    //    String name = sri.counter + " " + gene_name + "." + junction_counter.add(gene_name);
    // fail:
    // - spaces in name verboten
    // - double-quoting name no helpen

    //    String name = sri.counter + "_" + gene_name + "." + junction_counter.add(gene_name);
    // too cluttery

    String name = Integer.toString(sri.counter);
    // just the count

    //    fields.add(String.format("JUNC%08d", ++junction_counter));
    fields.add(name);
    // 4. name (e.g. JUNC00000001)

    fields.add(Integer.toString(sri.counter));
    // 5. score
    fields.add("+");
    // 6. strand
    fields.add(fields.get(1));
    fields.add(fields.get(2));
    // 7. thickStart
    // 8. thickEnd
    String color;
    if (rgs == null || rgs.size() == 0) {
      color = RGB_NOVEL;
    } else {
      color = RGB_KNOWN;
    }
    fields.add(color);
    // 9. itemRGB

    if (false) {
      // don't use (yet)
      fields.add("1");
      // 10. blockCount
      fields.add(Integer.toString(sri.segment_2_start - sri.segment_1_end));
      // 11. blockSizes: ???
      fields.add("0");
      // 12. blockStarts
    }

    ps.println(Str.join("\t", fields));
  }

  private boolean get_flank_qc (SplicedReadFlankingInfo fi, boolean is_read_edge) {

    return fi.count_aligned_bases >= (is_read_edge ?
				      QC_MIN_FLANKING_NT_READ_EDGE : QC_MIN_FLANKING_NT_INTERNAL);
    // if the sequence flanking the junction abuts a read edge,
    // use a more stringent check
  }


}
