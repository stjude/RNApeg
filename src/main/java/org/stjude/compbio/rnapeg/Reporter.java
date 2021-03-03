package org.stjude.compbio.rnapeg;
// utility

import java.util.*;
import java.io.*;

public class Reporter {
  private ArrayList<String> headers = new ArrayList<String>();
  private HashMap<String,Boolean> hl = new HashMap<String,Boolean>();

  private boolean headers_printed = false;
  private String delimiter = "\t";
  private HashMap<String,String> row;
  private BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(System.out));
  private WorkingFile wf = null;
  private boolean VERBOSE = true;
  private boolean VERBOSE_DF = true;

  public Reporter() {
    reset_row();
  }

  public void set_output_filename (String fn) throws IOException {
    if (fn.equals("-")) {
      wf = null;
      bw = new BufferedWriter(new OutputStreamWriter(System.out));
    } else {
      wf = new WorkingFile(fn);
      bw = new BufferedWriter(new OutputStreamWriter(wf.getPrintStream()));
    }
  }

  public void set_headers (ArrayList<String> headers) {
    this.headers = headers;
    hl = new HashMap<String,Boolean>();
    for (String h : headers) {
      hl.put(h, Boolean.valueOf(true));
    }
  }

  public void add_header (String h) {
    headers.add(h);
    hl.put(h, Boolean.valueOf(true));
  }

  public void reset_row() {
    row = new HashMap<String,String>();
  }

  public String get_value (String k) {
    return row.get(k);
  }

  public void set_value (String k, String v) {
    if (hl.containsKey(k)) {
      row.put(k,v);
    } else {
      System.err.println("ERROR, invalid column " + k);  // debug
    }
  }

  public void end_row() throws IOException {
    int count=0;
    if (headers_printed == false) {
      for(String header : headers) {
	if (count++ > 0) bw.write(delimiter);
	bw.write(header);
      }
      bw.newLine();
      headers_printed = true;
    }

    count = 0;
    for (String h : headers) {
      String v = row.get(h);
      if (v == null) v = "";
      if (count++ > 0) bw.write(delimiter);
      bw.write(v);
    }
    bw.newLine();
    reset_row();
  }

  public void close() throws IOException {
    if (VERBOSE) {
      // 6/2012: additional diagnostics for cluster troubleshooting
      // does PrintStream checkError() report anything useful?
      if (wf != null) {
	PrintStream ps = wf.getPrintStream();
	boolean error = ps.checkError();
	System.err.println("DEBUG: PrintStream checkError(): " + error);  // debug
      }
    }
    bw.close();
    
    if (wf != null) {
      wf.finish();
      if (VERBOSE) {
	// debug: report final file size
	debug_file(wf, "temp_output_file", false);
	debug_file(wf.get_file(), "final_output_file", true);
      }
    }
  }

  private void debug_file (File f, String label, boolean do_df) {
    try {
      System.err.print("DEBUG: " + label + ": " + 
		       "abs=" + f.getAbsoluteFile() + " " +
		       "canonical=" + f.getCanonicalFile() + " " +
		       "exists?: " + f.exists() + " ");  // debug
      if (f.exists()) {
	System.err.print("size=" + f.length());  // debug

      }
      System.err.println("");  // debug
      if (do_df && VERBOSE_DF) {
	String cmd = null;

	String df = find_binary("df");
	if (df == null) {
	  System.err.println("ERROR: can't find df!");  // debug
	} else {
	  if (f.exists()) {
	    cmd = df + " " + f.getAbsoluteFile();
	  } else {
	    cmd = df + " " + f;
	  }
	  exec_command(cmd);
	}

	String ls = find_binary("ls");
	if (ls == null) {
	  System.err.println("ERROR: can't find ls!");  // debug
	} else {
	  if (f.exists()) {
	    File dir = f.getAbsoluteFile().getParentFile();
	    if (dir == null) {
	      System.err.println("ERROR: no parent");  // debug
	    } else {
	      cmd = ls + " -la " + dir.getAbsoluteFile();
	      exec_command(cmd);
	    }
	  } else {
	    System.err.println("ERROR: can't run ls (need parent)");  // debug
	  }
	}
      }
    } catch (IOException e) {
      System.err.println("ERROR: " + e);  // debug
    }
  }

  private String find_binary (String name) {
    String result = null;
    ArrayList<String> locs = new ArrayList<String>();
    locs.add("/bin/" + name);
    locs.add("/usr/bin/" + name);
    locs.add("c:/cygwin/bin/" + name + ".exe");
    // cygwin
    for (String loc : locs) {
      File exe = new File(loc);
      if (exe.exists()) {
	result = exe.toString();
	break;
      }
    }
    return result;
  }

  private void exec_command (String cmd) {
    // FIX ME: replace w/Funk.Sys.exec_command()
    System.err.println("DEBUG: executing " + cmd);  // debug
    try {
      Runtime rt = Runtime.getRuntime();
      Process p = rt.exec(cmd);
      InputStream in = p.getInputStream();
      BufferedInputStream buf = new BufferedInputStream(in);
      InputStreamReader inread = new InputStreamReader(buf);
      BufferedReader bufferedreader = new BufferedReader(inread);

      String line;
      while ((line = bufferedreader.readLine()) != null) {
	System.err.println(line);
      }
    } catch (IOException e) {
      System.err.println("ERROR during exec: " + e);  // debug
    }
  }

  

}
