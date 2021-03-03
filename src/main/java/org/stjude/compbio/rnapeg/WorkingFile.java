package org.stjude.compbio.rnapeg;
// use a temporary filename until file is complete

import java.io.*;
import java.util.zip.*;

public class WorkingFile extends java.io.File {
  File final_file;
  private static final String TEMP_SUFFIX = ".tmp";
  private PrintStream ps = null;
  private boolean append_mode = false;
  private boolean scratch_dir_mode = false;

  public WorkingFile (String filename) {
    super(filename + TEMP_SUFFIX);
    final_file = new File(filename);
    append_mode = false;
    scratch_dir_mode = false;
  }

  public WorkingFile (String filename, File scratch_dir) {
    super(scratch_dir, filename + TEMP_SUFFIX);
    final_file = new File(filename);
    append_mode = false;
    scratch_dir_mode = true;
  }

  public WorkingFile (String filename, boolean append) {
    super(append ? filename : filename + TEMP_SUFFIX);
    final_file = new File(filename);
    append_mode = append;
    scratch_dir_mode = false;
  }

  private String temp_setup (String filename, boolean append) {
    String result = append ? filename : filename + TEMP_SUFFIX;
    return result;
  }

  public void finish() {
    // close stream (if opened from this class) and rename tempfile to final name
    if (ps != null) ps.close();
    if (!append_mode) {
      final_file.delete();
      if (scratch_dir_mode) {
	try {
	  NIOUtils.copyFile(this, final_file);
	  // copy from scratch to final destination
	  super.delete();
	  // delete scratch space copy
	} catch (IOException e) {
	  System.err.println("ERROR copying file: " + e);  // debug
	}
      } else {
	if (!renameTo(final_file)) {
	  System.err.println("error renaming to " + final_file.getAbsolutePath());  // debug
	}
      }
    }
  }

  public File get_file() {
    return final_file;
  }

  public boolean delete() {
    boolean v = final_file.delete();
    if (!v) System.err.println("can't delete " + final_file);  // debug
    return v;
  }

  public PrintStream getPrintStream() throws FileNotFoundException,IOException {
    if (ps == null) {
      OutputStream os = new BufferedOutputStream(new FileOutputStream(this, append_mode));
      String path = getAbsolutePath();
      //      System.err.println("path="+path);  // debug

      int gzi = path.lastIndexOf(".gz");

      if (gzi == (path.length() - (append_mode ? 0 : TEMP_SUFFIX.length()) - 3)) {
	ps = new PrintStream(new GZIPOutputStream(os));
      } else {
	ps = new PrintStream(os);
      }
    }    
    return ps;
  }

  public static void main (String[] argv) {
    try {
      if (true) {
	WorkingFile wf = new WorkingFile("some_output.txt", new File("c:/temp/"));
	PrintStream ps = wf.getPrintStream();
	ps.println("line1");
	ps.println("line2");
	ps.println("line3");
	ps.close();
	try {
	  System.err.println("sleeping...");  // debug
	  Thread.sleep(0);
	} catch (Exception e) {}
	wf.finish();
      } else {
	WorkingFile wf = new WorkingFile("some_output.txt");
	PrintStream ps = wf.getPrintStream();
	ps.println("whatever3");
	ps.close();
	wf.finish();
      }
    } catch (Exception e) {
      System.err.println("ERROR: "+e);  // debug
      e.printStackTrace();
    }
  }
  
}
