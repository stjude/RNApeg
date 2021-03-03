package org.stjude.compbio.rnapeg;

import java.io.*;
import java.util.zip.*;

public class FileUtils {

  public static BufferedReader getUniversalBufferedReader(String fn) throws FileNotFoundException,IOException {
    // transparently open regular or .gz files
    return getUniversalBufferedReader(new File(fn));
  }

  public static BufferedReader getUniversalBufferedReader(File f) throws FileNotFoundException,IOException {
    // transparently open regular or .gz files
    InputStream is = new FileInputStream(f);
    String fn = f.getName();
    if (fn.toLowerCase().indexOf(".gz") == fn.length() - 3) is = new GZIPInputStream(is);
    return new BufferedReader(new InputStreamReader(is));
  }

}
