// string and stream utility functions

package org.stjude.compbio.rnapeg;

import java.util.*;
import java.io.*;

public class Str {

  private static final String valid_number_characters = "0123456789";
  // hack for parseInt, MUST BE SORTED

  public static String basename (String s) {
    int len = s.length();
    int i;
    int last = -1;
    for (i=0; i < len; i++) {
      if (s.charAt(i) == '/') last = i;
    }
    if (last == -1) {
      return s;
    } else {
      return s.substring(last + 1);
    }
  }

  public static String join (Enumeration e) {
      return join(",", e);
  }

  public static String join (String delim, Collection c) {
    int count = 0;
    StringBuffer sb = new StringBuffer();
    for (Object o : c) {
      if (count++ > 0) sb.append(delim);
      //      sb.append((String) e.nextElement());
      sb.append(o.toString());
    }
    return sb.toString();
  }

  public static String join (String delim, Enumeration e) {
    int count = 0;
    StringBuffer sb = new StringBuffer();
    while (e.hasMoreElements()) {
      if (count++ > 0) sb.append(delim);
      //      sb.append((String) e.nextElement());
      sb.append(((Object) e.nextElement()).toString());
    }
    return sb.toString();
  }

  public static String join (String delim, Iterator i) {
    int count = 0;
    StringBuffer sb = new StringBuffer();
    while (i.hasNext()) {
      if (count++ > 0) sb.append(delim);
      sb.append(i.next().toString());
    }
    return sb.toString();
  }

  public static String join (String delim, String[] list) {
    return join(delim, list, 0);
  }

  public static String join (String delim, String[] list, int max_fields) {
    StringBuffer sb = new StringBuffer();
    int count = 0;
    int i;
    int max = list.length;
    if (max_fields > 0) max = max_fields;
    for (i=0; i<max; i++) {
      if (count++ > 0) sb.append(delim);
      sb.append(list[i]);
    }
    return sb.toString();
  }

  public static String double_print (double d, int places) {
    // double_print(3.14159, 2) => 3.14
    // like sprintf "%.2f" w/o rounding
    String result = (new Double(d)).toString();
    int max_len = result.indexOf(".") + places + 1;
    if (result.length() > max_len) {
      result = result.substring(0, max_len);
      // my kingdom for sprintf()  :/
    }
    return result;
  }

  public static boolean is_genbank_accession (String token) {
    // does the given string look like a GenBank accession number?
    int len = token.length();
    for (int i = 0; i < len; i++) {
      char c = token.charAt(i);
      if (Character.isLetter(c) && Character.isUpperCase(c) && i <= 1) {
	// OK if first and/or second character are uppercase letters
      } else if (Character.isDigit(c)) {
	// rest of string is digits
	if (i < 1) return false;
      } else if (c == '_') {
	if (i != 2) return false;
	// OK in 3rd position, eg NM_000927
      } else {
	return false;
      }
    }
    return true;
  }

  public static Vector read_lines (DataInputStream d) {
    Vector lines = new Vector();
    String line;

    InputStreamReader isr = new InputStreamReader(d);
    // horrible backwards-compatibility hack
    BufferedReader d2 = new BufferedReader(isr);

    try {
      while ((line = d2.readLine()) != null) {
	lines.addElement(line);
      }
    } catch (IOException e) {
      System.out.println("eek: " + e);  // debug
    }
    return lines;
  }

  public static String unpadded (String s) {
    StringBuffer sb = new StringBuffer(s);
    StringBuffer unpadded = new StringBuffer();
    for (int i=0; i < sb.length(); i++) {
      char c = sb.charAt(i);
      if (c != '*') {
	unpadded.append(c);
      }
    }
    return unpadded.toString();
  }

  public static String reverse_complement (String s) {
    StringBuffer rc = (new StringBuffer(s)).reverse();
    int len = rc.length();
    char c;
    for (int i=0; i < len; i++) {
      switch (rc.charAt(i)) {
      case 'a':
	c = 't'; break;
      case 'A':
	c = 'T'; break;
      case 'c':
	c = 'g'; break;
      case 'C':
	c = 'G'; break;
      case 'g':
	c = 'c'; break;
      case 'G':
	c = 'C'; break;
      case 't':
	c = 'a'; break;
      case 'T':
	c = 'A'; break;
      case '*':
	c = '*'; break;
      case 'n':
	c = 'n'; break;
      case 'N':
	c = 'N'; break;
      case 'x':
	c = 'x'; break;
      case 'X':
	c = 'X'; break;
      default:
	c = rc.charAt(i);
	System.out.println("reverse: warning, don't know how to rc " + c);
	break;
      }
      rc.setCharAt(i, c);
    }
    return rc.toString();
  }

    public static int count_instances_of_char (String s, char c) {
	// return how many times character c appears in string s.
	int count = 0;
	int last = 0;
	int index;

	while (true) {
	    index = s.indexOf(c, last);
	    if (index == -1) {
		break;
	    } else {
		count++;
		last = index + 1;
	    }
	}
	return(count);
    }

   public static void main (String argv[]) {
     // debug, hack
       //     System.out.println("enter:");  // debug
       //     Vector v = read_lines(new DataInputStream(System.in));
       //     System.out.println(" ");  // debug
       //     System.out.println(v.elementAt(0));  // debug

       //       System.out.println(count_instances_of_char("abcabcabz", 'q'));  // debug

       //       String s = "		a b c   ";
       //       System.out.println("->" + s + "<-");  // debug
       //       System.out.println("->" + trim_whitespace(s) + "<-");  // debug

       //       System.out.println(unquote("\"somestring\""));  // debug
       //       System.out.println(unquote("'someotherstring'"));  // debug

       // Vector v = FileToVector("Str.java");
       // System.err.println(v.size());  // debug

//        Enumeration e = tokenize_string("molid|100590|complex||node|27|AS|Calpain1", "|");
//        while (e.hasMoreElements()) {
// 	   String token = (String) e.nextElement();
// 	   System.err.println(token);  // debug
//        }

//     System.err.println(parseInt("-124521"));  // debug

     System.err.println(comma_format(900000));  // debug
     System.err.println(comma_format(1000000));  // debug
   }

  public static String comma_format (int v) {
    String vs = Integer.toString(v);
    int len = vs.length();
    int count = 0;
    StringBuffer sb = new StringBuffer();
    int i;
    boolean comma_needed = false;
    for (i=len - 1; i >= 0; i--) {
      if (comma_needed) {
	sb.append(',');
	comma_needed = false;
      }
      sb.append(vs.charAt(i));
      if (++count % 3 == 0) comma_needed=true;
    }
    return sb.reverse().toString();
  }

  public static int parseInt (String s) throws NumberFormatException {
    //
    //  alternative to Integer.parseInt(), which chews up a tremendous
    //  amount of memory during heavy flatfile parsing
    //
    int result = 0;
    int len = s.length();
    char c;
    int i;
    int index;
    boolean is_negative = false;
    if (s == null) {
      throw new NumberFormatException("null string");
    } else if (len == 0) {
      throw new NumberFormatException("empty string");
    } else {
      for(i=0; i < len; i++) {
	c = s.charAt(i);
	index = valid_number_characters.indexOf(c);
	if (index == -1) {
	  // not a number
	  if (i == 0) {
	    // allow + or - if first character
	    if (c == '-') {
	      is_negative = true;
	    } else {
	      if (c != '+') 
		throw new NumberFormatException("invalid character " + c);
	    }
	  } else {
	    throw new NumberFormatException("invalid character " + c);
	  }
	} else {
	  // a number
	  if (result > 0) result *= 10;
	  result += index;
	}
      }
    }

    if (is_negative) result *= -1;
    
    return result;
  }

    public static String trim_whitespace (String s) {
	// remove leading and trailing whitespace
	int len = s.length();
	int start_index, end_index;
	char c;

	for(start_index=0; start_index < len; start_index++) {
	    c = s.charAt(start_index);
	    if (c != ' ' && c != '\t') break;
	}

	for(end_index=len - 1; end_index >= 0; end_index--) {
	    c = s.charAt(end_index);
	    if (c != ' ' && c != '\t') break;
	}
	return s.substring(start_index,end_index + 1);
    }

    public static String unquote(String s) {
	int first = s.indexOf('"');
	int last = s.lastIndexOf('"');
	int len = s.length();

	if (first == 0 && last == len - 1) {
	    // match for doublequotes
	    return(s.substring(first + 1, last));
	} else {
	    first = s.indexOf('\'');
	    last = s.lastIndexOf('\'');
	    if (first == 0 && last == len - 1) {
		// match for singlequotes
		return(s.substring(first + 1, last));
	    }
	}
	
	return(s);
    }

    public static Hashtable delimited_string_to_hash(String s, String delimiter) {
	Hashtable result = new Hashtable();
	//	StringTokenizer st = new StringTokenizer(s, delimiter);
	Enumeration e = tokenize_string(s, delimiter);
	String key, value;
	//	while (st.hasMoreTokens()) {
	// key = st.nextToken();
	while (e.hasMoreElements()) {
	    key = (String) e.nextElement();
	    //	    System.err.println("ds2h:key = " + key);  // debug

	    //	    if (st.hasMoreTokens()) {
	    //		value = st.nextToken();
	    if (e.hasMoreElements()) {
		value = (String) e.nextElement();
		//		System.err.println("ds2h:value = " + value);  // debug
		result.put(key, value);
	    } else {
		System.err.println("WARNING: no value for key " + key);  // debug
		break;
	    }
	}
	return(result);
    }

    public static Vector FileToVector (String filename) {
	Vector lines = null;
	try {
	    BufferedReader br = new BufferedReader(new FileReader(filename));
	    lines = new Vector();
	    String line;
	    while (true) {
		line = br.readLine();
		if (line == null) break;
		lines.addElement(line);
	    }
	} catch (java.io.FileNotFoundException e) {
	    System.err.println("can't find file " + filename + " !");  // debug
	} catch (java.io.IOException e) {
	    System.err.println("i/o exception!");  // debug
	}
	return lines;
    }

    public static Enumeration tokenize_string (String thing, String token) {
	// kind of like StringTokenizer, but supports null tokens,
	// e.g. a key/value hash dump such as:
	//   molid|100590|complex||node|27|AS|Calpain1
	
	int delim_len = token.length();
	if (delim_len > 1) {
	    die("delim_len > 1 untested");
	}

	Vector tokens = new Vector();

	Vector stops = new Vector();
	stops.addElement(new Integer(0));
	int i = thing.indexOf(token);
	while (i > -1) {
	    //	    System.err.println(i);  // debug
	    stops.addElement(new Integer(i));
	    i = thing.indexOf(token, i + 1);
	}
	stops.addElement(new Integer(thing.length()));

	int slen = stops.size();
	int hither,yon;
	for (i=0; i < slen - 1; i++) {
	    hither = ((Integer) stops.elementAt(i)).intValue();
	    if (i > 0) hither += delim_len;
	    yon = ((Integer) stops.elementAt(i + 1)).intValue();
	    //	    System.err.println(hither + " -> " + yon);  // debug
	    //	    System.err.println("token: " + thing.substring(hither,yon));  // debug
	    tokens.addElement(thing.substring(hither, yon));
	}
	
	return tokens.elements();
    }

    private static void die (String msg) {
	System.err.println(msg);
	System.exit(1);
    }

  public static String url_query_string(Hashtable params) {
    // generate a simple "query" string of a URI/URL from the given parameters.
    // must be encoded with URI (URLEncoder is not appropriate)
    int count=0;
    StringBuffer sb = new StringBuffer();
    Enumeration e = params.keys();
    while (e.hasMoreElements()) {
      String key = (String) e.nextElement();
      String value = (String) params.get(key);
      if (count++ > 0) sb.append("&");
      sb.append(key + "=" + value);
    }
    return sb.toString();
  }

  public static String float_decimal_format (double d, int places) {
    return String.format("%." + places + "f", d);
  }

  public static String float_decimal_format (float f, int places) {
    return String.format("%." + places + "f", f);
  }




}
