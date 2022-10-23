package riemann;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.TreeSet;

public class GramInfo implements Comparable<GramInfo> {
   public final BigDecimal grampt;
   double zeta;
   double[] cos = new double[10];
   double[] sin = new double[9];
   double der;
   public final BigDecimal idx;
   Riemann.RZpoint point;

   public GramInfo(BigDecimal grampt, long idx) {
      this.grampt = grampt;
      this.idx = BigDecimal.valueOf(idx);
   }

   public GramInfo(BigDecimal grampt, BigDecimal idx) {
      this.grampt = grampt;
      this.idx = idx;
   }

   public static GramInfo[] RZGram(long n1, int count, MathContext mc) {
      return RZGram(BigDecimal.valueOf(n1), count, mc);
   }

   /**
    * find  sequential gram points beginning with the one just preceding the first parameter.
    * @param idx1
    * @param count
    * @param mc
    * @return
    */
   public static GramInfo[] RZGram(BigDecimal idx1, int count, MathContext mc) {
      BigDecimal idx = idx1.divideToIntegralValue(BigDecimal.ONE);
       GramInfo[] grampts = new GramInfo[count];
       BigDecimal target = Gram.pi.multiply(idx, mc);
      BigDecimal t2 = target.divide(Gram.bdTWO, mc);
      BigDecimal arg1 = Gram.sqrt(t2.divide(Gram.pi, mc), mc, 1.0E-15);
      Gram.initLogVals(arg1.intValue());
      BigDecimal ti = Gram.bdTWO.multiply(target.add(Gram.pi8,mc).divide(
            Gram.log(arg1,mc).subtract(BigDecimal.ONE, mc),mc),mc);
      for (int i = 0; i<count; i++) {
         for(int j = 1; j<15; j++) {
            t2 = ti.divide(Gram.bdTWO, mc);
            arg1 = Gram.sqrt(t2.divide(Gram.pi, mc), mc, 1.0E-15);
            BigDecimal logsqrtArg = Gram.log(arg1,mc);
            BigDecimal term = BigDecimal.ONE.divide(ti.multiply(BigDecimal.valueOf(48),mc), mc);
            BigDecimal thetai = ti.multiply(logsqrtArg,mc).subtract(t2,mc).
                  subtract(Gram.pi8 ,mc).add(term,mc);
            BigDecimal del = target.subtract(thetai,mc);
             if(Math.abs(del.doubleValue()) < 1.0E-15){break;}
            ti = ti.add(del.divide(logsqrtArg, mc),mc);
          }
         grampts[i] = new GramInfo(ti, idx);
         if(i>0) {
            ti= ti.add(grampts[i].grampt.subtract(grampts[i-1].grampt,mc),mc);
         } else {
            ti = ti.add(new BigDecimal(2*Math.PI/Math.log(ti.doubleValue()/(2*Math.PI))));
         }
         idx = idx.add(BigDecimal.ONE);
         target = Gram.pi.multiply(idx, mc);
}
      return grampts;
   }

   Riemann.RZpoint setRZpoint(BigDecimal bdOffset) {
      if (point == null) {
         point = new Riemann.RZpoint(this, bdOffset);
      }
      return point;
   }

   static TreeSet<Riemann.RZpoint> toRZPoint(Iterable<GramInfo> riemann,
                                             BigDecimal bdOffset) {
      TreeSet<Riemann.RZpoint> gram = new TreeSet<Riemann.RZpoint>();
      for (GramInfo gramInfo : riemann) {
         gram.add(gramInfo.setRZpoint(bdOffset));
      }
      return gram;
   }

   static GramInfo fromStringArray(String[] parsed, BigDecimal bdOffset) {
      BigDecimal grampt = BigDecimal.valueOf(
            Double.parseDouble(parsed[0])).add(bdOffset, Gram.mc);
      GramInfo gramInfo = new GramInfo(grampt, Long.parseLong(parsed[22]));
      gramInfo.zeta = Double.parseDouble(parsed[1]);
      gramInfo.der = Double.parseDouble(parsed[21]);
      int offset = 2;
      for (int j = 0; j < gramInfo.cos.length; j++) {
         gramInfo.cos[j] = Double.parseDouble(parsed[j + offset]);
      }
      offset += gramInfo.cos.length;
      for (int j = 0; j < gramInfo.sin.length; j++) {
         gramInfo.sin[j] = Double.parseDouble(parsed[j + offset]);
      }
      return gramInfo;
   }

   public String toString(BigDecimal offset) {
      String ret = Riemann.nf.format(offset == null ? grampt : grampt.subtract(offset)) + " \t" + zeta
            + " \t" + Arrays.toString(cos) + " \t"
            + Arrays.toString(sin) + " \t" + der + " \t" + idx;
      ret = ret.replaceAll("[\\[\\],]", "");
      return ret;
   }

   @Override
   public String toString() {
      return "GramInfo [grampt=" + grampt + ", idx=" + idx + "]";
   }

   public int compareTo(GramInfo arg0) {
      GramInfo other = (GramInfo) arg0;
      return idx.compareTo(other.idx);
   }
}
